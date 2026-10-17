#!/usr/bin/env python3
"""
00_quantify_real_osdr_rnaseq.py
Quantify every raw read of NASA OSDR study OSD-528 (GLDS-528) with kallisto, paired-end.

Steps
  1. Download the Mycobacterium marinum M RefSeq assembly GCF_000018345.1 (CDS FASTA and
     GenBank flat file) and check both against NCBI's md5checksums.txt.
  2. Build a kallisto index over every annotated CDS (chromosome NC_010612.1 and plasmid)
     plus the genome's rRNA, tRNA and other non-coding RNA genes (taken from the GenBank
     file). The non-coding genes act as decoys: 88-95% of read pairs in these libraries are
     rRNA (suggesting little or no rRNA depletion), and without decoys some rRNA reads could
     pseudoalign to coding genes. Only CDS rows are written to the count matrices.
  3. Download the 18 raw FASTQs listed by the OSDR files API and check each against the
     deposited GLDS-528_rna-seq_raw_md5sum.txt.
  4. Resolve the RFPRPM4 deposit defect (below).
  5. Infer library strandedness on a read subset, then quantify each sample paired-end.
  6. Write the estimated-count matrix, the TPM matrix and per-sample QC, including the
     share of read pairs assigned to rRNA, to tRNA/other non-coding RNA and to CDS.

The RFPRPM4 deposit defect
  The deposited RFPRPM4 R1 and R2 files are each the complete RFPRPM41 file followed by a
  second library with a different index. This script proves the prefix is RFPRPM41 by
  hashing it against RFPRPM41's deposited md5, then quantifies the remainder as
  "RFPRPM4". That remainder is *inferred* to be the true RFPRPM4 library: OSDR does not say
  so, and it should be confirmed with the depositors. Downstream scripts read
  config/samples.tsv to decide whether RFPRPM4 is included.

Read and reference files are large and git-ignored (data/raw/fastqs, data/reference,
data/raw/kallisto_cds_ncrna). The script is idempotent: verified downloads are reused.
"""

import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import urllib.request
from concurrent.futures import ThreadPoolExecutor

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.path.join(PROJECT_DIR, "data", "raw", "fastqs")
REF_DIR = os.path.join(PROJECT_DIR, "data", "reference")
QUANT_DIR = os.path.join(PROJECT_DIR, "data", "raw", "kallisto_cds_ncrna")
PROCESSED_DIR = os.path.join(PROJECT_DIR, "data", "processed")

USER_AGENT = "osd528-reanalysis/1.0"
KALLISTO = shutil.which("kallisto") or "/opt/homebrew/bin/kallisto"
THREADS = max(1, (os.cpu_count() or 4) - 2)

ASSEMBLY = "GCF_000018345.1_ASM1834v1"
NCBI_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/018/345/" + ASSEMBLY
REF_FILES = [f"{ASSEMBLY}_cds_from_genomic.fna.gz", f"{ASSEMBLY}_genomic.gbff.gz"]

OSDR_FILES_API = "https://osdr.nasa.gov/osdr/data/osd/files/528"
OSDR_BASE = "https://osdr.nasa.gov"
MD5_FILE = "GLDS-528_rna-seq_raw_md5sum.txt"

SAMPLES = ["RFP3D11", "RFP3D39", "RFP3D47", "RFPNG14", "RFPNG35", "RFPNG45",
           "RFPRPM4", "RFPRPM41", "RFPRPM6"]
INDEX = os.path.join(REF_DIR, f"{ASSEMBLY}_cds_ncrna.kallisto.idx")
NCRNA_FASTA = os.path.join(REF_DIR, f"{ASSEMBLY}_ncrna.fna")
NCRNA_TYPES = ("rRNA", "tRNA", "ncRNA", "tmRNA")


def fetch(url, dest):
    """Download url to dest via a .part file; raise on any failure."""
    tmp = dest + ".part"
    cmd = ["curl", "-sSfL", "--retry", "5", "--retry-delay", "5", "-A", USER_AGENT, "-o", tmp, url]
    subprocess.run(cmd, check=True)
    os.replace(tmp, dest)


def fetch_text(url):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read().decode()


def md5_of(path, limit=None):
    """md5 of a file, or of its first `limit` bytes."""
    h = hashlib.md5()
    remaining = limit
    with open(path, "rb") as f:
        while remaining is None or remaining > 0:
            chunk = f.read(1 << 22 if remaining is None else min(1 << 22, remaining))
            if not chunk:
                break
            h.update(chunk)
            if remaining is not None:
                remaining -= len(chunk)
    if remaining:
        raise ValueError(f"{path} is shorter than the {limit} bytes requested")
    return h.hexdigest()


def get_reference():
    os.makedirs(REF_DIR, exist_ok=True)
    sums = {}
    for line in fetch_text(NCBI_BASE + "/md5checksums.txt").splitlines():
        digest, name = line.split()
        sums[name.lstrip("./")] = digest
    for name in REF_FILES:
        dest = os.path.join(REF_DIR, name)
        if not (os.path.exists(dest) and md5_of(dest) == sums[name]):
            print(f"Downloading {name}")
            fetch(f"{NCBI_BASE}/{name}", dest)
        if md5_of(dest) != sums[name]:
            raise RuntimeError(f"md5 mismatch for {name}")
        print(f"  {name}: md5 verified against NCBI")
    return os.path.join(REF_DIR, REF_FILES[0])


def parse_cds_annotations(cds_gz):
    """CDS FASTA headers -> target_id: locus tag, gene symbol, protein, replicon."""
    ann = {}
    with gzip.open(cds_gz, "rt") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            tid = line[1:].split()[0]

            def field(k):
                m = re.search(rf"\[{k}=([^\]]+)\]", line)
                return m.group(1) if m else ""

            locus = field("locus_tag") or tid
            ann[tid] = {
                "gene_id": locus,
                "gene_symbol": field("gene") or locus,
                "protein": field("protein"),
                "protein_id": field("protein_id"),
                "pseudo": field("pseudo") == "true",
                "replicon": tid.split("|")[1].split("_cds")[0] if "|" in tid else "",
            }
    return ann


def write_ncrna_fasta(gbff_gz):
    """Extract rRNA, tRNA, ncRNA and tmRNA genes from the GenBank file as decoy targets.

    Target ids are "<type>|<locus_tag>|<product>" so they can never collide with CDS ids.
    Returns {target_id: type}.
    """
    comp = str.maketrans("acgtn", "tgcan")
    targets, records, cur, mode = {}, [], None, None
    with gzip.open(gbff_gz, "rt") as f:
        for line in f:
            if line.startswith("LOCUS"):
                cur = {"features": [], "seq": []}
                records.append(cur)
                mode = "head"
            elif line.startswith("FEATURES"):
                mode = "feat"
            elif line.startswith("ORIGIN"):
                mode = "seq"
            elif line.startswith("//"):
                mode = None
            elif mode == "feat":
                cur["features"].append(line.rstrip("\n"))
            elif mode == "seq":
                cur["seq"].append(re.sub(r"[^acgtn]", "", line.lower()))
    with open(NCRNA_FASTA, "w") as out:
        for rec in records:
            seq = "".join(rec["seq"])
            feats, feat = [], None
            for line in rec["features"]:
                m = re.match(r"^     (\S+)\s+(\S.*)$", line)
                if m:
                    feat = {"type": m.group(1), "loc": m.group(2), "qual": ""}
                    feats.append(feat)
                elif feat:
                    feat["qual"] += line.strip() + " "
            for feat in feats:
                if feat["type"] not in NCRNA_TYPES:
                    continue
                if "join" in feat["loc"] or "order" in feat["loc"]:
                    raise RuntimeError(f"Unexpected split location for a non-coding gene: {feat['loc']}")
                pos = [int(x) for x in re.findall(r"\d+", feat["loc"])]
                sub = seq[min(pos) - 1:max(pos)]
                if feat["loc"].startswith("complement"):
                    sub = sub.translate(comp)[::-1]
                tag = re.search(r'/locus_tag="([^"]+)"', feat["qual"])
                prod = re.search(r'/product="([^"]+)"', feat["qual"])
                tid = (f"{feat['type']}|{tag.group(1) if tag else len(targets)}|"
                       f"{(prod.group(1) if prod else '').replace(' ', '_')}")
                out.write(f">{tid}\n{sub}\n")
                targets[tid] = feat["type"]
    return targets


def build_index(cds_gz):
    """Index the CDS plus the non-coding RNA decoys; return {ncRNA target_id: type}."""
    ncrna = write_ncrna_fasta(os.path.join(REF_DIR, REF_FILES[1]))
    counts = {t: sum(v == t for v in ncrna.values()) for t in NCRNA_TYPES}
    print(f"  Non-coding RNA decoys from the GenBank file: {counts}")
    if not os.path.exists(INDEX):
        print("Building kallisto index (CDS + non-coding RNA decoys)")
        combined = INDEX + ".fa"
        with open(combined, "wb") as out:
            with gzip.open(cds_gz, "rb") as f:
                shutil.copyfileobj(f, out)
            with open(NCRNA_FASTA, "rb") as f:
                shutil.copyfileobj(f, out)
        subprocess.run([KALLISTO, "index", "-t", str(THREADS), "-i", INDEX, combined], check=True)
        os.remove(combined)
    return ncrna


def get_fastqs():
    os.makedirs(RAW_DIR, exist_ok=True)
    listing = json.loads(fetch_text(OSDR_FILES_API))
    files = {f["file_name"]: f for f in listing["studies"]["OSD-528"]["study_files"]}
    md5_path = os.path.join(RAW_DIR, MD5_FILE)
    if not os.path.exists(md5_path):
        fetch(OSDR_BASE + files[MD5_FILE]["remote_url"], md5_path)
    # The deposited md5 list names files without the "GLDS-528_rna-seq_" prefix.
    expected = {}
    for line in open(md5_path):
        parts = line.split()
        if len(parts) == 2:
            expected["GLDS-528_rna-seq_" + os.path.basename(parts[1])] = parts[0]

    names = [f"GLDS-528_rna-seq_{s}_{r}_raw.fastq.gz" for s in SAMPLES for r in ("R1", "R2")]
    missing = [n for n in names if n not in files or n not in expected]
    if missing:
        raise RuntimeError(f"Not listed by OSDR or missing from the md5 list: {missing}")

    def one(name):
        dest = os.path.join(RAW_DIR, name)
        size = files[name]["file_size"]
        if os.path.exists(dest) and os.path.getsize(dest) == size and md5_of(dest) == expected[name]:
            return f"  {name}: reused, md5 verified"
        fetch(OSDR_BASE + files[name]["remote_url"], dest)
        if os.path.getsize(dest) != size:
            raise RuntimeError(f"{name}: size {os.path.getsize(dest)} != OSDR {size}")
        if md5_of(dest) != expected[name]:
            raise RuntimeError(f"{name}: md5 does not match the deposited md5 list")
        return f"  {name}: downloaded, {size:,} bytes, md5 verified"

    with ThreadPoolExecutor(max_workers=4) as pool:
        for msg in pool.map(one, names):
            print(msg, flush=True)
    return {n: os.path.join(RAW_DIR, n) for n in names}, expected


def split_rpm4(paths, expected):
    """Prove RFPRPM4 = RFPRPM41 + another library, and write that library on its own."""
    out, report = {}, {}
    for r in ("R1", "R2"):
        rpm4 = paths[f"GLDS-528_rna-seq_RFPRPM4_{r}_raw.fastq.gz"]
        rpm41_name = f"GLDS-528_rna-seq_RFPRPM41_{r}_raw.fastq.gz"
        n41 = os.path.getsize(paths[rpm41_name])
        prefix_md5 = md5_of(rpm4, limit=n41)
        if prefix_md5 != expected[rpm41_name]:
            raise RuntimeError(f"RFPRPM4 {r}: first {n41:,} bytes are not RFPRPM41 "
                               f"({prefix_md5}); the deposit has changed, revisit this step")
        with open(rpm4, "rb") as f:
            f.seek(n41)
            magic = f.read(2)
        if magic != b"\x1f\x8b":
            raise RuntimeError(f"RFPRPM4 {r}: no gzip member starts at byte {n41:,}")
        tail = os.path.join(RAW_DIR, f"GLDS-528_rna-seq_RFPRPM4_{r}_second_library.fastq.gz")
        tail_size = os.path.getsize(rpm4) - n41
        if not (os.path.exists(tail) and os.path.getsize(tail) == tail_size):
            with open(rpm4, "rb") as src, open(tail + ".part", "wb") as dst:
                src.seek(n41)
                shutil.copyfileobj(src, dst, 1 << 22)
            os.replace(tail + ".part", tail)
        out[r] = tail
        report[r] = {"rpm41_prefix_bytes": n41, "rpm41_prefix_md5_matches": True,
                     "second_library_bytes": tail_size}
        print(f"  RFPRPM4 {r}: first {n41:,} bytes are byte-identical to RFPRPM41; "
              f"second library written ({tail_size:,} bytes)")
    return out, report


def first_read_index(path):
    with gzip.open(path, "rt") as f:
        return f.readline().strip().split(":")[-1]


def subsample(path, dest, n_reads):
    with gzip.open(path, "rt") as src, open(dest, "w") as dst:
        for i, line in enumerate(src):
            if i >= 4 * n_reads:
                break
            dst.write(line)


def kallisto_quant(r1, r2, out_dir, strand_flag):
    os.makedirs(out_dir, exist_ok=True)
    cmd = [KALLISTO, "quant", "-i", INDEX, "-o", out_dir, "-t", str(THREADS)]
    if strand_flag:
        cmd.append(strand_flag)
    subprocess.run(cmd + [r1, r2], check=True, capture_output=True, text=True)
    with open(os.path.join(out_dir, "run_info.json")) as f:
        return json.load(f)


def infer_strandedness(r1, r2):
    """Pseudoalign 200,000 pairs three ways and pick the library type they support."""
    tmp = os.path.join(QUANT_DIR, "_strand_test")
    os.makedirs(tmp, exist_ok=True)
    s1, s2 = os.path.join(tmp, "R1.fq"), os.path.join(tmp, "R2.fq")
    subsample(r1, s1, 200_000)
    subsample(r2, s2, 200_000)
    counts = {}
    for label, flag in (("unstranded", None), ("fr", "--fr-stranded"), ("rf", "--rf-stranded")):
        counts[label] = kallisto_quant(s1, s2, os.path.join(tmp, label), flag)["n_pseudoaligned"]
    shutil.rmtree(tmp)
    un = counts["unstranded"]
    if counts["rf"] >= 0.9 * un and counts["fr"] <= 0.3 * un:
        choice = ("rf", "--rf-stranded")
    elif counts["fr"] >= 0.9 * un and counts["rf"] <= 0.3 * un:
        choice = ("fr", "--fr-stranded")
    else:
        choice = ("unstranded", None)
    print(f"  Strandedness test (200,000 pairs): {counts} -> {choice[0]}")
    return choice, counts


def main():
    print("=== Phase 0: full paired-end quantification of NASA OSDR OSD-528 ===")
    cds_gz = get_reference()
    annotations = parse_cds_annotations(cds_gz)
    print(f"  {len(annotations):,} CDS records in {ASSEMBLY} "
          f"({sum(a['pseudo'] for a in annotations.values())} pseudogenes)")
    ncrna = build_index(cds_gz)

    paths, expected = get_fastqs()
    rpm4_tail, rpm4_report = split_rpm4(paths, expected)

    reads = {}
    for s in SAMPLES:
        if s == "RFPRPM4":
            reads[s] = (rpm4_tail["R1"], rpm4_tail["R2"])
        else:
            reads[s] = tuple(paths[f"GLDS-528_rna-seq_{s}_{r}_raw.fastq.gz"] for r in ("R1", "R2"))
    indices = {s: first_read_index(reads[s][0]) for s in SAMPLES}
    if len(set(indices.values())) != len(SAMPLES):
        raise RuntimeError(f"Library indices are not unique across samples: {indices}")
    print(f"  Library index of each sample's first read is unique: {indices}")

    (strand_label, strand_flag), strand_counts = infer_strandedness(*reads["RFPNG14"])

    version = subprocess.run([KALLISTO, "version"], capture_output=True, text=True).stdout.strip()
    qc = {"reference_assembly": ASSEMBLY, "kallisto_version": version,
          "strandedness": strand_label, "strandedness_test_pseudoaligned": strand_counts,
          "index_targets": {"cds": len(annotations), "ncrna_decoys": len(ncrna)},
          "rfprpm4_deposit_defect": rpm4_report, "samples": {}}
    abundances = {}
    for s in SAMPLES:
        out = os.path.join(QUANT_DIR, s)
        info_path = os.path.join(out, "run_info.json")
        if os.path.exists(info_path) and os.path.exists(os.path.join(out, "abundance.tsv")):
            info = json.load(open(info_path))
        else:
            print(f"  kallisto quant {s} ({strand_label})", flush=True)
            info = kallisto_quant(*reads[s], out, strand_flag)
        qc["samples"][s] = {
            "library_index": indices[s],
            "fastq_r1": os.path.basename(reads[s][0]),
            "fastq_r2": os.path.basename(reads[s][1]),
            "read_pairs_processed": info["n_processed"],
            "read_pairs_pseudoaligned": info["n_pseudoaligned"],
            "percent_pseudoaligned": info["p_pseudoaligned"],
        }
        with open(os.path.join(out, "abundance.tsv")) as f:
            abundances[s] = {r["target_id"]: r for r in csv.DictReader(f, delimiter="\t")}
        assigned = {"rRNA": 0.0, "other_ncRNA": 0.0, "CDS": 0.0}
        for t, row in abundances[s].items():
            kind = ncrna.get(t)
            key = "CDS" if kind is None else "rRNA" if kind == "rRNA" else "other_ncRNA"
            assigned[key] += float(row["est_counts"])
        qc["samples"][s]["percent_of_pairs_assigned"] = {
            k: round(100 * v / info["n_processed"], 2) for k, v in assigned.items()}
        print(f"  {s}: {info['n_processed']:,} pairs, {info['n_pseudoaligned']:,} "
              f"pseudoaligned ({info['p_pseudoaligned']}%); assigned "
              f"{qc['samples'][s]['percent_of_pairs_assigned']}", flush=True)

    targets = [t for t in abundances[SAMPLES[0]] if t in annotations]
    if len(targets) != len(annotations):
        raise RuntimeError(f"{len(annotations) - len(targets)} CDS missing from kallisto output")
    os.makedirs(PROCESSED_DIR, exist_ok=True)
    for fname, column, fmt in (("osd528_counts_raw.tsv", "est_counts", "{:.2f}"),
                               ("osd528_tpm.tsv", "tpm", "{:.4f}")):
        with open(os.path.join(PROCESSED_DIR, fname), "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["target_id", "gene_id", "gene_symbol", "replicon", "pseudo"]
                       + SAMPLES + ["protein"])
            for t in targets:
                a = annotations[t]
                w.writerow([t, a["gene_id"], a["gene_symbol"], a["replicon"],
                            "yes" if a["pseudo"] else "no"]
                           + [fmt.format(float(abundances[s][t][column])) for s in SAMPLES]
                           + [a["protein"]])
    with open(os.path.join(PROCESSED_DIR, "osd528_sequencing_qc_stats.json"), "w") as f:
        json.dump(qc, f, indent=2)
    print(f"Wrote counts and TPM for {len(targets):,} CDS x {len(SAMPLES)} samples, and QC.")


if __name__ == "__main__":
    main()
