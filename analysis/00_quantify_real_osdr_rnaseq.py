#!/usr/bin/env python3
"""
00_quantify_real_osdr_rnaseq.py
Direct empirical quantification of real NASA OSDR OSD-528 RNA-seq sequencing reads
using kallisto v0.52.0 against the Mycobacterium marinum NC_010612.1 reference transcriptome.
Uses pure Python standard library (csv, math, json, subprocess) for zero external dependencies.
"""

import os
import sys
import json
import re
import csv
import math
import subprocess

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF_FASTA = os.path.join(PROJECT_ROOT, "data", "reference", "M_marinum_NC_010612_cds.fasta")
KALLISTO_IDX = os.path.join(PROJECT_ROOT, "data", "reference", "m_marinum_kallisto.idx")
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")
TMP_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "fastq_slices")

os.makedirs(PROCESSED_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)

# Sample definitions and corresponding R1 FASTQ files in OSDR
SAMPLES = {
    "RFP3D11": "GLDS-528_rna-seq_RFP3D11_R1_raw.fastq.gz",
    "RFP3D39": "GLDS-528_rna-seq_RFP3D39_R1_raw.fastq.gz",
    "RFP3D47": "GLDS-528_rna-seq_RFP3D47_R1_raw.fastq.gz",
    "RFPNG14": "GLDS-528_rna-seq_RFPNG14_R1_raw.fastq.gz",
    "RFPNG35": "GLDS-528_rna-seq_RFPNG35_R1_raw.fastq.gz",
    "RFPNG45": "GLDS-528_rna-seq_RFPNG45_R1_raw.fastq.gz",
    "RFPRPM4": "GLDS-528_rna-seq_RFPRPM4_R1_raw.fastq.gz",
    "RFPRPM41": "GLDS-528_rna-seq_RFPRPM41_R1_raw.fastq.gz",
    "RFPRPM6": "GLDS-528_rna-seq_RFPRPM6_R1_raw.fastq.gz",
}

# Download slice size: 60 MB per sample gives ~2.5 million NextSeq reads per sample
SLICE_BYTES = 62914560  # 60 MB


def parse_reference_annotations(fasta_path):
    """Parse CDS FASTA headers to extract locus_tag, gene symbol, and protein description."""
    annotations = {}
    print(f"Parsing gene annotations from {fasta_path}...")
    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                tid = line[1:].split()[0]
                g_match = re.search(r"\[gene=([^\]]+)\]", line)
                l_match = re.search(r"\[locus_tag=([^\]]+)\]", line)
                p_match = re.search(r"\[protein=([^\]]+)\]", line)

                gene = g_match.group(1) if g_match else ""
                locus = l_match.group(1) if l_match else tid
                protein = p_match.group(1) if p_match else ""
                annotations[tid] = {
                    "gene_id": locus,
                    "gene_symbol": gene if gene else locus,
                    "protein": protein,
                }
    print(f"Loaded annotations for {len(annotations):,} coding sequences.")
    return annotations


def download_fastq_slice(sample_id, fastq_filename, out_path, slice_bytes=SLICE_BYTES):
    """Download a real slice of the FASTQ directly from OSDR S3 redirect with range request."""
    if os.path.exists(out_path) and os.path.getsize(out_path) >= slice_bytes * 0.8:
        print(f"  [{sample_id}] Reusing existing cached slice ({os.path.getsize(out_path):,} bytes).")
        return True

    osdr_url = f"https://osdr.nasa.gov/geode-py/ws/studies/OSD-528/download?source=datamanager&file={fastq_filename}"
    print(f"  [{sample_id}] Downloading {slice_bytes/(1024*1024):.0f} MB slice of real NextSeq reads from OSDR S3...")

    # Use curl with -L and -r to follow S3 redirect and grab byte slice
    cmd = [
        "curl",
        "-sSL",
        "-r",
        f"0-{slice_bytes}",
        osdr_url,
        "-o",
        out_path,
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not os.path.exists(out_path) or os.path.getsize(out_path) < 10000:
        print(f"  ERROR downloading {sample_id}: {res.stderr}")
        return False

    print(f"  [{sample_id}] Downloaded {os.path.getsize(out_path):,} bytes successfully.")
    return True


def run_kallisto_quant(sample_id, fastq_path, out_dir):
    """Run kallisto quant on the real sequencing reads."""
    os.makedirs(out_dir, exist_ok=True)
    kallisto_bin = "/opt/homebrew/bin/kallisto"
    cmd = [
        kallisto_bin,
        "quant",
        "-i",
        KALLISTO_IDX,
        "-o",
        out_dir,
        "--single",
        "-l",
        "99",
        "-s",
        "10",
        fastq_path,
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"  ERROR in kallisto quant for {sample_id}: {res.stderr}")
        return None

    # Parse stdout/stderr for processed and pseudoaligned read counts
    stats = {}
    for line in res.stderr.splitlines():
        if "processed" in line and "pseudoaligned" in line:
            m = re.search(r"processed ([\d,]+) reads, ([\d,]+) reads pseudoaligned", line)
            if m:
                stats["processed_reads"] = int(m.group(1).replace(",", ""))
                stats["aligned_reads"] = int(m.group(2).replace(",", ""))
    abundance_file = os.path.join(out_dir, "abundance.tsv")
    if not os.path.exists(abundance_file):
        return None
    return abundance_file, stats


def main():
    print("=== Phase 0: Real RNA-Seq Quantification for NASA OSDR OSD-528 ===")
    annotations = parse_reference_annotations(REF_FASTA)

    sample_abundances = {}
    sample_stats = {}
    all_targets = []

    sample_order = list(SAMPLES.keys())

    for sample_id in sample_order:
        fastq_name = SAMPLES[sample_id]
        slice_path = os.path.join(TMP_DIR, f"{sample_id}_R1.slice.fastq.gz")
        ok = download_fastq_slice(sample_id, fastq_name, slice_path)
        if not ok:
            print(f"Failed to obtain real reads for {sample_id}.")
            sys.exit(1)

        quant_out = os.path.join(TMP_DIR, f"kallisto_{sample_id}")
        res = run_kallisto_quant(sample_id, slice_path, quant_out)
        if res is None:
            print(f"Kallisto quantification failed for {sample_id}.")
            sys.exit(1)

        ab_file, stats = res
        sample_stats[sample_id] = stats
        print(f"  [{sample_id}] Processed {stats.get('processed_reads', 0):,} reads | {stats.get('aligned_reads', 0):,} pseudoaligned to M. marinum CDS.")

        # Read abundance.tsv
        counts = {}
        with open(ab_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                tid = row["target_id"]
                counts[tid] = float(row["est_counts"])
                if not all_targets:
                    # Capture order of targets from first sample
                    pass
        if not all_targets:
            all_targets = list(counts.keys())
        sample_abundances[sample_id] = counts

    # Calculate library sizes (sum of estimated counts)
    lib_sizes = {s: sum(sample_abundances[s].values()) for s in sample_order}
    print("\nLibrary Sizes (Mapped mRNA reads):")
    for s in sample_order:
        print(f"  {s}: {lib_sizes[s]:,.1f} mapped counts")

    # Save raw counts matrix
    raw_path = os.path.join(PROCESSED_DIR, "osd528_counts_raw.tsv")
    raw_header = ["target_id", "gene_id", "gene_symbol"] + sample_order + ["protein"]
    with open(raw_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(raw_header)
        for tid in all_targets:
            ann = annotations.get(tid, {})
            gid = ann.get("gene_id", tid)
            sym = ann.get("gene_symbol", gid)
            prot = ann.get("protein", "")
            counts_row = [f"{sample_abundances[s].get(tid, 0.0):.2f}" for s in sample_order]
            writer.writerow([tid, gid, sym] + counts_row + [prot])

    print(f"Saved raw empirical counts ({len(all_targets):,} genes x {len(sample_order)} samples) to {raw_path}")

    # Compute Normalized Expression: Log2(CPM + 1)
    norm_path = os.path.join(PROCESSED_DIR, "osd528_counts_normalized.tsv")
    norm_header = ["gene_id", "gene_symbol"] + sample_order + ["protein"]
    with open(norm_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(norm_header)
        for tid in all_targets:
            ann = annotations.get(tid, {})
            gid = ann.get("gene_id", tid)
            sym = ann.get("gene_symbol", gid)
            prot = ann.get("protein", "")

            cpm_vals = []
            for s in sample_order:
                cnt = sample_abundances[s].get(tid, 0.0)
                lsize = lib_sizes[s]
                cpm = (cnt / (lsize + 1e-9)) * 1e6
                log2_val = math.log2(cpm + 1.0)
                cpm_vals.append(f"{log2_val:.4f}")

            writer.writerow([gid, sym] + cpm_vals + [prot])

    print(f"Saved normalized expression matrix to {norm_path}")

    # Save QC statistics
    stats_path = os.path.join(PROCESSED_DIR, "osd528_sequencing_qc_stats.json")
    with open(stats_path, "w") as f:
        json.dump(sample_stats, f, indent=2)
    print(f"Saved sequencing QC stats to {stats_path}")
    print("Real empirical RNA-seq quantification complete.")


if __name__ == "__main__":
    main()
