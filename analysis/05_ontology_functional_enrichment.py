#!/usr/bin/env python3
"""
05_ontology_functional_enrichment.py

Gene Ontology (GO) over-representation analysis for the OSD-528 re-analysis
(Mycobacterium marinum strain M; 3D clinostat, RPM 2.0 and static 1g).

Method
------
1. Annotations. UniProtKB GO cross-references for the M. marinum strain M
   reference proteome (UP000001190; ATCC BAA-535; NCBI taxon 216594). UniProt
   names genes by the original ordered locus names (MMAR_xxxx). The DE tables
   use RefSeq locus tags (MMAR_RSxxxxx), so the two are joined through the
   /old_locus_tag qualifiers of the RefSeq GBFF (GCF_000018345.1) that the
   quantification step already uses. GO ids that the ontology release lists
   as alt_ids are mapped to their primary id. Obsolete ids are mapped to their
   replaced_by term when there is exactly one, and are otherwise dropped and
   counted.
2. Ontology. go-basic.obo. Each annotation is propagated to all of its is_a
   and part_of ancestors within the same GO namespace (true-path rule).
3. Test. One-sided hypergeometric test (over-representation) per GO term,
   run separately for BP, MF and CC. For each gene set, the universe is the
   genes in that set's background that carry at least one annotation in the
   ontology being tested. Terms with MIN_TERM_SIZE..MAX_TERM_SIZE annotated
   universe genes are tested. BH FDR is computed over every tested term within
   each gene set and ontology. Terms with no overlap stay in the BH family
   with p = 1 but are not written out.
4. Gene sets.
     - Per DESeq2 contrast: up and down DEGs (padj < DEG_PADJ, split by the
       sign of log2FoldChange). Background = genes with a non-missing padj in
       that contrast.
     - Shared simulated-microgravity response (deg_shared_simulated_microgravity.tsv),
       split by direction. Background = genes with a non-missing padj in both
       simulator-vs-static-1g contrasts.
     - Every WGCNA module except grey. Background = every gene in the WGCNA
       assignment table, grey included.

Reproducibility
---------------
The locus-to-GO snapshot (go_annotations_mmarinum.tsv), the GO subgraph needed
to propagate it (go_ontology_subset.tsv) and their provenance
(go_annotations_source.json) are written to data/processed/. By default the
script reuses that snapshot, so a fresh clone reproduces the results with no
network access. Raw downloads are cached in data/reference/go/ (git-ignored).
Outputs are sorted and contain no run timestamps, so re-running the script on
the same inputs gives byte-identical files.

Usage
-----
  python3 analysis/05_ontology_functional_enrichment.py
      Run the enrichment from the snapshot. Builds the snapshot first if it is missing.
  python3 analysis/05_ontology_functional_enrichment.py --rebuild-annotations
      Rebuild the snapshot from the cached raw downloads. Downloads them if absent.
  python3 analysis/05_ontology_functional_enrichment.py --refresh
      Re-download UniProt and go-basic.obo, then rebuild the snapshot.

Outputs (data/processed/)
-------------------------
  go_annotations_mmarinum.tsv   locus -> GO snapshot (direct annotations, with evidence codes)
  go_ontology_subset.tsv        GO terms reachable from the snapshot, with is_a / part_of parents
  go_annotations_source.json    source, query, releases, download dates, checksums, mapping counts
  go_enrichment_deg.tsv         results for the DEG and shared-response gene sets
  go_enrichment_modules.tsv     results for the WGCNA modules
  go_summary.json               method, coverage, per-gene-set significant counts and top terms, caveats
"""

from __future__ import annotations

import argparse
import collections
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import re
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PROCESSED = os.path.join(PROJECT_DIR, "data", "processed")
DATA_REFERENCE = os.path.join(PROJECT_DIR, "data", "reference")
GO_CACHE = os.path.join(DATA_REFERENCE, "go")

# ---------------------------------------------------------------- parameters
USER_AGENT = "osd528-reanalysis/1.0"
PROTEOME_ID = "UP000001190"
TAXON_ID = 216594
UNIPROT_URL = (
    "https://rest.uniprot.org/uniprotkb/stream?query=proteome:%s"
    "&fields=accession,gene_oln,gene_names,go_id&format=json" % PROTEOME_ID
)
GO_OBO_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"
GBFF_FILE = os.path.join(DATA_REFERENCE, "GCF_000018345.1_ASM1834v1_genomic.gbff.gz")

UNIPROT_RAW = os.path.join(GO_CACHE, "uniprot_%s_go.json" % PROTEOME_ID)
OBO_RAW = os.path.join(GO_CACHE, "go-basic.obo")

SNAP_ANNOT = os.path.join(DATA_PROCESSED, "go_annotations_mmarinum.tsv")
SNAP_ONTO = os.path.join(DATA_PROCESSED, "go_ontology_subset.tsv")
SNAP_SOURCE = os.path.join(DATA_PROCESSED, "go_annotations_source.json")
OUT_DEG = os.path.join(DATA_PROCESSED, "go_enrichment_deg.tsv")
OUT_MOD = os.path.join(DATA_PROCESSED, "go_enrichment_modules.tsv")
OUT_SUMMARY = os.path.join(DATA_PROCESSED, "go_summary.json")

CONTRASTS = [
    # (group name, file)
    ("clinostat_vs_static1g", "deg_3dclinostat_vs_static1g.tsv"),
    ("rpm2_vs_static1g", "deg_rpm2_vs_static1g.tsv"),
    ("clinostat_vs_rpm2", "deg_3dclinostat_vs_rpm2.tsv"),
]
SHARED_FILE = "deg_shared_simulated_microgravity.tsv"
SHARED_BACKGROUND = ("clinostat_vs_static1g", "rpm2_vs_static1g")
WGCNA_FILE = "wgcna_module_assignments.tsv"
WGCNA_UNASSIGNED = "grey"

DEG_PADJ = 0.05          # DEG threshold on DESeq2 BH-adjusted p
FDR_ALPHA = 0.05         # significance threshold for GO terms
MIN_TERM_SIZE = 5        # annotated universe genes per term, inclusive
MAX_TERM_SIZE = 500
TOP_N = 5                # terms listed per gene set in go_summary.json
SMALL_SET = 10           # gene sets with fewer annotated genes are flagged as low power

ONTOLOGIES = ["BP", "MF", "CC"]
NAMESPACE_CODE = {
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC",
}

RESULT_COLUMNS = [
    "gene_set", "group", "subset", "ontology", "go_id", "term_name",
    "term_size", "overlap", "set_size", "universe_size", "expected",
    "fold_enrichment", "pvalue", "fdr", "n_terms_tested",
    "fdr_bh_hit_terms_only", "overlap_gene_ids", "overlap_gene_symbols",
]


# ------------------------------------------------------------------ helpers
def log(msg: str) -> None:
    print(msg, flush=True)


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path: str) -> str:
    return os.path.relpath(path, PROJECT_DIR)


def sig(x: float, digits: int = 4) -> float:
    """Round to a fixed number of significant digits for JSON output."""
    return float("%.*e" % (digits - 1, x))


def write_json(path: str, obj) -> None:
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, sort_keys=False)
        f.write("\n")


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (monotone, capped at 1)."""
    p = np.asarray(p, dtype=float)
    m = p.size
    if m == 0:
        return p
    order = np.argsort(p, kind="mergesort")
    scaled = p[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(scaled[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1.0)
    return out


# ---------------------------------------------------------------- downloads
def download(url: str, dest: str, keep_headers: list[str]) -> dict:
    """Stream url to dest (atomic rename) and write a sidecar .meta.json."""
    import requests  # imported here so offline snapshot runs never need it

    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    last_err = None
    for attempt in range(1, 4):
        try:
            with requests.get(url, headers={"User-Agent": USER_AGENT},
                              stream=True, timeout=(30, 600)) as r:
                r.raise_for_status()
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
                headers = {k: r.headers.get(k) for k in keep_headers if r.headers.get(k)}
                final_url = r.url
            break
        except Exception as e:  # network errors: retry with a fixed backoff
            last_err = e
            log("  download attempt %d failed: %s" % (attempt, e))
            time.sleep(5 * attempt)
    else:
        raise SystemExit("Could not download %s: %s" % (url, last_err))
    os.replace(tmp, dest)
    meta = {
        "url": url,
        "final_url": final_url,
        "downloaded_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "headers": headers,
        "bytes": os.path.getsize(dest),
        "sha256": sha256(dest),
    }
    write_json(dest + ".meta.json", meta)
    log("  saved %s (%d bytes)" % (rel(dest), meta["bytes"]))
    return meta


def ensure_raw(refresh: bool) -> tuple[dict, dict]:
    if refresh or not (os.path.exists(UNIPROT_RAW) and os.path.exists(UNIPROT_RAW + ".meta.json")):
        log("Downloading UniProtKB GO annotations for %s ..." % PROTEOME_ID)
        download(UNIPROT_URL, UNIPROT_RAW,
                 ["X-UniProt-Release", "X-UniProt-Release-Date",
                  "X-API-Deployment-Date", "Date"])
    if refresh or not (os.path.exists(OBO_RAW) and os.path.exists(OBO_RAW + ".meta.json")):
        log("Downloading go-basic.obo ...")
        download(GO_OBO_URL, OBO_RAW, ["Last-Modified", "ETag", "Date"])
    with open(UNIPROT_RAW + ".meta.json") as f:
        up_meta = json.load(f)
    with open(OBO_RAW + ".meta.json") as f:
        obo_meta = json.load(f)
    if up_meta["sha256"] != sha256(UNIPROT_RAW) or obo_meta["sha256"] != sha256(OBO_RAW):
        raise SystemExit("Cached GO downloads do not match their .meta.json checksums; "
                         "re-run with --refresh.")
    return up_meta, obo_meta


# ------------------------------------------------------------------ parsers
def parse_obo(path: str) -> tuple[dict, dict]:
    """Parse go-basic.obo. Returns (header, terms)."""
    header: dict = {}
    terms: dict = {}
    cur = None
    in_header = True
    with open(path, encoding="utf-8") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith("["):
                in_header = False
                if cur is not None and cur.get("id"):
                    terms[cur["id"]] = cur
                cur = ({"is_a": [], "part_of": [], "alt_id": [], "replaced_by": [],
                        "obsolete": False} if line == "[Term]" else None)
                continue
            if not line or line.startswith("!"):
                continue
            key, _, val = line.partition(": ")
            if in_header:
                if key in ("format-version", "data-version", "date"):
                    header[key] = val
                continue
            if cur is None:
                continue
            val = val.split(" ! ")[0].strip()
            if key == "id":
                cur["id"] = val
            elif key == "name":
                cur["name"] = val
            elif key == "namespace":
                cur["namespace"] = val
            elif key == "is_a":
                cur["is_a"].append(val.split()[0])
            elif key == "relationship":
                parts = val.split()
                if parts[0] == "part_of":
                    cur["part_of"].append(parts[1])
            elif key == "alt_id":
                cur["alt_id"].append(val)
            elif key == "is_obsolete":
                cur["obsolete"] = val == "true"
            elif key == "replaced_by":
                cur["replaced_by"].append(val)
    if cur is not None and cur.get("id"):
        terms[cur["id"]] = cur
    return header, terms


def parse_gbff_old_tags(path: str) -> tuple[dict, dict]:
    """Map old locus tags (MMAR_xxxx) to RefSeq locus tags (MMAR_RSxxxxx)
    from the gene features of a GenBank flat file. Returns (old2new, stats)."""
    old2new: dict = collections.defaultdict(set)
    loci = set()
    loci_with_old = set()
    in_features = False
    cur_key = None
    cur_qual: list = []

    def flush():
        if cur_key != "gene":
            return
        q = "".join(cur_qual)
        m = re.search(r'/locus_tag="([^"]+)"', q)
        if not m:
            return
        lt = m.group(1)
        loci.add(lt)
        for old in re.findall(r'/old_locus_tag="([^"]+)"', q):
            old2new[old].add(lt)
            loci_with_old.add(lt)

    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("FEATURES"):
                in_features = True
                cur_key, cur_qual = None, []
                continue
            if not in_features:
                continue
            if not line.startswith(" "):  # ORIGIN, CONTIG, // ... end of feature table
                flush()
                in_features = False
                cur_key, cur_qual = None, []
                continue
            if len(line) > 5 and line[5] != " ":
                flush()
                cur_key = line[5:21].strip()
                cur_qual = []
            else:
                cur_qual.append(line[21:].rstrip("\n") + "\n")
    stats = {"gene_features": len(loci), "gene_features_with_old_locus_tag": len(loci_with_old),
             "old_tags_mapping_to_several_loci": sum(len(v) > 1 for v in old2new.values())}
    return {k: sorted(v) for k, v in old2new.items()}, stats


# ---------------------------------------------------------- snapshot builder
def build_snapshot(refresh: bool) -> None:
    up_meta, obo_meta = ensure_raw(refresh)
    if not os.path.exists(GBFF_FILE):
        raise SystemExit("Missing %s; run analysis/00_quantify_real_osdr_rnaseq.py first "
                         "(it downloads the RefSeq reference)." % rel(GBFF_FILE))

    log("Parsing go-basic.obo ...")
    obo_header, terms = parse_obo(OBO_RAW)
    alt2primary = {}
    for tid, t in terms.items():
        for a in t["alt_id"]:
            alt2primary[a] = tid

    def resolve(go_id: str):
        """Return (current id or None, how)."""
        if go_id in terms and not terms[go_id]["obsolete"]:
            return go_id, "current"
        if go_id in alt2primary and not terms[alt2primary[go_id]]["obsolete"]:
            return alt2primary[go_id], "alt_id"
        if go_id in terms and terms[go_id]["obsolete"]:
            rep = [r for r in terms[go_id]["replaced_by"]
                   if r in terms and not terms[r]["obsolete"]]
            if len(rep) == 1:
                return rep[0], "replaced_by"
            return None, "obsolete"
        return None, "unknown"

    log("Parsing RefSeq GBFF old_locus_tag map ...")
    old2new, gbff_stats = parse_gbff_old_tags(GBFF_FILE)

    log("Parsing UniProtKB entries ...")
    with open(UNIPROT_RAW) as f:
        entries = json.load(f)["results"]
    rows = set()
    n_with_go = 0
    entries_mapped = 0
    oln_unmapped = set()
    resolve_counts = collections.Counter()
    dropped_ids = collections.Counter()
    for e in entries:
        acc = e["primaryAccession"]
        olns = sorted({o["value"] for g in e.get("genes", [])
                       for o in g.get("orderedLocusNames", [])})
        gos = []
        for x in e.get("uniProtKBCrossReferences", []):
            if x.get("database") != "GO":
                continue
            props = {p["key"]: p["value"] for p in x.get("properties", [])}
            gos.append((x["id"], props.get("GoEvidenceType", "")))
        if gos:
            n_with_go += 1
        targets = [(o, lt) for o in olns for lt in old2new.get(o, [])]
        oln_unmapped.update(o for o in olns if o not in old2new)
        if targets:
            entries_mapped += 1
        for go_raw, evid in gos:
            go_id, how = resolve(go_raw)
            resolve_counts[how] += 1
            if go_id is None:
                dropped_ids[go_raw] += 1
                continue
            t = terms[go_id]
            for old, lt in targets:
                rows.add((lt, old, acc, NAMESPACE_CODE[t["namespace"]], go_id,
                          t["name"], go_raw, evid))

    rows = sorted(rows)
    # evidence codes over the locus-level annotations actually used (snapshot rows)
    evidence_counts = collections.Counter(r[7].split(":")[0] or "NA" for r in rows)
    with open(SNAP_ANNOT, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["gene_id", "old_locus_tag", "uniprot_accession", "ontology", "go_id",
                    "go_term_name", "go_id_in_uniprot", "evidence"])
        w.writerows(rows)
    log("  wrote %s (%d rows)" % (rel(SNAP_ANNOT), len(rows)))

    # GO subgraph: every term reachable from an annotated term via is_a/part_of
    # within its own namespace.
    seen = set()
    stack = sorted({r[4] for r in rows})
    while stack:
        tid = stack.pop()
        if tid in seen:
            continue
        seen.add(tid)
        t = terms[tid]
        for p in t["is_a"] + t["part_of"]:
            if p in terms and terms[p]["namespace"] == t["namespace"] and p not in seen:
                stack.append(p)
    with open(SNAP_ONTO, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["go_id", "ontology", "term_name", "is_a", "part_of"])
        for tid in sorted(seen):
            t = terms[tid]
            same = lambda ps: ";".join(sorted(p for p in ps if p in terms and
                                              terms[p]["namespace"] == t["namespace"]))
            w.writerow([tid, NAMESPACE_CODE[t["namespace"]], t["name"],
                        same(t["is_a"]), same(t["part_of"])])
    log("  wrote %s (%d terms)" % (rel(SNAP_ONTO), len(seen)))

    n_ann = sum(evidence_counts.values())
    source = {
        "annotation_source": {
            "database": "UniProtKB",
            "proteome": PROTEOME_ID,
            "organism": "Mycobacterium marinum strain M (ATCC BAA-535)",
            "ncbi_taxon": TAXON_ID,
            "query_url": up_meta["url"],
            "uniprot_release": up_meta["headers"].get("X-UniProt-Release"),
            "uniprot_release_date": up_meta["headers"].get("X-UniProt-Release-Date"),
            "api_deployment_date": up_meta["headers"].get("X-API-Deployment-Date"),
            "downloaded_utc": up_meta["downloaded_utc"],
            "raw_file": rel(UNIPROT_RAW),
            "raw_sha256": up_meta["sha256"],
            "entries": len(entries),
            "entries_with_go": n_with_go,
            "go_cross_references": sum(resolve_counts.values()),
            "go_id_resolution": dict(sorted(resolve_counts.items())),
            "go_ids_dropped_obsolete_or_unknown": dict(sorted(dropped_ids.items())),
            "evidence_code_counts_in_snapshot_rows": dict(sorted(evidence_counts.items(),
                                                key=lambda kv: (-kv[1], kv[0]))),
            "fraction_electronic_iea": sig(evidence_counts.get("IEA", 0) / n_ann) if n_ann else None,
        },
        "locus_mapping": {
            "method": ("UniProt ordered locus name (MMAR_xxxx) -> RefSeq /old_locus_tag on gene "
                       "features -> RefSeq /locus_tag (MMAR_RSxxxxx)"),
            "gbff_file": rel(GBFF_FILE),
            "gbff_sha256": sha256(GBFF_FILE),
            **gbff_stats,
            "uniprot_entries_mapped_to_a_locus": entries_mapped,
            "uniprot_ordered_locus_names_without_refseq_locus": len(oln_unmapped),
            "refseq_loci_with_go": len({r[0] for r in rows}),
        },
        "ontology_source": {
            "url": obo_meta["url"],
            "final_url": obo_meta.get("final_url"),
            "data_version": obo_header.get("data-version"),
            "format_version": obo_header.get("format-version"),
            "downloaded_utc": obo_meta["downloaded_utc"],
            "raw_file": rel(OBO_RAW),
            "raw_sha256": obo_meta["sha256"],
            "terms_in_file": len(terms),
            "relations_propagated": ["is_a", "part_of"],
            "subset_terms_saved": len(seen),
        },
        "snapshot_files": {
            "annotations": rel(SNAP_ANNOT),
            "annotations_sha256": sha256(SNAP_ANNOT),
            "ontology_subset": rel(SNAP_ONTO),
            "ontology_subset_sha256": sha256(SNAP_ONTO),
        },
    }
    write_json(SNAP_SOURCE, source)
    log("  wrote %s" % rel(SNAP_SOURCE))


# ------------------------------------------------------------ snapshot load
def load_snapshot():
    with open(SNAP_SOURCE) as f:
        source = json.load(f)
    for key, path in (("annotations_sha256", SNAP_ANNOT), ("ontology_subset_sha256", SNAP_ONTO)):
        if source["snapshot_files"][key] != sha256(path):
            raise SystemExit("%s does not match the checksum in %s; rebuild with "
                             "--rebuild-annotations." % (rel(path), rel(SNAP_SOURCE)))
    onto = pd.read_csv(SNAP_ONTO, sep="\t", dtype=str, keep_default_na=False)
    term_name = dict(zip(onto.go_id, onto.term_name))
    term_ont = dict(zip(onto.go_id, onto.ontology))
    parents = {r.go_id: [p for p in (r.is_a.split(";") + r.part_of.split(";")) if p]
               for r in onto.itertuples()}

    ancestors: dict = {}

    def anc(tid):
        if tid in ancestors:
            return ancestors[tid]
        # iterative DFS to avoid recursion limits on deep GO paths
        out = {tid}
        stack = list(parents.get(tid, []))
        while stack:
            p = stack.pop()
            if p in out:
                continue
            if p in ancestors:
                out |= ancestors[p]
                continue
            out.add(p)
            stack.extend(parents.get(p, []))
        ancestors[tid] = frozenset(out)
        return ancestors[tid]

    ann = pd.read_csv(SNAP_ANNOT, sep="\t", dtype=str, keep_default_na=False)
    direct = collections.defaultdict(set)
    for g, t in zip(ann.gene_id, ann.go_id):
        direct[g].add(t)
    # gene -> ontology -> set of GO terms (direct + all ancestors)
    gene_terms = {o: {} for o in ONTOLOGIES}
    for g in sorted(direct):
        full = set()
        for t in direct[g]:
            full |= anc(t)
        for o in ONTOLOGIES:
            s = frozenset(x for x in full if term_ont[x] == o)
            if s:
                gene_terms[o][g] = s
    return source, gene_terms, term_name


# ------------------------------------------------------------------- inputs
def read_tsv(name: str) -> pd.DataFrame:
    path = os.path.join(DATA_PROCESSED, name)
    if not os.path.exists(path):
        raise SystemExit("Missing input %s; run the upstream pipeline steps first." % rel(path))
    return pd.read_csv(path, sep="\t", dtype={"gene_id": str, "gene_symbol": str})


def define_gene_sets():
    """Return (gene_sets, backgrounds, symbols, input_info, checks)."""
    gene_sets = []          # dicts: id, group, subset, background, genes
    backgrounds = {}        # background id -> sorted gene list
    symbols = {}
    input_info = {}
    checks = []

    tested = {}
    for group, fname in CONTRASTS:
        df = read_tsv(fname)
        input_info[fname] = {"sha256": sha256(os.path.join(DATA_PROCESSED, fname)), "rows": int(len(df))}
        for g, s in zip(df.gene_id, df.gene_symbol):
            if isinstance(s, str) and s:
                symbols.setdefault(g, s)
        ok = df[df.padj.notna()]
        tested[group] = set(ok.gene_id)
        bg = "tested_" + group
        backgrounds[bg] = sorted(tested[group])
        deg = ok[ok.padj < DEG_PADJ]
        up = sorted(deg.gene_id[deg.log2FoldChange > 0])
        down = sorted(deg.gene_id[deg.log2FoldChange < 0])
        if "significant" in df.columns:
            n_flag = int((df.significant == "YES").sum())
            if n_flag != len(deg):
                checks.append("%s: %d genes flagged significant=YES but %d have padj < %g; "
                              "this step uses padj." % (fname, n_flag, len(deg), DEG_PADJ))
        input_info[fname].update({"genes_with_padj": int(len(ok)), "deg_up": len(up), "deg_down": len(down)})
        for subset, genes in (("up", up), ("down", down)):
            gene_sets.append({"id": "%s_%s" % (group, subset), "group": group, "subset": subset,
                              "background": bg, "genes": genes,
                              "description": "%s DEGs in %s (padj < %g, log2FoldChange %s 0)"
                                             % (subset, fname, DEG_PADJ, ">" if subset == "up" else "<")})

    # shared simulated-microgravity response
    sh = read_tsv(SHARED_FILE)
    input_info[SHARED_FILE] = {"sha256": sha256(os.path.join(DATA_PROCESSED, SHARED_FILE)), "rows": int(len(sh))}
    for g, s in zip(sh.gene_id, sh.gene_symbol):
        if isinstance(s, str) and s:
            symbols.setdefault(g, s)
    bg = "tested_both_simulators"
    backgrounds[bg] = sorted(tested[SHARED_BACKGROUND[0]] & tested[SHARED_BACKGROUND[1]])
    lc, lr = sh.log2FoldChange_clinostat, sh.log2FoldChange_rpm
    disagree = int(((lc > 0) != (lr > 0)).sum())
    if disagree:
        checks.append("%s: %d genes change direction between simulators; they are left out of "
                      "the shared up/down sets." % (SHARED_FILE, disagree))
    sh_up = sorted(sh.gene_id[(lc > 0) & (lr > 0)])
    sh_down = sorted(sh.gene_id[(lc < 0) & (lr < 0)])
    outside = sorted(set(sh.gene_id) - set(backgrounds[bg]))
    if outside:
        checks.append("%s: %d shared genes are not in the tested background of both contrasts."
                      % (SHARED_FILE, len(outside)))
    # cross-check against the contrast tables on disk (they may be regenerated separately)
    recomputed = {"up": set(), "down": set()}
    degs = {}
    for group, fname in CONTRASTS[:2]:
        df = read_tsv(fname)
        df = df[df.padj.notna() & (df.padj < DEG_PADJ)]
        degs[group] = {"up": set(df.gene_id[df.log2FoldChange > 0]),
                       "down": set(df.gene_id[df.log2FoldChange < 0])}
    for d in ("up", "down"):
        recomputed[d] = degs[SHARED_BACKGROUND[0]][d] & degs[SHARED_BACKGROUND[1]][d]
    if recomputed["up"] != set(sh_up) or recomputed["down"] != set(sh_down):
        checks.append("%s does not match the intersection of the current contrast tables "
                      "(file: %d up / %d down; recomputed: %d up / %d down). The file is used as "
                      "given; re-run the DE step if the inputs were regenerated."
                      % (SHARED_FILE, len(sh_up), len(sh_down),
                         len(recomputed["up"]), len(recomputed["down"])))
    input_info[SHARED_FILE].update({"shared_up": len(sh_up), "shared_down": len(sh_down)})
    for subset, genes in (("up", sh_up), ("down", sh_down)):
        gene_sets.append({"id": "shared_simulated_microgravity_%s" % subset,
                          "group": "shared_simulated_microgravity", "subset": subset,
                          "background": bg, "genes": genes,
                          "description": "%s in both 3D clinostat and RPM 2.0 vs static 1g (%s)"
                                         % (subset, SHARED_FILE)})

    # WGCNA modules
    wg = read_tsv(WGCNA_FILE)
    input_info[WGCNA_FILE] = {"sha256": sha256(os.path.join(DATA_PROCESSED, WGCNA_FILE)), "rows": int(len(wg))}
    for g, s in zip(wg.gene_id, wg.gene_symbol):
        if isinstance(s, str) and s:
            symbols.setdefault(g, s)
    bg = "wgcna_all_genes"
    backgrounds[bg] = sorted(set(wg.gene_id))
    modules = sorted(m for m in set(wg.module.astype(str)) if m != WGCNA_UNASSIGNED)
    input_info[WGCNA_FILE].update({"modules_tested": len(modules),
                                   "genes_unassigned_grey": int((wg.module == WGCNA_UNASSIGNED).sum())})
    for m in modules:
        gene_sets.append({"id": "module_%s" % m, "group": "wgcna", "subset": m, "background": bg,
                          "genes": sorted(wg.gene_id[wg.module == m]),
                          "description": "WGCNA module %s (%s)" % (m, WGCNA_FILE)})
    return gene_sets, backgrounds, symbols, input_info, checks


# --------------------------------------------------------------- enrichment
def term_index(universe: list, gene_terms_o: dict):
    """For one background and ontology: annotated universe and term -> genes."""
    uni = [g for g in universe if g in gene_terms_o]
    t2g = collections.defaultdict(set)
    for g in uni:
        for t in gene_terms_o[g]:
            t2g[t].add(g)
    eligible = sorted(t for t, gs in t2g.items() if MIN_TERM_SIZE <= len(gs) <= MAX_TERM_SIZE)
    return set(uni), {t: t2g[t] for t in eligible}


def enrich(gs: dict, ontology: str, uni: set, t2g: dict, term_name: dict, symbols: dict):
    query = sorted(set(gs["genes"]) & uni)
    N, M = len(query), len(uni)
    meta = {"set_size_annotated": N, "universe_size": M, "n_terms_tested": 0,
            "n_significant": 0, "n_significant_bh_hit_terms_only": 0}
    if N == 0 or M == 0 or not t2g:
        return [], meta
    terms = sorted(t2g)
    qs = set(query)
    n = np.array([len(t2g[t]) for t in terms])
    overlaps = [sorted(t2g[t] & qs) for t in terms]
    k = np.array([len(o) for o in overlaps])
    p = np.where(k > 0, hypergeom.sf(k - 1, M, n, N), 1.0)
    p = np.clip(p, 0.0, 1.0)
    q = bh_fdr(p)
    # Sensitivity: BH over terms with >= 1 overlapping gene only (the clusterProfiler
    # enricher convention). Less conservative; reported, not used for calls.
    hit = k > 0
    q_hit = np.ones(len(terms))
    q_hit[hit] = bh_fdr(p[hit])
    expected = n * N / M
    meta["n_terms_tested"] = len(terms)
    meta["n_significant"] = int(((q < FDR_ALPHA) & hit).sum())
    meta["n_significant_bh_hit_terms_only"] = int(((q_hit < FDR_ALPHA) & hit).sum())
    rows = []
    for i, t in enumerate(terms):
        if k[i] == 0:
            continue
        rows.append({
            "gene_set": gs["id"], "group": gs["group"], "subset": gs["subset"],
            "ontology": ontology, "go_id": t, "term_name": term_name[t],
            "term_size": int(n[i]), "overlap": int(k[i]), "set_size": N, "universe_size": M,
            "expected": float(expected[i]), "fold_enrichment": float(k[i] / expected[i]),
            "pvalue": float(p[i]), "fdr": float(q[i]), "n_terms_tested": len(terms),
            "fdr_bh_hit_terms_only": float(q_hit[i]),
            "overlap_gene_ids": ";".join(overlaps[i]),
            "overlap_gene_symbols": ";".join(symbols.get(g, g) for g in overlaps[i]),
        })
    rows.sort(key=lambda r: (r["fdr"], r["pvalue"], -r["fold_enrichment"], r["go_id"]))
    return rows, meta


def write_results(path: str, rows: list) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(RESULT_COLUMNS)
        for r in rows:
            w.writerow([
                r["gene_set"], r["group"], r["subset"], r["ontology"], r["go_id"], r["term_name"],
                r["term_size"], r["overlap"], r["set_size"], r["universe_size"],
                "%.4f" % r["expected"], "%.4f" % r["fold_enrichment"],
                "%.4e" % r["pvalue"], "%.4e" % r["fdr"], r["n_terms_tested"],
                "%.4e" % r["fdr_bh_hit_terms_only"],
                r["overlap_gene_ids"], r["overlap_gene_symbols"],
            ])
    log("  wrote %s (%d rows)" % (rel(path), len(rows)))


# --------------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rebuild-annotations", action="store_true",
                    help="rebuild the locus->GO snapshot from cached raw downloads")
    ap.add_argument("--refresh", action="store_true",
                    help="re-download UniProt and go-basic.obo, then rebuild the snapshot")
    args = ap.parse_args()

    log("=== Phase 5: GO over-representation analysis (UniProtKB GO, go-basic) ===")
    snapshot_ok = all(os.path.exists(p) for p in (SNAP_ANNOT, SNAP_ONTO, SNAP_SOURCE))
    if args.refresh or args.rebuild_annotations or not snapshot_ok:
        build_snapshot(refresh=args.refresh)
    else:
        log("Using annotation snapshot %s" % rel(SNAP_ANNOT))
    source, gene_terms, term_name = load_snapshot()

    gene_sets, backgrounds, symbols, input_info, checks = define_gene_sets()

    # coverage per background
    coverage = {}
    index = {}
    for bg, genes in backgrounds.items():
        any_ann = sum(1 for g in genes if any(g in gene_terms[o] for o in ONTOLOGIES))
        cov = {"genes": len(genes), "annotated_any": any_ann,
               "fraction_annotated_any": sig(any_ann / len(genes)) if genes else None}
        for o in ONTOLOGIES:
            uni, t2g = term_index(genes, gene_terms[o])
            index[(bg, o)] = (uni, t2g)
            cov["annotated_%s" % o] = len(uni)
            cov["eligible_terms_%s" % o] = len(t2g)
        coverage[bg] = cov

    deg_rows, mod_rows = [], []
    set_summaries = collections.OrderedDict()
    caveats = []
    for gs in gene_sets:
        all_rows = []
        per_ont = collections.OrderedDict()
        for o in ONTOLOGIES:
            uni, t2g = index[(gs["background"], o)]
            rows, meta = enrich(gs, o, uni, t2g, term_name, symbols)
            per_ont[o] = meta
            all_rows.extend(rows)
        (mod_rows if gs["group"] == "wgcna" else deg_rows).extend(all_rows)
        n_annot_any = sum(1 for g in gs["genes"]
                          if any(g in gene_terms[o] for o in ONTOLOGIES))
        top = sorted(all_rows, key=lambda r: (r["fdr"], r["pvalue"], -r["fold_enrichment"],
                                              ONTOLOGIES.index(r["ontology"]), r["go_id"]))[:TOP_N]
        notes = []
        if not gs["genes"]:
            notes.append("empty gene set: nothing tested")
        elif n_annot_any < SMALL_SET:
            notes.append("only %d genes with any GO annotation: very low power" % n_annot_any)
        set_summaries[gs["id"]] = {
            "description": gs["description"],
            "background": gs["background"],
            "n_genes": len(gs["genes"]),
            "n_genes_with_go": n_annot_any,
            "fraction_with_go": sig(n_annot_any / len(gs["genes"])) if gs["genes"] else None,
            "n_significant_fdr_lt_%g" % FDR_ALPHA: sum(m["n_significant"] for m in per_ont.values()),
            "per_ontology": per_ont,
            "top_terms_by_fdr": [
                {"ontology": r["ontology"], "go_id": r["go_id"], "term_name": r["term_name"],
                 "overlap": r["overlap"], "term_size": r["term_size"], "set_size": r["set_size"],
                 "fold_enrichment": sig(r["fold_enrichment"]), "pvalue": sig(r["pvalue"]),
                 "fdr": sig(r["fdr"]), "significant": bool(r["fdr"] < FDR_ALPHA)} for r in top],
            "notes": notes,
        }
        if notes:
            caveats.append("%s: %s." % (gs["id"], notes[0]))

    order = {gs["id"]: i for i, gs in enumerate(gene_sets)}
    sort_key = lambda r: (order[r["gene_set"]], ONTOLOGIES.index(r["ontology"]),
                          r["fdr"], r["pvalue"], -r["fold_enrichment"], r["go_id"])
    deg_rows.sort(key=sort_key)
    mod_rows.sort(key=sort_key)
    write_results(OUT_DEG, deg_rows)
    write_results(OUT_MOD, mod_rows)

    # global caveats, all computed from the data; per-gene-set notes follow them
    ann_src = source["annotation_source"]
    global_caveats = []
    iea = ann_src.get("fraction_electronic_iea")
    if iea is not None:
        global_caveats.append(
            "%.1f%% of the UniProtKB GO annotations used are electronic (IEA: InterPro2GO, "
            "UniRule, EC and similar pipelines), not experimental, so the results describe "
            "predicted function." % (100 * iea))
    same_cov = collections.OrderedDict()  # backgrounds with identical coverage share one line
    for bg, c in coverage.items():
        if c["genes"]:
            key = (c["annotated_any"], c["genes"], c["fraction_annotated_any"],
                   c["annotated_BP"], c["annotated_MF"], c["annotated_CC"])
            same_cov.setdefault(key, []).append(bg)
    for (n_any, n_all, frac, n_bp, n_mf, n_cc), bgs in same_cov.items():
        global_caveats.append(
            "%s: %d of %d genes (%.1f%%) have at least one GO annotation (BP %d, MF %d, CC %d). "
            "Unannotated genes, mostly hypothetical proteins, cannot contribute."
            % (", ".join(bgs), n_any, n_all, 100 * frac, n_bp, n_mf, n_cc))
    n_sig = sum(s["n_significant_fdr_lt_%g" % FDR_ALPHA] for s in set_summaries.values())
    n_sig_hit = sum(m["n_significant_bh_hit_terms_only"] for s in set_summaries.values()
                    for m in s["per_ontology"].values())
    global_caveats.append(
        "BH is applied over all size-eligible terms (conservative). With the clusterProfiler "
        "enricher convention (BH over terms with at least one hit, column fdr_bh_hit_terms_only) "
        "the total across gene sets is %d significant term-tests instead of %d." % (n_sig_hit, n_sig))
    global_caveats.append("GO terms overlap (parents and children share genes), so significant "
                          "terms within a gene set are not independent findings.")
    global_caveats.append("The hypergeometric test treats genes as independent. Co-transcribed "
                          "operons whose members share annotations (for example nuo or mce "
                          "clusters) can drive a term on their own; check overlap_gene_ids "
                          "before interpreting a term.")
    caveats = global_caveats + checks + caveats
    all_genes = sorted(set().union(*[set(v) for v in backgrounds.values()]))
    no_map = [g for g in all_genes if not any(g in gene_terms[o] for o in ONTOLOGIES)]

    summary = {
        "script": "analysis/05_ontology_functional_enrichment.py",
        "annotation_source": ann_src,
        "locus_mapping": source["locus_mapping"],
        "ontology_source": source["ontology_source"],
        "snapshot_files": source["snapshot_files"],
        "method": {
            "test": "one-sided hypergeometric (over-representation), scipy.stats.hypergeom.sf",
            "ontologies": "BP, MF and CC tested separately",
            "propagation": "true-path rule over is_a and part_of, within namespace",
            "universe": ("genes in the gene set's background with at least one annotation in "
                         "the ontology being tested"),
            "backgrounds": {
                "tested_<contrast>": "genes with non-missing DESeq2 padj in that contrast",
                "tested_both_simulators": ("genes with non-missing padj in both %s and %s"
                                           % SHARED_BACKGROUND),
                "wgcna_all_genes": "every gene in %s, grey included" % WGCNA_FILE,
            },
            "deg_definition": "DESeq2 padj < %g, direction by sign of log2FoldChange" % DEG_PADJ,
            "term_size_range": [MIN_TERM_SIZE, MAX_TERM_SIZE],
            "term_size_definition": "annotated universe genes in the term (after propagation)",
            "multiple_testing": ("Benjamini-Hochberg within each gene set and ontology, over all "
                                 "size-eligible terms (zero-overlap terms enter with p = 1)"),
            "multiple_testing_sensitivity": ("fdr_bh_hit_terms_only: BH over terms with at least "
                                             "one overlapping gene (clusterProfiler enricher "
                                             "convention); reported only, not used for calls"),
            "significance": "FDR < %g" % FDR_ALPHA,
            "evidence_codes": "all UniProtKB GO evidence codes kept, IEA included",
            "small_set_flag": "gene sets with fewer than %d GO-annotated genes" % SMALL_SET,
        },
        "inputs": input_info,
        "coverage": coverage,
        "genes_in_any_background_without_go": len(no_map),
        "gene_sets": set_summaries,
        "caveats": caveats,
        "outputs": [rel(OUT_DEG), rel(OUT_MOD), rel(OUT_SUMMARY)],
    }
    write_json(OUT_SUMMARY, summary)
    log("  wrote %s" % rel(OUT_SUMMARY))

    log("\nGene set                                   genes  GO-annot  sig(FDR<%g)" % FDR_ALPHA)
    for sid, s in set_summaries.items():
        log("  %-40s %5d  %8d  %5d" % (sid, s["n_genes"], s["n_genes_with_go"],
                                       s["n_significant_fdr_lt_%g" % FDR_ALPHA]))


if __name__ == "__main__":
    sys.exit(main())
