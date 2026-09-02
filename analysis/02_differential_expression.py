#!/usr/bin/env python3
"""
02_differential_expression.py
Empirical differential expression analysis across microgravity simulation modalities
calculated directly from real NASA OSDR OSD-528 RNA-seq normalized counts.

Contrasts evaluated:
  1. 3D Clinostat vs Static 1g
  2. RPM 2.0 vs Static 1g
  3. 3D Clinostat vs RPM 2.0
"""

import os
import sys
import math
import csv

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PROCESSED = os.path.join(PROJECT_DIR, "data", "processed")
NORM_MATRIX_FILE = os.path.join(DATA_PROCESSED, "osd528_counts_normalized.tsv")
SAMPLE_META_FILE = os.path.join(DATA_PROCESSED, "osd528_sample_metadata.tsv")


def t_test_two_sample(vals_a, vals_b):
    """Computes two-sample t-statistic and two-tailed p-value using standard normal/t approximation."""
    n_a = len(vals_a)
    n_b = len(vals_b)
    mean_a = sum(vals_a) / n_a
    mean_b = sum(vals_b) / n_b
    lfc = mean_a - mean_b

    var_a = sum((x - mean_a) ** 2 for x in vals_a) / max(1, n_a - 1)
    var_b = sum((x - mean_b) ** 2 for x in vals_b) / max(1, n_b - 1)

    # Pooled variance with moderation floor
    sp = math.sqrt(((n_a - 1) * var_a + (n_b - 1) * var_b) / max(1, n_a + n_b - 2)) + 0.05
    denom = sp * math.sqrt(1.0 / n_a + 1.0 / n_b)

    t_stat = lfc / (denom + 1e-9)
    z = abs(t_stat)
    # Approximation of complementary error function for two-tailed p-value
    pval = math.erfc(z / math.sqrt(2.0))
    pval = max(1e-15, min(1.0, pval))
    return lfc, pval


def calc_bh_fdr(records):
    """Computes Benjamini-Hochberg adjusted p-values (FDR)."""
    n = len(records)
    # records format: (gid, sym, desc, lfc, pval)
    sorted_idx = sorted(range(n), key=lambda i: records[i][4])
    fdrs = [1.0] * n
    min_fdr = 1.0

    for rank_rev, idx in enumerate(reversed(sorted_idx)):
        rank = n - rank_rev
        p = records[idx][4]
        q = (p * n) / rank
        min_fdr = min(min_fdr, q)
        fdrs[idx] = min(1.0, max(0.0, min_fdr))

    return [records[i] + (fdrs[i],) for i in range(n)]


def main():
    print("=== Phase 2: Empirical Differential Expression Analysis (OSD-528 Real RNA-Seq) ===")

    if not os.path.exists(NORM_MATRIX_FILE):
        print(f"Error: Normalized counts matrix not found at {NORM_MATRIX_FILE}")
        sys.exit(1)

    # Save canonical sample metadata
    samples_info = [
        ("RFP3D11", "3D_Clinostat", "Microgravity", "PDMS_silicone"),
        ("RFP3D39", "3D_Clinostat", "Microgravity", "PDMS_silicone"),
        ("RFP3D47", "3D_Clinostat", "Microgravity", "PDMS_silicone"),
        ("RFPNG14", "Static_1g", "NormalGravity", "PDMS_silicone"),
        ("RFPNG35", "Static_1g", "NormalGravity", "PDMS_silicone"),
        ("RFPNG45", "Static_1g", "NormalGravity", "PDMS_silicone"),
        ("RFPRPM4", "RPM_2.0", "Microgravity", "PDMS_silicone"),
        ("RFPRPM41", "RPM_2.0", "Microgravity", "PDMS_silicone"),
        ("RFPRPM6", "RPM_2.0", "Microgravity", "PDMS_silicone"),
    ]
    with open(SAMPLE_META_FILE, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id", "modality", "condition", "substrate"])
        for s in samples_info:
            writer.writerow(s)
    print(f"Saved sample metadata to {SAMPLE_META_FILE}")

    # Load normalized expression matrix
    genes_data = []
    with open(NORM_MATRIX_FILE, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gid = row["gene_id"]
            sym = row["gene_symbol"]
            protein = row.get("protein", "")

            # Extract sample values
            clin_vals = [float(row["RFP3D11"]), float(row["RFP3D39"]), float(row["RFP3D47"])]
            ng_vals = [float(row["RFPNG14"]), float(row["RFPNG35"]), float(row["RFPNG45"])]
            rpm_vals = [float(row["RFPRPM4"]), float(row["RFPRPM41"]), float(row["RFPRPM6"])]

            genes_data.append((gid, sym, protein, clin_vals, ng_vals, rpm_vals))

    print(f"Loaded {len(genes_data):,} genes from real expression matrix.")

    # Compute Contrasts
    res_clin_vs_ng = []
    res_rpm_vs_ng = []
    res_clin_vs_rpm = []

    for gid, sym, protein, clin_vals, ng_vals, rpm_vals in genes_data:
        # Filter completely unexpressed genes (all 0s)
        if sum(clin_vals + ng_vals + rpm_vals) == 0:
            continue

        lfc_c, pval_c = t_test_two_sample(clin_vals, ng_vals)
        lfc_r, pval_r = t_test_two_sample(rpm_vals, ng_vals)
        lfc_cr, pval_cr = t_test_two_sample(clin_vals, rpm_vals)

        res_clin_vs_ng.append((gid, sym, protein, lfc_c, pval_c))
        res_rpm_vs_ng.append((gid, sym, protein, lfc_r, pval_r))
        res_clin_vs_rpm.append((gid, sym, protein, lfc_cr, pval_cr))

    adj_clin_vs_ng = calc_bh_fdr(res_clin_vs_ng)
    adj_rpm_vs_ng = calc_bh_fdr(res_rpm_vs_ng)
    adj_clin_vs_rpm = calc_bh_fdr(res_clin_vs_rpm)

    def save_contrast_file(adj_data, filename, title):
        out_file = os.path.join(DATA_PROCESSED, filename)
        sig_count = 0
        with open(out_file, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["gene_id", "gene_symbol", "description", "category", "log2FoldChange", "pvalue", "padj", "significant"])
            for gid, sym, protein, lfc, pval, padj in sorted(adj_data, key=lambda x: x[5]):
                is_sig = "YES" if (padj < 0.05 and abs(lfc) >= 0.75) else "NO"
                if is_sig == "YES":
                    sig_count += 1
                cat = "General"
                writer.writerow([gid, sym, protein, cat, f"{lfc:.4f}", f"{pval:.4e}", f"{padj:.4e}", is_sig])
        print(f"Saved {title}: {out_file} (Significant DEGs: {sig_count:,} of {len(adj_data):,})")
        return sig_count

    sig_clin = save_contrast_file(adj_clin_vs_ng, "deg_3dclinostat_vs_static1g.tsv", "3D Clinostat vs Static 1g")
    sig_rpm = save_contrast_file(adj_rpm_vs_ng, "deg_rpm2_vs_static1g.tsv", "RPM 2.0 vs Static 1g")
    sig_cr = save_contrast_file(adj_clin_vs_rpm, "deg_3dclinostat_vs_rpm2.tsv", "3D Clinostat vs RPM 2.0")

    print("\nSummary of Empirical Contrasts (Real NASA OSDR OSD-528 RNA-Seq):")
    print(f"  3D Clinostat vs Static 1g : {sig_clin:,} significant DEGs (padj < 0.05, |log2FC| >= 0.75)")
    print(f"  RPM 2.0 vs Static 1g      : {sig_rpm:,} significant DEGs (padj < 0.05, |log2FC| >= 0.75)")
    print(f"  3D Clinostat vs RPM 2.0   : {sig_cr:,} significant DEGs (padj < 0.05, |log2FC| >= 0.75)")


if __name__ == "__main__":
    main()
