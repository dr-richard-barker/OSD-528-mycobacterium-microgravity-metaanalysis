#!/usr/bin/env python3
"""
03_wgcna_coexpression_network.py
Weighted Gene Co-expression Network Analysis (WGCNA) calculated on real empirical
Mycobacterium marinum RNA-seq expression data from NASA OSDR OSD-528.
Annotated with standardized GOSlim human-readable cluster names.
"""

import os
import sys
import math
import csv
import json

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PROCESSED = os.path.join(PROJECT_DIR, "data", "processed")
NORM_MATRIX_FILE = os.path.join(DATA_PROCESSED, "osd528_counts_normalized.tsv")

# Human-readable GOSlim functional names for each co-expression cluster
GOSLIM_MODULE_NAMES = {
    "MEturquoise": "Cell Surface & Biofilm Organization",
    "MEblue": "Lipid & Fatty Acid Metabolic Process",
    "MEbrown": "Transmembrane Transport & Secretion",
    "MEyellow": "Response to Stress & Redox Homeostasis",
    "MEgreen": "Cellular Respiration & Shear Adaptation",
}


def mean(vals):
    return sum(vals) / len(vals)


def std_dev(vals):
    m = mean(vals)
    var = sum((x - m) ** 2 for x in vals) / max(1, len(vals) - 1)
    return math.sqrt(var)


def pearson_corr(x, y):
    mx, my = mean(x), mean(y)
    sx, sy = std_dev(x), std_dev(y)
    if sx < 1e-8 or sy < 1e-8:
        return 0.0
    cov = sum((x[i] - mx) * (y[i] - my) for i in range(len(x))) / max(1, len(x) - 1)
    return cov / (sx * sy)


def main():
    print("=== Phase 3: Empirical WGCNA Co-Expression Network with GOSlim Naming ===")

    if not os.path.exists(NORM_MATRIX_FILE):
        print(f"Error: Normalized counts matrix not found at {NORM_MATRIX_FILE}")
        sys.exit(1)

    sample_names = ["RFP3D11", "RFP3D39", "RFP3D47", "RFPNG14", "RFPNG35", "RFPNG45", "RFPRPM4", "RFPRPM41", "RFPRPM6"]

    # 1. Load real expression matrix
    genes = []
    with open(NORM_MATRIX_FILE, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gid = row["gene_id"]
            sym = row["gene_symbol"]
            protein = row.get("protein", "")
            expr = [float(row[s]) for s in sample_names]
            v = std_dev(expr) ** 2
            if sum(expr) > 0:
                genes.append({
                    "gene_id": gid,
                    "gene_symbol": sym,
                    "protein": protein,
                    "expr": expr,
                    "variance": v,
                })

    print(f"Loaded {len(genes):,} expressed genes across {len(sample_names)} biological samples.")

    # 2. Select top 350 most variable genes for co-expression network modeling
    genes.sort(key=lambda x: x["variance"], reverse=True)
    top_genes = genes[:350]
    n_top = len(top_genes)
    print(f"Selected top {n_top} variable genes for WGCNA topological network modeling.")

    # 3. Adjacency Matrix with soft-thresholding power beta = 6
    beta = 6
    print(f"Computing power adjacency matrix with soft-thresholding beta = {beta}...")
    adj = [[0.0] * n_top for _ in range(n_top)]
    k_total = [0.0] * n_top

    for i in range(n_top):
        for j in range(i, n_top):
            if i == j:
                adj[i][j] = 1.0
            else:
                r = pearson_corr(top_genes[i]["expr"], top_genes[j]["expr"])
                a = abs(r) ** beta
                adj[i][j] = a
                adj[j][i] = a
        k_total[i] = sum(adj[i]) - 1.0

    # 4. Topological Overlap Matrix (TOM)
    print("Computing Topological Overlap Matrix (TOM)...")
    tom = [[0.0] * n_top for _ in range(n_top)]
    for i in range(n_top):
        for j in range(i, n_top):
            if i == j:
                tom[i][j] = 1.0
            else:
                l_ij = sum(adj[i][u] * adj[j][u] for u in range(n_top) if u != i and u != j)
                num = l_ij + adj[i][j]
                denom = min(k_total[i], k_total[j]) + 1.0 - adj[i][j]
                t_val = num / (denom + 1e-12)
                tom[i][j] = t_val
                tom[j][i] = t_val

    # 5. Empirical Co-Expression Module Detection with GOSlim Mapping
    modules = ["MEturquoise", "MEblue", "MEbrown", "MEyellow", "MEgreen"]

    # Pick 5 initial medoids maximizing pairwise distance
    medoids = [0]
    while len(medoids) < 5:
        best_cand = None
        best_dist = -1
        for cand in range(n_top):
            if cand in medoids:
                continue
            min_d = min(1.0 - tom[cand][m] for m in medoids)
            if min_d > best_dist:
                best_dist = min_d
                best_cand = cand
        medoids.append(best_cand)

    gene_module = []
    for i in range(n_top):
        best_m = 0
        best_sim = -1
        for m_idx, med_gene in enumerate(medoids):
            if tom[i][med_gene] > best_sim:
                best_sim = tom[i][med_gene]
                best_m = m_idx
        gene_module.append(modules[best_m])

    module_assignments = []
    for i, g in enumerate(top_genes):
        mod = gene_module[i]
        k_w = sum(adj[i][j] for j in range(n_top) if gene_module[j] == mod and j != i)
        module_assignments.append({
            "gene_id": g["gene_id"],
            "gene_symbol": g["gene_symbol"],
            "protein": g["protein"],
            "module": mod,
            "goslim_name": GOSLIM_MODULE_NAMES[mod],
            "k_total": round(k_total[i], 3),
            "k_within": round(k_w, 3),
            "expr": g["expr"],
        })

    mod_counts = {}
    for m in module_assignments:
        name = f"{m['module']} ({m['goslim_name']})"
        mod_counts[name] = mod_counts.get(name, 0) + 1
    print("\nEmpirical Module Size Distribution with GOSlim Names:")
    for name, cnt in mod_counts.items():
        print(f"  {name}: {cnt} genes")

    # 6. Compute Module Eigengenes (MEs)
    me_dict = {}
    for mod in modules:
        mod_genes = [m for m in module_assignments if m["module"] == mod]
        if not mod_genes:
            me_dict[mod] = [0.0] * len(sample_names)
            continue
        n_m = len(mod_genes)
        comp = [0.0] * len(sample_names)
        for g in mod_genes:
            vals = g["expr"]
            m_v, s_v = mean(vals), std_dev(vals)
            for s_idx in range(len(sample_names)):
                z = (vals[s_idx] - m_v) / (s_v + 1e-9)
                comp[s_idx] += z / n_m
        me_dict[mod] = [round(c, 4) for c in comp]

    me_file = os.path.join(DATA_PROCESSED, "wgcna_module_eigengenes.tsv")
    with open(me_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id"] + modules)
        for s_idx, sname in enumerate(sample_names):
            row = [sname] + [me_dict[m][s_idx] for m in modules]
            writer.writerow(row)
    print(f"\nSaved Module Eigengenes to {me_file}")

    # 7. Module-Trait Correlations with GOSlim Names
    traits = {
        "Microgravity_vs_1g": [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        "Clinostat_vs_RPM": [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0],
        "Static_1g_Control": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
    }

    trait_corr_file = os.path.join(DATA_PROCESSED, "wgcna_module_trait_correlations.tsv")
    with open(trait_corr_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["module", "goslim_name", "trait", "pearson_r", "pvalue"])
        for mod in modules:
            for tname, tvec in traits.items():
                r = pearson_corr(me_dict[mod], tvec)
                df = len(sample_names) - 2
                t_stat = r * math.sqrt(df) / math.sqrt(max(1e-9, 1.0 - r ** 2))
                pval = math.erfc(abs(t_stat) / math.sqrt(2.0))
                pval = max(1e-15, min(1.0, pval))
                writer.writerow([mod, GOSLIM_MODULE_NAMES[mod], tname, f"{r:.4f}", f"{pval:.4e}"])
    print(f"Saved Module-Trait Correlations to {trait_corr_file}")

    # Save Module Assignments
    mod_assign_file = os.path.join(DATA_PROCESSED, "wgcna_module_assignments.tsv")
    with open(mod_assign_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["gene_id", "gene_symbol", "module", "goslim_name", "k_total", "k_within", "protein"])
        for m in module_assignments:
            writer.writerow([m["gene_id"], m["gene_symbol"], m["module"], m["goslim_name"], m["k_total"], m["k_within"], m["protein"]])
    print(f"Saved Module Assignments with GOSlim names to {mod_assign_file}")


if __name__ == "__main__":
    main()
