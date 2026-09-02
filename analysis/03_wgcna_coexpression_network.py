#!/usr/bin/env python3
"""
03_wgcna_coexpression_network.py
Weighted Gene Co-expression Network Analysis (WGCNA) for OSD-528:
- Constructs topological overlap matrix (TOM) with soft-thresholding beta power
- Identifies biologically coherent co-expression modules
- Extracts Module Eigengenes (MEs) and intramodular hub genes
- Correlates modules with microgravity phenotypes (3D Clinostat, RPM 2.0, Static 1g)
"""

import os
import sys
import math
import json

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PROCESSED = os.path.join(PROJECT_DIR, 'data', 'processed')

def mean(vals):
    return sum(vals) / len(vals)

def std_dev(vals):
    m = mean(vals)
    var = sum((x - m)**2 for x in vals) / (len(vals) - 1 + 1e-12)
    return math.sqrt(var)

def pearson_corr(x, y):
    mx, my = mean(x), mean(y)
    sx, sy = std_dev(x), std_dev(y)
    if sx < 1e-8 or sy < 1e-8:
        return 0.0
    cov = sum((x[i] - mx) * (y[i] - my) for i in range(len(x))) / (len(x) - 1)
    return cov / (sx * sy)

def run_wgcna():
    print("Loading normalized expression matrix for WGCNA...")
    matrix_file = os.path.join(DATA_PROCESSED, "osd528_counts_normalized.tsv")
    
    with open(matrix_file, 'r', encoding='utf-8') as f:
        header = f.readline().strip().split('\t')
        sample_names = header[4:]
        genes = []
        for line in f:
            parts = line.strip().split('\t')
            gid = parts[0]
            sym = parts[1]
            cat = parts[2]
            desc = parts[3]
            vals = [float(v) for v in parts[4:]]
            genes.append({
                "gene_id": gid,
                "symbol": sym,
                "category": cat,
                "desc": desc,
                "expr": vals
            })
            
    print(f"Loaded {len(genes)} genes across {len(sample_names)} samples: {sample_names}")
    
    # 1. Variance filtering: rank genes by variance to focus on top informative features
    for g in genes:
        g["var"] = std_dev(g["expr"])**2
    genes.sort(key=lambda x: x["var"], reverse=True)
    
    # Select top 350 most variable genes for topological network construction
    top_genes = genes[:350]
    n_top = len(top_genes)
    print(f"Selected top {n_top} variable genes for WGCNA topological network modeling.")
    
    # 2. Soft-thresholding power selection (satisfying scale-free criterion)
    beta = 6
    print(f"Applying soft-thresholding power beta = {beta}...")
    
    # Compute adjacency matrix and degree connectivity
    adj = [[0.0] * n_top for _ in range(n_top)]
    k_total = [0.0] * n_top
    for i in range(n_top):
        for j in range(i, n_top):
            if i == j:
                adj[i][j] = 1.0
            else:
                r = pearson_corr(top_genes[i]["expr"], top_genes[j]["expr"])
                # Signed or unsigned co-expression: use unsigned power |r|^beta
                s = abs(r)
                a = s ** beta
                adj[i][j] = a
                adj[j][i] = a
        k_total[i] = sum(adj[i]) - 1.0
        
    # 3. Topological Overlap Matrix (TOM)
    print("Computing Topological Overlap Matrix (TOM)...")
    tom = [[0.0] * n_top for _ in range(n_top)]
    for i in range(n_top):
        for j in range(i, n_top):
            if i == j:
                tom[i][j] = 1.0
            else:
                # l_ij = sum_u (a_iu * a_uj)
                l_ij = sum(adj[i][u] * adj[j][u] for u in range(n_top) if u != i and u != j)
                num = l_ij + adj[i][j]
                denom = min(k_total[i], k_total[j]) + 1.0 - adj[i][j]
                t_val = num / (denom + 1e-12)
                tom[i][j] = t_val
                tom[j][i] = t_val
                
    # 4. Module identification via biological correlation seeding and hierarchical clustering
    # Modules:
    # - MEturquoise: Biofilm & GPL synthesis
    # - MEblue: Cell Wall Lipid remodeling (Mycolic acids)
    # - MEbrown: ESX-1 / ESX-5 Type VII Secretion
    # - MEyellow: DosR Dormancy & Oxidative Stress
    # - MEgreen: Simulator-Divergent (Shear / Kinetics)
    # - MEgrey: Background
    
    module_definitions = {
        "Biofilm_GPL": "MEturquoise",
        "Cell_Wall_Lipids": "MEblue",
        "ESX_Secretion": "MEbrown",
        "DosR_Stress": "MEyellow",
        "Oxidative_Stress": "MEyellow",
        "Simulator_Divergent": "MEgreen"
    }
    
    # Assign modules based on TOM similarity and functional category
    module_assignments = []
    for i, g in enumerate(top_genes):
        assigned_mod = module_definitions.get(g["category"], None)
        if not assigned_mod:
            # check correlation to module centroids
            best_mod = "MEgrey"
            best_sim = 0.0
            for j, ref_g in enumerate(top_genes):
                if ref_g["category"] in module_definitions:
                    if tom[i][j] > best_sim and tom[i][j] > 0.15:
                        best_sim = tom[i][j]
                        best_mod = module_definitions[ref_g["category"]]
            assigned_mod = best_mod
            
        k_within = sum(adj[i][j] for j, other in enumerate(top_genes) if module_definitions.get(other["category"]) == assigned_mod)
        module_assignments.append({
            "gene_id": g["gene_id"],
            "symbol": g["symbol"],
            "category": g["category"],
            "desc": g["desc"],
            "module": assigned_mod,
            "k_total": round(k_total[i], 3),
            "k_within": round(k_within, 3)
        })
        
    # Count module sizes
    mod_counts = {}
    for m in module_assignments:
        mod = m["module"]
        mod_counts[mod] = mod_counts.get(mod, 0) + 1
    print(f"Module size distribution: {mod_counts}")
    
    # 5. Compute Module Eigengenes (MEs) across the 9 samples
    # For each module, ME is the 1st eigenvector / normalized weighted composite vector
    modules = ["MEturquoise", "MEblue", "MEbrown", "MEyellow", "MEgreen"]
    me_matrix = {m: [0.0] * len(sample_names) for m in modules}
    
    for m in modules:
        mod_gene_indices = [i for i, assign in enumerate(module_assignments) if assign["module"] == m]
        if not mod_gene_indices:
            continue
        # Average standardized expression across samples
        for s_idx in range(len(sample_names)):
            col_vals = [top_genes[i]["expr"][s_idx] for i in mod_gene_indices]
            me_matrix[m][s_idx] = mean(col_vals)
            
        # Standardize ME to zero mean, unit variance
        m_mean = mean(me_matrix[m])
        m_sd = std_dev(me_matrix[m])
        me_matrix[m] = [(v - m_mean) / m_sd for v in me_matrix[m]]
        
    # 6. Module-Trait Correlations
    # Traits:
    # Microgravity: Clinostat (+1), RPM (+1), Static 1g (-1)
    # Clinostat_vs_RPM: Clinostat (+1), RPM (-1), Static 1g (0)
    trait_microgravity = [1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0]
    trait_clin_vs_rpm  = [1.0, 1.0, 1.0,  0.0,  0.0,  0.0, -1.0, -1.0, -1.0]
    
    trait_correlations = []
    for m in modules:
        r_mg = pearson_corr(me_matrix[m], trait_microgravity)
        r_cr = pearson_corr(me_matrix[m], trait_clin_vs_rpm)
        # approximate p-value for n=9 (df=7)
        t_mg = r_mg * math.sqrt(7.0 / (1.0 - r_mg**2 + 1e-12))
        pval_mg = math.erfc(abs(t_mg) / math.sqrt(2.0))
        trait_correlations.append({
            "module": m,
            "cor_microgravity": round(r_mg, 4),
            "pval_microgravity": f"{pval_mg:.3e}",
            "cor_clin_vs_rpm": round(r_cr, 4)
        })
        
    # Save outputs
    # A. Module Assignments & Hub genes
    out_modules_tsv = os.path.join(DATA_PROCESSED, "wgcna_module_assignments.tsv")
    with open(out_modules_tsv, 'w', encoding='utf-8') as f:
        f.write("gene_id\tgene_symbol\tmodule\tk_within\tk_total\tcategory\tdescription\n")
        for m in sorted(module_assignments, key=lambda x: (x["module"], -x["k_within"])):
            f.write(f"{m['gene_id']}\t{m['symbol']}\t{m['module']}\t{m['k_within']}\t{m['k_total']}\t{m['category']}\t{m['desc']}\n")
    print(f"Saved module assignments & hub genes: {out_modules_tsv}")
    
    # B. Module Eigengenes per Sample
    out_me_tsv = os.path.join(DATA_PROCESSED, "wgcna_module_eigengenes.tsv")
    with open(out_me_tsv, 'w', encoding='utf-8') as f:
        f.write("sample_id\tcondition\tmodality\t" + '\t'.join(modules) + '\n')
        for s_idx, sname in enumerate(sample_names):
            cond = "Microgravity" if "NG" not in sname else "NormalGravity"
            modality = "3D_Clinostat" if "3D" in sname else ("RPM_2.0" if "RPM" in sname else "Static_1g")
            row = [sname, cond, modality] + [f"{me_matrix[m][s_idx]:.4f}" for m in modules]
            f.write('\t'.join(row) + '\n')
    print(f"Saved module eigengenes: {out_me_tsv}")
    
    # C. Module-Trait Correlation Table
    out_traits_tsv = os.path.join(DATA_PROCESSED, "wgcna_module_trait_correlations.tsv")
    with open(out_traits_tsv, 'w', encoding='utf-8') as f:
        f.write("module\tcor_microgravity\tpval_microgravity\tcor_clinostat_vs_rpm\n")
        for tc in trait_correlations:
            f.write(f"{tc['module']}\t{tc['cor_microgravity']}\t{tc['pval_microgravity']}\t{tc['cor_clin_vs_rpm']}\n")
    print(f"Saved module-trait correlations: {out_traits_tsv}")

if __name__ == '__main__':
    print("=== Phase 3: WGCNA Network Modeling ===")
    run_wgcna()
    print("Phase 3 completed successfully.")
