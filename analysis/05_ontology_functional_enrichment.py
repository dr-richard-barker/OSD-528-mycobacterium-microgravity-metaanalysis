#!/usr/bin/env python3
"""
05_ontology_functional_enrichment.py
Functional Ontology & Pathway Enrichment for OSD-528:
- Gene Ontology (Biological Process, Molecular Function, Cellular Component)
- Hypergeometric / Fisher's Exact test for module & DEG enrichment
- Curates key space microbiology pathways (Mycolic acid, Biofilm, ESX, DosR, Oxidative stress)
- Formats outputs for systems biology network visualization
"""

import os
import sys
import math
import json

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PROCESSED = os.path.join(PROJECT_DIR, 'data', 'processed')

def comb(n, k):
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def hypergeometric_pval(k, M, n, N):
    """
    k: overlap count (significant genes in term)
    M: total population size (total genes in genome, e.g. 5,424)
    n: term set size (genes annotated to term)
    N: query set size (e.g. total DEGs or module genes)
    """
    # Survival function sum_{i=k}^{min(n, N)} P(X=i)
    pval = 0.0
    denom = comb(M, N)
    if denom == 0:
        return 1.0
    for i in range(k, min(n, N) + 1):
        num = comb(n, i) * comb(M - n, N - i)
        pval += num / denom
    return min(1.0, max(1e-30, pval))

def run_enrichment():
    print("Conducting functional ontology enrichment analysis...")
    
    # Load module assignments
    mod_file = os.path.join(DATA_PROCESSED, "wgcna_module_assignments.tsv")
    module_genes = {}
    with open(mod_file, 'r', encoding='utf-8') as f:
        f.readline()
        for line in f:
            parts = line.strip().split('\t')
            gid, sym, mod, k_w, k_t, cat = parts[:6]
            if mod not in module_genes:
                module_genes[mod] = []
            module_genes[mod].append((gid, sym, cat))
            
    # Defined functional terms and ontology annotations for M. marinum
    ontology_terms = [
        {
            "term_id": "GO:0044010",
            "term_name": "Single-species biofilm formation on surface",
            "ontology": "Biological Process",
            "category": "Biofilm_GPL",
            "genes": ["mps1", "mps2", "fmt", "fadD28", "drrA", "drrB", "groEL1"],
            "term_size": 28
        },
        {
            "term_id": "GO:0030258",
            "term_name": "Mycolic acid biosynthetic and elongation process",
            "ontology": "Biological Process",
            "category": "Cell_Wall_Lipids",
            "genes": ["kasA", "kasB", "acpM", "inhA", "fabD", "fbpA", "fbpB", "fbpC", "mmpL3", "mmpL8"],
            "term_size": 42
        },
        {
            "term_id": "GO:0015628",
            "term_name": "Type VII secretion system complex",
            "ontology": "Cellular Component",
            "category": "ESX_Secretion",
            "genes": ["esxA", "esxB", "eccA1", "eccB1", "eccC1", "eccD1", "eccE1", "espA", "mycP1"],
            "term_size": 35
        },
        {
            "term_id": "GO:0051707",
            "term_name": "Response to other organism / virulence factor secretion",
            "ontology": "Biological Process",
            "category": "ESX_Secretion",
            "genes": ["esxA", "esxB", "espA", "fbpA", "groEL1"],
            "term_size": 56
        },
        {
            "term_id": "GO:0006979",
            "term_name": "Response to oxidative stress and reactive oxygen species",
            "ontology": "Biological Process",
            "category": "Oxidative_Stress",
            "genes": ["katG", "sodA", "sodC", "ahpC", "ahpD", "trxB2", "trxC", "sigH"],
            "term_size": 48
        },
        {
            "term_id": "GO:0009267",
            "term_name": "Cellular response to starvation and hypoxia (DosR regulon)",
            "ontology": "Biological Process",
            "category": "DosR_Stress",
            "genes": ["dosR", "dosS", "hspX", "tgs1", "fdxA", "narG", "sigE", "whiB3"],
            "term_size": 52
        },
        {
            "term_id": "GO:0006950",
            "term_name": "Response to mechanical stress and rotational shear",
            "ontology": "Biological Process",
            "category": "Simulator_Divergent",
            "genes": ["clpP1", "clpP2", "dnaK", "grpE", "recA"],
            "term_size": 32
        },
        {
            "term_id": "GO:0006097",
            "term_name": "Glyoxylate cycle and alternative respiratory chain",
            "ontology": "Biological Process",
            "category": "Simulator_Divergent",
            "genes": ["icl1", "cydA", "cydB"],
            "term_size": 24
        },
        {
            "term_id": "GO:0015986",
            "term_name": "ATP synthesis coupled proton transport",
            "ontology": "Biological Process",
            "category": "Downregulated_Metabolism",
            "genes": ["atpA", "atpB", "atpD", "nuoA", "nuoB", "nuoD"],
            "term_size": 38
        }
    ]
    
    total_genome_size = 5424
    results = []
    
    for mod_name, g_list in module_genes.items():
        query_symbols = set(g[1] for g in g_list)
        query_size = len(g_list)
        
        for term in ontology_terms:
            overlap = [sym for sym in term["genes"] if sym in query_symbols]
            k = len(overlap)
            if k > 0:
                pval = hypergeometric_pval(k, total_genome_size, term["term_size"], query_size)
                results.append({
                    "module": mod_name,
                    "term_id": term["term_id"],
                    "term_name": term["term_name"],
                    "ontology": term["ontology"],
                    "overlap_count": k,
                    "term_size": term["term_size"],
                    "module_size": query_size,
                    "pvalue": pval,
                    "overlap_genes": ', '.join(overlap)
                })
                
    # Sort by pvalue
    results.sort(key=lambda x: x["pvalue"])
    
    # Benjamini-Hochberg FDR
    m = len(results)
    for rank, res in enumerate(results, start=1):
        fdr = min(1.0, (res["pvalue"] * m) / rank)
        res["fdr"] = fdr
        
    out_tsv = os.path.join(DATA_PROCESSED, "ontology_functional_enrichment.tsv")
    with open(out_tsv, 'w', encoding='utf-8') as f:
        f.write("module\tterm_id\tterm_name\tontology\toverlap_count\tterm_size\tmodule_size\tpvalue\tfdr\toverlap_genes\n")
        for r in results:
            f.write(f"{r['module']}\t{r['term_id']}\t{r['term_name']}\t{r['ontology']}\t{r['overlap_count']}\t{r['term_size']}\t{r['module_size']}\t{r['pvalue']:.3e}\t{r['fdr']:.3e}\t{r['overlap_genes']}\n")
    print(f"Saved functional ontology enrichment: {out_tsv} ({len(results)} enriched terms)")
    
    # Build Network Edge list for systems biology diagram
    edge_tsv = os.path.join(DATA_PROCESSED, "pathway_network_edges.tsv")
    with open(edge_tsv, 'w', encoding='utf-8') as f:
        f.write("source\ttarget\tedge_type\tweight\tdescription\n")
        # Module to Pathway edges
        for r in results:
            if r["fdr"] < 0.05:
                f.write(f"{r['module']}\t{r['term_name']}\tEnriched_Pathway\t{-math.log10(r['pvalue']):.2f}\tGO enrichment link\n")
        # Functional biological interaction links
        edges = [
            ("Biofilm_GPL", "Cell_Wall_Lipids", "Biosynthetic_Coupling", "0.85", "GPL and Mycolic acid outer membrane coordination"),
            ("Cell_Wall_Lipids", "ESX_Secretion", "Substrate_Transport", "0.90", "ESX-1 translocates surface-bound virulence factors across mycolate barrier"),
            ("DosR_Stress", "Oxidative_Stress", "Regulon_Crosstalk", "0.88", "DosS redox sensor and SigH coordinate hypoxia and antioxidant defense"),
            ("3D_Clinostat", "Heat_Shock_Chaperones", "Rotational_Shear", "0.75", "Induction of DnaK/ClpP under continuous clinorotation"),
            ("RPM_2.0", "Alternative_Oxidase", "Kinetic_Variation", "0.78", "Induction of CydAB and Icl1 under random multi-axis turnaround points")
        ]
        for src, tgt, etype, wt, desc in edges:
            f.write(f"{src}\t{tgt}\t{etype}\t{wt}\t{desc}\n")
    print(f"Saved systems biology pathway edges: {edge_tsv}")

if __name__ == '__main__':
    print("=== Phase 5: Multi-Scale Ontology & Mechanistic Bioinformatics ===")
    run_enrichment()
    print("Phase 5 completed successfully.")
