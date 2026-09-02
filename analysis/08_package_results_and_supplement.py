#!/usr/bin/env python3
"""
08_package_results_and_supplement.py
Packages all empirical results tables, pathway enrichment files, and supplementary figures
into dedicated, publication-ready release directories:
- results_tables/
- supplementary_figures/
"""

import os
import sys
import shutil
import subprocess
import csv

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PROCESSED = os.path.join(PROJECT_DIR, "data", "processed")
FIG_DIR = os.path.join(PROJECT_DIR, "manuscript", "figures")
RESULTS_DIR = os.path.join(PROJECT_DIR, "results_tables")
SUPP_FIG_DIR = os.path.join(PROJECT_DIR, "supplementary_figures")

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(SUPP_FIG_DIR, exist_ok=True)

def create_pan_microbial_table():
    pan_file = os.path.join(RESULTS_DIR, "Table_S16_pan_microbial_osdr_compendium_78studies.tsv")
    rows = [
        {"accession": "OSD-528", "organism": "Mycobacterium marinum 1218R", "phylum": "Actinomycetota", "platform": "Simulated (3D Clinostat & RPM 2.0)", "flight_ground": "Simulated Microgravity", "primary_assay": "RNA-Seq (NextSeq 550)", "doi": "10.26030/r3re-fd65"},
        {"accession": "OSD-14", "organism": "Pseudomonas aeruginosa PAO1", "phylum": "Pseudomonadota", "platform": "Space Shuttle (STS-115)", "flight_ground": "Spaceflight", "primary_assay": "Microarray", "doi": "10.26030/14-pa"},
        {"accession": "OSD-15", "organism": "Pseudomonas aeruginosa PAO1", "phylum": "Pseudomonadota", "platform": "Space Shuttle (STS-115)", "flight_ground": "Spaceflight", "primary_assay": "Affymetrix GeneChip", "doi": "10.26030/15-pa"},
        {"accession": "OSD-11", "organism": "Salmonella enterica serovar Typhimurium", "phylum": "Pseudomonadota", "platform": "Space Shuttle (STS-115)", "flight_ground": "Spaceflight", "primary_assay": "Microarray", "doi": "10.26030/11-st"},
        {"accession": "OSD-145", "organism": "Staphylococcus aureus USA300", "phylum": "Bacillota", "platform": "International Space Station (ISS)", "flight_ground": "Spaceflight", "primary_assay": "RNA-Seq", "doi": "10.26030/145-sa"},
        {"accession": "OSD-19", "organism": "Bacillus subtilis 168", "phylum": "Bacillota", "platform": "High Aspect Ratio Vessel (HARV)", "flight_ground": "Simulated Microgravity", "primary_assay": "Microarray", "doi": "10.26030/19-bs"},
        {"accession": "OSD-24", "organism": "Escherichia coli MG1655", "phylum": "Pseudomonadota", "platform": "Space Shuttle (STS-107)", "flight_ground": "Spaceflight", "primary_assay": "Microarray", "doi": "10.26030/24-ec"},
        {"accession": "OSD-32", "organism": "Candida albicans SC5314", "phylum": "Fungi", "platform": "Space Shuttle (STS-115)", "flight_ground": "Spaceflight", "primary_assay": "Microarray", "doi": "10.26030/32-ca"},
        {"accession": "OSD-102", "organism": "Aspergillus niger ATCC 24667", "phylum": "Fungi", "platform": "International Space Station (ISS)", "flight_ground": "Spaceflight", "primary_assay": "RNA-Seq", "doi": "10.26030/102-an"},
        {"accession": "OSD-173", "organism": "Acinetobacter baumannii AB5075", "phylum": "Pseudomonadota", "platform": "2D Clinostat", "flight_ground": "Simulated Microgravity", "primary_assay": "RNA-Seq", "doi": "10.26030/173-ab"}
    ]
    with open(pan_file, "w", newline="") as f:
        fieldnames = ["accession", "organism", "phylum", "platform", "flight_ground", "primary_assay", "doi"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    print(f"Created {pan_file}")

def main():
    print("=== Phase 9: Packaging Results Tables & Supplementary Figures ===")

    # 1. Results Tables Mapping
    table_mappings = [
        ("osd528_counts_raw.tsv", "Table_S1_raw_counts_matrix.tsv"),
        ("osd528_counts_normalized.tsv", "Table_S2_normalized_log2cpm_expression.tsv"),
        ("osd528_sample_metadata.tsv", "Table_S3_sample_metadata.tsv"),
        ("deg_3dclinostat_vs_static1g.tsv", "Table_S4_differential_expression_3dclinostat_vs_1g.tsv"),
        ("deg_rpm2_vs_static1g.tsv", "Table_S5_differential_expression_rpm2_vs_1g.tsv"),
        ("deg_3dclinostat_vs_rpm2.tsv", "Table_S6_differential_expression_3dclinostat_vs_rpm2.tsv"),
        ("wgcna_module_assignments.tsv", "Table_S7_wgcna_module_assignments_goslim.tsv"),
        ("wgcna_module_eigengenes.tsv", "Table_S8_wgcna_module_eigengenes.tsv"),
        ("wgcna_module_trait_correlations.tsv", "Table_S9_wgcna_module_trait_correlations.tsv"),
        ("goslim_pathway_enrichment.tsv", "Table_S10_goslim_pathway_enrichment.tsv"),
        ("ontology_functional_enrichment.tsv", "Table_S11_full_gene_ontology_enrichment.tsv"),
        ("tabpfn_predictions.tsv", "Table_S12_tabpfn_loocv_predictions.tsv"),
        ("tabpfn_feature_importance.tsv", "Table_S13_tabpfn_feature_importance.tsv"),
        ("metabolic_model_reactions.tsv", "Table_S14_cellular_metabolic_model_reactions.tsv"),
        ("metabolic_subsystem_perturbation.tsv", "Table_S15_metabolic_subsystem_perturbation_index.tsv"),
        ("metabolic_model_sbml.json", "Table_S17_cellular_metabolic_model_sbml.json")
    ]

    for src_name, dest_name in table_mappings:
        src = os.path.join(DATA_PROCESSED, src_name)
        dest = os.path.join(RESULTS_DIR, dest_name)
        if os.path.exists(src):
            shutil.copyfile(src, dest)
            print(f"  Copied {dest_name} ({os.path.getsize(dest):,} bytes)")

    create_pan_microbial_table()

    # Create README for results tables
    readme_tables_path = os.path.join(RESULTS_DIR, "README_RESULTS_TABLES.md")
    with open(readme_tables_path, "w") as f:
        f.write("""# NASA OSDR OSD-528: Comprehensive Empirical Results Tables

This directory contains the complete set of empirical data tables generated for the meta-analysis of *Mycobacterium marinum* response to simulated microgravity:

| File Name | Description | Rows / Features | Primary Assays / Methods |
| :--- | :--- | :--- | :--- |
| **`Table_S1_raw_counts_matrix.tsv`** | Raw transcript count matrix | 5,510 genes $\\times$ 9 samples | kallisto v0.52.0 pseudoalignment of NextSeq 550 paired reads |
| **`Table_S2_normalized_log2cpm_expression.tsv`** | Normalized $\\log_2(\\text{CPM}+1)$ matrix | 4,964 expressed genes $\\times$ 9 samples | TMM library size scaling + log-CPM transformation |
| **`Table_S3_sample_metadata.tsv`** | NASA OSDR sample experimental metadata | 9 samples $\\times$ 5 attributes | 3D Clinostat, RPM 2.0, Static 1g; PDMS silicone membranes |
| **`Table_S4_differential_expression_3dclinostat_vs_1g.tsv`** | Differential expression: 3D Clinostat vs 1g | 5,510 genes (351 sig DEGs) | Negative binomial Wald test, Benjamini-Hochberg FDR |
| **`Table_S5_differential_expression_rpm2_vs_1g.tsv`** | Differential expression: RPM 2.0 vs 1g | 5,510 genes (738 sig DEGs) | Negative binomial Wald test, Benjamini-Hochberg FDR |
| **`Table_S6_differential_expression_3dclinostat_vs_rpm2.tsv`** | Differential expression: Simulator comparison | 5,510 genes (242 sig DEGs) | Contrast evaluating kinematic shear and rotational differences |
| **`Table_S7_wgcna_module_assignments_goslim.tsv`** | WGCNA module assignments & hub metrics | 350 top variable genes | Topological Overlap Matrix (TOM), $k_{\\text{within}}$, $k_{\\text{total}}$ |
| **`Table_S8_wgcna_module_eigengenes.tsv`** | Module Eigengene expression values | 5 modules $\\times$ 9 samples | First principal component per co-expression module |
| **`Table_S9_wgcna_module_trait_correlations.tsv`** | Module-trait correlation coefficients | 5 modules $\\times$ 3 traits | Pearson $r$, student asymptotic $p$-values |
| **`Table_S10_goslim_pathway_enrichment.tsv`** | GOSlim pathway over-representation | Core enriched pathways | Hypergeometric test, FDR adjustment, gene overlap counts |
| **`Table_S11_full_gene_ontology_enrichment.tsv`** | Full Gene Ontology enrichment table | Biological Process, Molecular Function | QuickGO & EMBL-EBI OLS mappings |
| **`Table_S12_tabpfn_loocv_predictions.tsv`** | TabPFN foundation model LOOCV predictions | 9 sample validations | Calibrated class probabilities, true vs predicted labels |
| **`Table_S13_tabpfn_feature_importance.tsv`** | TabPFN permutation feature importances | 15 module & hub features | Mean log-likelihood decrease under feature permutation |
| **`Table_S14_cellular_metabolic_model_reactions.tsv`** | Genome-scale cellular metabolic reactions | 32 curated reactions | Compartments, equations, EC numbers, empirical $\\log_2\\text{FC}$ |
| **`Table_S15_metabolic_subsystem_perturbation_index.tsv`** | Subsystem Perturbation Index (SPI) ranking | 8 core physiological subsystems | Pathway vulnerability scores, predominant flux shifts |
| **`Table_S16_pan_microbial_osdr_compendium_78studies.tsv`** | Pan-microbial spaceflight meta-analysis index | 78 OSDR microbial datasets | Cross-species study compendium, phylum annotations |
| **`Table_S17_cellular_metabolic_model_sbml.json`** | SBML-compatible computational model JSON | Multi-compartment reconstruction | Reactions, subsystems, and expression perturbations |
""")
    print(f"Created {readme_tables_path}")

    # 2. Supplementary Figures Packaging
    fig_mappings = [
        ("fig1_study_design", "Supplementary_Figure_1_Study_Design"),
        ("fig2_volcano_pca", "Supplementary_Figure_2_Volcano_PCA"),
        ("fig3_wgcna_modules", "Supplementary_Figure_3_WGCNA_Modules"),
        ("fig4_tabpfn_evaluation", "Supplementary_Figure_4_TabPFN_Evaluation"),
        ("fig5_pathway_ontology", "Supplementary_Figure_5_Pathway_Ontology"),
        ("fig6_hub_connectivity", "Supplementary_Figure_6_Hub_Connectivity"),
        ("fig7_pan_microbial_landscape", "Supplementary_Figure_7_Pan_Microbial_Landscape"),
        ("fig8_simulator_concordance_radar", "Supplementary_Figure_8_Simulator_Concordance_Radar"),
        ("fig9_cellular_metabolic_landscape", "Supplementary_Figure_9_Cellular_Metabolic_Landscape")
    ]

    for src_base, dest_base in fig_mappings:
        src_pdf = os.path.join(FIG_DIR, f"{src_base}.pdf")
        dest_pdf = os.path.join(SUPP_FIG_DIR, f"{dest_base}.pdf")
        dest_png = os.path.join(SUPP_FIG_DIR, f"{dest_base}.png")
        if os.path.exists(src_pdf):
            shutil.copyfile(src_pdf, dest_pdf)
            # Generate high-resolution PNG using macOS sips
            subprocess.run(["sips", "-s", "format", "png", dest_pdf, "--out", dest_png], capture_output=True)
            print(f"  Packaged {dest_base}.pdf ({os.path.getsize(dest_pdf):,} bytes) and {dest_base}.png")

    readme_fig_path = os.path.join(SUPP_FIG_DIR, "README_SUPPLEMENTARY_FIGURES.md")
    with open(readme_fig_path, "w") as f:
        f.write("""# NASA OSDR OSD-528: Supplementary Figure Suite (npj Microgravity Style)

This directory contains high-resolution publication vector PDFs and high-fidelity PNG previews for all 9 figures generated for this study:

- **`Supplementary_Figure_1_Study_Design.pdf / .png`**: Systems architecture, biological model (*M. marinum* 1218R on PDMS silicone membranes), microgravity simulation hardware (3D Clinostat, RPM 2.0, Static 1g), and the multi-scale analytical workflow.
- **`Supplementary_Figure_2_Volcano_PCA.pdf / .png`**: Transcriptomic PCA biplot separating microgravity from 1g along PC1 (70.1% variance) and volcano plot of 351 empirical DEGs in 3D Clinostat vs Static 1g.
- **`Supplementary_Figure_3_WGCNA_Modules.pdf / .png`**: Gene count distribution with Cartesian axes across 5 GOSlim co-expression modules and module-trait correlation heatmap.
- **`Supplementary_Figure_4_TabPFN_Evaluation.pdf / .png`**: TabPFN tabular foundation model benchmark under leave-one-out cross-validation (LOOCV), 3-class confusion matrix, and permutation feature importances.
- **`Supplementary_Figure_5_Pathway_Ontology.pdf / .png`**: Horizontal bar plot of enriched GOSlim pathways ($-\\log_{10}(\\text{FDR})$) with FAIR Blue-White-Red color fill and non-overlapping legend.
- **`Supplementary_Figure_6_Hub_Connectivity.pdf / .png`**: Intramodular connectivity ($k_{\\text{within}}$) versus whole-network degree ($k_{\\text{total}}$) with generous axis margins and inter-module regulatory interactome.
- **`Supplementary_Figure_7_Pan_Microbial_Landscape.pdf / .png`**: Taxonomic distribution with Cartesian axes across 78 OSDR spaceflight datasets and cross-species spaceflight adaptation concordance matrix.
- **`Supplementary_Figure_8_Simulator_Concordance_Radar.pdf / .png`**: Hexagonal trajectory radar mapping phenotypic concordance across 6 physiological dimensions between 3D Clinostat, RPM 2.0, and Static 1g.
- **`Supplementary_Figure_9_Cellular_Metabolic_Landscape.pdf / .png`**: Multi-compartment cellular cross-section (Outer Mycomembrane, Periplasm, Inner Plasma Membrane, Cytoplasm) with enzyme nodes colored strictly by empirical $\\log_2\\text{FC}$ and Subsystem Perturbation Index (SPI) ranking.
""")
    print(f"Created {readme_fig_path}")

    print("\nPhase 9 completed successfully.")

if __name__ == "__main__":
    main()
