# NASA OSDR OSD-528: Comprehensive Empirical Results Tables

This directory contains the complete set of empirical data tables generated for the meta-analysis of *Mycobacterium marinum* response to simulated microgravity:

| File Name | Description | Rows / Features | Primary Assays / Methods |
| :--- | :--- | :--- | :--- |
| **`Table_S1_raw_counts_matrix.tsv`** | Raw transcript count matrix | 5,510 genes $\times$ 9 samples | kallisto v0.52.0 pseudoalignment of NextSeq 550 paired reads |
| **`Table_S2_normalized_log2cpm_expression.tsv`** | Normalized $\log_2(\text{CPM}+1)$ matrix | 4,964 expressed genes $\times$ 9 samples | TMM library size scaling + log-CPM transformation |
| **`Table_S3_sample_metadata.tsv`** | NASA OSDR sample experimental metadata | 9 samples $\times$ 5 attributes | 3D Clinostat, RPM 2.0, Static 1g; PDMS silicone membranes |
| **`Table_S4_differential_expression_3dclinostat_vs_1g.tsv`** | Differential expression: 3D Clinostat vs 1g | 5,510 genes (351 sig DEGs) | Negative binomial Wald test, Benjamini-Hochberg FDR |
| **`Table_S5_differential_expression_rpm2_vs_1g.tsv`** | Differential expression: RPM 2.0 vs 1g | 5,510 genes (738 sig DEGs) | Negative binomial Wald test, Benjamini-Hochberg FDR |
| **`Table_S6_differential_expression_3dclinostat_vs_rpm2.tsv`** | Differential expression: Simulator comparison | 5,510 genes (242 sig DEGs) | Contrast evaluating kinematic shear and rotational differences |
| **`Table_S7_wgcna_module_assignments_goslim.tsv`** | WGCNA module assignments & hub metrics | 350 top variable genes | Topological Overlap Matrix (TOM), $k_{\text{within}}$, $k_{\text{total}}$ |
| **`Table_S8_wgcna_module_eigengenes.tsv`** | Module Eigengene expression values | 5 modules $\times$ 9 samples | First principal component per co-expression module |
| **`Table_S9_wgcna_module_trait_correlations.tsv`** | Module-trait correlation coefficients | 5 modules $\times$ 3 traits | Pearson $r$, student asymptotic $p$-values |
| **`Table_S10_goslim_pathway_enrichment.tsv`** | GOSlim pathway over-representation | Core enriched pathways | Hypergeometric test, FDR adjustment, gene overlap counts |
| **`Table_S11_full_gene_ontology_enrichment.tsv`** | Full Gene Ontology enrichment table | Biological Process, Molecular Function | QuickGO & EMBL-EBI OLS mappings |
| **`Table_S12_tabpfn_loocv_predictions.tsv`** | TabPFN foundation model LOOCV predictions | 9 sample validations | Calibrated class probabilities, true vs predicted labels |
| **`Table_S13_tabpfn_feature_importance.tsv`** | TabPFN permutation feature importances | 15 module & hub features | Mean log-likelihood decrease under feature permutation |
| **`Table_S14_cellular_metabolic_model_reactions.tsv`** | Genome-scale cellular metabolic reactions | 32 curated reactions | Compartments, equations, EC numbers, empirical $\log_2\text{FC}$ |
| **`Table_S15_metabolic_subsystem_perturbation_index.tsv`** | Subsystem Perturbation Index (SPI) ranking | 8 core physiological subsystems | Pathway vulnerability scores, predominant flux shifts |
| **`Table_S16_pan_microbial_osdr_compendium_78studies.tsv`** | Pan-microbial spaceflight meta-analysis index | 78 OSDR microbial datasets | Cross-species study compendium, phylum annotations |
| **`Table_S17_cellular_metabolic_model_sbml.json`** | SBML-compatible computational model JSON | Multi-compartment reconstruction | Reactions, subsystems, and expression perturbations |
