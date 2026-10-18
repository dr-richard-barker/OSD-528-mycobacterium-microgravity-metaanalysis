# Supplementary tables

Written by `analysis/08_package_results_and_supplement.py` from the pipeline outputs in
`data/processed/`. Row counts are counted from the files.

| Table | Rows | Description | Source |
|---|---:|---|---|
| `Table_S1_kallisto_estimated_counts.tsv` | 5,539 rows | kallisto estimated read-pair counts per CDS, all nine deposited samples (RFPRPM4 = inferred second library) | `data/processed/osd528_counts_raw.tsv` |
| `Table_S2_kallisto_tpm.tsv` | 5,539 rows | Transcripts per million per CDS, all nine samples | `data/processed/osd528_tpm.tsv` |
| `Table_S3_samples_and_sequencing_qc.tsv` | 9 rows | Sample groups, primary/sensitivity use, library index, read pairs and read assignment (rRNA, other ncRNA, CDS) | `config/samples.tsv` + `osd528_sequencing_qc_stats.json` |
| `Table_S4_deseq2_clinostat_vs_static1g.tsv` | 3,968 rows | DESeq2: 3D clinostat vs static 1g (primary samples); MLE and apeglm-shrunk log2FC, Wald p, BH FDR | `data/processed/deg_3dclinostat_vs_static1g.tsv` |
| `Table_S5_deseq2_rpm2_vs_static1g.tsv` | 3,968 rows | DESeq2: RPM 2.0 vs static 1g (primary samples) | `data/processed/deg_rpm2_vs_static1g.tsv` |
| `Table_S6_deseq2_clinostat_vs_rpm2.tsv` | 3,968 rows | DESeq2: 3D clinostat vs RPM 2.0 (primary samples) | `data/processed/deg_3dclinostat_vs_rpm2.tsv` |
| `Table_S7_shared_simulated_microgravity_genes.tsv` | 533 rows | Genes differentially expressed under both simulators in the same direction | `data/processed/deg_shared_simulated_microgravity.tsv` |
| `Table_S8_deseq2_sensitivity_with_rfprpm4.tsv` | 11,910 rows | DESeq2 results for all three contrasts with the inferred RFPRPM4 library included (sensitivity analysis) | `data/processed/sensitivity_with_rfprpm4/deg_*.tsv` |
| `Table_S9_wgcna_module_assignments.tsv` | 3,968 rows | WGCNA module, module membership (kME) and connectivity (kTotal, kWithin) per gene | `data/processed/wgcna_module_assignments.tsv` |
| `Table_S10_wgcna_module_eigengenes.tsv` | 8 rows | Module eigengenes per primary sample | `data/processed/wgcna_module_eigengenes.tsv` |
| `Table_S11_wgcna_module_trait_correlations.tsv` | 63 rows | Module eigengene vs trait Pearson correlation, Student p and BH FDR | `data/processed/wgcna_module_trait_correlations.tsv` |
| `Table_S12_go_enrichment_deg_sets.tsv` | 3,201 rows | GO over-representation (hypergeometric, BH FDR) for up and down DEG sets and the shared response | `data/processed/go_enrichment_deg.tsv` |
| `Table_S13_go_enrichment_wgcna_modules.tsv` | 8,759 rows | GO over-representation for each WGCNA module | `data/processed/go_enrichment_modules.tsv` |
| `Table_S14_classification_loocv_predictions.tsv` | 272 rows | Leave-one-out held-out predictions and class probabilities (TabPFN, logistic regression, random forest) | `data/processed/ml_loocv_predictions.tsv` |
| `Table_S15_classification_metrics.tsv` | 32 rows | LOOCV accuracy, balanced accuracy, AUROC and label-permutation p-values per model, task and feature count | `data/processed/ml_metrics.tsv` |
| `Table_S16_pan_microbial_osdr_compendium_102studies.tsv` | 102 rows | Microbial spaceflight and simulated-microgravity studies in NASA OSDR, from the OSDR API | `analysis/01_fetch_osdr_microbial_data.py` |
| `Table_S17_kegg_mapped_reactions.tsv` | 28 rows | Curated reactions with their KEGG-supported M. marinum genes and per-contrast expression changes | `data/processed/metabolic_model_reactions.tsv` |
| `Table_S18_kegg_subsystem_summary.tsv` | 16 rows | Per-subsystem, per-contrast summaries of the KEGG-mapped reactions | `data/processed/metabolic_subsystem_perturbation.tsv` |
| `Table_S19_metabolic_model.json` | JSON | Reaction model with KEGG-supported gene links (JSON) | `data/processed/metabolic_model_sbml.json` |
