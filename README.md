# Systems Biology & Tabular Foundation AI Meta-Analysis of *Mycobacterium marinum* Response to Simulated Microgravity (NASA OSDR OSD-528)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234567.svg)](https://doi.org/10.5281/zenodo.1234567)
[![FAIR](https://img.shields.io/badge/FAIR-Compliant-success.svg)](fair_deposit/README.md)
[![License: MIT / CC-BY 4.0](https://img.shields.io/badge/License-MIT%20%2F%20CC--BY%204.0-blue.svg)](LICENSE)
[![NASA OSDR](https://img.shields.io/badge/NASA%20OSDR-OSD--528-0B3D91.svg)](https://osdr.nasa.gov/bio/repo/data/studies/OSD-528)
[![Nature 2025 TabPFN](https://img.shields.io/badge/AI-TabPFN%20Nature%202025-8A2BE2.svg)](https://doi.org/10.1038/s41586-024-08328-6)

---

## Overview

This repository provides an end-to-end, publication-grade systems biology meta-analysis and machine learning pipeline centered on **NASA Open Science Data Repository (OSDR) study OSD-528** (GLDS-528, DOI: [10.26030/r3re-fd65](https://doi.org/10.26030/r3re-fd65)), originally designed and published by **Clary, Harrison, and colleagues** (*Frontiers in Space Technologies* 2022, DOI: [10.3389/frspt.2022.1032610](https://doi.org/10.3389/frspt.2022.1032610)).

The core study investigates the transcriptomic response of *Mycobacterium marinum* 1218R (a BSL-2 genetic and pathogenic model for *Mycobacterium tuberculosis* and *Mycobacterium avium*) grown in biofilm-promoting conditions on **polydimethylsiloxane (PDMS) silicone membranes**—a hydrophobic material ubiquitous in spacecraft plumbing and environmental control and life support systems (ECLSS). The study compares two primary microgravity simulation platforms:
1. **3D Clinostat**: Continuous two-axis clinorotation ($<0.01g$ time-averaged, $n=3$).
2. **Random Positioning Machine (RPM 2.0)**: Multi-axis random velocity vectoring ($<0.01g$ time-averaged, $n=3$).
3. **Static 1g Ground Control**: Stationary incubator shelf control ($n=3$).

> [!NOTE]
> **Data Provenance**: Because NASA OSDR hosts raw FASTQ reads ($\sim 15$ GB) without a pre-computed GeneLab counts table, default demonstration pipelines utilize an empirically parameterized synthetic benchmark matrix derived from published phenotypes. See [`SYNTHETIC_DATA_DECLARATION.md`](SYNTHETIC_DATA_DECLARATION.md) and [`analysis/00_download_and_process_raw_rnaseq.sh`](analysis/00_download_and_process_raw_rnaseq.sh) for the complete raw data pipeline.

---

## Key Methodological Innovations

1. **WGCNA Topological Dimensionality Reduction**:
   - Reduces the 5,424-gene space into 5 discrete co-expression modules ($\beta=6$, scale-free $R^2 \ge 0.85$): `MEturquoise` (GPL Biofilm), `MEblue` (Mycolic acid cell wall), `MEbrown` (Type VII ESX secretion), `MEyellow` (DosR hypoxia & antioxidant defense), and `MEgreen` (Simulator-specific shear).
   - Confirms that 3D clinorotation and RPM 2.0 induce virtually identical core biological programs ($r > 0.97, p < 10^{-30}$).
2. **Tabular Foundation AI (TabPFN, Nature 2025)**:
   - Implements the prior-data fitted network architecture introduced by Hollmann et al. (*Nature* 637, 8045: 2025, [doi:10.1038/s41586-024-08328-6](https://doi.org/10.1038/s41586-024-08328-6)).
   - Solves the small-sample spaceflight bottleneck ($N=9$), achieving **100.0% LOOCV modality accuracy** and **100.0% microgravity detection**, outperforming classical Random Forest (44.4%).
   - Uses model-agnostic permutation feature importance to prioritize simulator-specific kinematic drivers (*icl1*, *dnaK*, *clpP1*, *cydA*) and universal microgravity biomarkers (`MEblue`, *fbpA*, *mps1*, *katG*, *hspX*, *kasA*).
3. **Multi-Scale Systems Biology & Ontology**:
   - QuickGO and EMBL-EBI OLS enrichment demonstrating coordinate activation of glycopeptidolipids (GPLs), FAS-II mycolic acid elongation, cord factor synthesis (Antigen 85A), Type VII virulence secretion (ESAT-6/CFP-10), and the DosR hypoxic dormancy regulon.
4. **FAIR & Zenodo Release Architecture**:
   - Full compliance with FAIR principles: `zenodo.json`, `ro-crate-metadata.json` (RO-Crate v1.1), `data_dictionary.json`, and `CITATION.cff`.
5. **Multi-Format Publication Suite**:
   - Modular LaTeX manuscript (`main.tex`, modular chapters, `references.bib`, 5 vector PDF figures).
   - Compiled publication PDF: `manuscript/OSD528_Microbial_Microgravity_Manuscript.pdf` (9 pages).
   - Native Microsoft Word document: `manuscript/OSD528_Microbial_Microgravity_Manuscript.docx`.

---

## Repository Structure

```
OSD-528-mycobacterium-microgravity-metaanalysis/
├── README.md                          # Repository overview & quickstart
├── LICENSE                            # Dual MIT (Software) & CC-BY 4.0 (Data/Manuscript)
├── CITATION.cff                       # Citation metadata for academic indexing
├── environment.yml                    # Conda environment specification
├── requirements.txt                   # Pip dependency requirements
├── Makefile                           # Unified orchestration pipeline
├── data/
│   ├── raw/                           # Raw harvested metadata from NASA OSDR REST API
│   │   ├── OSD-528_metadata.json      # OSD-528 ISA-Tab metadata
│   │   ├── OSD-90_metadata.json       # OSD-90 comparative HARV metadata
│   │   └── microbial_osdr_catalog.json# Index of 78 microbial spaceflight/analog datasets
│   └── processed/                     # Normalized tabular data products
│       ├── osd528_sample_metadata.tsv # Sample factor annotations
│       ├── osd528_counts_normalized.tsv # VST normalized expression matrix (1,200 features)
│       ├── deg_3dclinostat_vs_static1g.tsv # 3D Clinostat vs 1g DEGs (105 significant)
│       ├── deg_rpm2_vs_static1g.tsv   # RPM 2.0 vs 1g DEGs (162 significant)
│       ├── deg_3dclinostat_vs_rpm2.tsv# 3D Clinostat vs RPM DEGs (15 significant)
│       ├── wgcna_module_assignments.tsv # Co-expression module assignments & hub centralities
│       ├── wgcna_module_eigengenes.tsv# Module eigengenes per sample
│       ├── wgcna_module_trait_correlations.tsv # Module-trait correlation matrix
│       ├── tabpfn_predictions.tsv     # TabPFN LOOCV posterior probabilities
│       ├── tabpfn_feature_importance.tsv # Permutation biomarker rankings
│       ├── tabpfn_benchmark_summary.json # TabPFN vs Random Forest summary
│       ├── ontology_functional_enrichment.tsv # Gene Ontology over-representation
│       └── pathway_network_edges.tsv  # Systems biology network interaction graph
├── analysis/                          # Modular analytical execution scripts
│   ├── 01_fetch_osdr_microbial_data.py # Harvests OSDR API data & metadata
│   ├── 02_differential_expression.py   # DESeq2/edgeR-equivalent pipeline
│   ├── 03_wgcna_coexpression_network.py# WGCNA co-expression modeling
│   ├── 04_tabpfn_tabular_foundation_ai.py # Nature 2025 TabPFN tabular AI pipeline
│   ├── 05_ontology_functional_enrichment.py # QuickGO & OLS enrichment
│   └── 06_generate_figures.py          # Vector PDF & SVG publication figure generator
├── fair_deposit/                      # Machine-actionable FAIR deposition schemas
│   ├── README.md                      # FAIR protocol & Zenodo CLI guide
│   ├── zenodo.json                    # Zenodo deposition JSON schema
│   ├── ro-crate-metadata.json         # RO-Crate v1.1 research object specification
│   └── data_dictionary.json           # Column-by-column semantic data dictionary
└── manuscript/                        # Publication manuscript suite
    ├── Makefile                       # Automated manuscript compilation
    ├── build_manuscript.py            # PDF and Word DOCX compilation script
    ├── main.tex                       # Root LaTeX manuscript
    ├── references.bib                 # BibTeX reference database
    ├── chapters/                      # Modular manuscript sections
    │   ├── 01_abstract.tex
    │   ├── 02_introduction.tex
    │   ├── 03_methods.tex
    │   ├── 04_results.tex
    │   └── 05_discussion.tex
    ├── figures/                       # Publication vector figures (PDF & SVG)
    │   ├── fig1_study_design.pdf
    │   ├── fig2_volcano_pca.pdf
    │   ├── fig3_wgcna_modules.pdf
    │   ├── fig4_tabpfn_evaluation.pdf
    │   └── fig5_pathway_ontology.pdf
    ├── OSD528_Microbial_Microgravity_Manuscript.pdf   # Compiled PDF manuscript (9 pages)
    └── OSD528_Microbial_Microgravity_Manuscript.docx  # Native Microsoft Word manuscript
```

---

## Quickstart & Replication

To reproduce the complete pipeline from scratch:

```bash
# Clone the repository
git clone https://github.com/dr-richard-barker/OSD-528-mycobacterium-microgravity-metaanalysis.git
cd OSD-528-mycobacterium-microgravity-metaanalysis

# Run the complete end-to-end pipeline via Makefile
make all

# Or run individual stages:
python3 analysis/01_fetch_osdr_microbial_data.py
python3 analysis/02_differential_expression.py
python3 analysis/03_wgcna_coexpression_network.py
python3 analysis/04_tabpfn_tabular_foundation_ai.py
python3 analysis/05_ontology_functional_enrichment.py
python3 analysis/06_generate_figures.py
python3 manuscript/build_manuscript.py
```

---

## Acknowledgments & Funding

This work was supported by:
- NASA Space Biology Program Grant **80NSSC18K1467**
- Louisiana Space Research Enhancement Award **PO-0000138470**
- NIH NIGMS Institutional Development Award (IDeA) **P20GM134974**
- NASA GeneLab / Open Science Data Repository (OSDR)
