# Synthetic Data & Benchmark Provenance Declaration

## 1. Context and Rationale

NASA Open Science Data Repository (OSDR) study **OSD-528** contains raw paired-end Illumina RNA-sequencing reads (18 `.fastq.gz` files across 9 biological samples, totaling $\sim 15$ GB) and raw MultiQC reports. However, unlike certain other OSDR accessions, **no pre-computed GeneLab gene expression count matrix or differential expression table was deposited with the study files**.

To develop, validate, and demonstrate the integrated computational pipeline coupling:
- Variance-stabilizing transformation (VST) and median-of-ratios normalization,
- Weighted Gene Co-expression Network Analysis (WGCNA),
- Prior-data Fitted Tabular Foundation AI (TabPFN, *Nature* 2025), and
- Gene Ontology / pathway over-representation analysis,

the default files generated in `data/processed/osd528_counts_normalized.tsv` and downstream tables were constructed using an **empirically parameterized synthetic simulation matrix**.

---

## 2. Mathematical Formulation of the Simulation

The benchmark matrix was parameterized directly from the published biological findings of the primary OSD-528 study (*Clary et al., 2022, Frontiers in Space Technologies*, [doi:10.3389/frspt.2022.1032610](https://doi.org/10.3389/frspt.2022.1032610)) and canonical *Mycobacterium marinum* physiology:

- **Baseline Expression**: Each of the 1,200 simulated gene models drew a basal expression level $E_0 \sim \mathcal{U}(8.5, 12.5)$ (representing $\log_2$ VST-scaled counts).
- **Biological Effect Sizes ($\Delta_{fc}$)**: Literature-supported effect sizes were assigned to key functional regulons:
  - Glycopeptidolipid (GPL) biofilm synthesis (*mps1*, *mps2*, *fmt*, *fadD28*): $+1.7$ to $+2.6 \log_2\text{FC}$
  - Mycolic acid FAS-II cycle (*kasA*, *kasB*, *acpM*, *inhA*, *fbpA*): $+1.5$ to $+2.5 \log_2\text{FC}$
  - Type VII ESX-1 secretion (*esxA*, *esxB*, *eccA1*–*eccE1*): $+1.7$ to $+2.7 \log_2\text{FC}$
  - DosR dormancy & antioxidant defense (*dosR*, *dosS*, *hspX*, *tgs1*, *katG*, *sodA*): $+1.8$ to $+3.2 \log_2\text{FC}$
  - Simulator-divergent markers (*clpP1*, *dnaK* in clinostat vs. *cydA*, *icl1* in RPM): $+1.8$ to $+2.4 \log_2\text{FC}$
  - Downregulated central respiration (*nuoA*, *nuoB*, *gltA*, *atpA*): $-0.7$ to $-1.4 \log_2\text{FC}$
  - Background genome features: Gaussian noise $\Delta \sim \mathcal{N}(0, 0.35)$.
- **Replicate Biological Variance**: Individual replicate samples received independent stochastic variation $\epsilon \sim \mathcal{N}(0, 0.12)$.

---

## 3. Production Pipeline for Empirical FASTQ Quantification

For researchers wishing to process the raw $\sim 15$ GB sequencing reads deposited on NASA OSDR, we provide a complete, automated end-to-end bioinformatic script:
`analysis/00_download_and_process_raw_rnaseq.sh`

This script executes:
1. Programmatic downloading of all 18 FASTQ pairs from the OSDR S3/HTTP endpoints.
2. Adapter and quality trimming via `fastp`.
3. Alignment to the *Mycobacterium marinum* ATCC BAA-535 / M strain reference genome (`NC_010612.1`) using `hisat2`.
4. Feature quantification via `featureCounts` (Subread package).
5. Automatic replacement of `data/processed/osd528_counts_normalized.tsv` with the empirical read matrix.

All downstream analytical scripts (`02_differential_expression.py`, `03_wgcna_coexpression_network.py`, `04_tabpfn_tabular_foundation_ai.py`, `05_ontology_functional_enrichment.py`, and `06_generate_figures.py`) automatically detect and operate seamlessly on either the empirical count matrix or the synthetic benchmark.
