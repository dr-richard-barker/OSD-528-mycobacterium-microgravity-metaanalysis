# Data provenance: NASA OSDR OSD-528 re-analysis

This file records where every result in this repository comes from, and the repository's history,
including the parts of it that were wrong.

## Source data

- **Study:** NASA Open Science Data Repository OSD-528 (GLDS-528), *Development of an Inexpensive 3D
  Clinostat and Comparison with Other Microgravity Simulators using Mycobacterium marinum*.
  Dataset DOI: [10.26030/r3re-fd65](https://doi.org/10.26030/r3re-fd65).
- **Original study:** Clary et al. 2022 (see `manuscript/references.bib`, key `clary2022development`).
  The experiment, sequencing and OSDR deposit are theirs. This repository is an independent re-analysis
  of the deposited data by Richard Barker (Collaborative Science Environment Inc.).
- **Design (from the OSDR study description):** *M. marinum* expressing red fluorescent protein, grown in
  biofilm-promoting medium in flaskettes holding a piece of PDMS membrane, on an RPM 2.0, on a
  lab-designed 3D clinostat, or on the incubator shelf beside them (static 1g), for 4 days at 31 °C.
  RNA was prepared from the bacteria **in suspension**, so the transcriptomes describe planktonic cells
  grown alongside the membrane, not the membrane biofilm.
- **Raw reads:** 18 paired-end FASTQ files (9 samples) from the OSDR file API. Each was checked against the
  deposited `GLDS-528_rna-seq_raw_md5sum.txt` before use.

## The RFPRPM4 deposit defect

The deposited RFPRPM4 R1 and R2 files each consist of the complete RFPRPM41 file followed by a second
library with its own index. `analysis/00_quantify_real_osdr_rnaseq.py` proves this by hashing and splits
the files. The second library is treated as an **inferred** RFPRPM4: it is excluded from the primary
analysis (RPM n = 2) and used only in the sensitivity analysis (`config/samples.tsv`,
`data/processed/sensitivity_with_rfprpm4/`). A draft report to the OSDR curators is in
`docs/OSDR_RFPRPM4_deposit_report.md`; it has not been sent.

## How each result is produced

| Step | Script | Main outputs |
|---|---|---|
| Quantification: kallisto, paired-end, reverse-stranded; index = all CDS of RefSeq GCF_000018345.1 plus its rRNA/tRNA/ncRNA genes as decoys | `analysis/00_quantify_real_osdr_rnaseq.py` | `data/processed/osd528_counts_raw.tsv`, `osd528_tpm.tsv`, `osd528_sequencing_qc_stats.json` |
| Differential expression: DESeq2 + apeglm | `analysis/02_differential_expression.R` | `deg_*.tsv`, `de_summary.json`, `pca_scores.tsv`, `osd528_vst.tsv` |
| Co-expression network: WGCNA, signed | `analysis/03_wgcna_coexpression_network.R` | `wgcna_*.tsv`, `wgcna_summary.json` |
| Classification: TabPFN and baselines, leave-one-out CV | `analysis/04_tabpfn_tabular_foundation_ai.py` | `ml_*.tsv`, `ml_summary.json` |
| GO over-representation: UniProt annotations, go-basic | `analysis/05_ontology_functional_enrichment.py` | `go_*.tsv`, `go_summary.json` |
| Figures | `analysis/06_generate_figures.py` | `manuscript/figures/` |
| KEGG-mapped metabolic subsystems | `analysis/07_cellular_metabolic_modeling.py` | `metabolic_*.tsv`, `metabolic_model_sbml.json` |
| OSDR microbial study compendium | `analysis/01_fetch_osdr_microbial_data.py` | `results_tables/Table_S16_*.tsv` |
| Supplementary tables and figures | `analysis/08_package_results_and_supplement.py` | `results_tables/`, `supplementary_figures/` |

Every number in the manuscript, README and figures should trace to one of these output files. The
summary JSON files (`de_summary.json`, `wgcna_summary.json`, `go_summary.json`, `ml_summary.json`) hold
the headline values.

## Repository history

- **8972bd8 (2 Sep 2026):** the first version was generated from **synthetic data**. Expression values
  came from `random.gauss`/`random.uniform`, and the commit included a `SYNTHETIC_DATA_DECLARATION.md`.
- **4e265c1 (2 Sep 2026, 16 minutes later):** described as "100% empirical", but it quantified only the
  first 60 MB of each R1 file (single-end). It read RFPRPM4 as if it were a separate sample, so RFPRPM4
  duplicated RFPRPM41. The differential-expression and network steps used custom code. Figures were
  hand-placed TikZ with literal numbers, the GO "enrichment" p-values and term sizes were literals, the
  "TabPFN" model was a nearest-neighbour stand-in, and the metabolic model put 24 of its 32 reactions on
  genes KEGG does not support (`data/processed/metabolic_model_locus_audit.tsv`). Later commits to 6 Sep changed formatting, figures and metadata but kept these problems.
- **16–17 Oct 2026:** an ABAI evidence check found the problems above, and the analysis was rebuilt from
  the raw reads: every step in the table was re-implemented, and every figure and table is now
  regenerated from outputs. Fabricated DOIs (a placeholder Zenodo DOI and a fake journal DOI) were
  removed in pull requests #2 and #3.

## Known limitations

- Small samples: 3 ground, 3 clinostat and 2 RPM samples in the primary analysis.
- 88–95% of read pairs are rRNA (which suggests little or no rRNA depletion), leaving about 0.15–0.8
  million coding read pairs per sample (`osd528_sequencing_qc_stats.json`).
- GO annotations for *M. marinum* are almost all electronic (IEA); see `go_summary.json`.
- No Zenodo deposit exists yet, so there is no archive DOI for this repository.
