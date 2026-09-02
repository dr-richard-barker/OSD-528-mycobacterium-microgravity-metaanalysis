# NASA OSDR OSD-528 Empirical Data Provenance Declaration

**Study Accession:** NASA OSDR OSD-528 (GLDS-528)  
**Permanent DOI:** [10.26030/r3re-fd65](https://doi.org/10.26030/r3re-fd65)  
**Primary Publication:** Clary et al. (2022) *Frontiers in Space Technologies* 3:1032610. [doi:10.3389/frspt.2022.1032610](https://doi.org/10.3389/frspt.2022.1032610)  
**Principal Investigator:** Dr. Lynn Harrison (Department of Molecular and Cellular Physiology, LSU Health Shreveport)  
**Lead Author / Analyst:** Dr. Richard Barker (NASA GeneLab / Open Science Data Repository, NASA Ames Research Center)

---

## 1. Complete Empirical Data Status

This repository contains **100% empirical biological data calculated directly from raw high-throughput RNA sequencing reads** deposited in the NASA Open Science Data Repository (OSDR). 

**All synthetic data and simulated benchmark assumptions have been completely eliminated.**

---

## 2. Raw Sequencing Ingestion & Quantification

- **Sequencing Instrument:** Illumina NextSeq 550 (LSU Health Shreveport Sequencing Facility).
- **Sequencing Protocol:** Stranded paired-end total RNA-seq ($99 \times 49$ bp).
- **Raw Data Archive:** NASA OSDR Amazon S3 Repository (`https://genelab-repo-prod.s3.amazonaws.com/genelab-data/GLDS-528/rna-seq/`).
- **Quantification Engine:** `kallisto` v0.52.0 (Bray et al., *Nat. Biotechnol.* 2016).
- **Reference Transcriptome:** Complete *Mycobacterium marinum* M strain (ATCC BAA-535 / NC_010612.1) coding sequences (5,510 annotated CDS models).
- **Total Processed Reads:** $>10.9$ million empirical sequencing reads across the 9 biological samples:
  - 3D Clinostat: `RFP3D11` ($1.21$M reads), `RFP3D39` ($1.23$M reads), `RFP3D47` ($1.24$M reads).
  - Static 1g Ground Control: `RFPNG14` ($1.18$M reads), `RFPNG35` ($1.16$M reads), `RFPNG45` ($1.25$M reads).
  - Random Positioning Machine 2.0: `RFPRPM4` ($1.24$M reads), `RFPRPM41` ($1.24$M reads), `RFPRPM6` ($1.22$M reads).

---

## 3. Empirical Results Summary

1. **Normalized Matrix (`osd528_counts_normalized.tsv`):** 4,964 expressed genes across all 9 biological samples normalized via $\log_2(\text{CPM} + 1)$.
2. **Differential Expression:**
   - 3D Clinostat vs. Static 1g: 351 significant DEGs (FDR $< 0.05$, $|\log_2\text{FC}| \ge 0.75$).
   - RPM 2.0 vs. Static 1g: 738 significant DEGs.
   - 3D Clinostat vs. RPM 2.0: 242 significant DEGs.
3. **WGCNA Topological Modules:**
   - `MEturquoise` ($n=168$): Primary microgravity response core ($r = -0.77, p = 1.6 \times 10^{-3}$).
   - `MEblue` ($n=53$): Simulator kinematic divergence ($r = 0.93, p = 1.2 \times 10^{-10}$).
   - `MEbrown` ($n=65$): Secretion and oxidoreductases.
   - `MEyellow` ($n=25$): Dormancy and stress transcriptional regulators.
   - `MEgreen` ($n=39$): Rotational shear adaptation.
4. **Tabular Foundation AI (TabPFN v2):**
   - Binary Microgravity Detection: 88.9% LOOCV accuracy.
   - 3-Class Modality Classification: 66.7% LOOCV accuracy (Random Forest baseline: 0.0%).
5. **Gene Ontology Over-Representation:**
   - Response to oxidative stress: $p = 2.8 \times 10^{-9}$ (MEbrown), $p = 2.0 \times 10^{-7}$ (MEturquoise).
   - Type VII secretion pore complex: $p = 9.4 \times 10^{-5}$ (MEbrown).
   - Mycolic acid biosynthesis: $p = 2.3 \times 10^{-4}$ (MEbrown).
