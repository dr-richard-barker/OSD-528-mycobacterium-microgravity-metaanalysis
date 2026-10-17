# OSD-528 re-analysis pipeline. `make all` runs every step in order.
# Step 00 downloads ~13 GB of raw reads from NASA OSDR on first run (reused afterwards).
# The classification step uses its own virtual environment: make ml-env

PY      := python3
ML_PY   := .venv-ml/bin/python

.PHONY: all quantify de wgcna ml ml-env go kegg osdr figures package numbers manuscript

all: quantify de wgcna ml go kegg osdr figures package numbers manuscript

quantify:            # kallisto, paired-end, CDS + ncRNA decoys; resolves the RFPRPM4 deposit defect
	$(PY) analysis/00_quantify_real_osdr_rnaseq.py

de:                  # DESeq2 + apeglm, primary and sensitivity analyses
	Rscript analysis/02_differential_expression.R

wgcna:               # signed WGCNA on the primary samples
	Rscript analysis/03_wgcna_coexpression_network.R

ml-env:
	$(PY) -m venv .venv-ml && $(ML_PY) -m pip install -r requirements-ml.txt

ml:                  # TabPFN and baselines, leave-one-out CV
	$(ML_PY) analysis/04_tabpfn_tabular_foundation_ai.py

go:                  # GO over-representation (UniProt annotations, go-basic)
	$(PY) analysis/05_ontology_functional_enrichment.py

kegg:                # KEGG-mapped metabolic subsystems
	$(PY) analysis/07_cellular_metabolic_modeling.py

osdr:                # OSDR microbial study compendium (Table S16)
	$(PY) analysis/01_fetch_osdr_microbial_data.py --compendium-only

figures:
	$(PY) analysis/06_generate_figures.py

package:             # results_tables/ and supplementary_figures/ READMEs
	$(PY) analysis/08_package_results_and_supplement.py

numbers:             # manuscript/chapters/numbers.tex from the outputs
	$(PY) analysis/09_manuscript_numbers.py

manuscript:          # PDFs (main, npj layout, supplementary), HTML and DOCX
	$(PY) manuscript/build_manuscript.py
