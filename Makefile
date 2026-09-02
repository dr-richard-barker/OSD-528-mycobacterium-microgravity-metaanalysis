# Root Makefile for OSD-528 Meta-Analysis Pipeline

.PHONY: all ingest de wgcna tabpfn ontology figures manuscript clean

all: ingest de wgcna tabpfn ontology figures manuscript

ingest:
	python3 analysis/01_fetch_osdr_microbial_data.py

de:
	python3 analysis/02_differential_expression.py

wgcna:
	python3 analysis/03_wgcna_coexpression_network.py

tabpfn:
	python3 analysis/04_tabpfn_tabular_foundation_ai.py

ontology:
	python3 analysis/05_ontology_functional_enrichment.py

figures:
	python3 analysis/06_generate_figures.py

manuscript:
	python3 manuscript/build_manuscript.py

clean:
	rm -rf data/processed/*.tsv manuscript/figures/*.png manuscript/figures/*.svg manuscript/*.aux manuscript/*.log
