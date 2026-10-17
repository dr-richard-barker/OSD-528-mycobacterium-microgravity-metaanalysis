# FAIR metadata

Metadata for depositing this re-analysis in an archive such as Zenodo. **Nothing has been deposited yet**,
so there is no archive DOI.

| File | What it is | How it is made |
|---|---|---|
| `zenodo.json` | Zenodo deposition metadata: title, creator (with ORCID), description, licence, related identifiers (the OSD-528 dataset DOI and the original Clary et al. 2022 paper) | Maintained by hand |
| `data_dictionary.json` | Every supplementary table in `results_tables/` with its columns and their meaning | Written by `analysis/08_package_results_and_supplement.py` from the files' headers |
| `ro-crate-metadata.json` | RO-Crate 1.1 description of the supplementary tables, their author and licence, and the OSD-528 source dataset | Written by `analysis/08_package_results_and_supplement.py` |

Licences: code MIT (`LICENSE`); processed data, figures and text CC BY 4.0.

When depositing, upload the processed outputs (`results_tables/`, `data/processed/`), the code and the
manuscript. Do not upload `data/raw/` or `data/reference/`: they hold about 13 GB of raw reads and
reference files that anyone can re-download from OSDR and NCBI with `make quantify`.
