# FAIR Data Deposit & Zenodo Replication Protocol

This directory contains the machine-actionable metadata, schema definitions, and protocols for depositing the **OSD-528** meta-analysis into the **Zenodo Open Access Repository** in strict adherence to the **FAIR Principles** (Findable, Accessible, Interoperable, Reusable).

---

## 1. FAIR Compliance Matrix

| FAIR Principle | Implementation in this Repository |
| :--- | :--- |
| **Findable (F1–F4)** | Unique persistent identifiers (Zenodo DOI, OSDR accession linking `OSD-528` & `OSD-90`). Rich metadata registered in `zenodo.json` and `ro-crate-metadata.json`. Indexed by Google Dataset Search and OpenAIRE. |
| **Accessible (A1–A2)** | Unrestricted open access via standard HTTPS REST API. Long-term archival in CERN's Zenodo data center. Metadata remains permanently retrievable even if underlying files are versioned. |
| **Interoperable (I1–I3)** | Data formatted in standardized tabular TSV files with UTF-8 encoding. Formal ontology mapping (`NCBITaxon:1781`, `OBI:0001271`, `GO:0044010`, `GO:0030258`). Machine-actionable JSON-LD context via RO-Crate v1.1. |
| **Reusable (R1–R1.3)** | Complete column-level data dictionary (`data_dictionary.json`). Explicit permissive licensing: **CC-BY-4.0** for scientific data and text, **MIT** for code. Comprehensive provenance and computational build scripts. |

---

## 2. Directory Manifest

- `zenodo.json`: Zenodo deposition schema containing title, creators, ORCIDs, abstract, keywords, related OSDR DOIs, grant funding numbers (`NASA 80NSSC18K1467`), and license.
- `ro-crate-metadata.json`: Research Object Crate (RO-Crate v1.1) specification formalizing input data files, scripts, computational workflows, and resulting figures.
- `data_dictionary.json`: Comprehensive column-by-column semantic data dictionary describing all processed tabular files with units and ontology linkages.

---

## 3. Automated Zenodo Deposition Instructions

To publish or update this dataset package to Zenodo using the Zenodo REST API:

```bash
# Set your Zenodo access token (obtain from https://zenodo.org/account/settings/applications/tokens/)
export ZENODO_TOKEN="YOUR_ZENODO_API_TOKEN"

# Create a new deposition draft
DEPOSITION_ID=$(curl -s -H "Authorization: Bearer $ZENODO_TOKEN" \
  -H "Content-Type: application/json" \
  -X POST https://zenodo.org/api/deposit/depositions \
  -d @fair_deposit/zenodo.json | jq -r '.id')

echo "Created Zenodo Deposition ID: $DEPOSITION_ID"

# Upload the data archive
tar -czf osd528_fair_dataset.tar.gz data/ analysis/ manuscript/
curl -s -H "Authorization: Bearer $ZENODO_TOKEN" \
  -F "file=@osd528_fair_dataset.tar.gz" \
  https://zenodo.org/api/deposit/depositions/$DEPOSITION_ID/files

# Publish the deposition (mints official DOI)
# curl -s -H "Authorization: Bearer $ZENODO_TOKEN" -X POST https://zenodo.org/api/deposit/depositions/$DEPOSITION_ID/actions/publish
```
