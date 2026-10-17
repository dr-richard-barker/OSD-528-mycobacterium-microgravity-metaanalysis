#!/usr/bin/env python3
"""
08_package_results_and_supplement.py
Copy the pipeline outputs into results_tables/ under their supplementary-table names, and
write results_tables/README_RESULTS_TABLES.md and supplementary_figures/README.md.

Nothing here computes or types a result: every table is an output file of an earlier step
(or a join of output files), and the README row counts are counted from the files written.
Table S16 is produced directly by 01_fetch_osdr_microbial_data.py and only listed here.
Files in results_tables/ that are not in the current table list are removed, so a retired
table cannot linger.
"""

import csv
import glob
import json
import os
import shutil
import sys

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(PROJECT_DIR, "data", "processed")
RESULTS = os.path.join(PROJECT_DIR, "results_tables")
SUPP = os.path.join(PROJECT_DIR, "supplementary_figures")

# (table name, source under data/processed, description)
TABLES = [
    ("Table_S1_kallisto_estimated_counts.tsv", "osd528_counts_raw.tsv",
     "kallisto estimated read-pair counts per CDS, all nine deposited samples (RFPRPM4 = inferred second library)"),
    ("Table_S2_kallisto_tpm.tsv", "osd528_tpm.tsv", "Transcripts per million per CDS, all nine samples"),
    ("Table_S3_samples_and_sequencing_qc.tsv", None,
     "Sample groups, primary/sensitivity use, library index, read pairs and read assignment (rRNA, other ncRNA, CDS)"),
    ("Table_S4_deseq2_clinostat_vs_static1g.tsv", "deg_3dclinostat_vs_static1g.tsv",
     "DESeq2: 3D clinostat vs static 1g (primary samples); MLE and apeglm-shrunk log2FC, Wald p, BH FDR"),
    ("Table_S5_deseq2_rpm2_vs_static1g.tsv", "deg_rpm2_vs_static1g.tsv", "DESeq2: RPM 2.0 vs static 1g (primary samples)"),
    ("Table_S6_deseq2_clinostat_vs_rpm2.tsv", "deg_3dclinostat_vs_rpm2.tsv", "DESeq2: 3D clinostat vs RPM 2.0 (primary samples)"),
    ("Table_S7_shared_simulated_microgravity_genes.tsv", "deg_shared_simulated_microgravity.tsv",
     "Genes differentially expressed under both simulators in the same direction"),
    ("Table_S8_deseq2_sensitivity_with_rfprpm4.tsv", None,
     "DESeq2 results for all three contrasts with the inferred RFPRPM4 library included (sensitivity analysis)"),
    ("Table_S9_wgcna_module_assignments.tsv", "wgcna_module_assignments.tsv",
     "WGCNA module, module membership (kME) and connectivity (kTotal, kWithin) per gene"),
    ("Table_S10_wgcna_module_eigengenes.tsv", "wgcna_module_eigengenes.tsv", "Module eigengenes per primary sample"),
    ("Table_S11_wgcna_module_trait_correlations.tsv", "wgcna_module_trait_correlations.tsv",
     "Module eigengene vs trait Pearson correlation, Student p and BH FDR"),
    ("Table_S12_go_enrichment_deg_sets.tsv", "go_enrichment_deg.tsv",
     "GO over-representation (hypergeometric, BH FDR) for up and down DEG sets and the shared response"),
    ("Table_S13_go_enrichment_wgcna_modules.tsv", "go_enrichment_modules.tsv", "GO over-representation for each WGCNA module"),
    ("Table_S14_classification_loocv_predictions.tsv", "ml_loocv_predictions.tsv",
     "Leave-one-out held-out predictions and class probabilities (TabPFN, logistic regression, random forest)"),
    ("Table_S15_classification_metrics.tsv", "ml_metrics.tsv",
     "LOOCV accuracy, balanced accuracy, AUROC and label-permutation p-values per model, task and feature count"),
    ("Table_S17_kegg_mapped_reactions.tsv", "metabolic_model_reactions.tsv",
     "Curated reactions with their KEGG-supported M. marinum genes and per-contrast expression changes"),
    ("Table_S18_kegg_subsystem_summary.tsv", "metabolic_subsystem_perturbation.tsv",
     "Per-subsystem, per-contrast summaries of the KEGG-mapped reactions"),
    ("Table_S19_metabolic_model.json", "metabolic_model_sbml.json", "Reaction model with KEGG-supported gene links (JSON)"),
]
S16_GLOB = "Table_S16_*.tsv"
S16_DESC = "Microbial spaceflight and simulated-microgravity studies in NASA OSDR, from the OSDR API"


def tsv_rows(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def table_s3(dest):
    samples = tsv_rows(os.path.join(PROJECT_DIR, "config", "samples.tsv"))
    qc = json.load(open(os.path.join(PROC, "osd528_sequencing_qc_stats.json")))["samples"]
    out = []
    for s in samples:
        q = qc[s["sample_id"]]
        out.append({**s, "library_index": q["library_index"], "fastq_r1": q["fastq_r1"], "fastq_r2": q["fastq_r2"],
                    "read_pairs": q["read_pairs_processed"], "percent_pseudoaligned": q["percent_pseudoaligned"],
                    "percent_rRNA": q["percent_of_pairs_assigned"]["rRNA"],
                    "percent_other_ncRNA": q["percent_of_pairs_assigned"]["other_ncRNA"],
                    "percent_CDS": q["percent_of_pairs_assigned"]["CDS"]})
    write_tsv(dest, out, list(out[0].keys()))


def table_s8(dest):
    sens = os.path.join(PROC, "sensitivity_with_rfprpm4")
    out, fields = [], None
    for contrast, name in (("clinostat_vs_static1g", "deg_3dclinostat_vs_static1g.tsv"),
                           ("rpm2_vs_static1g", "deg_rpm2_vs_static1g.tsv"),
                           ("clinostat_vs_rpm2", "deg_3dclinostat_vs_rpm2.tsv")):
        for r in tsv_rows(os.path.join(sens, name)):
            out.append({"contrast": contrast, **r})
        fields = ["contrast"] + list(r.keys())
    write_tsv(dest, out, fields)


def count(path):
    if path.endswith(".json"):
        return "JSON"
    with open(path) as f:
        return f"{sum(1 for _ in f) - 1:,} rows"


COLUMN_DESCRIPTIONS = {
    "target_id": "kallisto target (RefSeq CDS record id)", "gene_id": "RefSeq locus tag (MMAR_RSxxxxx), GCF_000018345.1",
    "gene_symbol": "Gene symbol from the RefSeq annotation (locus tag if none)", "replicon": "Chromosome or plasmid accession",
    "pseudo": "yes if the CDS is annotated as a pseudogene", "protein": "RefSeq protein name",
    "baseMean": "DESeq2 mean of normalised counts", "log2FoldChange": "DESeq2 maximum-likelihood log2 fold change (numerator vs denominator)",
    "lfcSE": "Standard error of log2FoldChange", "log2FoldChange_shrunk": "apeglm-shrunk log2 fold change",
    "pvalue": "Wald test p-value", "padj": "Benjamini-Hochberg adjusted p-value", "significant": "YES if padj < 0.05",
    "contrast": "Comparison (numerator_vs_denominator)", "module": "WGCNA module colour (grey = unassigned)",
    "kME": "Module membership: correlation with the own module eigengene", "kTotal": "Whole-network connectivity",
    "kWithin": "Intramodular connectivity", "sample_id": "OSDR sample name", "trait": "Binary trait correlated with eigengenes",
    "pearson_r": "Pearson correlation coefficient", "group": "Condition (Ground, Clinostat, RPM) or gene-set group",
    "microgravity": "yes for clinostat and RPM samples", "primary": "yes if the sample enters the primary analysis",
    "note": "Free-text note", "gene_set": "Gene set tested", "subset": "Direction or module within the gene set",
    "ontology": "GO namespace: BP, MF or CC", "go_id": "GO term identifier", "term_name": "GO term name",
    "term_size": "Annotated universe genes in the term (after propagation)", "overlap": "Gene-set genes in the term",
    "set_size": "Annotated gene-set genes in the universe", "universe_size": "Annotated genes in the universe",
    "expected": "Expected overlap under the hypergeometric null", "fold_enrichment": "overlap / expected",
    "fdr": "Benjamini-Hochberg FDR within ontology and gene set", "n_terms_tested": "Terms tested in that family",
    "fdr_bh_hit_terms_only": "Sensitivity FDR over terms with at least one hit",
    "overlap_gene_ids": "Overlapping genes (locus tags)", "overlap_gene_symbols": "Overlapping genes (symbols)",
    "subsystem": "Metabolic or cellular subsystem", "bh_fdr": "Benjamini-Hochberg FDR", "hypergeom_p": "One-sided hypergeometric p-value",
    "n_deg": "DEGs (padj < 0.05)", "n_up": "Up-regulated DEGs", "n_down": "Down-regulated DEGs",
    "median_log2fc_shrunk": "Median apeglm-shrunk log2 fold change of tested genes",
    # samples and QC (Table S3)
    "library_index": "Index sequence in the first read header", "fastq_r1": "R1 FASTQ file used", "fastq_r2": "R2 FASTQ file used",
    "read_pairs": "Read pairs processed by kallisto", "percent_pseudoaligned": "Percent of pairs pseudoaligned",
    "percent_rRNA": "Percent of pairs assigned to rRNA decoys", "percent_other_ncRNA": "Percent of pairs assigned to tRNA/other ncRNA decoys",
    "percent_CDS": "Percent of pairs assigned to protein-coding genes",
    # classification (Tables S14, S15)
    "analysis": "primary (8 samples) or sensitivity_9 (adds the inferred RFPRPM4)", "task": "Classification task",
    "k": "Number of most variable genes selected within each training fold", "model": "Classifier",
    "n_samples": "Samples in the analysis", "class_counts": "Samples per class", "n_correct": "Held-out samples classified correctly",
    "accuracy": "LOOCV accuracy", "balanced_accuracy": "LOOCV balanced accuracy", "auroc": "AUROC of pooled held-out probabilities (binary task)",
    "permutation_method": "exact (all class-size-preserving relabellings) or Monte Carlo",
    "n_relabellings_total": "Distinct class-size-preserving relabellings", "n_relabellings_evaluated": "Relabellings evaluated",
    "p_accuracy": "Permutation p-value for accuracy", "p_balanced_accuracy": "Permutation p-value for balanced accuracy",
    "p_auroc": "Permutation p-value for AUROC", "null_mean_accuracy": "Mean accuracy over relabellings",
    "null_mean_balanced_accuracy": "Mean balanced accuracy over relabellings",
    "min_attainable_p": "Smallest p-value the permutation test can give", "true_label": "True class", "predicted_label": "Predicted class",
    "correct": "Whether the prediction was correct", "n_train": "Training samples in the fold",
    "prob_ground_1g": "Held-out probability of static 1g", "prob_simulated_microgravity": "Held-out probability of simulated microgravity",
    "prob_Clinostat": "Held-out probability of clinostat", "prob_Ground": "Held-out probability of ground", "prob_RPM": "Held-out probability of RPM",
    # OSDR compendium (Table S16)
    "accession": "OSDR study accession", "title": "Study title as given by OSDR", "organisms": "Microbial organisms listed by OSDR",
    "microbial_categories": "Categories from the NCBI Taxonomy lineage", "other_listed_organisms": "Other organisms listed",
    "host_organisms": "Host organisms listed", "assay_types": "Assay types as given by OSDR", "factors": "Study factors as given by OSDR",
    "project_type": "OSDR project type", "flight_program": "Flight programme", "mission": "Mission",
    "inclusion_basis": "Why the study meets the spaceflight/simulated-microgravity rule",
    "doi": "Study DOI provided by OSDR and resolved by DataCite (empty if none)",
    # KEGG-mapped reactions and subsystems (Tables S17, S18)
    "reaction_id": "Reaction identifier", "reaction_name": "Reaction name", "compartment": "Cellular compartment",
    "mapping_basis": "KEGG link type used (EC or KO)", "kegg_queries": "KEGG REST queries used", "ec_number": "EC number",
    "kegg_name": "KEGG entry name", "equation": "Reaction equation", "kegg_scope": "Scope of the KEGG link (gene family or specific)",
    "n_kegg_genes": "M. marinum genes linked by KEGG", "kegg_gene_ids": "KEGG gene ids (MMAR_xxxx)", "gene_ids": "RefSeq locus tags",
    "gene_symbols": "Gene symbols", "old_assignment": "Locus assigned in the earlier, uncorrected model",
    "old_locus_correct": "Whether that earlier locus is the KEGG-supported one", "n_reactions": "Reactions in the subsystem",
    "n_genes_mapped": "Genes mapped by KEGG", "n_genes_tested": "Mapped genes tested by DESeq2", "deg_fraction": "n_deg / n_genes_tested",
    "background_n_tested": "Genes tested in the contrast", "background_n_deg": "DEGs in the contrast",
    "lead_gene_id": "Gene with the smallest padj", "lead_gene_symbol": "Symbol of that gene",
    "lead_log2fc_shrunk": "Its apeglm-shrunk log2 fold change", "lead_padj": "Its adjusted p-value", "n_tested": "Genes tested",
}


def column_description(col, sample_ids):
    if col in COLUMN_DESCRIPTIONS:
        return COLUMN_DESCRIPTIONS[col]
    if col in sample_ids:
        return f"Value for sample {col}"
    if col.startswith("ME"):
        return f"Eigengene of WGCNA module {col[2:]}"
    for prefix, label in (("clinostat_", "clinostat vs ground"), ("rpm_", "RPM vs ground")):
        base = col[len(prefix):]
        if col.startswith(prefix) and base in COLUMN_DESCRIPTIONS:
            return f"{COLUMN_DESCRIPTIONS[base]} ({label})"
    for suffix, label in (("_clinostat", "clinostat vs ground"), ("_rpm", "RPM vs ground"),
                          ("_primary", "primary analysis"), ("_sensitivity", "sensitivity analysis")):
        if col.endswith(suffix) and col[:-len(suffix)] in COLUMN_DESCRIPTIONS:
            return f"{COLUMN_DESCRIPTIONS[col[:-len(suffix)]]} ({label})"
    return "Not described here; see the script that writes this file"


def write_fair(entries):
    """fair_deposit/data_dictionary.json and ro-crate-metadata.json from the files just written."""
    fair = os.path.join(PROJECT_DIR, "fair_deposit")
    sample_ids = {r["sample_id"] for r in tsv_rows(os.path.join(PROJECT_DIR, "config", "samples.tsv"))}
    files, undescribed = {}, set()
    for name, _, desc, source in entries:
        p = os.path.join(RESULTS, name)
        cols = {}
        if name.endswith(".tsv"):
            with open(p) as f:
                header = f.readline().rstrip("\n").split("\t")
            for c in header:
                d = column_description(c, sample_ids)
                cols[c] = d
                if d.startswith("Not described"):
                    undescribed.add(c)
        files[f"results_tables/{name}"] = {"description": desc, "source": source.strip("`"), "columns": cols}
    with open(os.path.join(fair, "data_dictionary.json"), "w") as f:
        json.dump({"title": "OSD-528 re-analysis: data dictionary",
                   "generated_by": "analysis/08_package_results_and_supplement.py", "files": files}, f, indent=2)
        f.write("\n")
    author = {"@id": "https://orcid.org/0000-0001-5681-9857", "@type": "Person", "name": "Richard Barker",
              "affiliation": "Collaborative Science Environment Inc., Madison, Wisconsin, USA"}
    graph = [
        {"@id": "ro-crate-metadata.json", "@type": "CreativeWork",
         "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"}, "about": {"@id": "./"}},
        {"@id": "./", "@type": "Dataset",
         "name": "A shared transcriptional response of Mycobacterium marinum to two microgravity simulators: "
                 "a re-analysis of NASA OSDR study OSD-528",
         "description": "Processed outputs, code and manuscript of an independent re-analysis of NASA OSDR OSD-528.",
         "author": {"@id": author["@id"]}, "license": {"@id": "https://creativecommons.org/licenses/by/4.0/"},
         "isBasedOn": {"@id": "https://doi.org/10.26030/r3re-fd65"},
         "hasPart": [{"@id": k} for k in files]},
        author,
        {"@id": "https://doi.org/10.26030/r3re-fd65", "@type": "Dataset", "name": "NASA OSDR OSD-528 (GLDS-528)"},
        {"@id": "https://creativecommons.org/licenses/by/4.0/", "@type": "CreativeWork", "name": "CC BY 4.0"},
    ]
    for k, v in files.items():
        graph.append({"@id": k, "@type": "File", "name": os.path.basename(k), "description": v["description"],
                      "encodingFormat": "application/json" if k.endswith(".json") else "text/tab-separated-values"})
    with open(os.path.join(fair, "ro-crate-metadata.json"), "w") as f:
        json.dump({"@context": "https://w3id.org/ro/crate/1.1/context", "@graph": graph}, f, indent=2)
        f.write("\n")
    print(f"Wrote fair_deposit/data_dictionary.json and ro-crate-metadata.json ({len(files)} files)")
    if undescribed:
        print(f"  columns without a description: {sorted(undescribed)}")


def main():
    os.makedirs(RESULTS, exist_ok=True)
    missing = [src for _, src, _ in TABLES if src and not os.path.exists(os.path.join(PROC, src))]
    s16 = sorted(glob.glob(os.path.join(RESULTS, S16_GLOB)))
    if missing or len(s16) != 1:
        sys.exit(f"Missing inputs: {missing}; Table S16 files found: {[os.path.basename(x) for x in s16]}")
    keep = {os.path.basename(s16[0]), "README_RESULTS_TABLES.md"}
    lines = ["# Supplementary tables", "",
             "Written by `analysis/08_package_results_and_supplement.py` from the pipeline outputs in",
             "`data/processed/`. Row counts are counted from the files.", "",
             "| Table | Rows | Description | Source |", "|---|---:|---|---|"]
    entries = []
    for name, src, desc in TABLES:
        dest = os.path.join(RESULTS, name)
        if name.startswith("Table_S3_"):
            table_s3(dest)
            source = "`config/samples.tsv` + `osd528_sequencing_qc_stats.json`"
        elif name.startswith("Table_S8_"):
            table_s8(dest)
            source = "`data/processed/sensitivity_with_rfprpm4/deg_*.tsv`"
        else:
            shutil.copyfile(os.path.join(PROC, src), dest)
            source = f"`data/processed/{src}`"
        keep.add(name)
        entries.append((name, count(dest), desc, source))
    entries.append((os.path.basename(s16[0]), count(s16[0]), S16_DESC, "`analysis/01_fetch_osdr_microbial_data.py`"))
    entries.sort(key=lambda e: int(e[0].split("_")[1][1:]))
    lines += [f"| `{n}` | {c} | {d} | {s} |" for n, c, d, s in entries]
    with open(os.path.join(RESULTS, "README_RESULTS_TABLES.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    # The same list for the supplementary information PDF (manuscript/npj_supplementary.tex).
    def tex(x):
        return x.replace("\\", "").replace("_", "\\_").replace("&", "\\&").replace("%", "\\%")
    tex_lines = ["% Generated by analysis/08_package_results_and_supplement.py. Do not edit by hand.",
                 "\\begin{description}"]
    for n, c, d, _ in entries:
        label = n.split("_")[1]
        tex_lines.append(f"  \\item[Table {label}] \\texttt{{{tex(n)}}} ({tex(c)}). {tex(d)}.")
    tex_lines.append("\\end{description}")
    with open(os.path.join(PROJECT_DIR, "manuscript", "chapters", "supplementary_tables.tex"), "w") as f:
        f.write("\n".join(tex_lines) + "\n")
    for old in os.listdir(RESULTS):
        if old not in keep:
            os.remove(os.path.join(RESULTS, old))
            print(f"  removed retired table {old}")
    print(f"Wrote {len(entries)} tables and README_RESULTS_TABLES.md")
    write_fair(entries)

    figs = sorted(f for f in os.listdir(SUPP) if f.endswith(".pdf")) if os.path.isdir(SUPP) else []
    with open(os.path.join(SUPP, "README.md"), "w") as f:
        f.write("# Supplementary figures\n\nDrawn by `analysis/06_generate_figures.py` (functions `figS*`) from the "
                "pipeline outputs. Each figure is provided as PDF (vector) and PNG.\n\n")
        for fig in figs:
            f.write(f"- `{fig[:-4]}`\n")


if __name__ == "__main__":
    main()
