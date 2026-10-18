#!/usr/bin/env python3
"""
09_manuscript_numbers.py
Write manuscript/chapters/numbers.tex: one LaTeX macro per number the manuscript quotes,
and the "Key results" block of README.md (between the numbers:start/end markers).

The manuscript text never types a result by hand. It uses these macros (for example
\\nDegClinoUp), and this script fills them from the pipeline's output files. Re-running
the pipeline and then this script updates every number in the text, the PDFs, the HTML
and the Word file at once. A missing input stops the script, so a stale or invented value
can't slip through.
"""

import csv
import json
import os
import sys

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(PROJECT_DIR, "data", "processed")
OUT = os.path.join(PROJECT_DIR, "manuscript", "chapters", "numbers.tex")

macros = {}
README_EXTRA = []


def load(name):
    p = os.path.join(PROC, name)
    if not os.path.exists(p):
        sys.exit(f"missing {os.path.relpath(p, PROJECT_DIR)}: run the pipeline first")
    with open(p) as f:
        return json.load(f)


def rows(name):
    p = os.path.join(PROC, name)
    if not os.path.exists(p):
        sys.exit(f"missing {os.path.relpath(p, PROJECT_DIR)}: run the pipeline first")
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def num(n):
    """Integer with a LaTeX thousands separator."""
    return f"{int(n):,}".replace(",", "{,}")


def sci(x, digits=2):
    """p-values and FDRs: plain decimal above 0.001, LaTeX scientific notation below."""
    x = float(x)
    if x >= 0.001:
        return f"{x:.{digits}g}"
    mant, exp = f"{x:.{digits - 1}e}".split("e")
    return f"{mant}\\times10^{{{int(exp)}}}"


def put(name, value):
    if not name.isalpha():
        raise ValueError(f"LaTeX macro names must be letters only: {name}")
    macros[name] = str(value)


# ------------------------------------------------------------------ quantification
qc = load("osd528_sequencing_qc_stats.json")
samples = {r["sample_id"]: r for r in csv.DictReader(open(os.path.join(PROJECT_DIR, "config", "samples.tsv")),
                                                     delimiter="\t")}
s = qc["samples"]
pairs = [v["read_pairs_processed"] for v in s.values()]
put("nSamplesDeposited", len(s))
put("nPrimary", sum(r["primary"] == "yes" for r in samples.values()))
for g, macro in (("Ground", "nGround"), ("Clinostat", "nClino"), ("RPM", "nRpm")):
    put(macro, sum(r["group"] == g and r["primary"] == "yes" for r in samples.values()))
put("pairsMinM", f"{min(pairs) / 1e6:.1f}")
put("pairsMaxM", f"{max(pairs) / 1e6:.1f}")
pa = [float(v["percent_pseudoaligned"]) for v in s.values()]
put("pctAlignedMin", f"{min(pa):.1f}")
put("pctAlignedMax", f"{max(pa):.1f}")
for key, macro in (("rRNA", "pctRrna"), ("CDS", "pctCds"), ("other_ncRNA", "pctNcrna")):
    v = [x["percent_of_pairs_assigned"][key] for x in s.values()]
    put(macro + "Min", f"{min(v):.1f}")
    put(macro + "Max", f"{max(v):.1f}")
put("kallistoVersion", qc["kallisto_version"].split()[-1])
put("strandedness", {"rf": "reverse", "fr": "forward", "unstranded": "un"}[qc["strandedness"]])
put("nCdsIndex", num(qc["index_targets"]["cds"]))
put("nNcrnaDecoys", num(qc["index_targets"]["ncrna_decoys"]))
d = qc["rfprpm4_deposit_defect"]
put("rpmFourPrefixBytesRone", num(d["R1"]["rpm41_prefix_bytes"]))
put("rpmFourPairs", num(s["RFPRPM4"]["read_pairs_processed"]))
put("rpmFourOnePairs", num(s["RFPRPM41"]["read_pairs_processed"]))

counts = rows("osd528_counts_raw.tsv")
totals = {smp: sum(float(r[smp]) for r in counts if r["pseudo"] == "no") for smp in samples}
prim_tot = [totals[k] for k, r in samples.items() if r["primary"] == "yes"]
put("cdsPairsPrimaryMinK", num(round(min(prim_tot) / 1e3)))
put("cdsPairsPrimaryMaxK", num(round(max(prim_tot) / 1e3)))

# ------------------------------------------------------------------ differential expression
de = load("de_summary.json")
p = de["primary"]
put("deseqVersion", de["method"]["deseq2"])
put("apeglmVersion", de["method"]["apeglm"])
put("nGenesInput", num(p["genes_input"]))
put("nGenesTested", num(p["genes_tested"]))
for key, tag in (("clinostat_vs_ground", "Clino"), ("rpm_vs_ground", "Rpm"), ("clinostat_vs_rpm", "ClinoRpm")):
    c = p["contrasts"][key]
    put(f"nDeg{tag}", num(c["deg"]))
    put(f"nDeg{tag}Up", num(c["up"]))
    put(f"nDeg{tag}Down", num(c["down"]))
    put(f"nDeg{tag}LfcOne", num(c["deg_abs_shrunk_lfc_ge_1"]))
ov = p["simulator_overlap"]
put("nDegShared", num(ov["deg_both"]))
put("nDegSharedSame", num(ov["deg_both_same_direction"]))
put("nDegEither", num(ov["deg_either"]))
put("jaccardShared", f"{ov['jaccard_same_direction']:.2f}")
put("rhoSimulators", f"{ov['spearman_lfc_all_tested']:.2f}")
shared = rows("deg_shared_simulated_microgravity.tsv")
put("nSharedUp", num(sum(float(r["log2FoldChange_clinostat"]) > 0 for r in shared)))
put("nSharedDown", num(sum(float(r["log2FoldChange_clinostat"]) < 0 for r in shared)))
ve = p["pca"]["variance_explained"]
put("pcaGenes", p["pca"]["genes"])
put("pcOnePct", f"{100 * ve[0]:.1f}")
put("pcTwoPct", f"{100 * ve[1]:.1f}")
pca = {r["sample_id"]: r for r in rows("pca_scores.tsv")}
put("pcOneThreeDElevenScore", f"{float(pca['RFP3D11']['PC1']):.1f}")
mg = [float(r["PC1"]) for k, r in pca.items() if samples[k]["microgravity"] == "yes" and k != "RFP3D11"]
gr = [float(r["PC1"]) for k, r in pca.items() if samples[k]["microgravity"] == "no"]
put("pcOneMgMin", f"{min(mg):.1f}")
put("pcOneMgMax", f"{max(mg):.1f}")
put("pcOneGroundMin", f"{min(gr):.1f}")
put("pcOneGroundMax", f"{max(gr):.1f}")
sens = de["sensitivity_with_rfprpm4"]
put("nDegClinoSens", num(sens["contrasts"]["clinostat_vs_ground"]["deg"]))
put("nDegRpmSens", num(sens["contrasts"]["rpm_vs_ground"]["deg"]))
put("nDegClinoRpmSens", num(sens["contrasts"]["clinostat_vs_rpm"]["deg"]))
pvs = de["primary_vs_sensitivity"]
put("nDegRpmBothSens", num(pvs["rpm_vs_ground"]["deg_both"]))
put("rLfcRpmSens", f"{pvs['rpm_vs_ground']['pearson_lfc']:.2f}")
put("nDegClinoBothSens", num(pvs["clinostat_vs_ground"]["deg_both"]))
put("rLfcClinoSens", f"{pvs['clinostat_vs_ground']['pearson_lfc']:.3f}")

# ------------------------------------------------------------------ WGCNA
w = load("wgcna_summary.json")
put("wgcnaVersion", w["wgcna_version"])
put("wgcnaPower", w["power"])
put("wgcnaRsq", f"{w['scale_free_R2_at_power']:.2f}")
put("wgcnaGenes", num(w["genes"]))
put("nModules", len(w["modules"]))
put("nGrey", w["unassigned_grey"])
put("nModuleTraitTests", w["module_trait_tests"])
put("nModuleTraitSig", w["module_trait_padj_lt_0.05"])
mt = [r for r in rows("wgcna_module_trait_correlations.tsv") if r["trait"] == "simulated_microgravity"]
for mod, tag in (("turquoise", "Turquoise"), ("blue", "Blue")):
    r = next(x for x in mt if x["module"] == "ME" + mod)
    put(f"size{tag}", num(w["modules"][mod]))
    put(f"r{tag}", f"{float(r['pearson_r']):.2f}")
    put(f"fdr{tag}", sci(r["padj"]))

# ------------------------------------------------------------------ GO
go = load("go_summary.json")
cov = go["coverage"]["tested_clinostat_vs_static1g"]
put("goTestedGenes", num(cov["genes"]))
put("goAnnotatedGenes", num(cov["annotated_any"]))
put("goAnnotatedPct", f"{100 * cov['fraction_annotated_any']:.0f}")
src = go["annotation_source"]
put("goIeaPct", f"{100 * src['fraction_electronic_iea']:.1f}")
put("uniprotRelease", src["uniprot_release"].replace("_", "\\_"))
put("goRelease", go["ontology_source"]["data_version"].split("/")[-1])
gd = rows("go_enrichment_deg.tsv")


def go_fdr(gset, term):
    hits = [r for r in gd if r["gene_set"] == gset and r["term_name"] == term]
    if not hits:
        sys.exit(f"GO term {term!r} not tested in {gset}; update the manuscript text and this script")
    r = min(hits, key=lambda x: float(x["fdr"]))
    return r


for gset, tag in (("shared_simulated_microgravity_up", "SharedUp"), ("shared_simulated_microgravity_down", "SharedDown"),
                  ("clinostat_vs_static1g_up", "ClinoUp"), ("clinostat_vs_static1g_down", "ClinoDown"),
                  ("rpm2_vs_static1g_up", "RpmUp"), ("rpm2_vs_static1g_down", "RpmDown")):
    put(f"nGoSig{tag}", sum(1 for r in gd if r["gene_set"] == gset and float(r["fdr"]) < 0.05))
for gset, term, tag in (
        ("shared_simulated_microgravity_up", "extracellular region", "ExtracellSharedUp"),
        ("shared_simulated_microgravity_up", "NADH dehydrogenase activity", "NadhSharedUp"),
        ("rpm2_vs_static1g_up", "NADH dehydrogenase activity", "NadhRpmUp"),
        ("clinostat_vs_static1g_up", "NADH dehydrogenase activity", "NadhClinoUp"),
        ("shared_simulated_microgravity_down", "DNA-binding transcription factor activity", "TfSharedDown"),
        ("shared_simulated_microgravity_down", "response to heat", "HeatSharedDown")):
    r = go_fdr(gset, term)
    put(f"goFdr{tag}", sci(r["fdr"]))
    put(f"goOverlap{tag}", f"{r['overlap']}/{r['term_size']}")



def symbol_list(gset, term):
    """Italic gene symbols behind a GO result; locus tags without a symbol are counted, not listed."""
    r = go_fdr(gset, term)
    syms = [x for x in r["overlap_gene_symbols"].split(";") if x]
    named = {}
    for x in syms:
        if not x.startswith("MMAR_"):
            named[x] = named.get(x, 0) + 1
    parts = [f"\\textit{{{k}}}" + (f" ({v} copies)" if v > 1 else "") for k, v in named.items()]
    unnamed = sum(x.startswith("MMAR_") for x in syms)
    return ", ".join(parts) + (f" and {unnamed} unnamed genes" if unnamed > 1 else " and 1 unnamed gene" if unnamed else "")


put("goGenesHeatSharedDown", symbol_list("shared_simulated_microgravity_down", "response to heat"))
put("goGenesNadhRpmUp", symbol_list("rpm2_vs_static1g_up", "NADH dehydrogenase activity"))
protein = {r["gene_id"]: r["protein"] for r in counts}
ext = go_fdr("shared_simulated_microgravity_up", "extracellular region")
ext_ids = [x for x in ext["overlap_gene_ids"].split(";") if x]
put("goMceSharedUp", sum("mce" in protein.get(i, "").lower() for i in ext_ids))

# ------------------------------------------------------------------ classification (TabPFN and baselines)
ml = rows("ml_metrics.tsv")
mls = load("ml_summary.json")
put("tabpfnVersion", mls["environment"]["tabpfn"] if isinstance(mls.get("environment"), dict) and "tabpfn" in mls["environment"]
    else mls["tabpfn_checkpoint"]["model_version"])
put("tabpfnCheckpoint", mls["tabpfn_checkpoint"]["checkpoint_file"])
TASK = {"binary_microgravity_vs_ground": "Bin", "three_class_modality": "Three"}
MODEL = {"TabPFN_v2": "Tabpfn", "LogReg_L2": "Logreg", "RandomForest_500": "Forest", "Majority_baseline": "Majority"}
KWORD = {"50": "Fifty", "200": "TwoHundred"}
ANALYSIS = {"primary": "", "sensitivity_9": "Sens"}
for r in ml:
    if r["analysis"] not in ANALYSIS or r["k"] not in KWORD:
        continue
    tag = f"ml{ANALYSIS[r['analysis']]}{TASK[r['task']]}{MODEL[r['model']]}{KWORD[r['k']]}"
    put(tag + "Correct", f"{r['n_correct']}/{r['n_samples']}")
    put(tag + "Bal", f"{float(r['balanced_accuracy']):.2f}")
    if r["auroc"] not in ("", "nan", "NaN"):
        put(tag + "Auroc", f"{float(r['auroc']):.2f}")
    if r["model"] != "Majority_baseline":
        if r["p_balanced_accuracy"] in ("", "nan", "NaN"):
            sys.exit(f"{tag}: permutation p-value missing; run the full classification step first")
        put(tag + "P", f"{float(r['p_balanced_accuracy']):.2g}")
        put(tag + "Pmin", f"{float(r['min_attainable_p']):.2g}")
        put(tag + "Nperm", num(r["n_relabellings_total"]))

pred = rows("ml_loocv_predictions.tsv")
wrong = {}
for r in pred:
    if r["analysis"] == "primary" and r["task"] == "binary_microgravity_vs_ground" and r["model"] != "Majority_baseline":
        wrong.setdefault((r["model"], r["k"]), set())
        if r["true_label"] != r["predicted_label"]:
            wrong[(r["model"], r["k"])].add(r["sample_id"])
sets = {frozenset(v) for v in wrong.values()}
if len(sets) != 1:
    sys.exit("models misclassify different samples in the binary task; update the manuscript text")
put("mlBinMisclassified", " and ".join(sorted(next(iter(sets)))))
bin_p = [float(r["p_balanced_accuracy"]) for r in ml if r["analysis"] == "primary"
         and r["task"] == "binary_microgravity_vs_ground" and r["model"] != "Majority_baseline"]
put("mlBinPMin", f"{min(bin_p):.2g}")
put("mlBinPMax", f"{max(bin_p):.2g}")
prim_bin = [r for r in ml if r["analysis"] == "primary" and r["task"] == "binary_microgravity_vs_ground"
            and r["model"] != "Majority_baseline"]
aurocs = [float(r["auroc"]) for r in prim_bin]
put("mlBinAurocMin", f"{min(aurocs):.2f}")
put("mlBinAurocMax", f"{max(aurocs):.2f}")
if not all(abs(float(r["p_auroc"]) - float(r["min_attainable_p"])) < 1e-9 for r in prim_bin):
    sys.exit("not every primary binary AUROC p equals the smallest attainable p; update the manuscript text")
put("mlBinPAuroc", f"{float(prim_bin[0]['p_auroc']):.2g}")
README_EXTRA.append(
    f"- **Classification** (leave-one-out, exact permutation tests): microgravity vs 1g, every model "
    f"{macros['mlBinTabpfnFiftyCorrect']} correct (balanced-accuracy p = {macros['mlBinPMin']}–{macros['mlBinPMax']}); "
    f"three classes, TabPFN {macros['mlThreeTabpfnFiftyCorrect']} (p = {macros['mlThreeTabpfnFiftyP']}) vs logistic "
    f"regression {macros['mlThreeLogregFiftyCorrect']}, but not robust to adding RFPRPM4 "
    f"({macros['mlSensThreeTabpfnFiftyCorrect']}, p = {macros['mlSensThreeTabpfnFiftyP']}).")

# ------------------------------------------------------------------ KEGG-mapped metabolic subsystems
audit = rows("metabolic_model_locus_audit.tsv")
rxn = rows("metabolic_model_reactions.tsv")
sub = rows("metabolic_subsystem_perturbation.tsv")
kegg = load("kegg_metabolic_snapshot.json")
put("keggOldReactions", len(audit))
put("keggOldCorrect", sum(r["old_locus_correct"] == "YES" for r in audit))
put("keggOldWrong", sum(r["old_locus_correct"] == "NO" for r in audit))
put("keggOldUnverifiable", sum(r["old_locus_correct"] == "UNVERIFIABLE" for r in audit))
put("keggDropped", sum(r["status"] == "dropped" for r in audit))
put("keggReactions", len(rxn))
put("keggGenes", len({g for r in rxn for g in r["gene_ids"].split(";") if g}))
put("keggSubsystems", len({r["subsystem"] for r in sub}))
put("keggGenesDate", kegg["source"]["kegg_database_dates"]["genes"].replace("/", "-"))


def subsystem(name, contrast):
    hit = [r for r in sub if r["subsystem"] == name and r["contrast"] == contrast]
    if len(hit) != 1:
        sys.exit(f"subsystem {name!r} / {contrast} not found once in metabolic_subsystem_perturbation.tsv")
    return hit[0]


for contrast, tag in (("rpm2_vs_static1g", "Rpm"), ("clinostat_vs_static1g", "Clino")):
    r = subsystem("Cellular Respiration & Energy", contrast)
    put(f"respDeg{tag}", f"{r['n_deg']}/{r['n_genes_tested']}")
    put(f"respUp{tag}", r["n_up"])
    put(f"respFdr{tag}", sci(r["bh_fdr"]))
    put(f"respMedian{tag}", f"{float(r['median_log2fc_shrunk']):+.2f}")
    m = subsystem("Mycolic Acid & FAS-II Envelope", contrast)
    put(f"mycolicDeg{tag}", f"{m['n_deg']}/{m['n_genes_tested']}")
    x = subsystem("Redox Homeostasis & Cofactors", contrast)
    put(f"redoxDeg{tag}", f"{x['n_deg']}/{x['n_genes_tested']}")
put("nSubsystemsSig", sum(float(r["bh_fdr"]) < 0.05 for r in sub))
argf = [r for r in rows("deg_3dclinostat_vs_static1g.tsv") if r["gene_symbol"] == "argF"]
if len(argf) != 1:
    sys.exit("argF not found once in the clinostat DESeq2 table")
put("argfLfcClino", f"{float(argf[0]['log2FoldChange_shrunk']):+.2f}")
put("argfFdrClino", sci(argf[0]["padj"]))
argf_r = [r for r in rows("deg_rpm2_vs_static1g.tsv") if r["gene_symbol"] == "argF"][0]
put("argfFdrRpm", sci(argf_r["padj"]))
README_EXTRA.append(
    f"- **KEGG-mapped metabolic subsystems:** {macros['keggReactions']} reactions with KEGG-supported genes "
    f"(of the old model's {macros['keggOldReactions']} reactions, {macros['keggOldWrong']} were on the wrong gene). "
    f"Only cellular respiration is enriched for DEGs, under the RPM ({macros['respDegRpm']} genes, "
    f"{'all' if macros['respDegRpm'].split('/')[0] == macros['respUpRpm'] else macros['respUpRpm']} up; "
    f"FDR {macros['respFdrRpm'].replace(chr(92) + 'times10^{', ' × 10^').replace('}', '')}); mycolic acid / FAS-II "
    f"shows {macros['mycolicDegRpm']} DEGs.")

# ------------------------------------------------------------------ OSDR microbial compendium (Table S16)
import glob  # noqa: E402

s16 = sorted(glob.glob(os.path.join(PROJECT_DIR, "results_tables", "Table_S16_*.tsv")))
if len(s16) != 1:
    sys.exit(f"expected one results_tables/Table_S16_*.tsv, found {len(s16)}")
with open(s16[0]) as f:
    comp = list(csv.DictReader(f, delimiter="\t"))
# Same rule as Figure 7: a study counts as spaceflight if its inclusion basis mentions spaceflight.
flight = [r for r in comp if "spaceflight" in r["inclusion_basis"].lower()]
put("nOsdrStudies", len(comp))
put("nOsdrSpaceflight", len(flight))
put("nOsdrGroundSim", len(comp) - len(flight))
put("nOsdrWithDoi", sum(bool(r["doi"].strip()) for r in comp))
cats = [c.strip() for r in comp for c in r["microbial_categories"].split(";") if c.strip()]
put("nOsdrBacteria", cats.count("bacteria"))
put("nOsdrCommunity", cats.count("microbial community"))
put("nOsdrFungi", cats.count("fungus"))
text = {r["accession"]: (r["inclusion_basis"] + " | " + r["factors"]).lower() for r in comp}
rpm_clino = sorted(a for a, t in text.items() if "random positioning" in t and "clinostat" in t)
put("osdrRpmClinostatStudies", ", ".join(rpm_clino))
put("nOsdrRpmClinostat", len(rpm_clino))
README_EXTRA.append(
    f"- **Context:** {macros['nOsdrStudies']} microbial spaceflight or simulated-microgravity studies in OSDR "
    f"({macros['nOsdrSpaceflight']} spaceflight, {macros['nOsdrGroundSim']} ground simulation; Table S16). "
    f"Studies pairing an RPM with a 3D clinostat: {macros['osdrRpmClinostatStudies']}.")

# ------------------------------------------------------------------ write
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    f.write("% Generated by analysis/09_manuscript_numbers.py from data/processed/. Do not edit by hand.\n")
    for k in sorted(macros):
        f.write(f"\\newcommand{{\\{k}}}{{{macros[k]}}}\n")
print(f"Wrote {len(macros)} macros to {os.path.relpath(OUT, PROJECT_DIR)}")

# ------------------------------------------------------------------ README block


def md(name):
    """A macro's value as Markdown: drop LaTeX thousands braces and italics."""
    import re
    v = macros[name].replace("{,}", ",")
    v = re.sub(r"\\textit\{([^}]*)\}", r"*\1*", v)
    return re.sub(r"\\times10\^\{(-?\d+)\}", r" × 10^\1", v)


readme_lines = [
    "Generated by `analysis/09_manuscript_numbers.py` from `data/processed/`; do not edit by hand.",
    "",
    f"- **Reads:** {md('pairsMinM')}–{md('pairsMaxM')} million read pairs per library, of which "
    f"{md('pctRrnaMin')}–{md('pctRrnaMax')}% are rRNA and {md('pctCdsMin')}–{md('pctCdsMax')}% protein-coding.",
    f"- **Samples:** primary analysis {md('nGround')} ground, {md('nClino')} clinostat and {md('nRpm')} RPM; "
    "the inferred RFPRPM4 library is used only in the sensitivity analysis.",
    f"- **Differential expression** (DESeq2, FDR < 0.05, {md('nGenesTested')} genes tested): "
    f"{md('nDegClino')} genes under the clinostat ({md('nDegClinoUp')} up, {md('nDegClinoDown')} down) and "
    f"{md('nDegRpm')} under the RPM ({md('nDegRpmUp')} up, {md('nDegRpmDown')} down).",
    f"- **Shared response:** {md('nDegShared')} DEGs under both simulators, {md('nDegSharedSame')} in the same "
    f"direction; fold changes correlate at Spearman ρ = {md('rhoSimulators')}. Clinostat vs RPM directly: "
    f"{md('nDegClinoRpm')} DEGs (low power, RPM n = {md('nRpm')}).",
    f"- **Sensitivity analysis** with the inferred RFPRPM4: {md('nDegRpmSens')} RPM DEGs; fold changes correlate "
    f"with the primary analysis at r = {md('rLfcRpmSens')}.",
    f"- **GO** (electronic annotations, {md('goAnnotatedPct')}% of tested genes annotated): shared up genes are "
    f"enriched for the extracellular region ({md('goMceSharedUp')} Mce-family genes) and NADH dehydrogenase "
    f"activity; shared down genes for transcription factor activity and the response to heat.",
    f"- **WGCNA:** {md('nModules')} modules; {md('nModuleTraitSig')} of {md('nModuleTraitTests')} module–trait "
    f"tests pass FDR < 0.05 (turquoise r = {md('rTurquoise')}, blue r = {md('rBlue')}, with simulated microgravity).",
]
for extra in README_EXTRA:
    readme_lines.append(extra)
readme_path = os.path.join(PROJECT_DIR, "README.md")
text = open(readme_path).read()
start, end = "<!-- numbers:start -->", "<!-- numbers:end -->"
if start in text and end in text:
    head, rest = text.split(start, 1)
    _, tail = rest.split(end, 1)
    open(readme_path, "w").write(head + start + "\n" + "\n".join(readme_lines) + "\n" + end + tail)
    print("Updated the Key results block of README.md")
else:
    sys.exit("README.md lacks the numbers:start/end markers")
