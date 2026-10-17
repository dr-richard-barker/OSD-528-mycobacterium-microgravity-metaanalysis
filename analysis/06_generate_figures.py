#!/usr/bin/env python3
"""
06_generate_figures.py
Draw the nine manuscript figures from the pipeline's output files.

Every number, point and label in a figure is read from data/processed/ (or config/), so a
figure can only change when its inputs change. Nothing is typed in by hand. A figure whose
inputs are missing is skipped with a message rather than drawn from placeholders.

Inputs (written by the numbered analysis scripts)
  config/samples.tsv                                  sample groups, primary/sensitivity use
  data/processed/osd528_sequencing_qc_stats.json      00: read pairs, rRNA share, deposit defect
  data/processed/de_summary.json, pca_scores.tsv,     02: DESeq2 results and PCA
    deg_*.tsv, deg_shared_simulated_microgravity.tsv
  data/processed/wgcna_*.tsv, wgcna_summary.json      03: co-expression network
  data/processed/ml_*.tsv, ml_summary.json            04: LOOCV classification
  data/processed/go_*.tsv, go_summary.json            05: GO enrichment
  data/processed/metabolic_*.tsv                      07: KEGG-mapped metabolic subsystems
  results_tables/Table_S16_*.tsv                      01: OSDR microbial study compendium

Outputs: manuscript/figures/figN_*.pdf (vector) and .png (200 dpi).
"""

import glob
import json
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402
from matplotlib.transforms import Bbox  # noqa: E402

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(PROJECT_DIR, "data", "processed")
FIG_DIR = os.path.join(PROJECT_DIR, "manuscript", "figures")
CONFIG = os.path.join(PROJECT_DIR, "config", "samples.tsv")

GROUP_COLOURS = {"Ground": "#2166ac", "Clinostat": "#ef8a62", "RPM": "#b2182b"}
GROUP_LABELS = {"Ground": "Static 1g", "Clinostat": "3D clinostat", "RPM": "RPM 2.0"}
UP, DOWN, NS = "#b2182b", "#2166ac", "#bdbdbd"
CONTRASTS = {
    "clinostat_vs_ground": ("deg_3dclinostat_vs_static1g.tsv", "3D clinostat vs static 1g"),
    "rpm_vs_ground": ("deg_rpm2_vs_static1g.tsv", "RPM 2.0 vs static 1g"),
}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9,
    "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "svg.hashsalt": "osd528",
})


def path(*parts):
    return os.path.join(PROC, *parts)


def have(*files):
    missing = [f for f in files if not glob.glob(f)]
    if missing:
        print(f"  skipped: missing {', '.join(os.path.relpath(m, PROJECT_DIR) for m in missing)}")
        return False
    return True


def panel_label(ax, letter):
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight="bold",
            va="bottom", ha="left")


def save(fig, name):
    os.makedirs(FIG_DIR, exist_ok=True)
    for ext, kw in (("pdf", {"metadata": {"CreationDate": None, "ModDate": None}}),
                    ("png", {"dpi": 200, "metadata": {"Software": None}})):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"), bbox_inches="tight", **kw)
    plt.close(fig)
    print(f"  wrote manuscript/figures/{name}.pdf/.png")


def module_colour(name):
    """Matplotlib colour for a WGCNA module name (R colour names, e.g. grey60)."""
    import re
    from matplotlib import colors as mcolors
    m = re.fullmatch(r"gr[ae]y(\d+)", name)
    if m:
        level = int(m.group(1)) / 100
        return (level, level, level)
    if name == "lightyellow":
        return "#e6e600"  # R's lightyellow is near-white; darken so the bar stays visible
    return name if mcolors.is_color_like(name) else "#9e9e9e"


def label_samples(ax, pca, bold=()):
    """Label PCA points without labels overlapping each other or other points.

    For each point, try four placements (right-above, left-above, right-below, left-below) and
    keep the first whose drawn box clears every earlier label and every other point.
    """
    ax.margins(x=0.15, y=0.12)
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    pts = ax.transData.transform(pca[["PC1", "PC2"]].values)
    taken = []
    for i, (_, r) in enumerate(pca.iterrows()):
        for dx, dy, ha, va in ((7, 7, "left", "bottom"), (-7, 7, "right", "bottom"),
                               (7, -7, "left", "top"), (-7, -7, "right", "top"),
                               (12, 14, "left", "bottom"), (-12, 14, "right", "bottom"),
                               (12, -14, "left", "top"), (-12, -14, "right", "top")):
            ann = ax.annotate(r.sample_id, (r.PC1, r.PC2), xytext=(dx, dy), textcoords="offset points",
                              ha=ha, va=va, fontsize=5.5,
                              fontweight="bold" if r.sample_id in bold else "normal",
                              arrowprops=dict(arrowstyle="-", lw=0.4, color="#636363", shrinkA=0, shrinkB=3))
            box = ann.get_window_extent(renderer).expanded(1.05, 1.15)  # the text box, not the leader line
            hits_point = any(box.overlaps(Bbox.from_bounds(x - 5, y - 5, 10, 10))
                             for j, (x, y) in enumerate(pts) if j != i)
            if not hits_point and not any(box.overlaps(t) for t in taken):
                taken.append(box)
                break
            ann.remove()
        else:  # no clear spot: fall back to the first placement
            ann = ax.annotate(r.sample_id, (r.PC1, r.PC2), xytext=(5, 3), textcoords="offset points",
                              fontsize=5.5, fontweight="bold" if r.sample_id in bold else "normal")
            taken.append(ann.get_window_extent(renderer))


def load_samples():
    return pd.read_csv(CONFIG, sep="\t")


def load_json(name):
    with open(path(name)) as f:
        return json.load(f)


# --------------------------------------------------------------------------- figure 1
def fig1_study_design():
    print("Figure 1: study design and workflow")
    if not have(CONFIG, path("osd528_sequencing_qc_stats.json"), path("de_summary.json")):
        return
    samples = load_samples()
    qc = load_json("osd528_sequencing_qc_stats.json")
    de = load_json("de_summary.json")
    pct = pd.DataFrame({s: v["percent_of_pairs_assigned"] for s, v in qc["samples"].items()}).T
    pairs = [v["read_pairs_processed"] for v in qc["samples"].values()]

    fig = plt.figure(figsize=(7.2, 4.6))
    ax = fig.add_axes([0.0, 0.52, 1.0, 0.46])
    ax.axis("off")
    panel_label(ax, "a")
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 1)
    for i, group in enumerate(["Ground", "Clinostat", "RPM"]):
        rows = samples[samples.group == group]
        prim = rows[rows.primary == "yes"].sample_id.tolist()
        sens = rows[rows.primary != "yes"].sample_id.tolist()
        ax.add_patch(FancyBboxPatch((i + 0.06, 0.08), 0.88, 0.8, boxstyle="round,pad=0.01",
                                    fc="white", ec=GROUP_COLOURS[group], lw=1.5))
        ax.text(i + 0.5, 0.80, GROUP_LABELS[group], ha="center", va="center",
                fontweight="bold", color=GROUP_COLOURS[group], fontsize=9)
        ax.text(i + 0.5, 0.67, f"primary n = {len(prim)}", ha="center", fontsize=8)
        ax.text(i + 0.5, 0.42, "\n".join(prim), ha="center", va="center", fontsize=7.5)
        if sens:
            ax.text(i + 0.5, 0.15, "+ " + ", ".join(sens) + "\n(inferred library;\nsensitivity analysis only)",
                    ha="center", va="center", fontsize=6.5, color="#555555", style="italic")
    ax.set_title("OSD-528: Mycobacterium marinum RNA-seq under two microgravity simulators "
                 "and static 1g", fontsize=9, loc="left", x=0.02)

    axb = fig.add_axes([0.07, 0.06, 0.40, 0.36])
    panel_label(axb, "b")
    order = list(qc["samples"].keys())
    bottom = np.zeros(len(order))
    for col, colour, lab in (("rRNA", "#9e9e9e", "rRNA"), ("other_ncRNA", "#d9d9d9", "tRNA / other ncRNA"),
                             ("CDS", "#2166ac", "protein-coding (CDS)")):
        vals = pct.loc[order, col].values
        axb.bar(range(len(order)), vals, bottom=bottom, color=colour, label=lab, width=0.75)
        bottom += vals
    axb.set_xticks(range(len(order)))
    axb.set_xticklabels(order, rotation=60, ha="right", fontsize=6.5)
    axb.set_ylabel("% of read pairs")
    axb.set_ylim(0, 100)
    axb.legend(loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=3, frameon=False, fontsize=6)
    axb.set_title(f"Read assignment ({min(pairs)/1e6:.1f}-{max(pairs)/1e6:.1f} M pairs per sample)",
                  fontsize=8, loc="left")

    axc = fig.add_axes([0.53, 0.0, 0.46, 0.44])
    axc.axis("off")
    axc.set_xlim(0, 1)
    axc.set_ylim(0, 1)
    panel_label(axc, "c")
    steps = [
        "18 raw FASTQs from OSDR, md5-verified\nRFPRPM4 files = RFPRPM41 + a second library",
        f"kallisto {qc['kallisto_version'].split()[-1]}, paired-end, {qc['strandedness']}-stranded\n"
        f"index: {qc['index_targets']['cds']:,} CDS + {qc['index_targets']['ncrna_decoys']} non-coding RNA decoys",
        f"DESeq2 {de['method']['deseq2']} + apeglm shrinkage\n"
        f"{de['primary']['genes_tested']:,} genes tested, FDR < {de['method']['alpha']}",
        "WGCNA  |  GO enrichment  |  KEGG subsystems\nTabPFN and baselines (leave-one-out CV)",
    ]
    for k, text in enumerate(steps):
        y = 0.86 - k * 0.25
        axc.add_patch(FancyBboxPatch((0.03, y - 0.085), 0.94, 0.17, boxstyle="round,pad=0.008",
                                     fc="#f7f7f7", ec="#636363", lw=0.8))
        axc.text(0.5, y, text, ha="center", va="center", fontsize=6.3, linespacing=1.3)
        if k < len(steps) - 1:
            axc.annotate("", xy=(0.5, y - 0.165), xytext=(0.5, y - 0.09),
                         arrowprops=dict(arrowstyle="-|>", color="#636363", lw=0.8))
    axc.set_title("Analysis workflow", fontsize=8, loc="left")
    save(fig, "fig1_study_design")


# --------------------------------------------------------------------------- figure 2
def volcano(ax, table, title, alpha):
    d = pd.read_csv(path(table), sep="\t").dropna(subset=["padj"])
    y = -np.log10(d.padj.clip(lower=1e-300))
    x = d.log2FoldChange_shrunk
    sig = d.padj < alpha
    up, down = sig & (d.log2FoldChange > 0), sig & (d.log2FoldChange < 0)
    ax.scatter(x[~sig], y[~sig], s=3, c=NS, lw=0, rasterized=True)
    ax.scatter(x[up], y[up], s=4, c=UP, lw=0, rasterized=True, label=f"up ({up.sum()})")
    ax.scatter(x[down], y[down], s=4, c=DOWN, lw=0, rasterized=True, label=f"down ({down.sum()})")
    ax.axhline(-np.log10(alpha), color="#636363", lw=0.6, ls="--")
    ax.set_xlabel("log$_2$ fold change (apeglm-shrunk)")
    ax.set_ylabel("$-$log$_{10}$ FDR")
    ax.set_title(title, loc="left")
    ax.legend(loc="upper left", markerscale=2, framealpha=0.85, edgecolor="none")


def fig2_volcano_pca():
    print("Figure 2: PCA and volcano plots")
    files = [path("pca_scores.tsv"), path("de_summary.json")] + [path(t) for t, _ in CONTRASTS.values()]
    if not have(*files):
        return
    de = load_json("de_summary.json")
    ve = de["primary"]["pca"]["variance_explained"]
    pca = pd.read_csv(path("pca_scores.tsv"), sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.0), constrained_layout=True,
                             gridspec_kw={"width_ratios": [1.25, 1, 1]})
    ax = axes[0]
    for group, sub in pca.groupby("group"):
        ax.scatter(sub.PC1, sub.PC2, s=28, c=GROUP_COLOURS[group], label=GROUP_LABELS[group],
                   edgecolor="black", lw=0.4, zorder=3)
    label_samples(ax, pca)
    ax.set_xlabel(f"PC1 ({100 * ve[0]:.1f}% of variance)")
    ax.set_ylabel(f"PC2 ({100 * ve[1]:.1f}% of variance)")
    ax.set_title(f"PCA, {de['primary']['pca']['genes']} most variable genes (VST)", loc="left")
    ax.legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3,
              handletextpad=0.2, columnspacing=0.8)
    panel_label(ax, "a")
    for ax, (letter, (table, title)) in zip(axes[1:], zip("bc", CONTRASTS.values())):
        volcano(ax, table, title, de["method"]["alpha"])
        panel_label(ax, letter)
    save(fig, "fig2_volcano_pca")


# --------------------------------------------------------------------------- figure 3
def fig3_wgcna_modules():
    print("Figure 3: WGCNA network construction, modules and module-trait associations")
    files = [path("wgcna_soft_threshold.tsv"), path("wgcna_summary.json"),
             path("wgcna_module_trait_correlations.tsv")]
    if not have(*files):
        return
    fit = pd.read_csv(path("wgcna_soft_threshold.tsv"), sep="\t")
    summ = load_json("wgcna_summary.json")
    mt = pd.read_csv(path("wgcna_module_trait_correlations.tsv"), sep="\t")
    fig = plt.figure(figsize=(7.2, 4.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.35])
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(fit.Power, fit.scale_free_R2, "o-", color="#636363", ms=3, lw=0.8)
    ax.axhline(0.8, color=UP, lw=0.6, ls="--")
    ax.axvline(summ["power"], color=DOWN, lw=0.6, ls=":")
    ax.set_xlabel("soft-thresholding power")
    ax.set_ylabel("scale-free fit (signed $R^2$)")
    ax.set_title(f"Power {summ['power']} chosen ($R^2$ = {summ['scale_free_R2_at_power']:.2f})", loc="left")
    panel_label(ax, "a")

    ax = fig.add_subplot(gs[1, 0])
    sizes = pd.Series(summ["modules"]).sort_values(ascending=False)
    ax.bar(range(len(sizes)), sizes.values, color=[module_colour(c) for c in sizes.index],
           edgecolor="black", lw=0.3)
    ax.set_xticks(range(len(sizes)))
    ax.set_xticklabels(sizes.index, rotation=70, ha="right", fontsize=6)
    ax.set_ylabel("genes")
    ax.set_title(f"{len(sizes)} modules from {summ['genes']:,} genes\n"
                 f"({summ['unassigned_grey']} unassigned)", loc="left")
    panel_label(ax, "b")

    ax = fig.add_subplot(gs[:, 1])
    traits = ["simulated_microgravity", "clinostat", "rpm"]
    r = mt.pivot(index="module", columns="trait", values="pearson_r")[traits]
    q = mt.pivot(index="module", columns="trait", values="padj")[traits]
    order = r["simulated_microgravity"].sort_values().index
    r, q = r.loc[order], q.loc[order]
    im = ax.imshow(r.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    for i in range(r.shape[0]):
        for j in range(r.shape[1]):
            star = "*" if q.values[i, j] < 0.05 else ""
            ax.text(j, i, f"{r.values[i, j]:.2f}{star}", ha="center", va="center", fontsize=5.5,
                    color="white" if abs(r.values[i, j]) > 0.7 else "black")
    ax.set_xticks(range(len(traits)))
    ax.set_xticklabels(["simulated\nmicrogravity", "clinostat", "RPM"], fontsize=6.5)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([m.replace("ME", "") for m in order], fontsize=6)
    ax.set_title(f"Module-trait correlation (Pearson r; * FDR < 0.05, "
                 f"{int((mt.padj < 0.05).sum())} of {len(mt)})", loc="left")
    fig.colorbar(im, ax=ax, shrink=0.6, label="r")
    panel_label(ax, "c")
    save(fig, "fig3_wgcna_modules")


# --------------------------------------------------------------------------- figure 4
ML_LABELS = {"TabPFN_v2": "TabPFN v2", "LogReg_L2": "logistic regression (L2)",
             "RandomForest_500": "random forest (500 trees)", "Majority_baseline": "majority class"}
ML_COLOURS = {"TabPFN_v2": "#b2182b", "LogReg_L2": "#2166ac", "RandomForest_500": "#636363",
              "Majority_baseline": "#d9d9d9"}


def fig4_classification():
    print("Figure 4: leave-one-out classification (TabPFN and baselines)")
    if not have(path("ml_metrics.tsv")):
        return
    m = pd.read_csv(path("ml_metrics.tsv"), sep="\t")
    need = {"analysis", "task", "k", "model", "n_samples", "balanced_accuracy", "p_balanced_accuracy", "class_counts"}
    if not need.issubset(m.columns):
        print(f"  skipped: ml_metrics.tsv lacks {sorted(need - set(m.columns))}")
        return
    prim = m[m.analysis == "primary"]
    tasks = [("binary_microgravity_vs_ground", "Simulated microgravity vs static 1g"),
             ("three_class_modality", "Clinostat vs RPM vs static 1g")]
    models = [x for x in ML_LABELS if x in set(prim.model)]
    ks = sorted(prim.k.unique())
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), constrained_layout=True)
    width = 0.8 / len(models)
    for ax, (task, title), letter in zip(axes, tasks, "ab"):
        sub = prim[prim.task == task]
        for i, model in enumerate(models):
            for j, k in enumerate(ks):
                r = sub[(sub.model == model) & (sub.k == k)]
                if r.empty:
                    continue
                r = r.iloc[0]
                x = j + (i - (len(models) - 1) / 2) * width
                ax.bar(x, r.balanced_accuracy, width, color=ML_COLOURS[model], edgecolor="black", lw=0.3,
                       label=ML_LABELS[model] if j == 0 else None)
                if model != "Majority_baseline" and pd.notna(r.p_balanced_accuracy):
                    ax.text(x, r.balanced_accuracy + 0.02, f"p={r.p_balanced_accuracy:.2g}", ha="center",
                            va="bottom", fontsize=4.8, rotation=90)
        n_classes = len(sub.class_counts.iloc[0].split(";"))
        ax.axhline(1 / n_classes, color="#636363", lw=0.6, ls="--")
        ax.text(len(ks) - 0.55, 1 / n_classes + 0.01, "chance", fontsize=5.5, color="#636363", ha="right")
        ax.set_xticks(range(len(ks)))
        ax.set_xticklabels([f"top {k} genes" for k in ks])
        ax.set_ylim(0, 1.3)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_ylabel("LOOCV balanced accuracy")
        counts = sub.class_counts.iloc[0].replace("=", " ").replace(";", ", ").replace("_", " ")
        ax.set_title(f"{title}\nn = {int(sub.n_samples.iloc[0])} ({counts})", loc="left", fontsize=7)
        panel_label(ax, letter)
    axes[0].legend(frameon=False, fontsize=5.5, loc="upper left", ncol=2)
    save(fig, "fig4_classification")


# --------------------------------------------------------------------------- figure 5
GO_SET_LABELS = {
    "shared_simulated_microgravity_up": "shared up", "shared_simulated_microgravity_down": "shared down",
    "clinostat_vs_static1g_up": "clinostat up", "clinostat_vs_static1g_down": "clinostat down",
    "rpm2_vs_static1g_up": "RPM up", "rpm2_vs_static1g_down": "RPM down",
}


def fig5_go_enrichment():
    print("Figure 5: GO enrichment of the shared response")
    if not have(path("go_enrichment_deg.tsv"), path("deg_shared_simulated_microgravity.tsv")):
        return
    g = pd.read_csv(path("go_enrichment_deg.tsv"), sep="\t")
    need = {"gene_set", "ontology", "go_id", "term_name", "overlap", "term_size", "fold_enrichment", "fdr"}
    if not need.issubset(g.columns):
        print(f"  skipped: go_enrichment_deg.tsv lacks {sorted(need - set(g.columns))}")
        return
    shared = pd.read_csv(path("deg_shared_simulated_microgravity.tsv"), sep="\t")
    n_up = int((shared.log2FoldChange_clinostat > 0).sum())
    n_down = int((shared.log2FoldChange_clinostat < 0).sum())
    sets = ["clinostat_vs_static1g_up", "rpm2_vs_static1g_up", "shared_simulated_microgravity_up",
            "clinostat_vs_static1g_down", "rpm2_vs_static1g_down", "shared_simulated_microgravity_down"]
    picked = []
    for gset in ("shared_simulated_microgravity_up", "shared_simulated_microgravity_down"):
        picked += g[(g.gene_set == gset) & (g.fdr < 0.05)].sort_values("fdr").head(8).go_id.tolist()
    picked = list(dict.fromkeys(picked))
    if not picked:
        print("  no GO term at FDR < 0.05 in the shared sets; drawing a note instead")
        fig, ax = plt.subplots(figsize=(7.2, 1.2))
        ax.axis("off")
        ax.text(0.5, 0.5, "No GO term reached FDR < 0.05 in the shared up or down genes.", ha="center")
        save(fig, "fig5_go_enrichment")
        return
    sub = g[g.go_id.isin(picked) & g.gene_set.isin(sets)].set_index(["go_id", "gene_set"])
    names = g.drop_duplicates("go_id").set_index("go_id").loc[picked]
    fig, ax = plt.subplots(figsize=(7.2, 0.9 + 0.24 * len(picked)), constrained_layout=True)
    fdrs = sub.fdr.clip(lower=1e-300)
    vmax = float(-np.log10(fdrs.min()))
    for i, go in enumerate(picked):
        for k, gset in enumerate(sets):
            if (go, gset) not in sub.index:
                ax.text(k, i, "n.t.", ha="center", va="center", fontsize=5, color="#9e9e9e")
                continue
            r = sub.loc[(go, gset)]
            size = 12 + 18 * min(float(r.fold_enrichment), 8)
            colour = plt.cm.viridis(min(-np.log10(max(r.fdr, 1e-300)) / vmax, 1.0))
            ax.scatter(k, i, s=size, color=colour, edgecolor="black" if r.fdr < 0.05 else "none",
                       lw=0.7, zorder=3)
            if r.fdr < 0.05:
                ax.text(k + 0.28, i, f"{int(r.overlap)}/{int(r.term_size)}", va="center", fontsize=4.8)
    ax.set_yticks(range(len(picked)))
    ax.set_yticklabels([f"{n if len(n) <= 60 else n[:58] + '...'} [{o}]"
                        for n, o in zip(names.term_name, names.ontology)], fontsize=6)
    ax.invert_yaxis()
    labels = {"clinostat_vs_static1g_up": "clinostat\nup", "rpm2_vs_static1g_up": "RPM\nup",
              "shared_simulated_microgravity_up": f"shared up\n(n = {n_up})",
              "clinostat_vs_static1g_down": "clinostat\ndown", "rpm2_vs_static1g_down": "RPM\ndown",
              "shared_simulated_microgravity_down": f"shared down\n(n = {n_down})"}
    ax.set_xticks(range(len(sets)))
    ax.set_xticklabels([labels[x] for x in sets], fontsize=6.5)
    ax.set_xlim(-0.5, len(sets) - 0.2)
    ax.axvline(2.5, color="#bdbdbd", lw=0.6)
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(0, vmax))
    fig.colorbar(sm, ax=ax, shrink=0.5, label="$-$log$_{10}$ FDR")
    for fe in (2, 4, 8):
        ax.scatter([], [], s=12 + 18 * fe, color="#bdbdbd", label=f"{fe}x")
    ax.legend(title="fold enrichment", frameon=False, fontsize=5.5, title_fontsize=6,
              loc="upper left", bbox_to_anchor=(1.0, 0.25))
    ax.set_title("GO terms enriched in the shared response (top 8 per direction); outlined: FDR < 0.05",
                 loc="left", fontsize=7)
    save(fig, "fig5_go_enrichment")


# --------------------------------------------------------------------------- figure 6
def fig6_hub_connectivity():
    print("Figure 6: hub genes of the microgravity-associated modules")
    if not have(path("wgcna_module_assignments.tsv"), path("wgcna_module_trait_correlations.tsv")):
        return
    a = pd.read_csv(path("wgcna_module_assignments.tsv"), sep="\t")
    mt = pd.read_csv(path("wgcna_module_trait_correlations.tsv"), sep="\t")
    hits = mt[(mt.trait == "simulated_microgravity") & (mt.padj < 0.05)].sort_values("padj")
    modules = [m.replace("ME", "") for m in hits.module]
    if not modules:
        modules = [mt[mt.trait == "simulated_microgravity"].sort_values("padj").module.iloc[0].replace("ME", "")]
    fig, axes = plt.subplots(len(modules), 2, figsize=(7.2, 2.6 * len(modules)), constrained_layout=True,
                             squeeze=False, gridspec_kw={"width_ratios": [1, 1.7]})
    letters = iter("abcdefgh")
    for j, mod in enumerate(modules):
        sub = a[a.module == mod]
        row = mt[(mt.module == "ME" + mod) & (mt.trait == "simulated_microgravity")].iloc[0]
        colour = UP if row.pearson_r > 0 else DOWN
        top = sub.sort_values("kWithin", ascending=False).head(10)
        ax = axes[j][0]
        ax.scatter(sub.kME, sub.kWithin, s=5, c="#636363", lw=0, alpha=0.6, rasterized=True)
        ax.scatter(top.kME, top.kWithin, s=10, c=colour, lw=0)
        ax.set_xlabel("module membership (kME)")
        ax.set_ylabel("intramodular connectivity")
        ax.set_title(f"{mod}: {len(sub):,} genes\nr = {row.pearson_r:.2f} with simulated "
                     f"microgravity (FDR {row.padj:.2g})", loc="left", fontsize=7)
        panel_label(ax, next(letters))
        ax = axes[j][1]
        t = top[::-1]
        labels = []
        for _, r in t.iterrows():
            prot = str(r.protein) if pd.notna(r.protein) else ""
            prot = prot if len(prot) <= 48 else prot[:46] + "..."
            name = r.gene_symbol if r.gene_symbol != r.gene_id else r.gene_id
            labels.append(f"{name}  {prot}")
        ax.barh(range(len(t)), t.kWithin, color=colour, edgecolor="black", lw=0.3)
        ax.set_yticks(range(len(t)))
        ax.set_yticklabels(labels, fontsize=6)
        ax.set_xlabel("intramodular connectivity (kWithin)")
        ax.set_title(f"Top 10 hub genes ({mod})", loc="left", fontsize=7)
        panel_label(ax, next(letters))
    save(fig, "fig6_hub_connectivity")


# --------------------------------------------------------------------------- figure 7
def fig7_osdr_microbial_compendium():
    print("Figure 7: OSDR microbial spaceflight and simulated-microgravity studies")
    tables = sorted(glob.glob(os.path.join(PROJECT_DIR, "results_tables", "Table_S16_*.tsv")))
    if len(tables) != 1:
        print(f"  skipped: expected one results_tables/Table_S16_*.tsv, found {len(tables)}")
        return
    t = pd.read_csv(tables[0], sep="\t")
    need = {"accession", "organisms", "microbial_categories", "assay_types", "inclusion_basis"}
    if not need.issubset(t.columns):
        print(f"  skipped: Table S16 lacks {sorted(need - set(t.columns))}")
        return
    t["setting"] = np.where(t.inclusion_basis.str.contains("spaceflight", case=False), "spaceflight",
                            "ground simulation")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9), constrained_layout=True,
                             gridspec_kw={"width_ratios": [1, 1.1, 1.3]})
    ax = axes[0]
    cats = t.assign(cat=t.microbial_categories.str.split(r"\s*;\s*")).explode("cat")
    tab = pd.crosstab(cats.cat, cats.setting).reindex(columns=["spaceflight", "ground simulation"], fill_value=0)
    tab = tab.loc[tab.sum(axis=1).sort_values().index]
    left = np.zeros(len(tab))
    for col, colour in (("spaceflight", "#2166ac"), ("ground simulation", "#ef8a62")):
        ax.barh(range(len(tab)), tab[col], left=left, color=colour, label=col, edgecolor="black", lw=0.3)
        left += tab[col].values
    ax.set_yticks(range(len(tab)))
    ax.set_yticklabels(tab.index, fontsize=6.5)
    ax.set_xlabel("studies")
    ax.set_title(f"{len(t)} microbial OSDR studies", loc="left", fontsize=7)
    ax.legend(frameon=False, fontsize=5.5, loc="upper center", bbox_to_anchor=(0.4, -0.2), ncol=1)
    panel_label(ax, "a")
    ax = axes[1]
    single = t[t.microbial_categories != "microbial community"]
    genus = single.organisms.str.split(r"\s*;\s*").explode().str.split().str[0].value_counts().head(12)[::-1]
    ax.barh(range(len(genus)), genus.values, color="#636363", edgecolor="black", lw=0.3)
    ax.set_yticks(range(len(genus)))
    ax.set_yticklabels(genus.index, style="italic", fontsize=6.5)
    ax.set_xlabel("studies")
    ax.set_title("Most-studied genera\n(single-organism studies)", loc="left", fontsize=7)
    panel_label(ax, "b")
    ax = axes[2]
    # One count per study per assay type; drop the platform detail in parentheses.
    assays = (t.assign(a=t.assay_types.str.split(r"\s*;\s*")).explode("a")
              .assign(a=lambda d: d.a.str.split(" (", regex=False).str[0].str.strip().str.lower())
              .drop_duplicates(["accession", "a"]).a.value_counts().head(10)[::-1])
    ax.barh(range(len(assays)), assays.values, color="#9e9e9e", edgecolor="black", lw=0.3)
    ax.set_yticks(range(len(assays)))
    def tidy(a):  # OSDR assay names vary in case; restore the usual capitals after lower-casing
        for low, up in (("16s", "16S"), ("its", "ITS"), ("rna", "RNA"), ("dna", "DNA")):
            a = re.sub(rf"\b{low}\b", up, a)
        return a if len(a) <= 40 else a[:38] + "..."
    ax.set_yticklabels([tidy(a) for a in assays.index], fontsize=5.8)
    ax.set_xlabel("studies")
    ax.set_title("Assay types (studies using each)", loc="left", fontsize=7)
    panel_label(ax, "c")
    save(fig, "fig7_osdr_microbial_compendium")


# --------------------------------------------------------------------------- figure 8
def fig8_simulator_concordance():
    print("Figure 8: concordance between the two simulators")
    files = [path(t) for t, _ in CONTRASTS.values()] + [path("de_summary.json")]
    if not have(*files):
        return
    de = load_json("de_summary.json")
    ov = de["primary"]["simulator_overlap"]
    a = pd.read_csv(path(CONTRASTS["clinostat_vs_ground"][0]), sep="\t")
    b = pd.read_csv(path(CONTRASTS["rpm_vs_ground"][0]), sep="\t")
    m = a.merge(b, on="gene_id", suffixes=("_c", "_r"))
    both = (m.significant_c == "YES") & (m.significant_r == "YES")
    only_c = (m.significant_c == "YES") & ~both
    only_r = (m.significant_r == "YES") & ~both
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), constrained_layout=True,
                             gridspec_kw={"width_ratios": [1.4, 1]})
    ax = axes[0]
    rest = ~(both | only_c | only_r)
    ax.scatter(m.log2FoldChange_shrunk_c[rest], m.log2FoldChange_shrunk_r[rest], s=3, c=NS, lw=0, rasterized=True,
               label=f"neither ({rest.sum()})")
    ax.scatter(m.log2FoldChange_shrunk_c[only_c], m.log2FoldChange_shrunk_r[only_c], s=4, c=GROUP_COLOURS["Clinostat"],
               lw=0, rasterized=True, label=f"clinostat only ({only_c.sum()})")
    ax.scatter(m.log2FoldChange_shrunk_c[only_r], m.log2FoldChange_shrunk_r[only_r], s=4, c=GROUP_COLOURS["RPM"],
               lw=0, rasterized=True, label=f"RPM only ({only_r.sum()})")
    ax.scatter(m.log2FoldChange_shrunk_c[both], m.log2FoldChange_shrunk_r[both], s=4, c="black", lw=0,
               rasterized=True, label=f"both ({both.sum()})")
    lo = np.nanmin(np.r_[m.log2FoldChange_shrunk_c, m.log2FoldChange_shrunk_r]) * 1.05
    hi = np.nanmax(np.r_[m.log2FoldChange_shrunk_c, m.log2FoldChange_shrunk_r]) * 1.05
    ax.plot([lo, hi], [lo, hi], color="#636363", lw=0.5, ls="--")
    ax.set_xlim(np.nanmin(m.log2FoldChange_shrunk_c) * 1.1, np.nanmax(m.log2FoldChange_shrunk_c) * 1.1)
    ax.set_ylim(np.nanmin(m.log2FoldChange_shrunk_r) * 1.05, np.nanmax(m.log2FoldChange_shrunk_r) * 1.1)
    ax.set_xlabel("log$_2$FC, 3D clinostat vs 1g (shrunk)")
    ax.set_ylabel("log$_2$FC, RPM 2.0 vs 1g (shrunk)")
    ax.set_title(f"All {len(m):,} tested genes; Spearman $\\rho$ = {ov['spearman_lfc_all_tested']:.2f} "
                 f"(MLE fold changes)", loc="left", fontsize=7)
    ax.legend(frameon=False, fontsize=6, loc="upper left", markerscale=2)
    panel_label(ax, "a")

    ax = axes[1]
    cont = de["primary"]["contrasts"]
    labels = ["clinostat vs 1g", "RPM vs 1g", "clinostat vs RPM"]
    keys = ["clinostat_vs_ground", "rpm_vs_ground", "clinostat_vs_rpm"]
    ups = [cont[k]["up"] for k in keys]
    downs = [cont[k]["down"] for k in keys]
    x = np.arange(len(keys))
    ax.bar(x, ups, color=UP, label="up", edgecolor="black", lw=0.3)
    ax.bar(x, [-d for d in downs], color=DOWN, label="down", edgecolor="black", lw=0.3)
    for i in range(len(keys)):
        ax.text(i, ups[i] + 8, str(ups[i]), ha="center", fontsize=6)
        ax.text(i, -downs[i] - 8, str(downs[i]), ha="center", va="top", fontsize=6)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=6.5)
    ax.set_ylabel(f"DEGs (FDR < {de['method']['alpha']})")
    ax.set_ylim(-max(downs) * 1.25 - 20, max(ups) * 1.35 + 20)
    ax.set_title(f"{ov['deg_both_same_direction']} shared DEGs, all same direction "
                 f"(Jaccard {ov['jaccard_same_direction']:.2f})"
                 if ov["deg_both_same_direction"] == ov["deg_both"] else
                 f"{ov['deg_both']} shared DEGs, {ov['deg_both_same_direction']} same direction",
                 loc="left", fontsize=7)
    ax.legend(frameon=False, fontsize=6, loc="upper right")
    panel_label(ax, "b")
    save(fig, "fig8_simulator_concordance")


# --------------------------------------------------------------------------- figure 9
def fig9_metabolic_subsystems():
    print("Figure 9: KEGG-mapped metabolic subsystems")
    if not have(path("metabolic_subsystem_perturbation.tsv")):
        return
    s = pd.read_csv(path("metabolic_subsystem_perturbation.tsv"), sep="\t")
    need = {"subsystem", "contrast", "median_log2fc_shrunk", "n_deg", "n_genes_tested", "bh_fdr", "n_up", "n_down"}
    if not need.issubset(s.columns):
        print(f"  skipped: metabolic_subsystem_perturbation.tsv lacks {sorted(need - set(s.columns))}")
        return
    labels = {"clinostat_vs_static1g": "3D clinostat vs 1g", "rpm2_vs_static1g": "RPM 2.0 vs 1g"}
    contrasts = [c for c in labels if c in set(s.contrast)]
    piv = s.pivot(index="subsystem", columns="contrast", values="median_log2fc_shrunk")[contrasts]
    order = (s.groupby("subsystem").bh_fdr.min().sort_values(ascending=False)).index
    piv = piv.loc[order]
    fig, ax = plt.subplots(figsize=(5.6, 0.7 + 0.42 * len(piv)), constrained_layout=True)
    lim = float(np.nanmax(np.abs(piv.values))) or 1.0
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    for i, sub in enumerate(piv.index):
        for j, c in enumerate(contrasts):
            r = s[(s.subsystem == sub) & (s.contrast == c)].iloc[0]
            star = "*" if r.bh_fdr < 0.05 else ""
            ax.text(j, i, f"{r.median_log2fc_shrunk:+.2f}{star}\n{int(r.n_deg)}/{int(r.n_genes_tested)} DEGs "
                          f"({int(r.n_up)} up, {int(r.n_down)} down)", ha="center", va="center", fontsize=5.6,
                    color="white" if abs(r.median_log2fc_shrunk) > 0.6 * lim else "black")
    ax.set_xticks(range(len(contrasts)))
    ax.set_xticklabels([labels[c] for c in contrasts], fontsize=7)
    ax.set_yticks(range(len(piv.index)))
    ax.set_yticklabels(piv.index, fontsize=6.5)
    ax.set_title("KEGG-mapped subsystems: median shrunk log$_2$FC, DEGs/tested genes\n"
                 "(* DEG enrichment BH FDR < 0.05)", loc="left", fontsize=7)
    fig.colorbar(im, ax=ax, shrink=0.7, label="median shrunk log$_2$FC")
    save(fig, "fig9_metabolic_subsystems")


# --------------------------------------------------------------------------- supplementary
SUPP_DIR = os.path.join(PROJECT_DIR, "supplementary_figures")


def save_supp(fig, name):
    os.makedirs(SUPP_DIR, exist_ok=True)
    for ext, kw in (("pdf", {"metadata": {"CreationDate": None, "ModDate": None}}),
                    ("png", {"dpi": 200, "metadata": {"Software": None}})):
        fig.savefig(os.path.join(SUPP_DIR, f"{name}.{ext}"), bbox_inches="tight", **kw)
    plt.close(fig)
    print(f"  wrote supplementary_figures/{name}.pdf/.png")


def figS1_sensitivity():
    print("Supplementary Figure S1: sensitivity analysis with the inferred RFPRPM4 library")
    sens = os.path.join(PROC, "sensitivity_with_rfprpm4")
    files = [os.path.join(sens, "pca_scores.tsv"), os.path.join(sens, "deg_rpm2_vs_static1g.tsv"),
             path("deg_rpm2_vs_static1g.tsv"), path("de_summary.json")]
    if not have(*files):
        return
    de = load_json("de_summary.json")
    ve = de["sensitivity_with_rfprpm4"]["pca"]["variance_explained"]
    pca = pd.read_csv(files[0], sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), constrained_layout=True)
    ax = axes[0]
    for group, sub in pca.groupby("group"):
        ax.scatter(sub.PC1, sub.PC2, s=28, c=GROUP_COLOURS[group], label=GROUP_LABELS[group],
                   edgecolor="black", lw=0.4, zorder=3)
    label_samples(ax, pca, bold=("RFPRPM4",))
    ax.set_xlabel(f"PC1 ({100 * ve[0]:.1f}% of variance)")
    ax.set_ylabel(f"PC2 ({100 * ve[1]:.1f}% of variance)")
    ax.set_title("PCA with all nine samples (inferred RFPRPM4 in bold)", loc="left", fontsize=7)
    ax.legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3)
    panel_label(ax, "a")
    ax = axes[1]
    a = pd.read_csv(path("deg_rpm2_vs_static1g.tsv"), sep="\t")
    b = pd.read_csv(files[1], sep="\t")
    m = a.merge(b, on="gene_id", suffixes=("_p", "_s"))
    cmp = de["primary_vs_sensitivity"]["rpm_vs_ground"]
    sig_p, sig_s = m.significant_p == "YES", m.significant_s == "YES"
    for mask, colour, label in ((~sig_p & ~sig_s, NS, "neither"), (sig_p & ~sig_s, GROUP_COLOURS["Clinostat"], "primary only"),
                                (~sig_p & sig_s, DOWN, "with RFPRPM4 only"), (sig_p & sig_s, "black", "both")):
        ax.scatter(m.log2FoldChange_p[mask], m.log2FoldChange_s[mask], s=3, c=colour, lw=0, rasterized=True,
                   label=f"{label} ({int(mask.sum())})")
    lo = np.nanmin(np.r_[m.log2FoldChange_p, m.log2FoldChange_s])
    hi = np.nanmax(np.r_[m.log2FoldChange_p, m.log2FoldChange_s])
    ax.plot([lo, hi], [lo, hi], color="#636363", lw=0.5, ls="--")
    ax.set_xlabel("log$_2$FC RPM vs 1g, primary (RPM n = 2)")
    ax.set_ylabel("log$_2$FC RPM vs 1g, with RFPRPM4 (n = 3)")
    ax.set_title(f"RPM vs 1g fold changes; Pearson r = {cmp['pearson_lfc']:.2f}", loc="left", fontsize=7)
    ax.legend(frameon=False, fontsize=5.5, loc="upper left", markerscale=2)
    panel_label(ax, "b")
    save_supp(fig, "Supplementary_Figure_S1_sensitivity_RFPRPM4")


def figS2_go_modules():
    print("Supplementary Figure S2: GO enrichment of WGCNA modules")
    if not have(path("go_enrichment_modules.tsv")):
        return
    g = pd.read_csv(path("go_enrichment_modules.tsv"), sep="\t")
    sig = g[g.fdr < 0.05].sort_values("fdr")
    if sig.empty:
        print("  no module GO term at FDR < 0.05; skipped")
        return
    top = sig.groupby("gene_set").head(3)
    order = list(dict.fromkeys(top.sort_values(["gene_set", "fdr"]).gene_set))
    rows_ = []
    for gset in order:
        for _, r in top[top.gene_set == gset].iterrows():
            rows_.append((gset.replace("module_", ""), r))
    fig, ax = plt.subplots(figsize=(7.2, 0.6 + 0.2 * len(rows_)), constrained_layout=True)
    y = np.arange(len(rows_))[::-1]
    vals = [-np.log10(r.fdr) for _, r in rows_]
    ax.barh(y, vals, color=[module_colour(mod) for mod, _ in rows_], edgecolor="black", lw=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{mod}: {r.term_name if len(r.term_name) <= 55 else r.term_name[:53] + '...'} [{r.ontology}]"
                        for mod, r in rows_], fontsize=5.8)
    for yy, (_, r), v in zip(y, rows_, vals):
        ax.text(v + 0.03, yy, f"{int(r.overlap)}/{int(r.term_size)}", va="center", fontsize=5)
    ax.set_xlabel("$-$log$_{10}$ FDR")
    ax.set_title("WGCNA modules: top 3 GO terms at FDR < 0.05 (overlap/term size)", loc="left", fontsize=7)
    save_supp(fig, "Supplementary_Figure_S2_go_wgcna_modules")


FIGURES = [fig1_study_design, fig2_volcano_pca, fig3_wgcna_modules, fig4_classification,
           fig5_go_enrichment, fig6_hub_connectivity, fig7_osdr_microbial_compendium,
           fig8_simulator_concordance, fig9_metabolic_subsystems, figS1_sensitivity, figS2_go_modules]


def main(selected):
    for f in FIGURES:
        if not selected or any(f.__name__.startswith(s) for s in selected):
            f()


if __name__ == "__main__":
    main(sys.argv[1:])
