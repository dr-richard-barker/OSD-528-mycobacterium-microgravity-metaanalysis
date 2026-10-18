#!/usr/bin/env python3
"""
04_tabpfn_tabular_foundation_ai.py

Leave-one-out classification of the OSD-528 samples with the real TabPFN
(Hollmann et al., Nature 2025; TabPFN v2 weights from Prior Labs) and two
conventional baselines, with exact label-permutation tests.

The file name is kept for continuity with the rest of the pipeline. An earlier
version of this script did not run TabPFN at all (its "TabPFNClassifier" was a
Gaussian-kernel nearest-neighbour rule) and reported an OSD-90 transfer that was
never run. Its outputs are the old ``tabpfn_*`` files. This script writes
``data/processed/ml_*`` files only.

Design
------
* Features: DESeq2 variance-stabilised expression (``osd528_vst.tsv``), samples
  as rows and genes as columns.
* Tasks: (a) binary, simulated microgravity (clinostat + RPM) vs 1 g ground;
  (b) three-class, clinostat vs RPM vs ground.
* Evaluation: leave-one-out cross-validation (LOOCV) over samples. Inside each
  training fold, and using the training samples only, the k most variable genes
  are selected (k = 50 and k = 200, both reported) and standardised with the
  training mean and SD. The held-out sample is transformed with those
  statistics. Every model sees the same folds and the same features.
* Models: TabPFN v2 (``TabPFNClassifier``), L2 logistic regression, a 500-tree
  random forest, and a majority-class baseline that shows what LOOCV gives with
  no information at all.
* Metrics: LOOCV accuracy and balanced accuracy for both tasks, and AUROC from
  the pooled held-out probabilities for the binary task.
* Permutation test: the whole LOOCV pipeline is rerun under every distinct
  relabelling of the samples that keeps the class sizes (exact test). If there
  are more relabellings than ``--max-exact``, a Monte-Carlo test with
  ``--n-mc`` relabellings and a fixed seed is used instead.
* Feature importance (descriptive only): L2 logistic regression refitted on all
  samples for the binary task; standardised coefficients of the top genes,
  with how often each gene was selected, and kept the same sign, across the
  LOOCV training folds. Not validated.
* OSD-90 cross-study transfer: not performed.

Analyses: ``primary`` (the 8 samples with primary == yes in config/samples.tsv)
and ``sensitivity_9`` (adds the inferred RFPRPM4 library).

Outputs (data/processed/): ml_loocv_predictions.tsv, ml_metrics.tsv,
ml_feature_importance.tsv, ml_permutation_null.tsv, ml_summary.json.

Run with the dedicated environment (see requirements-ml.txt):

    .venv-ml/bin/python analysis/04_tabpfn_tabular_foundation_ai.py --n-jobs 12

TabPFN runs on the CPU. The script never sends a Hugging Face or Prior Labs
token, never opens the Prior Labs licence browser flow, and turns off TabPFN
telemetry. The v2 weights are an ungated download.
"""

import os

# Set before tabpfn / huggingface_hub are imported; worker processes inherit them.
os.environ["TABPFN_DISABLE_TELEMETRY"] = "1"
os.environ["TABPFN_NO_BROWSER"] = "1"
os.environ["HF_HUB_DISABLE_IMPLICIT_TOKEN"] = "1"
os.environ["HF_HUB_DISABLE_TELEMETRY"] = "1"
os.environ.pop("HF_TOKEN", None)
os.environ.pop("TABPFN_TOKEN", None)
for _var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS"):
    os.environ[_var] = "1"

import argparse
import hashlib
import importlib.metadata as importlib_metadata
import json
import math
import platform
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_PROCESSED = PROJECT_DIR / "data" / "processed"
SAMPLES_FILE = PROJECT_DIR / "config" / "samples.tsv"
VST_FILES = {
    "primary": DATA_PROCESSED / "osd528_vst.tsv",
    "sensitivity_9": DATA_PROCESSED / "sensitivity_with_rfprpm4" / "osd528_vst.tsv",
}
ANNOTATION_FILE = DATA_PROCESSED / "osd528_tpm.tsv"  # optional, for gene symbols

SEED = 528
K_VALUES = (50, 200)
MODELS = ("TabPFN_v2", "LogReg_L2", "RandomForest_500", "Majority_baseline")
TASKS = {
    "binary_microgravity_vs_ground": ("ground_1g", "simulated_microgravity"),
    "three_class_modality": ("Clinostat", "Ground", "RPM"),
}
POSITIVE_CLASS = "simulated_microgravity"
N_TOP_IMPORTANCE = 25
LOGREG_PARAMS = dict(penalty="l2", C=1.0, solver="lbfgs", max_iter=10000)
RF_PARAMS = dict(n_estimators=500, random_state=SEED, n_jobs=1)
TOL = 1e-9


# --------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_samples():
    s = pd.read_csv(SAMPLES_FILE, sep="\t", dtype=str)
    for col in ("sample_id", "group", "microgravity", "primary"):
        if col not in s.columns:
            raise SystemExit(f"{SAMPLES_FILE} lacks column {col!r}")
    s["microgravity"] = s["microgravity"].str.strip().str.lower()
    s["primary"] = s["primary"].str.strip().str.lower()
    bad = set(s["group"]) - set(TASKS["three_class_modality"])
    if bad:
        raise SystemExit(f"Unexpected group values in {SAMPLES_FILE}: {sorted(bad)}")
    return s


def load_analysis(name, samples):
    """Return (sample table in samples.tsv order, genes x samples VST frame)."""
    path = VST_FILES[name]
    vst = pd.read_csv(path, sep="\t", index_col="gene_id")
    if name == "primary":
        expected = list(samples.loc[samples["primary"] == "yes", "sample_id"])
    else:
        expected = list(samples["sample_id"])
    if sorted(vst.columns) != sorted(expected):
        raise SystemExit(
            f"{path} has samples {sorted(vst.columns)}, expected {sorted(expected)} "
            f"for the {name} analysis"
        )
    s = samples.set_index("sample_id").loc[expected].reset_index()
    vst = vst[expected].apply(pd.to_numeric, errors="coerce")
    vst = vst[np.isfinite(vst).all(axis=1)]
    return s, vst


def labels_for(task, s):
    classes = TASKS[task]
    if task == "binary_microgravity_vs_ground":
        names = np.where(s["microgravity"] == "yes", "simulated_microgravity", "ground_1g")
    else:
        names = s["group"].to_numpy()
    y = np.array([classes.index(v) for v in names], dtype=int)
    return y, list(names)


def _top_k(X, genes, k):
    var = X.var(axis=0, ddof=1)
    sel = np.lexsort((genes, -var))[:k]  # highest variance first, ties by gene id
    if var[sel].min() <= 0:
        raise SystemExit(f"Fewer than k={k} non-constant genes")
    return sel


def fold_features(vst, k):
    """Label-free LOOCV features: top-k variance genes and scaling, training only."""
    X_all = vst.to_numpy().T  # samples x genes
    genes = vst.index.to_numpy()
    n = X_all.shape[0]
    folds = []
    for i in range(n):
        train = np.array([j for j in range(n) if j != i])
        sel = _top_k(X_all[train], genes, k)
        mu = X_all[train][:, sel].mean(axis=0)
        sd = X_all[train][:, sel].std(axis=0, ddof=0)
        folds.append({
            "test": i,
            "train": train,
            "genes": genes[sel],
            "X_train": (X_all[train][:, sel] - mu) / sd,
            "X_test": (X_all[[i]][:, sel] - mu) / sd,
        })
    return folds


def full_data_features(vst, k):
    X_all = vst.to_numpy().T
    genes = vst.index.to_numpy()
    sel = _top_k(X_all, genes, k)
    mu = X_all[:, sel].mean(axis=0)
    sd = X_all[:, sel].std(axis=0, ddof=0)
    return genes[sel], (X_all[:, sel] - mu) / sd


# --------------------------------------------------------------------------
# Models
# --------------------------------------------------------------------------

def _quiet():
    warnings.filterwarnings("ignore")
    import torch
    torch.set_num_threads(1)


def make_model(name):
    if name == "TabPFN_v2":
        from tabpfn import TabPFNClassifier
        from tabpfn.constants import ModelVersion
        return TabPFNClassifier.create_default_for_version(
            ModelVersion.V2, device="cpu", random_state=SEED, n_preprocessing_jobs=1)
    if name == "LogReg_L2":
        from sklearn.linear_model import LogisticRegression
        return LogisticRegression(random_state=SEED, **LOGREG_PARAMS)
    if name == "RandomForest_500":
        from sklearn.ensemble import RandomForestClassifier
        return RandomForestClassifier(**RF_PARAMS)
    if name == "Majority_baseline":
        from sklearn.dummy import DummyClassifier
        return DummyClassifier(strategy="prior")
    raise ValueError(name)


def loocv(folds, y, n_classes):
    """Held-out class probabilities, one (n, n_classes) array per model."""
    out = {m: np.zeros((len(y), n_classes)) for m in MODELS}
    for f in folds:
        y_train = y[f["train"]]
        for m in MODELS:
            clf = make_model(m)
            clf.fit(f["X_train"], y_train)
            p = np.asarray(clf.predict_proba(f["X_test"]), dtype=float)[0]
            out[m][f["test"], np.asarray(clf.classes_, dtype=int)] = p
    return out


def score(y, prob, n_classes):
    from sklearn.metrics import balanced_accuracy_score, roc_auc_score
    pred = prob.argmax(axis=1)  # ties go to the first class in TASKS order
    acc = float(np.mean(pred == y))
    bacc = float(balanced_accuracy_score(y, pred))
    auc = float(roc_auc_score(y, prob[:, 1])) if n_classes == 2 else float("nan")
    return acc, bacc, auc


def run_batch(folds, labelings, n_classes):
    """Worker: LOOCV for each labelling; per-model (metrics, held-out probs)."""
    _quiet()
    results = []
    for y in labelings:
        y = np.asarray(y, dtype=int)
        probs = loocv(folds, y, n_classes)
        results.append({m: (score(y, probs[m], n_classes), probs[m]) for m in MODELS})
    return results


# --------------------------------------------------------------------------
# Permutations
# --------------------------------------------------------------------------

def distinct_labelings(y):
    """All distinct orderings of the multiset y, in lexicographic order."""
    counts = list(np.bincount(y))
    n = len(y)
    out, cur = [], []

    def rec():
        if len(cur) == n:
            out.append(tuple(cur))
            return
        for c in range(len(counts)):
            if counts[c]:
                counts[c] -= 1
                cur.append(c)
                rec()
                cur.pop()
                counts[c] += 1

    rec()
    return out


def n_distinct_labelings(y):
    total = math.factorial(len(y))
    for c in np.bincount(y):
        total //= math.factorial(int(c))
    return total


def partition_of(y):
    y = np.asarray(y)
    return frozenset(frozenset(np.flatnonzero(y == c).tolist()) for c in set(y.tolist()))


def permutation_plan(y, max_exact, n_mc, rng):
    total = n_distinct_labelings(y)
    if total <= max_exact:
        return distinct_labelings(y), "exact", total
    labs = [tuple(y)] + [tuple(rng.permutation(y)) for _ in range(n_mc)]
    return labs, f"monte_carlo_{n_mc}", total


def p_values(stats, obs_idx, method):
    if method == "not_run":
        return [float("nan")] * 3
    out = []
    for j in range(stats.shape[1]):
        obs = stats[obs_idx, j]
        if np.isnan(obs):
            out.append(float("nan"))
        elif method == "exact":  # observed labelling is one of the enumerated ones
            out.append(float(np.mean(stats[:, j] >= obs - TOL)))
        else:  # Monte-Carlo: (1 + #random >= observed) / (1 + #random)
            others = np.delete(stats[:, j], obs_idx)
            out.append(float((1 + np.sum(others >= obs - TOL)) / (1 + len(others))))
    return out


# --------------------------------------------------------------------------
# Provenance
# --------------------------------------------------------------------------

def versions():
    pkgs = ["tabpfn", "torch", "scikit-learn", "numpy", "pandas", "scipy",
            "joblib", "huggingface-hub"]
    v = {"python": platform.python_version(), "platform": platform.platform()}
    for p in pkgs:
        try:
            v[p] = importlib_metadata.version(p)
        except importlib_metadata.PackageNotFoundError:
            v[p] = None
    return v


def checkpoint_info(clf):
    from tabpfn.model_loading import ModelSource
    path = Path(str(clf.model_path))
    meta = path.parent / ".cache" / "huggingface" / "download" / f"{path.name}.metadata"
    revision = None
    if meta.exists():
        first = meta.read_text().splitlines()[:1]
        revision = first[0].strip() if first else None
    return {
        "model_version": "v2",
        "checkpoint_file": path.name,
        "checkpoint_sha256": sha256(path),
        "huggingface_repo": ModelSource.get_classifier_v2().repo_id,
        "huggingface_revision": revision,
        "weights_licence": "Prior Labs License 1.1 (Apache 2.0 plus an attribution "
                           "clause; Hugging Face license_name priorlabs-1-1); ungated "
                           "download, no account or licence acceptance needed",
        "constructor": "TabPFNClassifier.create_default_for_version(ModelVersion.V2, "
                       f"device='cpu', random_state={SEED}, n_preprocessing_jobs=1)",
        "n_estimators": int(clf.n_estimators),
        "softmax_temperature": float(clf.softmax_temperature),
        "device": "cpu",
        "why_not_default_v3": "tabpfn 8.x defaults to the TabPFN-3 weights, whose "
                              "download requires logging in to a Prior Labs account "
                              "and accepting the tabpfn-3-license-v1.0 licence (v2.5 "
                              "and v2.6 likewise). That was not done. TabPFN v2 is the "
                              "model described in Hollmann et al., Nature 2025.",
    }


def annotation():
    if not ANNOTATION_FILE.exists():
        return {}
    try:
        a = pd.read_csv(ANNOTATION_FILE, sep="\t", dtype=str,
                        usecols=lambda c: c in ("gene_id", "gene_symbol", "protein"))
    except ValueError:
        return {}
    if "gene_id" not in a.columns:
        return {}
    a = a.drop_duplicates("gene_id").set_index("gene_id")
    return {g: (r.get("gene_symbol"), r.get("protein")) for g, r in a.iterrows()}


def feature_importance(analysis, s, vst, ann):
    """Descriptive L2 logistic-regression coefficients for the binary task."""
    from sklearn.linear_model import LogisticRegression
    task = "binary_microgravity_vs_ground"
    y, _ = labels_for(task, s)
    rows = []
    for k in K_VALUES:
        genes, X = full_data_features(vst, k)
        coef = LogisticRegression(random_state=SEED, **LOGREG_PARAMS).fit(X, y).coef_[0]
        fold_coefs = []
        for f in fold_features(vst, k):
            c = LogisticRegression(random_state=SEED, **LOGREG_PARAMS).fit(
                f["X_train"], y[f["train"]]).coef_[0]
            fold_coefs.append(dict(zip(f["genes"], c)))
        order = np.argsort(-np.abs(coef), kind="stable")[:N_TOP_IMPORTANCE]
        for rank, j in enumerate(order, 1):
            g = genes[j]
            sel = [fc[g] for fc in fold_coefs if g in fc]
            sym, prot = ann.get(g, (None, None))
            rows.append({
                "analysis": analysis, "task": task, "k": k, "rank": rank,
                "gene_id": g, "gene_symbol": sym, "product": prot,
                "std_coef_toward_microgravity": coef[j],
                "abs_std_coef": abs(coef[j]),
                "loocv_folds_selected": len(sel),
                "loocv_folds_same_sign": int(sum(np.sign(c) == np.sign(coef[j]) for c in sel)),
                "n_loocv_folds": len(fold_coefs),
                "status": "descriptive_not_validated",
            })
    return rows


def _clean(v):
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return None if np.isnan(v) else round(float(v), 6)
    return v


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="OSD-528 LOOCV classification with TabPFN v2")
    ap.add_argument("--n-jobs", type=int, default=max(1, min(12, (os.cpu_count() or 2) - 2)))
    ap.add_argument("--max-exact", type=int, default=20000,
                    help="enumerate all relabellings up to this many, else Monte-Carlo")
    ap.add_argument("--n-mc", type=int, default=2000, help="Monte-Carlo relabellings")
    ap.add_argument("--no-permutation", action="store_true",
                    help="skip the permutation test (quick checks only; p-values blank)")
    args = ap.parse_args()

    from joblib import Parallel, delayed

    _quiet()
    t0 = time.time()
    samples = load_samples()
    rng = np.random.default_rng(SEED)

    # Load (and, the first time, download) the TabPFN v2 weights once in this
    # process before the worker pool starts.
    warm = make_model("TabPFN_v2")
    warm.fit(np.random.default_rng(SEED).normal(size=(6, 4)), np.array([0, 0, 0, 1, 1, 1]))
    ckpt = checkpoint_info(warm)
    print(f"tabpfn {importlib_metadata.version('tabpfn')} | {ckpt['checkpoint_file']} "
          f"| sha256 {ckpt['checkpoint_sha256'][:12]}")

    pred_rows, metric_rows, null_rows, imp_rows, inputs = [], [], [], [], {}
    ann = annotation()

    for analysis in VST_FILES:
        s, vst = load_analysis(analysis, samples)
        inputs[analysis] = {
            "vst_file": str(VST_FILES[analysis].relative_to(PROJECT_DIR)),
            "vst_sha256": sha256(VST_FILES[analysis]),
            "samples": list(s["sample_id"]),
            "n_genes_finite": int(vst.shape[0]),
        }
        print(f"\n[{analysis}] {len(s)} samples x {vst.shape[0]} genes")
        imp_rows += feature_importance(analysis, s, vst, ann)

        for k in K_VALUES:
            folds = fold_features(vst, k)
            for task, classes in TASKS.items():
                y, names = labels_for(task, s)
                n_classes = len(classes)
                if args.no_permutation:
                    labs, method, total = [tuple(y)], "not_run", n_distinct_labelings(y)
                else:
                    labs, method, total = permutation_plan(y, args.max_exact, args.n_mc, rng)
                obs_idx = labs.index(tuple(y))
                size = max(1, math.ceil(len(labs) / (4 * args.n_jobs)))
                batches = [labs[i:i + size] for i in range(0, len(labs), size)]
                t1 = time.time()
                res = Parallel(n_jobs=args.n_jobs, backend="loky")(
                    delayed(run_batch)(folds, b, n_classes) for b in batches)
                res = [r for b in res for r in b]
                print(f"  k={k:<3} {task:<31} {len(labs):>5} labellings ({method}) "
                      f"{time.time() - t1:6.0f} s", flush=True)

                obs_part = partition_of(y)
                same_part = sum(partition_of(l) == obs_part for l in labs)
                for m in MODELS:
                    stats = np.array([r[m][0] for r in res])  # (n_labs, 3)
                    acc, bacc, auc = stats[obs_idx]
                    p = p_values(stats, obs_idx, method)
                    ran = method != "not_run"
                    metric_rows.append({
                        "analysis": analysis, "task": task, "k": k, "model": m,
                        "n_samples": len(y),
                        "class_counts": ";".join(f"{c}={int(np.sum(y == i))}"
                                                 for i, c in enumerate(classes)),
                        "n_correct": int(round(acc * len(y))),
                        "accuracy": acc, "balanced_accuracy": bacc, "auroc": auc,
                        "permutation_method": method,
                        "n_relabellings_total": total,
                        "n_relabellings_evaluated": len(labs) if ran else 0,
                        "p_accuracy": p[0], "p_balanced_accuracy": p[1], "p_auroc": p[2],
                        "null_mean_accuracy": float(stats[:, 0].mean()) if ran else float("nan"),
                        "null_mean_balanced_accuracy":
                            float(stats[:, 1].mean()) if ran else float("nan"),
                        "min_attainable_p":
                            same_part / len(labs) if method == "exact" else float("nan"),
                    })
                    if ran:
                        for li, l in enumerate(labs):
                            null_rows.append({
                                "analysis": analysis, "task": task, "k": k, "model": m,
                                "relabelling_index": li, "is_observed": li == obs_idx,
                                "labels": "".join(str(c) for c in l),
                                "accuracy": stats[li, 0],
                                "balanced_accuracy": stats[li, 1],
                                "auroc": stats[li, 2],
                            })
                    prob = res[obs_idx][m][1]
                    for i, sid in enumerate(s["sample_id"]):
                        predicted = classes[int(prob[i].argmax())]
                        row = {
                            "analysis": analysis, "task": task, "k": k, "model": m,
                            "sample_id": sid, "group": s.loc[i, "group"],
                            "true_label": names[i], "predicted_label": predicted,
                            "correct": predicted == names[i], "n_train": len(y) - 1,
                        }
                        for task_classes in TASKS.values():
                            for c in task_classes:
                                row[f"prob_{c}"] = (prob[i, classes.index(c)]
                                                    if c in classes else float("nan"))
                        pred_rows.append(row)

    fmt = "%.6g"
    pd.DataFrame(pred_rows).to_csv(DATA_PROCESSED / "ml_loocv_predictions.tsv", sep="\t",
                                   index=False, float_format=fmt)
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(DATA_PROCESSED / "ml_metrics.tsv", sep="\t", index=False, float_format=fmt)
    pd.DataFrame(imp_rows).to_csv(DATA_PROCESSED / "ml_feature_importance.tsv", sep="\t",
                                  index=False, float_format=fmt)
    outputs = ["data/processed/ml_loocv_predictions.tsv", "data/processed/ml_metrics.tsv",
               "data/processed/ml_feature_importance.tsv"]
    if null_rows:
        pd.DataFrame(null_rows).to_csv(DATA_PROCESSED / "ml_permutation_null.tsv", sep="\t",
                                       index=False, float_format=fmt)
        outputs.append("data/processed/ml_permutation_null.tsv")
    outputs.append("data/processed/ml_summary.json")

    summary = {
        "dataset": "NASA OSDR OSD-528, Mycobacterium marinum under simulated microgravity",
        "script": "analysis/04_tabpfn_tabular_foundation_ai.py",
        "environment": versions(),
        "requirements_file": "requirements-ml.txt",
        "tabpfn_checkpoint": ckpt,
        "credentials_and_telemetry": "No Hugging Face or Prior Labs token used "
            "(HF_HUB_DISABLE_IMPLICIT_TOKEN=1; HF_TOKEN and TABPFN_TOKEN unset); "
            "licence browser flow disabled (TABPFN_NO_BROWSER=1); TabPFN telemetry "
            "disabled (TABPFN_DISABLE_TELEMETRY=1).",
        "inputs": {"samples_file": "config/samples.tsv",
                   "samples_sha256": sha256(SAMPLES_FILE), **inputs},
        "design": {
            "analyses": {
                "primary": "the samples with primary == yes in config/samples.tsv",
                "sensitivity_9": "all samples in config/samples.tsv (primary plus the "
                                 "inferred RFPRPM4 library)",
            },
            "tasks": {t: list(c) for t, c in TASKS.items()},
            "positive_class_for_auroc": POSITIVE_CLASS,
            "evaluation": "leave-one-out cross-validation over samples",
            "features": "DESeq2 VST expression; in each training fold the k genes "
                        "with the highest variance across the training samples only "
                        "(ties broken by gene id), standardised with training mean "
                        "and SD; the held-out sample is transformed with those "
                        "statistics",
            "k_values": list(K_VALUES),
            "models": {
                "TabPFN_v2": ckpt["constructor"],
                "LogReg_L2": f"sklearn LogisticRegression({LOGREG_PARAMS}, "
                             f"random_state={SEED})",
                "RandomForest_500": f"sklearn RandomForestClassifier({RF_PARAMS})",
                "Majority_baseline": "sklearn DummyClassifier(strategy='prior'): "
                                     "predicts the training-fold majority class",
            },
            "prediction_rule": "predicted label = argmax of held-out class "
                               "probabilities; ties go to the first class in the task's "
                               "class list",
            "permutation_test": "the full LOOCV pipeline (feature selection, scaling, "
                                "fitting) rerun under every distinct relabelling that "
                                "keeps the class sizes (exact); p = fraction of "
                                "relabellings, observed included, scoring >= the "
                                "observed value. Monte-Carlo with seed "
                                f"{SEED} above {args.max_exact} relabellings.",
            "permutation_test_run": not args.no_permutation,
            "min_attainable_p": "fraction of relabellings that induce the same "
                                "partition of samples as the observed labels; a model "
                                "that ignores class names cannot get a smaller p",
            "seed": SEED,
        },
        "metrics": [{k: _clean(v) for k, v in r.items()} for r in metric_rows],
        "feature_importance": {
            "file": "data/processed/ml_feature_importance.tsv",
            "method": "L2 logistic regression refitted on all samples of each analysis "
                      "(binary task), top-k variance genes computed on all samples, "
                      "standardised; coefficients > 0 point toward simulated "
                      f"microgravity. Top {N_TOP_IMPORTANCE} genes by |coefficient| per "
                      "analysis and k, with the number of LOOCV training folds in which "
                      "each gene was selected and kept the same coefficient sign.",
            "status": "descriptive only, not validated. With n = 8 or 9 and k >> n the "
                      "coefficients are heavily shrunk and spread over correlated "
                      "genes, so the ranking is unstable and is not a biomarker list. "
                      "No importances are reported for TabPFN or the random forest.",
        },
        "osd90_cross_study_transfer": "not performed",
        "caveats": [
            "n is 8 (primary) or 9 (sensitivity) samples; each LOOCV fold trains on 7 "
            "or 8 samples. Accuracy moves in steps of 1/8 or 1/9.",
            "The permutation p-value cannot go below min_attainable_p, so even a "
            "perfect LOOCV score is weak evidence on its own.",
            "LOOCV with tiny classes penalises classifiers that lean on class "
            "frequencies: holding out a sample shrinks its class in the training "
            "fold. The majority baseline shows the resulting no-information level "
            "(it can sit below chance), and the permutation null includes this effect.",
            "AUROC is computed from held-out probabilities pooled across LOOCV folds, "
            "which come from different fitted models; treat it as descriptive.",
            "The three-class task has 2 (primary) or 3 (sensitivity) RPM samples.",
            "The VST was fitted by DESeq2 with blind = TRUE on all samples of each "
            "analysis, the held-out one included. It does not use labels, so it "
            "cannot leak class information, but it is not refitted per fold.",
            "Top-k genes are chosen by unsupervised variance, not by association with "
            "the labels; k = 50 and k = 200 are both reported and neither was tuned.",
            "The sensitivity analysis adds RFPRPM4, a library inferred from the "
            "deposited RFPRPM4 files (which begin with the complete RFPRPM41 files); "
            "it may not be an independent biological replicate.",
            "TabPFN v2 is used because the package default (TabPFN-3) requires a Prior "
            "Labs account and licence acceptance; later TabPFN versions may give "
            "different results.",
            "OSD-90 cross-study transfer was not performed; this analysis supports no "
            "claim about transfer to other studies.",
        ],
        "outputs": outputs,
    }
    with open(DATA_PROCESSED / "ml_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, allow_nan=False)
        fh.write("\n")

    cols = ["analysis", "task", "k", "model", "accuracy", "balanced_accuracy", "auroc",
            "p_accuracy", "p_balanced_accuracy", "p_auroc", "min_attainable_p"]
    with pd.option_context("display.width", 220, "display.max_rows", 200):
        print("\n" + metrics[cols].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\nDone in {time.time() - t0:.0f} s. Outputs: data/processed/ml_*")


if __name__ == "__main__":
    sys.exit(main())
