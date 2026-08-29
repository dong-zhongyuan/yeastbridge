#!/usr/bin/env python
"""Norman program-retrieval benchmark (registered; retrieval_registration.md).

Direction-retrieval exam for foundation-model gene embeddings on unseen
perturbation genes. Run under yeastbridge_vs_models python.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from sklearn.cluster import KMeans
from sklearn.linear_model import Ridge

ROOT = Path("/public/home/mengxl/dzy/yeastbridge_re")
SRC = ROOT / "product/norman_assets"
FEATS = ROOT / "feasibility/norman_foundation/features"
OUT = ROOT / "feasibility/norman_foundation/retrieval"
SEED = 20260821
N_PROGRAMS = 6
ALPHA_GRID = np.logspace(-3, 6, 10)

oof = np.load(SRC / "raw_oof_full.npz", allow_pickle=True)
conditions = [str(c) for c in oof["condition_ids"]]
y = oof["y_true"].astype(np.float64)
genes_of = [tuple(g for g in c.split("+") if g and g != "ctrl") for c in conditions]
single = np.array([len(g) == 1 for g in genes_of])

km = KMeans(n_clusters=N_PROGRAMS, random_state=SEED, n_init=10)
programs = km.fit_predict(y)
print("program sizes:", np.bincount(programs).tolist(), flush=True)

single_genes = sorted({genes_of[i][0] for i in np.where(single)[0]})
rng = np.random.default_rng(SEED)
fold_of = {g: int(rng.integers(0, 5)) for g in single_genes}

models = {}
for m in ("scgpt", "geneformer", "scfoundation"):
    d = np.load(FEATS / f"{m}.npz", allow_pickle=True)
    vocab = {str(s): i for i, s in enumerate(d["gene_symbols"])}
    models[m] = (vocab, d["embeddings"].astype(np.float64))


def embed_condition(gene_tuple, kind):
    if kind == "zero":
        return None
    vocab, emb = models[kind]
    idx = [vocab.get(g) for g in gene_tuple]
    if any(i is None for i in idx):
        return None
    return emb[idx].mean(axis=0)


def build_X(rows_idx, kind):
    if kind == "zero":
        return np.zeros((len(rows_idx), 1))
    vectors = [embed_condition(genes_of[i], kind) for i in rows_idx]
    return np.array(vectors)


def ridge_predict(kind, train_idx, test_idx):
    Xtr = build_X(list(train_idx), kind)
    Xte = build_X(list(test_idx), kind)
    if kind == "zero":
        return np.tile(y[train_idx].mean(axis=0), (len(test_idx), 1))
    train_ok = ~np.isnan(Xtr).any(axis=1)
    if train_ok.sum() < 5:
        return None
    best_alpha, best_err = None, np.inf
    Xtr_ok, y_tr_ok = Xtr[train_ok], y[np.array(train_idx)[train_ok]]
    n = Xtr_ok.shape[0]
    for alpha in ALPHA_GRID:
        errs = []
        for f in range(4):
            va = np.arange(n) % 4 == f
            model = Ridge(alpha=alpha).fit(Xtr_ok[~va], y_tr_ok[~va])
            errs.append(((model.predict(Xtr_ok[va]) - y_tr_ok[va]) ** 2).mean())
        if np.mean(errs) < best_err:
            best_err, best_alpha = np.mean(errs), alpha
    model = Ridge(alpha=best_alpha).fit(Xtr[train_ok], y[np.array(train_idx)[train_ok]])
    preds = np.full((len(test_idx), y.shape[1]), np.nan)
    mask = ~np.isnan(Xte).any(axis=1)
    if mask.any():
        preds[mask] = model.predict(Xte[mask])
    return preds


def cosine_matrix(a, b):
    a = a / (np.linalg.norm(a, axis=1, keepdims=True) + 1e-12)
    b = b / (np.linalg.norm(b, axis=1, keepdims=True) + 1e-12)
    return a @ b.T


anchors_cache = {}


def program_anchors(exclude_gene_conditions, exclude_gene):
    key = (exclude_gene,)
    if key in anchors_cache:
        return anchors_cache[key]
    anchors = np.zeros((N_PROGRAMS, y.shape[1]))
    for p in range(N_PROGRAMS):
        rows = [
            i
            for i in range(len(conditions))
            if programs[i] == p and i not in exclude_gene_conditions and genes_of[i][0] != exclude_gene
        ]
        anchors[p] = y[rows].mean(axis=0) if rows else 0.0
    anchors_cache[key] = anchors
    return anchors


results = {}
per_gene_truth = {}
for kind in ("zero", "scgpt", "geneformer", "scfoundation"):
    correct, total, per_gene = 0, 0, {}
    for fold in range(5):
        test_genes_f = [g for g in single_genes if fold_of[g] == fold]
        test_idx = [i for i in np.where(single)[0] if genes_of[i][0] in test_genes_f]
        if not test_idx:
            continue
        train_idx = [i for i in range(len(conditions)) if genes_of[i][0] not in test_genes_f]
        preds = ridge_predict(kind, train_idx, test_idx)
        if preds is None:
            continue
        for pos, i in enumerate(test_idx):
            gene = genes_of[i][0]
            anchors = program_anchors(set(test_idx), gene)
            row = preds[pos]
            if np.isnan(row).any():
                continue
            sims = cosine_matrix(row[None, :], anchors)[0]
            assigned = int(np.argmax(sims))
            if kind == "zero":
                assigned = int(rng.integers(0, N_PROGRAMS))
            hit = int(assigned == programs[i])
            correct += hit
            total += 1
            per_gene.setdefault(gene, []).append(hit)
    accuracy = round(correct / max(total, 1), 4)
    results[kind] = {
        "accuracy": accuracy,
        "n_evaluated": total,
        "n_genes": len(per_gene),
        "mean_gene_accuracy": round(float(np.mean([np.mean(v) for v in per_gene.values()])), 4) if per_gene else None,
    }
    if kind != "zero":
        per_gene_truth[kind] = {g: float(np.mean(v)) for g, v in per_gene.items()}
    print(kind, results[kind], flush=True)

zero_acc = results["zero"]["accuracy"]
perm_p = {}
for kind in ("scgpt", "geneformer", "scfoundation"):
    acc = results[kind]["accuracy"]
    n = results[kind]["n_evaluated"]
    hits = np.random.default_rng(SEED).binomial(n, zero_acc, 10000)
    perm_p[kind] = round(float((hits >= acc * n).mean()), 5)

best = max(("scgpt", "geneformer", "scfoundation"), key=lambda k: results[k]["accuracy"])
gene_scores = per_gene_truth[best]
genes_list = sorted(gene_scores)
values = np.array([gene_scores[g] for g in genes_list])
boot = np.random.default_rng(SEED)
zeros = np.random.default_rng(SEED).binomial(1, zero_acc, size=(10000, len(values)))
boot_means = zeros.mean(axis=1)
ci = [round(float(np.percentile(boot_means, 2.5)), 4), round(float(np.percentile(boot_means, 97.5)), 4)]
gate_pass = bool(results[best]["accuracy"] > zero_acc and perm_p[best] < 0.05)

# registered consistency check: three-model fusion (cosine-average and vote)
fusion_results = {}
for mode in ("avg", "vote"):
    correct, total = 0, 0
    rng_f = np.random.default_rng(SEED)
    for fold in range(5):
        test_genes_f = [g for g in single_genes if fold_of[g] == fold]
        test_idx = [i for i in np.where(single)[0] if genes_of[i][0] in test_genes_f]
        if not test_idx:
            continue
        train_idx = [i for i in range(len(conditions)) if genes_of[i][0] not in test_genes_f]
        preds_per_model = {}
        for kind in ("scgpt", "geneformer", "scfoundation"):
            preds_per_model[kind] = ridge_predict(kind, train_idx, test_idx)
        for pos, i in enumerate(test_idx):
            gene = genes_of[i][0]
            anchors = program_anchors(set(test_idx), gene)
            sims_per_model = {}
            for kind, preds in preds_per_model.items():
                if preds is None or np.isnan(preds[pos]).any():
                    continue
                sims_per_model[kind] = cosine_matrix(preds[pos][None, :], anchors)[0]
            if not sims_per_model:
                continue
            keys = sorted(sims_per_model)
            if mode == "avg":
                assigned = int(np.argmax(np.mean([sims_per_model[k] for k in keys], axis=0)))
            else:
                votes = np.bincount([int(np.argmax(sims_per_model[k])) for k in keys], minlength=N_PROGRAMS)
                top = np.flatnonzero(votes == votes.max())
                assigned = int(rng_f.choice(top))
            correct += int(assigned == programs[i])
            total += 1
    fusion_results[f"fusion_{mode}"] = {"accuracy": round(correct / max(total, 1), 4), "n_evaluated": total}
    print(f"fusion_{mode}", fusion_results[f"fusion_{mode}"], flush=True)

report = {
    "fusion_consistency": fusion_results,
    "registration": "reports/human_foundation_norman/retrieval/retrieval_registration.md",
    "seed": SEED,
    "n_programs": N_PROGRAMS,
    "program_sizes": np.bincount(programs).tolist(),
    "results": results,
    "permutation_p_vs_zero": perm_p,
    "best_model": best,
    "zero_accuracy": zero_acc,
    "chance": round(1 / N_PROGRAMS, 4),
    "gene_level_bootstrap_ci_of_best": ci,
    "gate_pass": gate_pass,
    "claim": "pretrained gene representations carry program-direction information for unseen perturbations"
    if gate_pass
    else "NO_GO: embeddings do not beat the zero-info control on direction retrieval",
    "claim_boundary": "one dataset's program structure; not a response-regression claim",
}
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "retrieval_result.json").write_text(json.dumps(report, indent=2, sort_keys=True))
print(json.dumps(report, indent=2, sort_keys=True))
