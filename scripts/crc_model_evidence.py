#!/usr/bin/env python
"""CRC model-evidence layer (yeastbridge_re product step 1), scFoundation-only.

Copied verbatim from yeastbridge_vs/scripts/crc_model_evidence.py with the
minimal edits required by the reproduction project: paths point into
yeastbridge_re product directories, and the engine is the single registered
best retrieval reader (scFoundation; see feasibility/norman_foundation/
retrieval/retrieval_result.json) instead of the three-model fusion. Scores
every target-universe gene: does perturbing it push cells along the
registered Becker disease->healthy recovery direction? Run under
yeastbridge_vs_models python.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import Ridge

ROOT = Path("/public/home/mengxl/dzy/yeastbridge_re")
SRC = ROOT / "product/norman_assets"
SCAN = ROOT / "product/target_scan/inputs"
OUT = ROOT / "product/target_scan/model_evidence_scfoundation_v1"
SEED = 20260821
ALPHA_GRID = np.logspace(-3, 6, 10)

oof = np.load(SRC / "raw_oof_full.npz", allow_pickle=True)
conditions = [str(c) for c in oof["condition_ids"]]
y = oof["y_true"].astype(np.float64)
response_symbols = [str(s) for s in oof["response_gene_symbols"]]


def condition_embedding(gene_tuple, vocab, emb):
    idx = [vocab.get(g) for g in gene_tuple]
    if any(i is None for i in idx):
        return None
    return emb[idx].mean(axis=0)


def build_model(name, vocab, emb):
    rows, keep = [], []
    for i, c in enumerate(conditions):
        v = condition_embedding(tuple(g for g in c.split("+") if g and g != "ctrl"), vocab, emb)
        if v is not None:
            rows.append(v)
            keep.append(i)
    X = np.array(rows)
    Y = y[keep]
    best_alpha, best_err = None, np.inf
    n = len(X)
    for alpha in ALPHA_GRID:
        errs = []
        for f in range(4):
            va = np.arange(n) % 4 == f
            model = Ridge(alpha=alpha).fit(X[~va], Y[~va])
            errs.append(((model.predict(X[va]) - Y[va]) ** 2).mean())
        if np.mean(errs) < best_err:
            best_err, best_alpha = np.mean(errs), alpha
    model = Ridge(alpha=best_alpha).fit(X, Y)
    return {"name": name, "vocab": vocab, "emb": emb, "model": model, "alpha": best_alpha, "n_train": n}


# --- single registered reader: scFoundation full gene table, loaded via the
# bridge_reuse module copied verbatim from the old project; only the
# scfoundation entry is consumed (registered best retrieval reader). ---
import sys as _sys

_sys.path.insert(0, str(ROOT / "src"))
from yeastbridge_re.bridge_reuse import load_three_model_embeddings  # noqa: E402

_tables = load_three_model_embeddings()
_vocab, _emb = _tables["scfoundation"]
models = [build_model("scfoundation", _vocab, _emb)]
print("scfoundation loaded via bridge_reuse:", _emb.shape, flush=True)
# --- Becker recovery signature (step-1 frozen asset) ---
sig = pd.read_csv(SCAN / "state_signature.tsv", sep="\t")
sig_by_symbol = dict(zip(sig["gene"].astype(str), sig["state_effect_disease_minus_desired"].astype(float)))
recovery = {g: -v for g, v in sig_by_symbol.items()}  # healthy minus disease
common = [g for g in response_symbols if g in recovery]
print(f"signature/response overlap genes: {len(common)}", flush=True)
common_idx = np.array([response_symbols.index(g) for g in common])
recovery_vec = np.array([recovery[g] for g in common])

# --- score every universe gene ---
def score_gene(gene):
    per_model_pred = []
    used = []
    for m in models:
        v = condition_embedding((gene,), m["vocab"], m["emb"])
        if v is None:
            continue
        per_model_pred.append(m["model"].predict(v[None, :])[0])
        used.append(m["name"])
    if not per_model_pred:
        return {"target_id": gene, "model_direction_score": None, "models_used": ""}
    pred = np.mean(per_model_pred, axis=0)[common_idx]
    if pred.std() < 1e-12:
        score = 0.0
    else:
        score = float(np.corrcoef(pred, recovery_vec)[0, 1])
    return {"target_id": gene, "model_direction_score": round(score, 4), "models_used": ";".join(used)}


targets = []
with (SCAN / "target_universe.tsv").open() as h:
    reader = csv.DictReader(h, delimiter="\t")
    for row in reader:
        targets.append(row["target_id"])

rows = [score_gene(gene) for gene in targets]

scores = pd.DataFrame(rows)
OUT.mkdir(parents=True, exist_ok=True)
scores.to_csv(OUT / "model_direction_scores.tsv", sep="\t", index=False)

# --- dual-track signature: model direction scores for every signature gene ---
sig_rows = []
for gene, stat_effect in sig_by_symbol.items():
    entry = score_gene(gene)
    entry["gene"] = gene
    entry["statistical_state_effect"] = stat_effect
    entry["direction_agrees"] = (
        None if entry["model_direction_score"] is None
        else bool(np.sign(entry["model_direction_score"]) == np.sign(stat_effect))
    )
    sig_rows.append(entry)
sig_scores = pd.DataFrame(sig_rows)
sig_scores.to_csv(OUT / "signature_model_scores.tsv", sep="\t", index=False)
n_agree = int(sig_scores["direction_agrees"].fillna(False).sum())
n_scored = int(sig_scores["model_direction_score"].notna().sum())
print(f"signature genes scored by scfoundation: {n_scored}/{len(sig_rows)}, direction agrees: {n_agree}", flush=True)

# --- dual-evidence ranking of the registered 47 ---
base = pd.read_csv(SCAN / "candidate_baseline.tsv", sep="\t")
base["model_direction_score"] = base["target_id"].map(
    dict(zip(scores["target_id"], scores["model_direction_score"]))
)
eligible = base[base["baseline_eligible"].astype(str).str.lower() == "true"].copy()
ranked = eligible.sort_values("model_direction_score", ascending=False)[
    ["target_id", "target_family", "intended_direction", "model_direction_score",
     "state_signal_to_noise", "directional_xatlas_reversal", "xatlas_batch_direction_fraction"]
]
ranked.to_csv(OUT / "dual_evidence_ranking.tsv", sep="\t", index=False)

report = {
    "engine": "registered best retrieval reader: frozen scFoundation embeddings + full-task ridge adapter (single model)",
    "models": [{"name": m["name"], "alpha": m["alpha"], "n_train_conditions": m["n_train"]} for m in models],
    "universe_scored": int(scores["model_direction_score"].notna().sum()),
    "signature_overlap_genes": len(common),
    "top15_dual_evidence": ranked.head(15).to_dict(orient="records"),
    "model_score_agreement_with_statistical_direction": {
        "median_score_eligible": round(float(eligible["model_direction_score"].median()), 4),
        "median_score_universe": round(float(scores["model_direction_score"].median()), 4),
    },
    "claim_boundary": "model direction scores are retrieval-based priors on unseen genes, not perturbation evidence; they rank, they do not validate",
}
(OUT / "model_evidence_report.json").write_text(json.dumps(report, indent=2, sort_keys=True, default=float))
print(json.dumps({k: report[k] for k in ("universe_scored", "signature_overlap_genes")}, indent=1))
print(ranked.head(15).to_string(index=False), flush=True)
