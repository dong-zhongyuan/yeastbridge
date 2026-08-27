#!/usr/bin/env python3
"""Product step 3 (endpoint v2): execute the transferred yeast tasks against
the public HIP/HOP chemogenomic screen.

v2 statistic (SENSITIVITY.md follow-up): full-profile Spearman correlation
between each target's COMPLETE step-2 task ranking (rank vector over
measured strains) and each compound condition's |z| sensitivity profile.
Cutoff-free - no arbitrary gene-set constant; the v1 K-sweep showed the
signal is diffuse across the whole ranking, which is exactly the shape a
global rank correlation measures. Null = strain-label permutation
(1,000 draws, seed 42, same permutation applied to every target's rank
vector so target margins are preserved); empirical p -> best dose per
InChIKey -> BH-FDR across all (target, InChIKey) pairs.
"""
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / "configs/product_execute.json").read_text())
RESULTS = ROOT / CFG["results_dir"]
TASK_DIR = ROOT / CFG["task_dir"]
N_PERM = CFG["n_perm"]
SEED = CFG["seed"]


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def rank_cols(A):
    """Column-wise ranks (average ties) for a 2-D float array."""
    n = A.shape[0]
    out = np.empty(A.shape, dtype=np.float64)
    for j in range(A.shape[1]):
        col = A[:, j]
        o = np.argsort(col, kind="stable")
        r = np.empty(n)
        r[o] = np.arange(1, n + 1)
        s = col[o]
        i = 0
        while i < n:
            k = i
            while k + 1 < n and s[k + 1] == s[i]:
                k += 1
            if k > i:
                r[o[i:k + 1]] = (i + 1 + k + 1) / 2.0
            i = k + 1
        out[:, j] = r
    return out


def standardize_cols(A):
    A = A - A.mean(axis=0, keepdims=True)
    sd = A.std(axis=0, keepdims=True)
    return A / np.maximum(sd, 1e-12)


def main() -> None:
    z = np.load(CFG["response_npz"], allow_pickle=True)
    orfs = z["strain_orfs"].astype(str)
    inchi = z["compound_inchikeys"].astype(str)
    doses = z["doses"].astype(str)
    is_vehicle = z["is_vehicle"].astype(bool)
    M = z["z_score"]
    strain_pos = {g: i for i, g in enumerate(orfs)}

    keep = (~is_vehicle) & np.array([len(k) > 0 for k in inchi])
    cols = np.nonzero(keep)[0]
    n_str, n_cond = len(orfs), len(cols)
    print(f"strains={n_str} non-vehicle conditions={n_cond}", flush=True)

    # compound condition matrix, |z|, rank-transformed and standardized
    Y = standardize_cols(rank_cols(np.abs(M[:, cols]).astype(np.float64)))  # (n_str, n_cond)

    # target task rank vectors over the same strains (better rank = lower value)
    task_files = sorted(TASK_DIR.glob("yeast_task_*.tsv"))
    targets, X, coverage = [], [], []
    for tf in task_files:
        target = tf.stem.replace("yeast_task_", "")
        genes = [r["yeast_gene"] for r in csv.DictReader(open(tf), delimiter="\t")]
        pos = [strain_pos[g] for g in genes if g in strain_pos]
        if len(pos) < 100:
            continue
        v = np.full(n_str, np.nan)
        v[np.asarray(pos)] = np.arange(1, len(pos) + 1)
        v = np.where(np.isnan(v), len(pos) + (n_str - len(pos)) / 2.0, v)
        targets.append(target)
        X.append(v)
        coverage.append(len(pos))
    X = np.stack(X, axis=0)
    X = standardize_cols(rank_cols(X.T).T)
    print(f"targets={len(targets)} task-strain coverage median={int(np.median(coverage))}", flush=True)

    rng = np.random.default_rng(SEED)
    rho = (X @ Y) / n_str  # (n_targets, n_cond) - Pearson on standardized ranks = Spearman
    null_ge = np.zeros((len(targets), n_cond))
    for b in range(N_PERM):
        perm = rng.permutation(n_str)
        null_ge += ((X[:, perm] @ Y) / n_str) >= rho
        if (b + 1) % 100 == 0:
            print(f"perm {b+1}/{N_PERM}", flush=True)
    emp_p = null_ge / N_PERM

    comp_tab = pd.read_csv(CFG["compound_table"], sep="\t", compression="gzip").fillna("")
    ik2smiles = dict(zip(comp_tab["inchikey"], comp_tab["smiles"]))
    col_inchi = inchi[cols]
    col_dose = doses[cols]
    rows = []
    for ti, target in enumerate(targets):
        for j in range(n_cond):
            rows.append({"target_id": target, "inchikey": col_inchi[j], "dose": col_dose[j],
                         "spearman_rho": float(rho[ti, j]), "emp_p": float(emp_p[ti, j])})
    D = pd.DataFrame(rows)
    best = (D.sort_values("spearman_rho", ascending=False)
              .groupby(["target_id", "inchikey"], as_index=False).first())
    best["q"] = bh_fdr(best["emp_p"].to_numpy())
    best["smiles"] = best["inchikey"].map(ik2smiles)
    best = best.sort_values(["target_id", "spearman_rho"], ascending=[True, False])

    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "exec_matrix_v2.tsv").write_text(best.to_csv(sep="\t", index=False))
    top_dir = RESULTS / "per_target_top_v2"
    top_dir.mkdir(exist_ok=True)
    for target, grp in best.groupby("target_id"):
        grp.head(CFG["top_report_per_target"]).to_csv(top_dir / f"{target}.tsv", sep="\t", index=False)
    sig = best[best["q"] < CFG["fdr_alpha"]]
    summary = {
        "endpoint": "v2 full-profile Spearman (cutoff-free)",
        "n_targets": int(best["target_id"].nunique()),
        "n_compounds": int(best["inchikey"].nunique()),
        "n_pairs": int(len(best)),
        "n_pairs_q_below_alpha": int(len(sig)),
        "max_rho": float(best["spearman_rho"].max()),
        "n_perm": N_PERM, "seed": SEED,
    }
    (RESULTS / "execute_summary_v2.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)
    print("[done] ->", RESULTS, flush=True)


if __name__ == "__main__":
    main()
