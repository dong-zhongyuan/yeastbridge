#!/usr/bin/env python3
"""Product step 3: execute the transferred yeast tasks against the public
HIP/HOP chemogenomic screen (design: product/execute_hiphop/DESIGN.md;
config: configs/product_execute.json).

Endpoint v1-topk (registered 2026-08-27, superseded): top-K task-gene
enrichment; the K-sweep (SENSITIVITY.md) showed individual pairs are
K-fragile.

Endpoint v2-spearman (current primary): cutoff-free full-profile Spearman
between each target's COMPLETE step-2 ranking and each compound
condition's |z| profile; strain-label permutation null (n_perm, seed),
best dose per InChIKey, BH-FDR. Selected by config "endpoint" or
--endpoint; v1 kept runnable for reproducibility.
"""
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def _load_cfg():
    import argparse
    ap0 = argparse.ArgumentParser(add_help=False)
    ap0.add_argument("--config", default="configs/product_execute.json")
    known, _ = ap0.parse_known_args()
    return known.config, json.loads((ROOT / known.config).read_text())


CFG_NAME, CFG = _load_cfg()
RESULTS = ROOT / CFG["results_dir"]
TASK_DIR = ROOT / CFG["task_dir"]
TOP_K = CFG["task_top_k"]
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
    """Column-wise ranks with average ties (Spearman 的秩变换)."""
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
    return A / np.maximum(A.std(axis=0, keepdims=True), 1e-12)


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=CFG_NAME)  # 已在模块级消费,此处仅为放行
    # 诊断用途:覆盖配置常数做 K 敏感性扫描;默认取配置值=注册行为
    ap.add_argument("--task-top-k", type=int, default=CFG["task_top_k"])
    ap.add_argument("--results-suffix", default="")
    ap.add_argument("--endpoint", default=CFG.get("endpoint", "v1-topk"),
                    choices=["v1-topk", "v2-spearman"])
    args = ap.parse_args()
    global TOP_K, RESULTS, ENDPOINT
    TOP_K = args.task_top_k
    RESULTS = ROOT / (CFG["results_dir"] + args.results_suffix)
    ENDPOINT = args.endpoint
    print(f"[config] endpoint={ENDPOINT} task_top_k={TOP_K} results={RESULTS}", flush=True)
    z = np.load(CFG["response_npz"], allow_pickle=True)
    orfs = z["strain_orfs"].astype(str)
    inchi = z["compound_inchikeys"].astype(str)
    doses = z["doses"].astype(str)
    is_vehicle = z["is_vehicle"].astype(bool)
    M = z["z_score"]  # (5668, 3850)
    strain_pos = {g: i for i, g in enumerate(orfs)}

    keep = (~is_vehicle) & np.array([len(k) > 0 for k in inchi])
    cols = np.nonzero(keep)[0]
    print(f"strains={len(orfs)} non-vehicle conditions={len(cols)}", flush=True)

    absM = np.abs(M[:, cols])

    comp_tab = pd.read_csv(CFG["compound_table"], sep="\t", compression="gzip").fillna("")
    ik2smiles = dict(zip(comp_tab["inchikey"], comp_tab["smiles"]))
    col_inchi = inchi[cols]
    col_dose = doses[cols]

    task_files = sorted(TASK_DIR.glob("yeast_task_*.tsv"))
    rng = np.random.default_rng(SEED)
    rows = []
    if ENDPOINT == "v2-spearman":
        Y = standardize_cols(rank_cols(absM.astype(np.float64)))
        targets, X, coverage = [], [], []
        for tf in task_files:
            target = tf.stem.replace("yeast_task_", "")
            genes = [r["yeast_gene"] for r in csv.DictReader(open(tf), delimiter="\t")]
            pos = [strain_pos[g] for g in genes if g in strain_pos]
            if len(pos) < 100:
                continue
            v = np.full(len(orfs), len(pos) + (len(orfs) - len(pos)) / 2.0)
            v[np.asarray(pos)] = np.arange(1, len(pos) + 1)
            targets.append(target)
            X.append(v)
            coverage.append(len(pos))
        X = standardize_cols(rank_cols(np.stack(X).T).T)
        print(f"[v2] targets={len(targets)} task-strain coverage median={int(np.median(coverage))}", flush=True)
        rho = (X @ Y) / len(orfs)  # 标准化秩的 Pearson = Spearman
        null_ge = np.zeros(rho.shape)
        for b in range(N_PERM):
            perm = rng.permutation(len(orfs))
            null_ge += ((X[:, perm] @ Y) / len(orfs)) >= rho
            if (b + 1) % 200 == 0:
                print(f"[v2] perm {b+1}/{N_PERM}", flush=True)
        emp_p = null_ge / N_PERM
        for ti, target in enumerate(targets):
            for j in range(len(cols)):
                rows.append({"target_id": target, "inchikey": col_inchi[j], "dose": col_dose[j],
                             "spearman_rho": float(rho[ti, j]), "emp_p": float(emp_p[ti, j])})
    else:
        for tf in task_files:
            target = tf.stem.replace("yeast_task_", "")
            genes = [r["yeast_gene"] for i, r in enumerate(csv.DictReader(open(tf), delimiter="\t")) if i < TOP_K]
            hits = [strain_pos[g] for g in genes if g in strain_pos]
            if not hits:
                continue
            obs = absM[hits, :].mean(axis=0)  # (n_conditions,)
            null_draws = np.zeros(N_PERM)
            draw_sets = [rng.choice(len(orfs), size=len(hits), replace=False) for _ in range(N_PERM)]
            idx_sets = np.stack(draw_sets)
            null_mat = absM[idx_sets, :].mean(axis=1)  # (N_PERM, n_conditions)
            mu = null_mat.mean(axis=0)
            sd = np.maximum(null_mat.std(axis=0), 1e-9)
            z_exec = (obs - mu) / sd
            emp_p = (null_mat >= obs[None, :]).mean(axis=0)
            for j in range(len(cols)):
                rows.append({"target_id": target, "inchikey": col_inchi[j], "dose": col_dose[j],
                             "mean_absz": float(obs[j]), "z_exec": float(z_exec[j]),
                             "emp_p": float(emp_p[j]), "n_task_strains": len(hits)})
            print(f"{target}: task strains={len(hits)}/{TOP_K}", flush=True)

    stat_col = "spearman_rho" if ENDPOINT == "v2-spearman" else "z_exec"
    D = pd.DataFrame(rows)
    best = (D.sort_values(stat_col, ascending=False)
              .groupby(["target_id", "inchikey"], as_index=False).first())
    best["q"] = bh_fdr(best["emp_p"].to_numpy())
    best["smiles"] = best["inchikey"].map(ik2smiles)
    best = best.sort_values(["target_id", stat_col], ascending=[True, False])

    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "exec_matrix.tsv").write_text(
        best.to_csv(sep="\t", index=False))
    top_dir = RESULTS / "per_target_top"
    top_dir.mkdir(exist_ok=True)
    for target, grp in best.groupby("target_id"):
        grp.head(CFG["top_report_per_target"]).to_csv(
            top_dir / f"{target}.tsv", sep="\t", index=False)
    n_sig = int((best["q"] < CFG["fdr_alpha"]).sum())
    summary = {
        "config": CFG_NAME,
        "endpoint": ENDPOINT,
        "n_targets": int(best["target_id"].nunique()),
        "n_compounds": int(best["inchikey"].nunique()),
        "n_pairs": int(len(best)),
        "n_pairs_q_below_alpha": n_sig,
        "task_top_k": TOP_K if ENDPOINT != "v2-spearman" else None,
        "n_perm": N_PERM, "seed": SEED,
        "statistic": ("full-profile Spearman, strain-label permutation null"
                      if ENDPOINT == "v2-spearman" else
                      "standardized mean |z| of task strains vs per-target empirical null"),
    }
    (RESULTS / "execute_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)
    print("[done] ->", RESULTS, flush=True)


if __name__ == "__main__":
    main()
