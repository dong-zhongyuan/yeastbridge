#!/usr/bin/env python3
"""Product step 3: execute the transferred yeast tasks against the public
HIP/HOP chemogenomic screen (design: product/execute_hiphop/DESIGN.md;
config: configs/product_execute.json).

For each target's top-K task genes and each non-vehicle compound
condition, computes the standardized mean |z| of the task strains against
a per-target empirical null (random matched-size gene sets), aggregates to
best-dose per InChIKey, and applies BH-FDR across all pairs.
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


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    # 诊断用途:覆盖配置常数做 K 敏感性扫描;默认取配置值=注册行为
    ap.add_argument("--task-top-k", type=int, default=CFG["task_top_k"])
    ap.add_argument("--results-suffix", default="")
    args = ap.parse_args()
    global TOP_K, RESULTS
    TOP_K = args.task_top_k
    RESULTS = ROOT / (CFG["results_dir"] + args.results_suffix)
    print(f"[config] task_top_k={TOP_K} results={RESULTS}", flush=True)
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

    D = pd.DataFrame(rows)
    best = (D.sort_values("z_exec", ascending=False)
              .groupby(["target_id", "inchikey"], as_index=False).first())
    best["q"] = bh_fdr(best["emp_p"].to_numpy())
    best["smiles"] = best["inchikey"].map(ik2smiles)
    best = best.sort_values(["target_id", "z_exec"], ascending=[True, False])

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
        "n_targets": int(best["target_id"].nunique()),
        "n_compounds": int(best["inchikey"].nunique()),
        "n_pairs": int(len(best)),
        "n_pairs_q_below_alpha": n_sig,
        "task_top_k": TOP_K, "n_perm": N_PERM, "seed": SEED,
        "statistic": "standardized mean |z| of task strains vs per-target empirical null",
    }
    (RESULTS / "execute_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)
    print("[done] ->", RESULTS, flush=True)


if __name__ == "__main__":
    main()
