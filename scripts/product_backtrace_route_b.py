#!/usr/bin/env python3
"""Product step 4: back-trace yeast results to human targets (route B,
transposed) + convergence layer. Design: product/backtrace_route_b/DESIGN.md;
config: configs/product_backtrace.json."""
import csv
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / "configs/product_backtrace.json").read_text())
YB = Path(CFG["yeastbridge_root"])
RESULTS = ROOT / CFG["results_dir"]
N_PERM = CFG["n_perm_compound"]
SEED = CFG["seed"]


def norm_rows(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-12)


def load_universe_projections(table, order_pos):
    idx = pd.read_csv(ROOT / CFG["universe_embeddings_dir"] / "index.tsv", sep="\t", dtype=str).fillna("")
    X = np.load(ROOT / CFG["universe_embeddings_dir"] / "esm2_mean_fp32.npy")
    keep = idx[idx["common"] != ""]
    genes = keep["common"].tolist()
    ck = torch.load(ROOT / CFG["route_b_model"], map_location="cpu", weights_only=False)["state_dict"]
    W = ck["pos_emb.proj.weight"].float().numpy()
    b = ck["pos_emb.proj.bias"].float().numpy()
    H = norm_rows(X[keep.index.to_numpy()] @ W.T + b)  # (n_uni, 768)
    return genes, H


def main() -> None:
    table = np.load(ROOT / CFG["route_b_table"]).astype(np.float64)
    order = pd.read_csv(ROOT / CFG["gene_order"], sep="\t", dtype=str)["systematic"].tolist()
    tn = norm_rows(table)
    uni_genes, H = load_universe_projections(table, order)
    uni_pos = {g: i for i, g in enumerate(uni_genes)}
    gene_pos = {g: i for i, g in enumerate(order)}

    # ---------- Part 1: precision calibration (OrthoDB as exam) ----------
    pairs, multi = {}, set()
    for r in csv.DictReader(open(ROOT / CFG["ortho_table"]), delimiter="\t"):
        y, h = r["yeast_sgd"].strip(), r["human_symbol"].strip()
        if y in gene_pos and h in uni_pos:
            if y in pairs:
                multi.add(y)
            pairs.setdefault(y, set()).add(h)
    gold = {y: list(hs)[0] for y, hs in pairs.items() if y not in multi and len(hs) == 1}
    ranks = []
    for y, h in gold.items():
        scores = H @ tn[gene_pos[y]]
        r = 1 + int(np.count_nonzero(scores > scores[uni_pos[h]]))
        ranks.append(r)
    ranks = np.array(ranks)
    calibration = {
        "n_gold_1to1": len(ranks),
        "top1": round(float((ranks == 1).mean()), 4),
        "top5": round(float((ranks <= 5).mean()), 4),
        "mrr": round(float((1.0 / ranks).mean()), 4),
        "median_rank": float(np.median(ranks)),
        "random_median": (len(uni_genes) + 1) / 2,
    }
    print("[calibration]", json.dumps(calibration), flush=True)

    # ---------- Part 2: closure + specificity + knowledge ----------
    z = np.load(CFG["response_npz"], allow_pickle=True)
    orfs = z["strain_orfs"].astype(str)
    inchi = z["compound_inchikeys"].astype(str)
    doses = z["doses"].astype(str)
    is_vehicle = z["is_vehicle"].astype(bool)
    M = z["z_score"]
    strain_rows = {g: i for i, g in enumerate(orfs)}
    # 酵母表行在菌株上的映射(缺失菌株权重置零)
    tab_rows = np.zeros((len(orfs), table.shape[1]))
    have = []
    for g, i in strain_rows.items():
        if g in gene_pos:
            tab_rows[i] = table[gene_pos[g]]
            have.append(i)
    have = np.array(have)
    tn_rows = tab_rows / np.maximum(np.linalg.norm(tab_rows, axis=1, keepdims=True), 1e-12)

    exec_df = pd.read_csv(ROOT / CFG["exec_matrix"], sep="\t")
    sig = exec_df[exec_df["q"] < CFG["fdr_alpha"]].copy()

    rng = np.random.default_rng(SEED)
    closure = {}
    for ik in sig["inchikey"].unique():
        cols = np.nonzero((inchi == ik) & ~is_vehicle)[0]
        if not len(cols):
            continue
        w = np.abs(M[:, cols].mean(axis=1))  # 化合物条件平均 |z| 谱(靶点中立,不偷看某靶点的最优剂量)
        wvec = np.zeros(len(orfs))
        for g, i in strain_rows.items():
            if g in gene_pos:
                wvec[i] = w[i]
        q = wvec @ tn_rows
        nq = np.linalg.norm(q)
        if nq == 0:
            continue
        q = q / nq
        scores = H @ q
        null_scores = np.zeros((N_PERM, len(uni_genes)))
        for b in range(N_PERM):
            wperm = np.zeros(len(orfs))
            perm = rng.permutation(len(have))
            wperm[have[perm]] = wvec[have]
            qp = wperm @ tn_rows
            nqp = np.linalg.norm(qp)
            if nqp == 0:
                continue
            null_scores[b] = H @ (qp / nqp)
        mu = null_scores.mean(axis=0)
        sd = np.maximum(null_scores.std(axis=0), 1e-9)
        closure[ik] = (scores, mu, sd)

    # 知识通道: moa 已知靶点(蛋白名规范化匹配宇宙)
    uni_meta = {}
    for fn in CFG["universe_files"]:
        for r in csv.DictReader(open(Path(CFG["universe_snapshots"]) / fn), delimiter="\t"):
            g = r["Gene Names"].split()[0].rstrip(";") if r.get("Gene Names") else ""
            pn = r.get("Protein Names", "") or ""
            if g:
                uni_meta[g] = re.sub(r"[^a-z0-9]", "", pn.lower())[:60]
    moa = {}
    for r in csv.DictReader(open(CFG["moa_labels"]), delimiter="\t"):
        key = re.sub(r"[^a-z0-9]", "", (r["target_pref_name"] or "").lower())[:60]
        moa.setdefault(r["inchikey"], set()).add(key)
    gene_by_meta = {}
    for g, m in uni_meta.items():
        gene_by_meta.setdefault(m, []).append(g)

    scf = {r["target_id"]: r["model_direction_score"] for r in
           csv.DictReader(open(ROOT / CFG["universe_targets"]), delimiter="\t")}

    rows = []
    for _, r in sig.iterrows():
        ik, t = r["inchikey"], r["target_id"]
        if ik not in closure or t not in uni_pos:
            continue
        scores, mu, sd = closure[ik]
        pi = uni_pos[t]
        rank = 1 + int(np.count_nonzero(scores > scores[pi]))
        zc = float((scores[pi] - mu[pi]) / sd[pi])
        known = any(k in gene_by_meta for k in moa.get(ik, ()))
        known_here = known and any(t in gene_by_meta[k] for k in moa.get(ik, ()) if k in gene_by_meta)
        rows.append({
            "target_id": t, "inchikey": ik, "dose": r["dose"],
            "exec_rho": r["spearman_rho"], "exec_q": r["q"],
            "closure_rank": rank, "closure_pct": round(rank / len(uni_genes), 4),
            "z_closure": round(zc, 3),
            "chembl_known_target": known, "chembl_known_is_target": known_here,
            "scf_direction_score": scf.get(t, ""),
        })
    D = pd.DataFrame(rows)

    def tier(r):
        n = 0
        n += r["exec_q"] < CFG["fdr_alpha"]
        n += r["z_closure"] >= 2.0
        n += bool(r["chembl_known_is_target"])
        n += str(r["scf_direction_score"]) not in ("", "nan")
        return "A" if n == 4 else ("B" if n == 3 else "C")

    D["tier"] = D.apply(tier, axis=1)
    D = D.sort_values(["tier", "z_closure", "exec_rho"], ascending=[True, False, False])
    RESULTS.mkdir(parents=True, exist_ok=True)
    D.to_csv(RESULTS / "backtrace_matrix.tsv", sep="\t", index=False)
    (RESULTS / "calibration.json").write_text(json.dumps(calibration, indent=2))
    D.head(100).to_csv(RESULTS / "convergence_top.tsv", sep="\t", index=False)
    n_ka = int(D["chembl_known_is_target"].sum())
    print(f"[done] pairs={len(D)} tierA={int((D['tier']=='A').sum())} "
          f"tierB={int((D['tier']=='B').sum())} known-target-hits={n_ka}", flush=True)


if __name__ == "__main__":
    main()
