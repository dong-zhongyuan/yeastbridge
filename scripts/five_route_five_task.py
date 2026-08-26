#!/usr/bin/env python
"""Five-route five-task evaluation driver (v3 final form; see
feasibility/transfer_routes/REGISTRATION_v3.md).

The evaluation is the ORIGINAL framework's own tasks (yeastbridge/eval,
frozen protocols, seed 42) run on the new scFoundation-backbone gene
tables. The three scF tables are registered at runtime into
eval.data._FEATURE_LOADERS (read-only reuse of the old project). Route E
(SGA propagation) enters by registered graph-native protocols on T2/T4
and is N/A on T1/T3/T5. No anchor construction, no intent queries.

Run with:
  cd /public/home/mengxl/dzy/yeastbridge_re && env \
  PYTHONPATH=/public/home/mengxl/dzy/yeastbridge_re/src:/public/home/mengxl/dzy/yeastbridge_re/scripts:/public/home/mengxl/dzy/yeastbridge \
  /public/home/mengxl/dzy/envs/yeastbridge/bin/python scripts/five_route_five_task.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/public/home/mengxl/dzy/yeastbridge_re")
YB = Path("/public/home/mengxl/dzy/yeastbridge")
for p in (str(ROOT / "scripts"), str(ROOT / "src"), str(YB)):
    if p not in sys.path:
        sys.path.insert(0, p)

from transfer_method_selection import ALPHA, arm_incumbent  # noqa: E402
from yeastbridge_re.second_round import _load_sga_neighbors  # noqa: E402

import eval.data as ed  # noqa: E402
import eval.tasks as et  # noqa: E402

ASSETS = ROOT / "feasibility/transfer_routes/assets"
SCF_ROUTES = ROOT / "feasibility/transfer_routes/scf_routes"
OUT = ROOT / "feasibility/transfer_routes/results"
SEED = 42
FEATURES = ["A2_scf_ortholog", "B2_scf_esm2inject", "C2_scf_scratch",
            "routed_scyeast", "esm2_mean"]
REFERENCE = "esm2_mean"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def register_scf_features():
    genes = pd.read_csv(ASSETS / "scf_yeast/A2_init.tsv", sep="\t", dtype=str)["systematic"].tolist()
    idx_df = pd.DataFrame({"systematic": genes})
    for name in ("A2_scf_ortholog", "B2_scf_esm2inject", "C2_scf_scratch"):
        X = np.load(SCF_ROUTES / name.split("_")[0] / "gene_table_final.npy")
        assert X.shape[0] == len(genes), (name, X.shape, len(genes))
        ed._FEATURE_LOADERS[name] = (lambda df=idx_df, X=X: (df, X))
        print(f"[feature] {name}: {X.shape}", flush=True)


def run_original_tasks() -> dict:
    results = {}
    for f in FEATURES:
        row = {}
        for task, fn in (("T2", et.run_t2), ("T3", et.run_t3), ("T4", et.run_t4), ("T5", et.run_t5)):
            metrics, _detail = fn(feature=f, seed=SEED)
            row[task] = {k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                         for k, v in metrics.items() if isinstance(v, (int, float, bool, np.floating, np.integer))}
            print(f"[task] {f} {task}: "
                  + json.dumps({k: row[task][k] for k in list(row[task])[:4]}), flush=True)
        results[f] = row
    return results


# ---------------- Route E: registered graph-native protocols ----------------
def e_t2_essentiality(nodes, index_of, neighbors) -> dict:
    lab = ed.essentiality()
    lab = lab[lab["essentiality"].isin(["essential", "nonessential"])]
    ess = [g for g in lab.loc[lab["essentiality"] == "essential", "systematic"] if g in index_of]
    non = [g for g in lab.loc[lab["essentiality"] == "nonessential", "systematic"] if g in index_of]
    print(f"[E-T2] essential on graph={len(ess)} nonessential pool={len(non)}", flush=True)
    non_idx = np.array([index_of[g] for g in non])
    ranks = []
    for k, o in enumerate(ess):
        seeds = {g: 1.0 for g in ess if g != o}
        vec = arm_incumbent(seeds, nodes, neighbors, alpha=ALPHA)
        s_pool = vec[non_idx]
        s_o = vec[index_of[o]]
        r = 1.0 + np.count_nonzero(s_pool > s_o) + np.count_nonzero(s_pool == s_o) / 2.0
        ranks.append(r / len(non))
        if (k + 1) % 200 == 0:
            print(f"[E-T2] {k+1}/{len(ess)} median={np.median(ranks):.4f}", flush=True)
    ranks = np.array(ranks)
    return {"median_normalized_rank": round(float(np.median(ranks)), 5),
            "pooled_auroc": round(float(np.mean(1.0 - ranks)), 5),
            "n_essential_loo": len(ess), "n_pool": len(non)}


def e_t4_engineering(nodes, index_of, neighbors) -> dict:
    rec = ed.engineering_records()
    pw = ed.pathway_genes()
    hits3 = hits5 = dir_ok = n = 0
    for _, r in rec.iterrows():
        targets = [g.strip() for g in str(r["systematic"]).split("/") if g.strip()]
        pid = r["pathway_id"]
        pw_genes = pw.loc[pw["pathway_id"] == pid, "systematic"].tolist()
        universe = sorted(set(pw_genes) | set(targets))
        if r["direction"] not in ("overexpress", "knockdown"):
            continue
        for g in targets:
            seeds = {m: 1.0 for m in universe if m != g and m in index_of}
            if not seeds or g not in index_of:
                continue
            vec = arm_incumbent(seeds, nodes, neighbors, alpha=ALPHA)
            sc = {m: float(vec[index_of[m]]) for m in universe if m in index_of}
            order_oe = sorted(sc, key=lambda x: sc[x])
            order_kd = sorted(sc, key=lambda x: -sc[x])
            order = order_oe if r["direction"] == "overexpress" else order_kd
            rank = order.index(g) + 1
            n += 1
            hits3 += rank <= 3
            hits5 += rank <= 5
            dir_ok += (sc[g] < 0) == (r["direction"] == "overexpress")
    return {"n_usable": n, "hit_at_3": round(hits3 / max(n, 1), 4),
            "hit_at_5": round(hits5 / max(n, 1), 4),
            "direction_consistency": round(dir_ok / max(n, 1), 4)}


def main() -> None:
    register_scf_features()
    results = run_original_tasks()

    neighbors = _load_sga_neighbors(ROOT, "feasibility/transfer/assets/sga_significant_p005_absge008.corrected.tsv.gz")
    all_nodes = set()
    for orf, conn in neighbors.items():
        all_nodes.add(orf)
        all_nodes.update(conn)
    nodes = sorted(all_nodes)
    index_of = {n: i for i, n in enumerate(nodes)}
    print(f"[E] graph {len(nodes)} nodes", flush=True)

    results["E_sga_propagation"] = {
        "T2": e_t2_essentiality(nodes, index_of, neighbors),
        "T4": e_t4_engineering(nodes, index_of, neighbors),
        "T1": "N/A (no gene-embedding product)",
        "T3": "N/A (no gene-embedding product)",
        "T5": "N/A (no gene-embedding product)",
    }

    # pre-declared selection: best T3 spearman among scF/scYeast routes with
    # T2 AUROC >= esm2_mean - 0.01; ties within 0.002 go to higher T2.
    ref_auroc = results[REFERENCE]["T2"]["auroc"]
    cand = {f: r for f, r in results.items() if f != REFERENCE and isinstance(r, dict) and "T3" in r}
    best_t3 = max(c[f]["T3"]["spearman_mean"] for f in cand)
    qualified = [f for f in cand if cand[f]["T2"]["auroc"] >= ref_auroc - 0.01]
    sel = None
    if qualified:
        top = [f for f in qualified if cand[f]["T3"]["spearman_mean"] >= best_t3 - 0.002]
        sel = max(top, key=lambda f: cand[f]["T2"]["auroc"])
    selection = (f"{sel} (selected: highest T3 spearman among T2-qualified routes)" if sel
                 else "no route selected; E' remains incumbent graph mechanism "
                      "(retained without new positive evidence)")

    out = {
        "registration": "feasibility/transfer_routes/REGISTRATION_v3.md",
        "form": "original five-task harness (eval/tasks.py verbatim) on scF-backed tables",
        "selection": selection,
        "reference": REFERENCE,
        "tasks": results,
        "seed": SEED,
        "asset_hashes": {f: sha256(SCF_ROUTES / f.split("_")[0] / "gene_table_final.npy")
                         for f in ("A2_scf_ortholog", "B2_scf_esm2inject", "C2_scf_scratch")},
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "five_route_results.json").write_text(json.dumps(out, indent=2, sort_keys=True))
    print(f"[done] selection: {selection}", flush=True)
    print(f"[done] wrote {OUT}/five_route_results.json", flush=True)


if __name__ == "__main__":
    main()
