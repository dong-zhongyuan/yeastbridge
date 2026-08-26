#!/usr/bin/env python
"""Five-route five-task transfer evaluation driver (registered;
feasibility/transfer_routes/REGISTRATION.md).

Adaptation of the original YeastBridge-Eval framework to the
functional-target transfer step. Five routes compete under five tasks;
route E' reuses the registered SGA propagation verbatim, dense routes
query the original project's frozen embedding tables read-only, and the
human-side knowledge source is the scFoundation intent of product step 1.

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

from transfer_method_selection import (  # noqa: E402
    ALPHA,
    arm_incumbent,
    build_anchor_weights,
    load_intent_weights,
)
from yeastbridge_re.second_round import _load_orthology, _load_sga_neighbors  # noqa: E402
from eval.data import esm2_embeddings, routea_gene_embeddings, scyeast_gene_embeddings  # noqa: E402

ASSETS = ROOT / "feasibility/transfer_routes/assets"
SHARED = ROOT / "feasibility/transfer/assets"
OUT = ROOT / "feasibility/transfer_routes/results"
KEMM = YB / "data/ko_kemmeren/kemmeren_2014"
PATHWAY = YB / "data/pathway"
BOOT_SEED = 20260827
PERM_SEED = 20260827
N_BOOT = 10000
T1_BOOT = 2000
MATERIAL_DELTA = 0.02
INCUMBENT = "E_sga_propagation"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def row_normalized(X: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(X, axis=1, keepdims=True)
    return X / np.maximum(n, 1e-12)


def make_dense_route(table_rows: dict, Xn: np.ndarray, query_fn, index_of, n_nodes):
    """Dense cosine route: score = cos(embedding, query(reduced anchors)),
    assigned to graph nodes; nodes without an embedding score 0."""
    genes = [g for g in table_rows if g in index_of]
    rows = np.array([table_rows[g] for g in genes])
    node_idx = np.array([index_of[g] for g in genes])

    def produce(reduced):
        q = query_fn(reduced)
        vec = np.zeros(n_nodes)
        if q is not None:
            vec[node_idx] = Xn[rows] @ q
        return vec

    return produce


def anchor_space_query(table_rows: dict, X: np.ndarray):
    """Query = signed-weighted mean of the reduced anchors' own embeddings."""
    def q(reduced):
        idxs, ws = [], []
        for orf, w in reduced.items():
            r = table_rows.get(orf)
            if r is not None:
                idxs.append(r)
                ws.append(w)
        if not idxs:
            return None
        v = np.asarray(ws, dtype=np.float64) @ X[np.asarray(idxs)]
        n = np.linalg.norm(v)
        return v / n if n > 0 else None

    return q


def human_ortholog_query(orthology, human_rows: dict, Xh: np.ndarray):
    """Route B' query: signed-weighted mean of ESM2 embeddings of the HUMAN
    ortholog proteins of the reduced anchors (mean over each anchor's human
    genes present in the human table)."""
    per_anchor = {}
    for orf, genes in orthology.items():
        rs = [human_rows[g] for g in genes if g in human_rows]
        if rs:
            per_anchor[orf] = Xh[np.asarray(rs)].mean(axis=0)

    def q(reduced):
        ws, vs = [], []
        for orf, w in reduced.items():
            v = per_anchor.get(orf)
            if v is not None:
                ws.append(w)
                vs.append(v)
        if not vs:
            return None
        v = np.asarray(ws, dtype=np.float64) @ np.asarray(vs)
        n = np.linalg.norm(v)
        return v / n if n > 0 else None

    return q


def load_human_esm2():
    d = ASSETS / "esm2_human_anchors"
    idx = pd.read_csv(d / "index.tsv", sep="\t", dtype=str).fillna("")
    X = np.load(d / "esm2_mean_fp32.npy")
    rows: dict = {}
    for i, sym in enumerate(idx["common"]):
        if sym:
            rows.setdefault(sym, []).append(i)
    merged = {s: k for k, s in enumerate(rows)}
    Xm = np.stack([X[np.asarray(rs)].mean(axis=0) for s, rs in rows.items()])
    return merged, row_normalized(Xm)


def run_loo(produce, anchor_source, eligible, index_of, n_nodes):
    """LOO anchor recovery on the corrected v2 endpoint (non-anchor pool).
    Returns per-anchor normalized rank, AUROC, hit@1%, hit@5% and the score
    vectors (the latter reused by T3')."""
    mr = np.zeros(len(eligible))
    au = np.zeros(len(eligible))
    h1 = np.zeros(len(eligible))
    h5 = np.zeros(len(eligible))
    vecs = []
    for k, o in enumerate(eligible):
        reduced = {x: w for x, w in anchor_source.items() if x != o}
        vec = produce(reduced)
        vecs.append(vec)
        score = vec * np.sign(anchor_source[o])
        pool_mask = np.ones(n_nodes, dtype=bool)
        for x in reduced:
            xi = index_of.get(x)
            if xi is not None:
                pool_mask[xi] = False
        pool_scores = score[pool_mask]
        s_o = score[index_of[o]]
        n_greater = int(np.count_nonzero(pool_scores > s_o))
        n_tied = int(np.count_nonzero(pool_scores == s_o)) + 1
        rank = 1.0 + n_greater + (n_tied - 1) / 2.0
        pool_size = len(pool_scores)
        mr[k] = rank / pool_size
        au[k] = 1.0 - (rank - 1.0) / pool_size
        h1[k] = 1.0 if rank <= max(1, int(np.ceil(0.01 * pool_size))) else 0.0
        h5[k] = 1.0 if rank <= max(1, int(np.ceil(0.05 * pool_size))) else 0.0
    return mr, au, h1, h5, vecs


def summarize(name, mr, au, h1, h5):
    return {
        "route": name,
        "median_normalized_rank": round(float(np.median(mr)), 5),
        "mean_normalized_rank": round(float(mr.mean()), 5),
        "pooled_auroc": round(float(au.mean()), 5),
        "hit_at_1pct": round(float(h1.mean()), 4),
        "hit_at_5pct": round(float(h5.mean()), 4),
    }


def boot_median_ci(diff: np.ndarray):
    rng = np.random.default_rng(BOOT_SEED)
    idx = rng.integers(0, len(diff), size=(N_BOOT, len(diff)))
    boot_meds = np.median(diff[idx], axis=1)
    return float(np.median(diff)), float(np.percentile(boot_meds, 2.5)), float(np.percentile(boot_meds, 97.5))


def main() -> None:
    orthology = _load_orthology(ROOT, "feasibility/transfer/assets/orthodb_yeast_human_s288c.tsv")
    neighbors = _load_sga_neighbors(ROOT, "feasibility/transfer/assets/sga_significant_p005_absge008.corrected.tsv.gz")
    weights = load_intent_weights(apply_gate=True)
    anchors = build_anchor_weights(orthology, weights)

    all_nodes = set()
    for orf, connected in neighbors.items():
        all_nodes.add(orf)
        all_nodes.update(connected)
    all_nodes.update(anchors)
    nodes = sorted(all_nodes)
    index_of = {n: i for i, n in enumerate(nodes)}
    n_nodes = len(nodes)
    print(f"graph: {n_nodes} nodes; anchors: {len(anchors)}", flush=True)

    # frozen embedding tables (read-only reuse of the original project)
    tables = {}
    for name, kind in (("A_homolog_repr", "ft"), ("C_scratch_control", "scratch")):
        idx, X = routea_gene_embeddings(kind)
        tables[name] = ({g: i for i, g in enumerate(idx["systematic"])}, X)
    idx, X = scyeast_gene_embeddings()
    tables["D_native_scyeast"] = ({g: i for i, g in enumerate(idx["systematic"])}, X)
    idx, X = esm2_embeddings()
    tables["B_yeast_esm2"] = ({g: i for i, g in enumerate(idx["systematic"])}, X)
    Xn = {k: row_normalized(v[1]) for k, v in tables.items()}
    human_rows, Xh_n = load_human_esm2()
    print("tables: " + ", ".join(f"{k}={v[1].shape}" for k, v in tables.items())
          + f"; human={Xh_n.shape}", flush=True)

    routes = {
        "A_homolog_repr": make_dense_route(
            tables["A_homolog_repr"][0], Xn["A_homolog_repr"],
            anchor_space_query(tables["A_homolog_repr"][0], tables["A_homolog_repr"][1]),
            index_of, n_nodes),
        "B_protein_bridge": make_dense_route(
            tables["B_yeast_esm2"][0], Xn["B_yeast_esm2"],
            human_ortholog_query(orthology, human_rows, Xh_n), index_of, n_nodes),
        "C_scratch_control": make_dense_route(
            tables["C_scratch_control"][0], Xn["C_scratch_control"],
            anchor_space_query(tables["C_scratch_control"][0], tables["C_scratch_control"][1]),
            index_of, n_nodes),
        "D_native_scyeast": make_dense_route(
            tables["D_native_scyeast"][0], Xn["D_native_scyeast"],
            anchor_space_query(tables["D_native_scyeast"][0], tables["D_native_scyeast"][1]),
            index_of, n_nodes),
        "E_sga_propagation": lambda reduced: arm_incumbent(reduced, nodes, neighbors, alpha=ALPHA),
    }
    t4_table_of = {
        "A_homolog_repr": "A_homolog_repr",
        "B_protein_bridge": "B_yeast_esm2",
        "C_scratch_control": "C_scratch_control",
        "D_native_scyeast": "D_native_scyeast",
    }

    perm_values = list(anchors.values())
    rng = np.random.default_rng(PERM_SEED)
    pp = rng.permutation(len(perm_values))
    anchors_perm = {orf: perm_values[pp[k]] for k, orf in enumerate(anchors)}

    eligible = sorted(o for o, w in anchors.items() if w != 0.0 and o in index_of)
    print(f"eligible LOO anchors: {len(eligible)}", flush=True)

    # ---------------- T2' primary + cached LOO vectors ----------------
    t2 = {}
    vec_cache = {}
    perm_mr = {}
    for name, produce in routes.items():
        mr, au, h1, h5, vecs = run_loo(produce, anchors, eligible, index_of, n_nodes)
        vec_cache[name] = vecs
        pmr, pau, ph1, ph5, _ = run_loo(produce, anchors_perm, eligible, index_of, n_nodes)
        perm_mr[name] = pmr
        t2[name] = {"summary": summarize(name, mr, au, h1, h5), "mr": mr,
                    "permuted_summary": summarize(name + "_permuted", pmr, pau, ph1, ph5)}
        print("T2'", t2[name]["summary"], flush=True)

    # G1: route vs incumbent on T2' median rank (route must be BETTER by
    # material delta and CI excluding zero; lower rank is better)
    g1 = {}
    for name in routes:
        if name == INCUMBENT:
            continue
        delta, lo, hi = boot_median_ci(t2[name]["mr"] - t2[INCUMBENT]["mr"])
        g1[name] = {"delta_median_rank_route_minus_incumbent": round(delta, 5),
                    "ci95": [round(lo, 5), round(hi, 5)],
                    "route_beats_incumbent": bool(delta <= -MATERIAL_DELTA and hi < 0)}
    # G2 water line: route vs its own permuted-anchor control
    # (CI of median(route - permuted) upper bound < 0; no material delta)
    g2 = {}
    for name in routes:
        delta, lo, hi = boot_median_ci(t2[name]["mr"] - perm_mr[name])
        g2[name] = {"delta_median_rank_route_minus_permuted": round(delta, 5),
                    "ci95": [round(lo, 5), round(hi, 5)],
                    "route_beats_permuted": bool(hi < 0)}
    selection = f"{INCUMBENT} (retained: no route passed G1 and G2)"
    for name, c in g1.items():
        if c["route_beats_incumbent"] and g2[name]["route_beats_permuted"]:
            selection = f"{name} (passed G1 vs incumbent and G2 water line)"
            break
    print(f"T2' G1: {json.dumps({k: v['route_beats_incumbent'] for k, v in g1.items()})}", flush=True)
    print(f"T2' G2: {json.dumps({k: v['route_beats_permuted'] for k, v in g2.items()})}", flush=True)
    print(f"T2' selection: {selection}", flush=True)

    # ---------------- T1' state-direction geometry ----------------
    kem = pd.read_parquet(KEMM / "kemmeren_t3_expr_log2fc.parquet")
    meta = pd.read_csv(KEMM / "kemmeren_t3_sample_metadata.tsv", sep="\t", dtype=str).fillna("")
    m2s = dict(zip(meta["mutant_name"], meta["mutant_systematic"]))
    kem_genes = [c for c in kem.columns if c != "mutant"]
    mutant_names = kem["mutant"].tolist()
    prof = kem[kem_genes].to_numpy(dtype=np.float64)
    sys_of = {m: m2s.get(m, "") for m in mutant_names}
    measured = [k for k, m in enumerate(mutant_names) if sys_of[m] in index_of]
    target_vec = np.array([anchors.get(g, 0.0) for g in kem_genes])
    t_cols = np.nonzero(target_vec != 0.0)[0]
    print(f"T1' kemmeren: {len(measured)} mutants on graph; intent-defined columns={len(t_cols)}", flush=True)

    t1 = {}
    for name, produce in routes.items():
        vec = produce(anchors)
        wts = np.array([vec[index_of[sys_of[mutant_names[k]]]] for k in measured])
        predicted = wts @ prof[measured]
        pr = pd.Series(predicted[t_cols]).rank().to_numpy()
        tr = pd.Series(target_vec[t_cols]).rank().to_numpy()
        rho = float(np.corrcoef(pr, tr)[0, 1])
        rb = np.random.default_rng(BOOT_SEED + 1)
        boots = np.array([np.corrcoef(pr[i], tr[i])[0, 1]
                          for i in rb.integers(0, len(t_cols), size=(T1_BOOT, len(t_cols)))])
        t1[name] = {"spearman": round(rho, 5),
                    "ci95": [round(float(np.nanpercentile(boots, 2.5)), 5),
                             round(float(np.nanpercentile(boots, 97.5)), 5)]}
        print("T1'", name, t1[name], flush=True)

    # ---------------- T3' perturbation matching ----------------
    mutant_by_sys: dict = {}
    for k, m in enumerate(mutant_names):
        s = sys_of[m]
        if s:
            mutant_by_sys.setdefault(s, []).append(k)
    eligible3 = [o for o in eligible if o in mutant_by_sys]
    shared_cols = np.array([i for i, g in enumerate(kem_genes) if g in index_of])
    P = prof[:, shared_cols]
    Pn = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-12)
    shared_node_idx = [index_of[kem_genes[i]] for i in shared_cols]
    t3 = {}
    for name in routes:
        ranks, hit5 = [], 0.0
        for o in eligible3:
            vec = vec_cache[name][eligible.index(o)]
            v = vec[shared_node_idx]
            nv = np.linalg.norm(v)
            if nv == 0:
                for _ in mutant_by_sys[o]:
                    ranks.append(1.0)
                continue
            a = Pn @ (v / nv)
            for kk in mutant_by_sys[o]:
                r = 1.0 + float(np.count_nonzero(a > a[kk]))
                ranks.append(r / len(mutant_names))
                hit5 += r <= max(1, int(np.ceil(0.05 * len(mutant_names))))
        status = "descriptive (eligible < 15, registered downgrade)" if len(eligible3) < 15 else "eligible"
        t3[name] = {"eligible_anchors": len(eligible3),
                    "n_replicate_ranks": len(ranks),
                    "median_normalized_rank": round(float(np.median(ranks)), 5) if ranks else None,
                    "hit_at_5pct": round(float(hit5 / max(len(ranks), 1)), 4) if ranks else None,
                    "status": status}
        print("T3'", name, t3[name], flush=True)

    # ---------------- T4' engineering-record ranking ----------------
    records = pd.read_csv(PATHWAY / "engineering_records.tsv", sep="\t", dtype=str).fillna("")
    pathway_genes = {}
    for pid in sorted(set(records["pathway_id"])):
        pathway_genes[pid] = set(
            pd.read_csv(PATHWAY / f"pathway_{pid}_genes.tsv", sep="\t", dtype=str)["systematic"])
    t4_rows = []
    for name, produce in routes.items():
        hits3 = hits5 = dir_ok = n = 0
        for _, rec in records.iterrows():
            targets = [s for s in str(rec["systematic"]).split("/") if s]
            pid = rec["pathway_id"]
            members = pathway_genes.get(pid, set()) | set(targets)
            if name == "E_sga_propagation":
                query_members = [g for g in members if g not in targets and g in index_of]
                if not query_members or not any(g in index_of for g in targets):
                    continue
                vec = produce({g: 1.0 for g in query_members})
                sc = {g: float(vec[index_of[g]]) for g in members if g in index_of}
            else:
                tkey = t4_table_of[name]
                trows, tXn = tables[tkey][0], Xn[tkey]
                query_members = [g for g in members if g not in targets and g in trows]
                if not query_members or not any(g in trows for g in targets):
                    continue
                qv = tXn[np.array([trows[g] for g in query_members])].mean(axis=0)
                nq = np.linalg.norm(qv)
                if nq == 0:
                    continue
                qv = qv / nq
                sc = {g: float(tXn[trows[g]] @ qv) for g in members if g in trows}
            present = [g for g in targets if g in sc]
            if not present:
                continue
            best = max(present, key=lambda g: sc[g])
            r = 1 + sum(1 for v in sc.values() if v > sc[best])
            n += 1
            hits3 += r <= 3
            hits5 += r <= 5
            want = 1.0 if str(rec["direction"]) == "overexpress" else -1.0
            dir_ok += np.sign(sc[best]) == want
        t4_rows.append({"route": name, "n_records": n,
                        "hit_at_3": round(hits3 / max(n, 1), 4),
                        "hit_at_5": round(hits5 / max(n, 1), 4),
                        "direction_consistency": round(dir_ok / max(n, 1), 4)})
        print("T4'", t4_rows[-1], flush=True)

    # ---------------- T5' anchor-efficiency curve ----------------
    t5 = {name: {} for name in routes}
    for frac in (0.25, 0.5, 0.75, 1.0):
        seeds = (1,) if frac == 1.0 else (1, 2, 3)
        for seed in seeds:
            rb = np.random.default_rng(10_000 * seed + int(frac * 100))
            size = max(2, int(round(len(eligible) * frac)))
            sub = [str(x) for x in rb.choice(np.array(eligible), size=size, replace=False)]
            sub_src = {o: anchors[o] for o in sub}
            key = f"frac{frac}_s{seed}"
            for name, produce in routes.items():
                mr5, _, _, _, _ = run_loo(produce, sub_src, sub, index_of, n_nodes)
                t5[name][key] = round(float(np.median(mr5)), 5)
        print(f"T5' frac={frac} done", flush=True)

    results = {
        "registration": "feasibility/transfer_routes/REGISTRATION.md",
        "selection": selection,
        "T2_anchor_recovery": {
            "arms": {k: {**v["summary"], "permuted": v["permuted_summary"]} for k, v in t2.items()},
            "G1_route_vs_incumbent": g1,
            "G2_water_lines": g2,
            "n_loo_anchors": len(eligible),
            "n_nodes": n_nodes,
        },
        "T1_state_geometry": t1,
        "T3_perturbation_matching": t3,
        "T4_engineering_ranking": t4_rows,
        "T5_anchor_efficiency": t5,
        "statistics": {"n_boot": N_BOOT, "t1_boot": T1_BOOT, "seed": BOOT_SEED,
                       "material_delta": MATERIAL_DELTA},
        "asset_hashes": {
            "human_anchor_proteins.fasta": sha256(ASSETS / "human_anchor_proteins.fasta"),
            "esm2_human_anchors/esm2_mean_fp32.npy": sha256(ASSETS / "esm2_human_anchors/esm2_mean_fp32.npy"),
            "orthodb": sha256(SHARED / "orthodb_yeast_human_s288c.tsv"),
            "sga_corrected": sha256(SHARED / "sga_significant_p005_absge008.corrected.tsv.gz"),
            "signature_model_scores": sha256(SHARED / "signature_model_scores.tsv"),
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "five_route_results.json").write_text(json.dumps(results, indent=2, sort_keys=True))
    print(f"[done] wrote {OUT}/five_route_results.json", flush=True)


if __name__ == "__main__":
    main()
