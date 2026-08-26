#!/usr/bin/env python
"""Transfer-method selection driver (registered; feasibility/transfer/REGISTRATION.md).

F1 method comparison + F2 component ablation of the human-intent -> yeast
transfer step, evaluated by leave-one-out anchor recovery on the SGA graph.

All core logic is reused from verbatim copies of the legacy project:
`yeastbridge_re.second_round` (propagation, loaders; copied from
yeastbridge_vs discovery_release/second_round.py) and the intent
fusion-gating rule copied from its `_load_intent`. Only the experiment
wiring (arms, LOO loop, paired bootstrap) is new. Run with:
PYTHONPATH=src /public/home/mengxl/dzy/envs/yeastbridge_vs/bin/python
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from yeastbridge_re.second_round import (
    _load_orthology,
    _load_sga_neighbors,
    _propagate_personalized_pagerank,
)

ROOT = Path("/public/home/mengxl/dzy/yeastbridge_re")
ASSETS = ROOT / "feasibility/transfer/assets"
OUT = ROOT / "feasibility/transfer/results"
ALPHA = 0.30  # production design value
BOOT_SEED = 20260826
REWIRE_SEED = 20260826
PERM_SEED = 20260826
N_BOOT = 10000
MATERIAL_DELTA = 0.02
TOL = 1e-10
MAX_ITERS = 300

# --- intent weights: loading + predeclared fusion-gating rule copied from
# legacy second_round.py::_load_intent, pointed at the scFoundation signature
# scores from product step 1 ---
def load_intent_weights(apply_gate: bool) -> dict[str, float]:
    with (ASSETS / "state_signature.tsv").open(encoding="utf-8-sig", newline="") as h:
        weights = {
            str(row["gene"]): float(row["state_effect_disease_minus_desired"])
            for row in csv.DictReader(h, delimiter="\t")
        }
    if not apply_gate:
        return weights
    fused = {}
    with (ASSETS / "signature_model_scores.tsv").open(encoding="utf-8-sig", newline="") as h:
        for row in csv.DictReader(h, delimiter="\t"):
            fused[str(row["gene"])] = row
    for gene, stat in list(weights.items()):
        entry = fused.get(gene)
        if entry is None or entry.get("model_direction_score") in (None, "", "None"):
            weights[gene] = stat * 0.75
            continue
        m = float(entry["model_direction_score"])
        if np.sign(m) != np.sign(stat):
            weights[gene] = stat * 0.5
    return weights


# --- anchor construction copied verbatim from legacy _feature_weights ---
def build_anchor_weights(orthology, weights):
    anchor_weights: dict[str, float] = {}
    for orf, genes in orthology.items():
        human_weights = [weights[gene] for gene in genes if gene in weights]
        if human_weights:
            anchor_weights[orf] = -float(np.mean(human_weights))
    return anchor_weights


# --- walk-matrix construction copied verbatim from legacy
# _propagate_personalized_pagerank ---
def build_walk(nodes, neighbors):
    index_of = {node: i for i, node in enumerate(nodes)}
    size = len(nodes)
    rows, cols = [], []
    seen_pairs: set[tuple[int, int]] = set()
    for orf, connected in neighbors.items():
        i = index_of.get(orf)
        if i is None:
            continue
        for other in connected:
            j = index_of.get(other)
            if j is None or (i, j) in seen_pairs:
                continue
            seen_pairs.add((i, j))
            rows.append(i)
            cols.append(j)
    data = np.ones(len(rows), dtype=np.float64)
    adjacency = sp.coo_matrix((data, (rows, cols)), shape=(size, size)).tocsr()
    row_sums = np.asarray(adjacency.sum(axis=1)).ravel()
    walk = sp.diags(1.0 / np.maximum(row_sums, 1.0)) @ adjacency
    return walk, index_of


def signed_mass(anchor_weights, index_of, size):
    mass = np.zeros(size)
    total = sum(abs(w) for w in anchor_weights.values())
    if total <= 0.0:
        return mass
    for orf, w in anchor_weights.items():
        i = index_of.get(orf)
        if i is not None:
            mass[i] = w / total
    return mass


# --- arm definitions (incumbent calls the verbatim legacy propagation) ---
def arm_incumbent(anchors, nodes, neighbors, **kw):
    vec, _ = _propagate_personalized_pagerank(nodes, neighbors, anchors, alpha=kw.get("alpha", ALPHA))
    return vec


def arm_no_propagation(anchors, nodes, neighbors, **kw):
    walk, index_of = kw["walk"], kw["index_of"]
    return signed_mass(anchors, index_of, len(nodes))


def arm_one_hop(anchors, nodes, neighbors, **kw):
    walk, index_of = kw["walk"], kw["index_of"]
    return walk.T @ signed_mass(anchors, index_of, len(nodes))


def arm_full_diffusion(anchors, nodes, neighbors, **kw):
    walk, index_of = kw["walk"], kw["index_of"]
    value = np.abs(signed_mass(anchors, index_of, len(nodes)))
    for _ in range(MAX_ITERS):
        nxt = walk.T @ value
        if float(np.abs(nxt - value).max()) < TOL:
            value = nxt
            break
        value = nxt
    return value


def arm_unsigned(anchors, nodes, neighbors, **kw):
    vec, _ = _propagate_personalized_pagerank(
        nodes, neighbors, {o: abs(w) for o, w in anchors.items()}, alpha=ALPHA
    )
    return vec


def rewire_degree_preserving(neighbors, seed):
    """Degree-preserving edge-swap rewiring of the undirected SGA graph."""
    edges = set()
    for a, conn in neighbors.items():
        for b in conn:
            if a != b:
                edges.add((a, b) if a < b else (b, a))
    edges = sorted(edges)
    rng = np.random.default_rng(seed)
    edge_set = set(edges)
    edges = list(edges)
    n_swaps = 10 * len(edges)
    for _ in range(n_swaps):
        i, j = rng.integers(0, len(edges), size=2)
        (a, b), (c, d) = edges[i], edges[j]
        if len({a, b, c, d}) < 4:
            continue
        e1 = (a, d) if a < d else (d, a)
        e2 = (c, b) if c < b else (b, c)
        if e1 in edge_set or e2 in edge_set or e1 == e2:
            continue
        edge_set.discard((a, b))
        edge_set.discard((c, d))
        edge_set.add(e1)
        edge_set.add(e2)
        edges[i], edges[j] = e1, e2
    rewired: dict[str, set[str]] = {}
    for a, b in edge_set:
        rewired.setdefault(a, set()).add(b)
        rewired.setdefault(b, set()).add(a)
    return rewired


def main() -> None:
    orthology = _load_orthology(ROOT, "feasibility/transfer/assets/orthodb_yeast_human_s288c.tsv")
    neighbors = _load_sga_neighbors(ROOT, "feasibility/transfer/assets/sga_significant_p005_absge008.corrected.tsv.gz")
    weights_gated = load_intent_weights(apply_gate=True)
    weights_raw = load_intent_weights(apply_gate=False)
    anchors_gated = build_anchor_weights(orthology, weights_gated)
    anchors_raw = build_anchor_weights(orthology, weights_raw)
    print(f"anchors: gated={len(anchors_gated)} raw={len(anchors_raw)}", flush=True)

    all_nodes: set[str] = set()
    for orf, connected in neighbors.items():
        all_nodes.add(orf)
        all_nodes.update(connected)
    all_nodes.update(anchors_gated)
    nodes = sorted(all_nodes)
    walk, index_of = build_walk(nodes, neighbors)
    n_nodes = len(nodes)
    print(f"graph: {n_nodes} nodes", flush=True)

    # controls prepared once (registered seeds)
    rewired = rewire_degree_preserving(neighbors, REWIRE_SEED)
    perm_values = list(anchors_gated.values())
    rng_perm = np.random.default_rng(PERM_SEED)
    perm_idx = rng_perm.permutation(len(perm_values))
    anchors_permuted = {
        orf: perm_values[perm_idx[k]] for k, orf in enumerate(anchors_gated)
    }

    arms_f1 = {
        "incumbent": lambda a: arm_incumbent(a, nodes, neighbors),
        "no_propagation": lambda a: arm_no_propagation(a, nodes, neighbors, walk=walk, index_of=index_of),
        "one_hop": lambda a: arm_one_hop(a, nodes, neighbors, walk=walk, index_of=index_of),
        "full_diffusion": lambda a: arm_full_diffusion(a, nodes, neighbors, walk=walk, index_of=index_of),
        "rewired_graph_control": lambda a: _propagate_personalized_pagerank(nodes, rewired, a, alpha=ALPHA)[0],
        "permuted_anchor_control": lambda a: arm_incumbent(
            {o: anchors_permuted[o] for o in a}, nodes, neighbors
        ),
    }
    arms_f2 = {
        "full": arms_f1["incumbent"],
        "no_fusion_gate": lambda a: arm_incumbent(a, nodes, neighbors),  # anchors built from raw weights
        "unsigned": lambda a: arm_unsigned(a, nodes, neighbors),
        "no_restart": arms_f1["full_diffusion"],
        "alpha_0.10": lambda a: arm_incumbent(a, nodes, neighbors, alpha=0.10),
        "alpha_0.50": lambda a: arm_incumbent(a, nodes, neighbors, alpha=0.50),
        "alpha_0.70": lambda a: arm_incumbent(a, nodes, neighbors, alpha=0.70),
        "uniform_weights": lambda a: arm_incumbent({o: np.sign(w) * 1.0 for o, w in a.items()}, nodes, neighbors),
    }

    eligible = sorted(o for o, w in anchors_gated.items() if w != 0.0 and o in index_of)
    top1 = max(1, int(np.ceil(0.01 * n_nodes)))
    top5 = max(1, int(np.ceil(0.05 * n_nodes)))
    print(f"eligible LOO anchors: {len(eligible)}; top1%={top1} top5%={top5} nodes", flush=True)

    def run_arm(arm_fn, anchor_source):
        med_ranks = np.zeros(len(eligible))
        hit1 = np.zeros(len(eligible))
        hit5 = np.zeros(len(eligible))
        for k, o in enumerate(eligible):
            reduced = {x: w for x, w in anchor_source.items() if x != o}
            vec = arm_fn(reduced)
            sign_o = np.sign(anchor_source[o])
            score = vec * sign_o
            # v2 endpoint: rank o within the non-anchor pool only (remaining
            # true anchors occupy the top of a full-graph ranking by
            # construction under restart methods and are excluded)
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
            med_ranks[k] = rank / pool_size
            hit1[k] = 1.0 if rank <= max(1, int(np.ceil(0.01 * pool_size))) else 0.0
            hit5[k] = 1.0 if rank <= max(1, int(np.ceil(0.05 * pool_size))) else 0.0
        return med_ranks, hit1, hit5

    def summarize(name, med_ranks, hit1, hit5):
        return {
            "arm": name,
            "median_normalized_rank": round(float(np.median(med_ranks)), 5),
            "mean_normalized_rank": round(float(med_ranks.mean()), 5),
            "hit_at_1pct": round(float(hit1.mean()), 4),
            "hit_at_5pct": round(float(hit5.mean()), 4),
        }

    def paired_boot(inc_mr, arm_mr):
        # primary statistic: per-ORF normalized-rank difference
        # (arm minus incumbent; positive = incumbent ranks better)
        diff = arm_mr - inc_mr
        rng = np.random.default_rng(BOOT_SEED)
        idx = rng.integers(0, len(diff), size=(N_BOOT, len(diff)))
        boot_meds = np.median(diff[idx], axis=1)
        ci = [float(np.percentile(boot_meds, 2.5)), float(np.percentile(boot_meds, 97.5))]
        delta = float(np.median(diff))
        return {
            "delta_median_rank_arm_minus_incumbent": round(delta, 5),
            "ci95": [round(c, 5) for c in ci],
            "incumbent_beats_arm": bool(delta >= MATERIAL_DELTA and ci[0] > 0),
        }

    results = {}
    for label, arms, anchor_source in (("F1", arms_f1, anchors_gated), ("F2", arms_f2, anchors_gated)):
        per_arm = {}
        for name, fn in arms.items():
            src = anchors_raw if (label == "F2" and name == "no_fusion_gate") else anchor_source
            mr, h1, h5 = run_arm(fn, src)
            per_arm[name] = {"summary": summarize(name, mr, h1, h5), "hit1": h1, "hit5": h5, "mr": mr}
            print(label, per_arm[name]["summary"], flush=True)
        incumbent_key = "incumbent" if label == "F1" else "full"
        comparisons = {}
        for name in per_arm:
            if name == incumbent_key:
                continue
            comparisons[name] = paired_boot(per_arm[incumbent_key]["mr"], per_arm[name]["mr"])
        results[label] = {
            "n_loo_anchors": len(eligible),
            "n_nodes": n_nodes,
            "arms": {k: v["summary"] for k, v in per_arm.items()},
            "comparisons_vs_incumbent": comparisons,
        }

    # F1 gates (registered): incumbent must beat both controls
    f1 = results["F1"]
    f1["gates"] = {
        "structure_informed": bool(
            f1["comparisons_vs_incumbent"]["rewired_graph_control"]["incumbent_beats_arm"]
            and f1["comparisons_vs_incumbent"]["permuted_anchor_control"]["incumbent_beats_arm"]
        ),
        "network_step_justified": bool(
            f1["comparisons_vs_incumbent"]["no_propagation"]["incumbent_beats_arm"]
            and f1["comparisons_vs_incumbent"]["one_hop"]["incumbent_beats_arm"]
        ),
    }
    # F2: arms not beaten by full flag unjustified components
    f2 = results["F2"]
    f2["flags"] = {
        name: ("full NOT significantly better -> component unjustified"
               if not c["incumbent_beats_arm"] else "full significantly better")
        for name, c in f2["comparisons_vs_incumbent"].items()
    }

    results["registration"] = "feasibility/transfer/REGISTRATION.md"
    results["endpoint"] = "leave-one-out anchor recovery on the SGA graph"
    results["statistics"] = {"n_boot": N_BOOT, "seed": BOOT_SEED, "material_delta": MATERIAL_DELTA}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "transfer_method_selection_results.json").write_text(json.dumps(results, indent=2, sort_keys=True))
    print(json.dumps({"F1_gates": f1["gates"], "F2_flags": f2["flags"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
