#!/usr/bin/env python3
"""Product step 2: transfer human functional targets to yeast-executable
tasks via route B's mechanism, verbatim (design: product/transfer_route_b/
DESIGN.md; config: configs/product_transfer.json).

Mechanism (route B, unchanged): protein space is shared across species.
Query = ESM2 embedding of the human target protein passed through the
TRAINED injection projection (from the finetuned scF route-B model);
score every yeast gene by cosine against the trained route-B gene table;
the ranking is the yeast-executable task.

This works for membrane targets with no yeast ortholog by construction -
no orthology table participates anywhere in this step.
"""
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / "configs/product_transfer.json").read_text())
YB = Path(CFG["yeastbridge_root"])

INPUTS = ROOT / CFG["inputs_dir"]
RESULTS = ROOT / CFG["results_dir"]
TOPK = CFG["top_k_report"]


def load_projection():
    """Trained injection-layer projection from the route-B finetuned model."""
    ck = torch.load(ROOT / CFG["route_b_model"], map_location="cpu", weights_only=False)
    sd = ck["state_dict"]
    W = sd["pos_emb.proj.weight"].float().numpy()
    b = sd["pos_emb.proj.bias"].float().numpy()
    return W, b


def load_candidate_embeddings():
    idx = pd.read_csv(INPUTS / "esm2_candidates/index.tsv", sep="\t", dtype=str).fillna("")
    idx["seq_len"] = idx["seq_len"].astype(int)
    X = np.load(INPUTS / "esm2_candidates/esm2_mean_fp32.npy")
    # 提取器的表头解析不认 ">GENE|ACC" 双段式头,uniprot 列为空;改用
    # seq_len 对回 fasta 记录(先验证长度全唯一)。
    lens = {}
    name, n = None, 0
    for line in open(INPUTS / "candidate_proteins.fasta"):
        if line.startswith(">"):
            if name:
                lens[name] = n
            name, n = line[1:].strip().split("|")[0], 0
        else:
            n += len(line.strip())
    if name:
        lens[name] = n
    assert len(set(lens.values())) == len(lens), "seq_len 非唯一,需重写 fasta 头重提取"
    len2row = {int(r.seq_len): i for i, r in idx.iterrows()}
    genes = list(lens)
    rows = [len2row[lens[g]] for g in genes]
    return genes, X[np.asarray(rows)]


def main() -> None:
    table = np.load(ROOT / CFG["route_b_table"]).astype(np.float64)
    order = pd.read_csv(ROOT / CFG["gene_order"], sep="\t", dtype=str)["systematic"].tolist()
    assert table.shape[0] == len(order)
    tn = table / np.maximum(np.linalg.norm(table, axis=1, keepdims=True), 1e-12)

    W, b = load_projection()
    cgenes, Xc = load_candidate_embeddings()
    q = Xc @ W.T + b
    q = q / np.maximum(np.linalg.norm(q, axis=1, keepdims=True), 1e-12)

    # sanity module genes (e.g. sce04011 pheromone/MAPK) from KEGG membership
    kegg = {}
    for line in open(YB / CFG["sanity_module_kegg"]):
        g, p = line.strip().split("\t")
        kegg.setdefault(p.replace("path:", ""), []).append(g.replace("sce:", ""))
    module = [g for g in kegg.get(CFG["sanity_module_id"], []) if g in set(order)]

    cands = {r["target_id"]: r for r in
             csv.DictReader(open(ROOT / CFG["candidate_ranking"]), delimiter="\t")}

    RESULTS.mkdir(parents=True, exist_ok=True)
    summary = []
    for k, g in enumerate(cgenes):
        scores = tn @ q[k]
        rank = np.argsort(-scores)
        genes_sorted = [order[i] for i in rank]
        top = [(order[i], float(scores[i])) for i in rank[:TOPK]]
        out = RESULTS / f"yeast_task_{g}.tsv"
        with open(out, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["rank", "yeast_gene", "cosine"])
            for i, (gene, s) in enumerate([(order[i], float(scores[i])) for i in rank], 1):
                w.writerow([i, gene, round(s, 5)])
        pos = {gene: i for i, gene in enumerate(genes_sorted)}
        mod_ranks = np.array([pos[m] for m in module]) + 1
        all_ranks = np.arange(1, len(order) + 1)
        summary.append({
            "target_id": g,
            "family": cands[g]["target_family"],
            "intended_direction": cands[g]["intended_direction"],
            "top5": [t[0] for t in top[:5]],
            "top1pct": [t[0] for t in top[:67]],
            "module_median_rank": float(np.median(mod_ranks)),
            "module_random_median_rank": (len(order) + 1) / 2.0,
            "module_top50_count": int(sum(1 for r in mod_ranks if r <= 50)),
        })
        print(f"{g} ({cands[g]['target_family']}, {cands[g]['intended_direction']}): "
              f"top5={summary[-1]['top5']} sce04011 median rank={summary[-1]['module_median_rank']:.0f} "
              f"(random~{summary[-1]['module_random_median_rank']:.0f})", flush=True)

    (RESULTS / "transfer_summary.json").write_text(
        json.dumps({"config": "configs/product_transfer.json",
                    "mechanism": "route B verbatim: ESM2(protein) -> trained proj -> cosine vs trained route-B gene table",
                    "n_targets": len(summary),
                    "sanity_module": CFG["sanity_module_id"],
                    "targets": summary}, indent=2))
    print(f"[done] {len(summary)} targets -> {RESULTS}", flush=True)


if __name__ == "__main__":
    main()
