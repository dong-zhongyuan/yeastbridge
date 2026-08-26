#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scf_build_yeast_assets.py — scFoundation 酵母适配资产构建(五路线 v3)。

复用两个原脚本的规则,只换嵌入源:
- scripts/routeA/build_routeA_init.py 的同源初始化规则(OrthoDB+OMA+InParanoid 并集,
  命中取均值,未命中按逐维 mean/std 高斯采样 seed=42),源从 scGPT human vocab 换成
  scFoundation models.ckpt 的 model.pos_emb.weight(19267x768)+OS 基因索引。
- scripts/routeB/protein_inject.py 的 build_esm2_matrix 规则(按 vocab 行序组装 1280d
  蛋白矩阵,缺失基因留零行,由注入层投影偏置充当共享行)。

输出(feasibility/transfer_routes/assets/scf_yeast/):
  A2_init.npy         (6739,768) 同源移植初始化(基因行 0..6735 + 2 分辨率位 + 1 pad 位)
  C2_init.npy         (6739,768) 全随机初始化(尾 3 行拷贝预训练,保持输入机制一致)
  A2_init.tsv                    每基因初始化元数据
  B2_esm2_matrix.npy  (6736,1280) 路线B''注入层蛋白矩阵
  corpus_counts.npy   (38225,6736) float32 原始 count,行=细胞,列=gene_master 序
  corpus_cells.tsv               细胞条码(行序对齐)
  build_stats.json
"""
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / "configs/scf_training.json").read_text())
YB = Path(CFG["yeastbridge_root"])
OUT = ROOT / CFG["out_assets_dir"]
SCF_CKPT = YB / CFG["scf_ckpt"]
SCF_INDEX = YB / CFG["scf_gene_index"]
M = YB / "data/mappings"
RAW_COUNTS = YB / CFG["raw_counts"]
SEED = CFG["seed"]

SOURCES = [
    ("OrthoDB", M / "orthodb_yeast_human_s288c.tsv", "yeast_sgd", "human_symbol", None),
    ("OMA", M / "oma_yeast_human.tsv", "yeast_systematic", "human_symbol", "ortholog_type"),
    ("InParanoid", M / "inparanoid_yeast_human.tsv", "yeast_systematic", "human_symbol", "ortholog_type"),
]


def load_master():
    genes = []
    with open(M / "gene_master.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            g = row["systematic"].strip()
            if g:
                genes.append((g, row["common"].strip()))
    return genes


def load_orthologs():
    table = {}
    for name, path, yc, hc, tc in SOURCES:
        with open(path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                y, h = row[yc].strip(), row[hc].strip()
                if not y or not h:
                    continue
                rec = table.setdefault(y, {"symbols": set(), "sources": set(), "one_to_one": False})
                rec["symbols"].add(h)
                rec["sources"].add(name)
                if tc and row.get(tc, "").strip() == "1:1":
                    rec["one_to_one"] = True
    return table


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    master = load_master()
    ortho = load_orthologs()
    n_genes = len(master)
    vocab_genes = [g for g, _ in master]

    ck = torch.load(SCF_CKPT, map_location="cpu", weights_only=False)
    sd = ck["gene"]["state_dict"]
    pos = sd["model.pos_emb.weight"].numpy().astype(np.float32)
    symbols = pd.read_csv(SCF_INDEX, sep="\t")["gene_name"].astype(str).tolist()
    assert len(symbols) == 19264 and pos.shape == (19267, 768), (len(symbols), pos.shape)
    hpos = {s: i for i, s in enumerate(symbols)}
    print(f"[load] master={n_genes} scf_pos={pos.shape}", flush=True)

    dim_mean = pos[:19264].mean(axis=0)
    dim_std = pos[:19264].std(axis=0)
    rng = np.random.default_rng(SEED)
    tail = pos[19264:].copy()  # 2 分辨率位 + 1 pad 位

    a2 = np.empty((n_genes, 768), dtype=np.float32)
    c2 = rng.normal(dim_mean, dim_std, size=(n_genes, 768)).astype(np.float32)
    meta_rows = []
    n_ortholog = n_random = 0
    for i, (sysname, common) in enumerate(master):
        rec = ortho.get(sysname)
        used = sorted(s for s in rec["symbols"] if s in hpos) if rec else []
        if used:
            a2[i] = pos[[hpos[s] for s in used]].mean(axis=0)
            n_ortholog += 1
        else:
            a2[i] = rng.normal(dim_mean, dim_std)
            n_random += 1
        meta_rows.append({
            "systematic": sysname, "common": common,
            "init_type": "ortholog" if used else "random",
            "n_symbols_in_scf": len(used),
            "human_symbols_used": ",".join(used),
            "sources": ",".join(sorted(rec["sources"])) if rec else "",
        })

    np.save(OUT / "A2_init.npy", np.vstack([a2, tail]))
    np.save(OUT / "C2_init.npy", np.vstack([c2, tail]))
    with open(OUT / "A2_init.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(meta_rows[0].keys()))
        w.writeheader()
        w.writerows(meta_rows)

    # B'' ESM2 蛋白矩阵(规则复制自 protein_inject.build_esm2_matrix)
    idx = pd.read_csv(YB / "data/processed/esm2_650m/index.tsv", sep="\t", dtype=str).fillna("")
    X = np.load(YB / "data/processed/esm2_650m/esm2_mean_fp32.npy")
    epos = {g: i for i, g in enumerate(idx["systematic"])}
    b2 = np.zeros((n_genes, X.shape[1]), dtype=np.float32)
    n_esm = 0
    for i, g in enumerate(vocab_genes):
        if g in epos:
            b2[i] = X[epos[g]]
            n_esm += 1
    np.save(OUT / "B2_esm2_matrix.npy", b2)

    # 原始 count 语料(行=细胞, 列=gene_master 序, 缺失基因=0)
    df = pd.read_csv(RAW_COUNTS, sep="\t", index_col=0)
    print(f"[corpus] raw tsv {df.shape}", flush=True)
    colpos = {c: i for i, c in enumerate(df.columns)}
    counts = np.zeros((df.shape[0], n_genes), dtype=np.float32)
    hit = 0
    for j, g in enumerate(vocab_genes):
        if g in colpos:
            counts[:, j] = df.iloc[:, colpos[g]].to_numpy(dtype=np.float32)
            hit += 1
    np.save(OUT / "corpus_counts.npy", counts)
    pd.DataFrame({"cell": df.index}).to_csv(OUT / "corpus_cells.tsv", sep="\t", index=False)

    stats = {
        "n_genes": n_genes,
        "A2_ortholog_init": n_ortholog, "A2_random_init": n_random,
        "B2_esm2_rows": n_esm, "B2_esm2_missing": n_genes - n_esm,
        "corpus_cells": int(counts.shape[0]),
        "corpus_genes_hit": hit,
        "scf_pos_shape": list(pos.shape), "seed": SEED,
    }
    (OUT / "build_stats.json").write_text(json.dumps(stats, indent=2))
    print("[done]", json.dumps(stats), flush=True)


if __name__ == "__main__":
    main()
