"""Shared assets for YeastBridge-VS product layer, extracted verbatim from
the YeastBridge-VS evidence repository so every new component reuses one
implementation instead of re-writing its own copy.

Sources (VS commit anchors):
- SGA edge parsing: scripts/paired_four_arm.py::_sga_modules (edge block)
- Norman program clustering: scripts/norman_retrieval.py (KMeans block)
- Foundation embedding extraction: scripts/crc_model_evidence.py
  (product three-model merged-vocab version, commit 1137ac2)
"""

from __future__ import annotations

import gzip
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path("/public/home/mengxl/dzy/yeastbridge_vs")
SGA_TSV = next(
    (ROOT / "data/curated/yeast/sga").glob(
        "*/sga_significant_p005_absge008.corrected.tsv.gz"
    )
)
NORMAN = ROOT / "reports/human_foundation_norman/full"
FEATURES = ROOT / "reports/human_foundation_norman/features"
SCGPT_VOCAB = "/public/home/mengxl/dzy/yeastbridge/models/scgpt/scGPT_human/vocab.json"
SCGPT_CKPT = "/public/home/mengxl/dzy/yeastbridge/models/scgpt/scGPT_human/best_model.pt"
SCF_CKPT = "/public/home/mengxl/dzy/yeastbridge/models/scfoundation/models.ckpt"
SCF_INDEX = (
    "/public/home/mengxl/dzy/yeastbridge/src/external/scfoundation/model/"
    "OS_scRNA_gene_index.19264.tsv"
)
SEED = 20260821


def sga_edges():
    """Repaired-SGA undirected adjacency (verbatim logic from paired_four_arm.py)."""
    edges = defaultdict(set)
    with gzip.open(SGA_TSV, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        si = {name: i for i, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            q, a = fields[si["query_orf"]], fields[si["array_orf"]]
            if q != a:
                edges[q].add(a)
                edges[a].add(q)
    return edges


def norman_programs(n_programs: int = 6):
    """KMeans programs over the frozen Norman task (verbatim from norman_retrieval.py).

    Returns (program label per condition, anchors, response gene symbols)."""
    from sklearn.cluster import KMeans

    oof = np.load(NORMAN / "raw_oof_full.npz", allow_pickle=True)
    y = oof["y_true"].astype(np.float64)
    km = KMeans(n_clusters=n_programs, random_state=SEED, n_init=10)
    programs = km.fit_predict(y)
    return programs, km.cluster_centers_, [str(s) for s in oof["response_gene_symbols"]]


def load_three_model_embeddings():
    """Three-model gene embedding tables (the registered winning fusion form).

    Returns a dict model_name -> (vocab: symbol->row, matrix (n, d)), exactly the
    merged-vocab construction from crc_model_evidence.py (product 1137ac2)."""
    import torch

    tables = {}

    vocab = json.load(open(SCGPT_VOCAB))
    ckpt = torch.load(SCGPT_CKPT, map_location="cpu", weights_only=False)
    tables["scgpt"] = (vocab, ckpt["encoder.embedding.weight"].numpy().astype(np.float64))

    import pandas as pd

    symbols = pd.read_csv(SCF_INDEX, sep="\t").iloc[:, 0].astype(str).tolist()
    scf = torch.load(SCF_CKPT, map_location="cpu", weights_only=False)
    pos = (
        scf["gene"]["state_dict"]["model.pos_emb.weight"]
        .numpy()
        .astype(np.float64)[: len(symbols)]
    )
    tables["scfoundation"] = ({s: i for i, s in enumerate(symbols)}, pos)

    gf_targets = np.load(FEATURES / "geneformer_target_embeddings.npz", allow_pickle=True)
    gf_norman = np.load(FEATURES / "geneformer.npz", allow_pickle=True)
    gvocab, rows = {}, []
    for s, e in zip(gf_norman["gene_symbols"], gf_norman["embeddings"], strict=True):
        gvocab[str(s)] = len(rows)
        rows.append(np.asarray(e, dtype=np.float64))
    for s, e in zip(gf_targets["symbols"], gf_targets["embeddings"], strict=True):
        s = str(s)
        if s not in gvocab:
            gvocab[s] = len(rows)
            rows.append(np.asarray(e, dtype=np.float64))
    tables["geneformer"] = (gvocab, np.array(rows))

    return tables


def fused_ortholog_embedding(gene_symbols, tables=None):
    """Per-gene fused embedding = concatenation of every model's embedding of
    its (possibly multiple) human orthologs; genes absent from a model's vocab
    contribute zeros plus availability is tracked by the caller.

    Input: iterable of human gene symbols; multiple orthologs are averaged
    per model first (same rule as crc_model_evidence.py).
    Returns (matrix (n, sum_dims), per-model availability dict of boolean arrays).
    """
    if tables is None:
        tables = load_three_model_embeddings()
    blocks, availability = [], {}
    for name, (vocab, matrix) in tables.items():
        dim = matrix.shape[1]
        block = np.zeros((len(gene_symbols), dim), dtype=np.float64)
        avail = np.zeros(len(gene_symbols), dtype=bool)
        for i, symbol in enumerate(gene_symbols):
            if symbol in vocab:
                block[i] = matrix[vocab[symbol]]
                avail[i] = True
        blocks.append(block)
        availability[name] = avail
    return np.hstack(blocks), availability
