#!/usr/bin/env python
"""Extract ESM2 per-protein mean-pooled embeddings for the HUMAN anchor
proteins of yeastbridge_re route B' (feasibility/transfer_routes).

Copied from yeastbridge/scripts/data/extract_esm2_embeddings.py (W1); the
only change is the --device option (auto: prefer GPU with enough free VRAM,
fall back to CPU when both GPUs are occupied by root services; GPU is
retried-preferred on later runs). Batch budgeting logic unchanged.
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np
import pandas as pd
import torch
import esm

PROJECT = "/public/home/mengxl/dzy/yeastbridge"


def read_fasta(path):
    records, name, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:], []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def parse_header(header):
    """>sp|A0A0B7P3V8|YP41B_YEAST Desc ... GN=TY4B-P ... -> (accession, gene_name)"""
    accession, gene_name = "", ""
    parts = header.split("|")
    if len(parts) >= 3:
        accession = parts[1]
    if " GN=" in header:
        gene_name = header.split(" GN=")[1].split()[0]
    return accession, gene_name


def make_batches(items, tok_per_batch, attn_budget):
    """items: [(idx, label, seq)] sorted by len desc. Greedy pack with two caps:
    - padded tokens: max_len * batch_size <= tok_per_batch
    - attention elements: batch_size * max_len^2 <= attn_budget
      (attn scores are O(B * heads * L^2) and get a fp32 upcast in fair-esm,
      so long-sequence batches blow up without the quadratic cap)"""
    batches, cur, cur_max = [], [], 0
    for it in items:
        need = len(it[2]) + 2  # BOS + EOS
        new_max = max(cur_max, need)
        n = len(cur) + 1
        if cur and (new_max * n > tok_per_batch or new_max * new_max * n > attn_budget):
            batches.append(cur)
            cur, cur_max = [it], need
        else:
            cur.append(it)
            cur_max = new_max
    if cur:
        batches.append(cur)
    return batches


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", default=f"{PROJECT}/data/proteome/teacher_sa.fasta")
    ap.add_argument("--gene-master", default=f"{PROJECT}/data/mappings/gene_master.tsv")
    ap.add_argument("--model", default=f"{PROJECT}/models/esm2/esm2_t33_650M_UR50D.pt")
    ap.add_argument("--outdir", default=f"{PROJECT}/data/processed/esm2_650m")
    ap.add_argument("--repr-layer", type=int, default=33)
    ap.add_argument("--tok-per-batch", type=int, default=65536)
    ap.add_argument("--attn-budget", type=float, default=8e7,
                    help="max B*L^2 attention elements per batch")
    ap.add_argument("--limit", type=int, default=0, help="debug: only first N proteins")
    ap.add_argument("--device", default="auto",
                    help="auto: prefer GPU when CUDA is available with >=8GB free VRAM, else CPU; or force cuda/cpu")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    if args.device == "auto":
        if torch.cuda.is_available() and torch.cuda.mem_get_info()[0] / 1024**3 >= 8.0:
            device = torch.device("cuda")
        else:
            device = torch.device("cpu")
    else:
        device = torch.device(args.device)
    assert device.type == "cpu" or torch.cuda.is_available(), "CUDA not available"
    print(f"[setup] device={device} (requested={args.device}) model={args.model} layer={args.repr_layer}", flush=True)

    t0 = time.time()
    model, alphabet = esm.pretrained.load_model_and_alphabet_local(args.model)
    if device.type == "cuda":
        model = model.eval().half().to(device)
    else:
        model = model.eval().to(device)
    for p in model.parameters():
        p.requires_grad_(False)
    batch_converter = alphabet.get_batch_converter()
    dim = model.embed_dim
    print(f"[setup] model loaded in {time.time()-t0:.1f}s, embed_dim={dim}, "
          f"dtype={'fp16' if device.type == 'cuda' else 'fp32'}", flush=True)

    records = read_fasta(args.fasta)
    if args.limit:
        records = records[: args.limit]
    print(f"[data] {len(records)} proteins from {args.fasta}", flush=True)

    # uniprot accession -> (systematic, common)
    gm = pd.read_csv(args.gene_master, sep="\t", dtype=str).fillna("")
    acc2sys = dict(zip(gm["uniprot"], gm["systematic"]))
    acc2common = dict(zip(gm["uniprot"], gm["common"]))

    items = []  # (orig_idx, label, seq)
    meta = []
    for i, (header, seq) in enumerate(records):
        acc, gname = parse_header(header)
        items.append((i, acc or f"prot{i}", seq))
        meta.append((acc2sys.get(acc, ""), acc2common.get(acc, "") or gname, acc, len(seq)))
    items.sort(key=lambda x: -len(x[2]))
    batches = make_batches(items, args.tok_per_batch, args.attn_budget)
    total_res = sum(len(s) for _, _, s in items)
    print(f"[data] {len(batches)} batches (budget {args.tok_per_batch} tok), "
          f"total_residues={total_res}, mapped_to_gene_master={sum(1 for m in meta if m[0])}", flush=True)

    N = len(items)
    embs = np.zeros((N, dim), dtype=np.float32)
    nan_count = 0
    t1, done_res = time.time(), 0
    for bi, batch in enumerate(batches):
        labels = [lab for _, lab, _ in batch]
        seqs = [s for _, _, s in batch]
        with torch.inference_mode():
            _, _, toks = batch_converter(list(zip(labels, seqs)))
            toks = toks.to(device, non_blocking=True)
            out = model(toks, repr_layers=[args.repr_layer])
            reps = out["representations"][args.repr_layer]  # (B, T, C) fp16
            for j, (orig_idx, _, seq) in enumerate(batch):
                n = len(seq)
                v = reps[j, 1 : n + 1].float().mean(0).cpu().numpy()
                if np.isnan(v).any():
                    nan_count += 1
                embs[orig_idx] = v
        done_res += sum(len(s) for s in seqs)
        if bi % 5 == 0 or bi == len(batches) - 1:
            el = time.time() - t1
            print(f"[run] batch {bi+1}/{len(batches)} seqs={len(seqs)} "
                  f"res={done_res}/{total_res} ({done_res/el:.0f} res/s) elapsed={el:.0f}s", flush=True)
        del toks, out, reps

    order = np.array([orig for batch in batches for (orig, _, _) in batch], dtype=np.int64)
    rows = [meta[k] for k in order]
    idx = pd.DataFrame(rows, columns=["systematic", "common", "uniprot", "seq_len"])
    X = np.stack([embs[k] for k in order]).astype(np.float32)

    np.save(os.path.join(args.outdir, "esm2_mean_fp32.npy"), X)
    idx.to_csv(os.path.join(args.outdir, "index.tsv"), sep="\t", index=False)
    pq = idx.copy()
    pq["embedding"] = list(X)
    pq.to_parquet(os.path.join(args.outdir, "esm2_mean.parquet"), index=False)
    info = dict(model=os.path.basename(args.model), repr_layer=args.repr_layer,
                n_proteins=int(N), dim=int(dim), total_residues=int(total_res),
                nan_embeddings=int(nan_count), seconds=round(time.time() - t1, 1),
                dtype="float32 mean-pooled over residues (excl BOS/EOS), model run in fp16")
    with open(os.path.join(args.outdir, "run_info.json"), "w") as fh:
        json.dump(info, fh, indent=2, ensure_ascii=False)
    print(f"[done] {json.dumps(info)}", flush=True)
    print(f"[done] wrote {args.outdir}/{{esm2_mean_fp32.npy,index.tsv,esm2_mean.parquet,run_info.json}}", flush=True)


if __name__ == "__main__":
    main()
