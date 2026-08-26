#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scf_finetune_yeast.py — scFoundation 酵母微调(五路线 v3, REGISTRATION_v3.md)。

用法:
  cd /public/home/mengxl/dzy/yeastbridge_re
  /public/home/mengxl/dzy/envs/yeastbridge/bin/python scripts/scf_finetune_yeast.py --route A2

结构与官方 demo 对齐:
- 模型/配置加载: vendored scfoundation load.load_model_frommmf(key='gene')
- 输入配方: get_embedding.py singlecell 't4' 逐字(log1p CPM-10k + [4.0, log10(total)])
- 数据构造: load.getEncoerDecoderData
- 训练粒度: 官方 finetune_model.py LinearProbingClassifier 同款 —— 冻结 encoder
  第 0..9 层,训练末 2 层 + token_emb + pos_emb(或注入层) + decoder
- 遮蔽: 表达基因 p=0.30 / 零值基因 p=0.03,纯置 mask_token_id(注册偏差:未复刻
  replace/random 腐蚀);损失 = 掩码位置连续值 MSE
- 优化: AdamW lr 1e-4, clip 1.0, batch 8 x accum 4, epochs 6, seed 42, bf16 autocast

路线 A2(同源移植)/B2(ESM2 注入, ProteinEmbeddingInjector 规则复制自
scripts/routeB/protein_inject.py)/C2(全随机)仅 pos_emb 来源不同。
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

YB = Path("/public/home/mengxl/dzy/yeastbridge")
ASSETS = Path("/public/home/mengxl/dzy/yeastbridge_re/feasibility/transfer_routes/assets/scf_yeast")
OUTDIR = Path("/public/home/mengxl/dzy/yeastbridge_re/feasibility/transfer_routes/scf_routes")
MODEL_DIR = YB / "src/external/scfoundation/model"
sys.path.insert(0, str(MODEL_DIR))

from load import gatherData, getEncoerDecoderData, load_model_frommmf  # noqa: E402

N_GENES = 6733  # gene_master 真实基因数(旧 routeA_vocab 6736 = 6733 基因 + 3 特殊符号)
SEQ_LEN = N_GENES + 2  # + 2 分辨率位; pos_emb 行数 = SEQ_LEN + 1(pad 位)
MAX_SEQ = 1200  # 原 finetune_scgpt.py train_args.max_seq_len=1200, trunc_by_sample 官方行为
EPOCHS = 6
BATCH = 32
ACCUM = 1
LR = 1e-4
CLIP = 1.0
SEED = 42
MASK_P = 0.30
ZERO_MASK_P = 0.03


class ProteinEmbeddingInjector(nn.Module):
    """逐字复制自 scripts/routeB/protein_inject.py,仅 d_model 由调用方给定。"""

    def __init__(self, prot_mat: np.ndarray, n_special: int, d_model: int):
        super().__init__()
        self.n_genes = prot_mat.shape[0]
        self.register_buffer("prot", torch.from_numpy(prot_mat.astype(np.float32)))
        self.proj = nn.Linear(prot_mat.shape[1], d_model)
        self.special = nn.Embedding(n_special, d_model)

    def forward(self, idx: torch.Tensor) -> torch.Tensor:
        is_special = idx >= self.n_genes
        gene_idx = idx.clamp(max=self.n_genes - 1)
        out = self.proj(self.prot[gene_idx])
        if bool(is_special.any()):
            sp = self.special((idx - self.n_genes).clamp(min=0))
            out = torch.where(is_special.unsqueeze(-1), sp, out)
        return out


class Checkpointed(nn.Module):
    """梯度检查点包装:训练时重算换显存(否则 token_emb 可训练会把 12 层注意力
    矩阵全部驻留,~47GB)。use_reentrant=False 下末 2 层可训练参数照常收梯度。"""

    def __init__(self, m):
        super().__init__()
        self.m = m

    def forward(self, x, padding_mask=None):
        if self.training and x.requires_grad:
            return torch.utils.checkpoint.checkpoint(
                self.m, x, padding_mask=padding_mask, use_reentrant=False)
        return self.m(x, padding_mask=padding_mask)


def encdec_subset(data, data_raw, gene_ids, config):
    """load.getEncoerDecoderData 的逐字逻辑,唯一改动:data_gene_ids 由调用方给出
    (子集行的真实基因列号),不再假设 arange(全词表宽)。"""
    decoder_data = data.clone().detach()
    decoder_data_padding = torch.full_like(data, False, dtype=torch.bool).to(data.device)
    encoder_data_labels = data_raw > 0
    encoder_data, encoder_data_padding = gatherData(decoder_data, encoder_data_labels,
                                                    config["pad_token_id"])
    new_data_raw = data_raw
    encoder_position_gene_ids, _ = gatherData(gene_ids, encoder_data_labels,
                                              config["pad_token_id"])
    decoder_position_gene_ids = gene_ids.clone()
    data_mask_labels = None
    encoder_position_gene_ids[encoder_data_padding] = config["seq_len"]
    decoder_position_gene_ids[decoder_data_padding] = config["seq_len"]
    return (encoder_data, encoder_position_gene_ids, encoder_data_padding,
            encoder_data_labels, decoder_data, decoder_data_padding,
            new_data_raw, data_mask_labels, decoder_position_gene_ids)


def pick_device(requested):
    if requested == "auto":
        if torch.cuda.is_available() and torch.cuda.mem_get_info()[0] / 1024**3 >= 8.0:
            return torch.device("cuda")
        return torch.device("cpu")
    return torch.device(requested)


def build_model(route, device):
    model, cfg = load_model_frommmf(str(YB / "models/scfoundation/models.ckpt"), "gene")
    cfg = dict(cfg)
    cfg["seq_len"] = SEQ_LEN
    cfg["gene_num"] = SEQ_LEN
    for sub in ("encoder", "decoder"):
        if isinstance(cfg.get(sub), dict):
            cfg[sub]["seq_len"] = SEQ_LEN

    from pretrainmodels import select_model  # noqa: E402
    model = select_model(cfg)
    sd = torch.load(YB / "models/scfoundation/models.ckpt", map_location="cpu", weights_only=False)
    sd = sd["gene"]["state_dict"]
    sd = {k[len("model."):] if k.startswith("model.") else k: v for k, v in sd.items()}
    sd.pop("pos_emb.weight", None)  # 酵母词表形状不同,加载后整体替换
    missing, unexpected = model.load_state_dict(sd, strict=False)
    bad = [k for k in missing if "pos_emb" not in k]
    assert not bad and not unexpected, f"权重对齐失败 missing={bad[:5]} unexpected={unexpected[:5]}"

    if route in ("A2", "C2"):
        init = np.load(ASSETS / f"{route}_init.npy")
        assert init.shape == (SEQ_LEN + 1, 768), init.shape
        emb = nn.Embedding(SEQ_LEN + 1, 768)
        with torch.no_grad():
            emb.weight.copy_(torch.from_numpy(init))
        model.pos_emb = emb
    else:
        prot = np.load(ASSETS / "B2_esm2_matrix.npy")
        model.pos_emb = ProteinEmbeddingInjector(prot, n_special=3, d_model=768)

    enc_layers = model.encoder.transformer_encoder
    enc_layers = enc_layers.layers if hasattr(enc_layers, "layers") else enc_layers
    for layer in enc_layers[:-2]:
        for p in layer.parameters():
            p.requires_grad = False
    model.encoder = Checkpointed(model.encoder)
    model.decoder = Checkpointed(model.decoder)
    model = model.to(device)
    return model, cfg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--route", required=True, choices=["A2", "B2", "C2"])
    ap.add_argument("--device", default="auto")
    ap.add_argument("--smoke", action="store_true", help="single train step + table save, then exit")
    args = ap.parse_args()

    torch.manual_seed(SEED)
    np.random.seed(SEED)
    device = pick_device(args.device)
    print(f"[setup] route={args.route} device={device}", flush=True)

    counts = np.load(ASSETS / "corpus_counts.npy", mmap_mode="r")
    n_cells = counts.shape[0]
    model, cfg = build_model(args.route, device)
    trainable = [p for p in model.parameters() if p.requires_grad]
    print(f"[setup] cells={n_cells} trainable_params={sum(p.numel() for p in trainable)/1e6:.1f}M", flush=True)
    opt = torch.optim.AdamW(trainable, lr=LR)

    rng = np.random.default_rng(SEED)
    perm = rng.permutation(n_cells)
    n_val = int(0.05 * n_cells)
    val_idx, train_idx = perm[:n_val], perm[n_val:]
    outdir = OUTDIR / args.route
    outdir.mkdir(parents=True, exist_ok=True)

    def build_subset(batch_raw):
        """每细胞: 表达基因中均匀随机抽 MAX_SEQ 个(trunc_by_sample 官方行为),
        附 2 个分辨率位; 返回 (x_sub, gene_ids_sub), 位置 id = 真实基因列号。"""
        raw = torch.from_numpy(np.asarray(batch_raw, dtype=np.float32))
        total = raw.sum(dim=1, keepdim=True).clamp(min=1.0)
        vals = torch.log1p(raw / total * 1e4)
        res_v = torch.cat([torch.full((raw.shape[0], 1), 4.0), torch.log10(total)], dim=1)
        xs, gs = [], []
        for b in range(vals.shape[0]):
            idx = torch.nonzero(vals[b] > 0, as_tuple=True)[0]
            if idx.numel() > MAX_SEQ:
                perm = torch.randperm(idx.numel())[:MAX_SEQ]
                idx = idx[perm]
            xs.append(torch.cat([vals[b, idx], res_v[b]]))
            gs.append(torch.cat([idx, torch.tensor([N_GENES, N_GENES + 1])]))
        K = max(x.numel() for x in xs)
        x_sub = torch.zeros((len(xs), K), dtype=torch.float32)
        g_sub = torch.full((len(gs), K), SEQ_LEN, dtype=torch.long)  # 尾部填充位 → pad 行
        for i, (xv, gv) in enumerate(zip(xs, gs)):
            x_sub[i, :xv.numel()] = xv
            g_sub[i, :gv.numel()] = gv
        return x_sub, g_sub

    def step(batch_raw, train):
        x, gene_ids = build_subset(batch_raw)
        x, gene_ids = x.to(device), gene_ids.to(device)
        expressed = x[:, :-2] > 0
        r = torch.rand_like(x[:, :-2])
        mask = (expressed & (r < MASK_P)) | (~expressed & (r < ZERO_MASK_P))
        masked = x.clone()
        masked[:, :-2][mask] = float(cfg["mask_token_id"])
        enc = encdec_subset(masked, x, gene_ids, cfg)
        (enc_data, enc_pos, enc_pad, enc_labels, dec_data, dec_pad, _, _, dec_pos) = enc
        with torch.autocast("cuda", dtype=torch.bfloat16, enabled=device.type == "cuda"):
            out = model(x=enc_data, padding_label=enc_pad,
                        encoder_position_gene_ids=enc_pos, encoder_labels=enc_labels,
                        decoder_data=dec_data, mask_gene_name=False, mask_labels=None,
                        decoder_position_gene_ids=dec_pos,
                        decoder_data_padding_labels=dec_pad)
            loss = nn.functional.mse_loss(out[:, :-2][mask], x[:, :-2][mask])
        if train:
            (loss / ACCUM).backward()
        return float(loss)

    def save_table(tag):
        model.eval()
        with torch.no_grad():
            if args.route == "B2":
                inj = model.pos_emb
                table = inj.proj(inj.prot).cpu().numpy().astype(np.float32)
            else:
                table = model.pos_emb.weight[:N_GENES].detach().cpu().numpy().astype(np.float32)
        np.save(outdir / f"gene_table_{tag}.npy", table)
        model.train()

    t0 = time.time()
    if args.smoke:
        batch = counts[np.sort(train_idx[:BATCH])]
        loss = step(batch, train=True)
        nn.utils.clip_grad_norm_(trainable, CLIP)
        opt.step()
        save_table("smoke")
        print(f"[smoke] one step OK, loss={loss:.4f}; table saved", flush=True)
        return
    for epoch in range(1, EPOCHS + 1):
        order = rng.permutation(len(train_idx))
        opt.zero_grad(set_to_none=True)
        running, seen = 0.0, 0
        for b in range(0, len(order), BATCH):
            batch = counts[np.sort(train_idx[order[b:b + BATCH]])]
            loss = step(batch, train=True)
            running += loss
            seen += 1
            if seen % ACCUM == 0:
                nn.utils.clip_grad_norm_(trainable, CLIP)
                opt.step()
                opt.zero_grad(set_to_none=True)
            if seen % 200 == 0:
                el = time.time() - t0
                print(f"[ep{epoch}] step {seen}/{len(order)//BATCH} loss={running/seen:.4f} "
                      f"elapsed={el/60:.1f}min", flush=True)
        vloss, vb = 0.0, 0
        for b in range(0, len(val_idx), BATCH):
            vloss += step(counts[np.sort(val_idx[b:b + BATCH])], train=False)
            vb += 1
        print(f"[ep{epoch}] train={running/max(seen,1):.4f} val={vloss/max(vb,1):.4f}", flush=True)
        save_table(f"ep{epoch}")

    torch.save({"state_dict": {k: v.cpu() for k, v in model.state_dict().items()},
                "route": args.route, "config": cfg},
               outdir / "final_model.pt")
    save_table("final")
    meta = {"route": args.route, "epochs": EPOCHS, "batch": BATCH, "accum": ACCUM,
            "lr": LR, "mask_p": MASK_P, "zero_mask_p": ZERO_MASK_P, "seed": SEED,
            "n_train": len(train_idx), "n_val": len(val_idx)}
    (outdir / "train_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[done] {outdir} ({(time.time()-t0)/60:.1f} min)", flush=True)


if __name__ == "__main__":
    main()
