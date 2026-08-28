#!/usr/bin/env python
"""pd_validation application: evaluate the branch (ChEMBL-annotated) docked
pairs with the DeepPurpose pretrained ensemble and apply the two gates
(docking affinity <= dock_gate_kcal AND ensemble mean pIC50 >= pic50_gate).
Input pairs: product/chembl_branch/results/gpu_dock_pairs.tsv (best pocket
per pair). Sequences from the universe fasta; neutral SMILES from the
original exec matrix. Registered in product/pd_validation/DESIGN.md;
parameters in configs/pd_validation.json.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MODELS = [
    ("MPNN", "CNN", "mpnn_cnn_bindingdb_ic50"),
    ("CNN", "CNN", "cnn_cnn_bindingdb_ic50"),
    ("Morgan", "CNN", "morgan_cnn_bindingdb_ic50"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/pd_validation.json")
    args = ap.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())

    import pandas as pd
    from DeepPurpose import utils, DTI

    pairs = pd.read_csv(
        ROOT / "product/chembl_branch/results/gpu_dock_pairs.tsv", sep="\t")

    seq_of = {}
    name = None
    for line in (ROOT / cfg["universe_fasta"]).read_text().splitlines():
        if line.startswith(">"):
            name = line.split("|")[1]
        elif name:
            seq_of.setdefault(name, "")
            seq_of[name] += line.strip()

    em = pd.read_csv(
        ROOT / "product/execute_hiphop/results/exec_matrix.tsv", sep="\t")
    smi_of = dict(zip(em["inchikey"], em["smiles"]))

    rows = []
    for r in pairs.itertuples():
        seq = seq_of.get(r.acc, "")
        smi = smi_of.get(r.inchikey, "")
        if seq and smi:
            rows.append(dict(inchikey=r.inchikey, target_gene=r.target_gene,
                             acc=r.acc, docking_affinity=r.affinity,
                             pchembl_measured=r.pchembl,
                             smiles=smi, seq=seq))
    df = pd.DataFrame(rows)
    print(f"pairs to evaluate: {len(df)}", flush=True)

    preds = {}
    for drug_enc, tgt_enc, name_m in MODELS:
        path = utils.download_pretrained_model(name_m)
        model = DTI.model_pretrained(path)
        X = utils.data_process(df["smiles"].tolist(), df["seq"].tolist(),
                               [0] * len(df), drug_enc, tgt_enc,
                               split_method="no_split")
        p = model.predict(X)
        preds[name_m] = [float(v) for v in p]
        print(f"{name_m}: done", flush=True)
    for k, v in preds.items():
        df[k] = v
    df["pd_pic50_mean"] = df[[m[2] for m in MODELS]].mean(axis=1)

    df["dock_gate"] = df["docking_affinity"] <= cfg["dock_gate_kcal"]
    df["pd_gate"] = df["pd_pic50_mean"] >= cfg["pic50_gate"]
    df["final_pass"] = df["dock_gate"] & df["pd_gate"]

    res = ROOT / cfg["results_dir"]
    res.mkdir(parents=True, exist_ok=True)
    out_cols = ["inchikey", "target_gene", "acc", "docking_affinity",
                "pd_pic50_mean"] + [m[2] for m in MODELS] + \
        ["pchembl_measured", "dock_gate", "pd_gate", "final_pass"]
    df[out_cols].sort_values(
        ["final_pass", "docking_affinity"],
        ascending=[False, True]).to_csv(
        res / "final_pairs.tsv", sep="\t", index=False)
    summary = dict(n_pairs=len(df),
                   n_dock_gate=int(df["dock_gate"].sum()),
                   n_pd_gate=int(df["pd_gate"].sum()),
                   n_final=int(df["final_pass"].sum()))
    (res / "apply_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
