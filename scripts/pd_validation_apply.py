#!/usr/bin/env python
"""pd_validation: calibration (wide-range holdout) and application of the
pretrained DeepPurpose ensemble to the branch docked pairs.

Stages:
  --stage calibrate : fetch a wide-range ChEMBL holdout (quantitative
      activities for the branch targets), evaluate all pretrained models,
      fit per-model linear rescaling onto the pChEMBL scale (required to
      mix BindingDB-IC50 / KIBA / DAVIS output units), persist
      results/calibration.json. Metric note (registered): narrow-range
      Pearson r on the delivered pairs is attenuation-dominated; engine
      performance is reported as wide-range r/RMSE, delivered pairs as MAE.
  --stage apply     : load calibration.json, rescale + average the ensemble
      for the delivered pairs, recompute gates, write
      results/final_pairs_calibrated.tsv.

Input pairs: product/chembl_branch/results/gpu_dock_pairs.tsv. Sequences
from the universe fasta; neutral SMILES from the exec matrix. Registered in
product/pd_validation/DESIGN.md; parameters in configs/pd_validation.json.
"""
from __future__ import annotations

import argparse
import json
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MODELS = [
    ("MPNN", "CNN", "mpnn_cnn_bindingdb_ic50"),
    ("CNN", "CNN", "cnn_cnn_bindingdb_ic50"),
    ("Morgan", "CNN", "morgan_cnn_bindingdb_ic50"),
    ("MPNN", "CNN", "mpnn_cnn_kiba"),
    ("Morgan", "AAC", "morgan_aac_kiba"),
    ("Morgan", "CNN", "morgan_cnn_kiba"),
    ("MPNN", "CNN", "mpnn_cnn_davis"),
    ("CNN", "CNN", "cnn_cnn_davis"),
    ("Morgan", "AAC", "morgan_aac_bindingdb_ic50"),
]


def get_json(url, tries=6):
    for t in range(tries):
        try:
            with urllib.request.urlopen(
                    urllib.request.Request(url), timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception:  # noqa: BLE001
            time.sleep(min(120, 5 * 2 ** t))
    return None


def seq_map(cfg):
    seq_of = {}
    name = None
    for line in (ROOT / cfg["universe_fasta"]).read_text().splitlines():
        if line.startswith(">"):
            name = line.split("|")[1]
        elif name:
            seq_of.setdefault(name, "")
            seq_of[name] += line.strip()
    return seq_of


def predict_all(smiles, seqs):
    """Per-model raw predictions for aligned lists."""
    import numpy as np
    from DeepPurpose import utils, DTI

    out = {}
    for drug_enc, tgt_enc, name_m in MODELS:
        path = utils.download_pretrained_model(name_m)
        model = DTI.model_pretrained(path)
        X = utils.data_process(smiles, seqs, [0] * len(smiles),
                               drug_enc, tgt_enc, split_method="no_split")
        out[name_m] = np.array([float(v) for v in model.predict(X)])
        print(f"  {name_m}: done", flush=True)
    return out


def calibrate(cfg, per_target_limit):
    import numpy as np
    import pandas as pd

    res = ROOT / cfg["results_dir"]
    inp = ROOT / cfg["inputs_dir"]
    inp.mkdir(parents=True, exist_ok=True)
    base = cfg["chembl_base"]

    bp = pd.read_csv(
        ROOT / "product/chembl_branch/results/branch_pairs.tsv", sep="\t")
    acc_of = dict(zip(bp["target_chembl_id"], bp["acc"]))
    rows = []
    for j, tc in enumerate(sorted(acc_of), 1):
        a = get_json(f"{base}/activity.json?target_chembl_id={tc}"
                     f"&pchembl_value__isnull=false&limit={per_target_limit}")
        for act in (a or {}).get("activities", []):
            rows.append(dict(target_chembl_id=tc,
                             molecule=act.get("molecule_chembl_id"),
                             pchembl=float(act["pchembl_value"])))
        time.sleep(0.2)
        if j % 10 == 0:
            print(f"  activities {j}/{len(acc_of)} targets", flush=True)
    ho = (pd.DataFrame(rows).dropna()
            .sort_values("pchembl", ascending=False)
            .drop_duplicates(["target_chembl_id", "molecule"]))
    mols = sorted(ho["molecule"].unique())
    smi_of = {}
    for s in range(0, len(mols), 50):
        q = ",".join(mols[s:s + 50])
        rec = get_json(f"{base}/molecule.json?molecule_chembl_id__in={q}"
                       f"&limit=50")
        for m in (rec or {}).get("molecules", []):
            ms = m.get("molecule_structures") or {}
            if ms.get("canonical_smiles"):
                smi_of[m["molecule_chembl_id"]] = ms["canonical_smiles"]
        time.sleep(0.2)
    ho = ho[ho["molecule"].isin(smi_of)].reset_index(drop=True)
    ho["smiles"] = ho["molecule"].map(smi_of)
    ho["acc"] = ho["target_chembl_id"].map(acc_of)
    ho.to_csv(inp / "calibration_holdout.tsv", sep="\t", index=False)

    seq_of = seq_map(cfg)
    P = predict_all(ho["smiles"].tolist(),
                    [seq_of[a] for a in ho["acc"]])
    y = ho["pchembl"].to_numpy()
    report = dict(
        n_holdout=int(len(ho)),
        pchembl_range=[float(y.min()), float(y.max())],
        per_model={},
        caveat="holdout shares source DBs with DeepPurpose pretrained data; "
               "metrics may be mildly optimistic; used for rescaling")
    calib = {}
    ens = None
    for name_m, p in P.items():
        a1, b1 = np.polyfit(p, y, 1)
        calib[name_m] = [round(float(a1), 4), round(float(b1), 4)]
        rs = a1 * p + b1
        ens = rs if ens is None else ens + rs
        report["per_model"][name_m] = dict(
            r=round(float(np.corrcoef(p, y)[0, 1]), 3),
            rmse_raw=round(float(np.sqrt(np.mean((p - y) ** 2))), 2),
            rmse_rescaled=round(float(np.sqrt(np.mean((rs - y) ** 2))), 2))
        print(f"  {name_m}: r={report['per_model'][name_m]['r']}", flush=True)
    ens /= len(P)
    report["ensemble"] = dict(
        r=round(float(np.corrcoef(ens, y)[0, 1]), 3),
        rmse=round(float(np.sqrt(np.mean((ens - y) ** 2))), 2),
        mae=round(float(np.mean(np.abs(ens - y))), 2))
    report["rescale"] = calib
    (res / "calibration.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report["ensemble"], indent=2), flush=True)


def apply_stage(cfg):
    import numpy as np
    import pandas as pd

    res = ROOT / cfg["results_dir"]
    cal = json.loads((res / "calibration.json").read_text())
    df = pd.read_csv(res / "final_pairs.tsv", sep="\t")
    seq_of = seq_map(cfg)
    em = pd.read_csv(
        ROOT / "product/execute_hiphop/results/exec_matrix.tsv", sep="\t")
    smi_of = dict(zip(em["inchikey"], em["smiles"]))
    smis = [smi_of.get(i, "") for i in df["inchikey"]]
    seqs = [seq_of.get(a, "") for a in df["acc"]]
    ok = [bool(s) and bool(q) for s, q in zip(smis, seqs)]
    P = predict_all([s for s, o in zip(smis, ok) if o],
                    [q for q, o in zip(seqs, ok) if o])
    ens = None
    k = 0
    for name_m, p in P.items():
        a1, b1 = cal["rescale"][name_m]
        full = np.full(len(df), np.nan)
        full[np.array(ok)] = a1 * p + b1
        df[f"pd_raw_{name_m}"] = full
        ens = full if ens is None else ens + np.nan_to_num(full)
        k += 1
    df["pd_pic50_calibrated"] = ens / k
    m = df["pd_pic50_calibrated"].notna()
    yy = df.loc[m, "pchembl_measured"].astype(float)
    pp = df.loc[m, "pd_pic50_calibrated"]
    df["pd_gate_calibrated"] = df["pd_pic50_calibrated"] >= cfg["pic50_gate"]
    df["final_pass_calibrated"] = df["dock_gate"] & df["pd_gate_calibrated"]
    df.to_csv(res / "final_pairs_calibrated.tsv", sep="\t", index=False)
    summary = dict(
        n=int(m.sum()),
        mae=round(float(np.mean(np.abs(pp - yy))), 2),
        r_internal_only=round(float(np.corrcoef(pp, yy)[0, 1]), 3),
        n_final=int(df["final_pass_calibrated"].sum()))
    (res / "apply_calibrated_summary.json").write_text(
        json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/pd_validation.json")
    ap.add_argument("--stage", choices=["calibrate", "apply"],
                    default="apply")
    ap.add_argument("--per-target-limit", type=int, default=200)
    args = ap.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())
    if args.stage == "calibrate":
        calibrate(cfg, args.per_target_limit)
    else:
        apply_stage(cfg)


if __name__ == "__main__":
    main()
