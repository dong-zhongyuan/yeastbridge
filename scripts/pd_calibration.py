#!/usr/bin/env python
"""pd_validation calibration: leakage-aware holdout, per-model performance,
weighted ensemble + linear recalibration, and re-scoring of the delivered
pairs.

Holdout: ChEMBL quantitative activities for the branch targets (GPCR/ion
channel, rich data), pChEMBL not null, dedup per (molecule, target) by max.
Caveat registered in DESIGN: DeepPurpose was pretrained on BindingDB which
shares compounds with ChEMBL, so holdout metrics may be mildly optimistic;
they calibrate the ensemble, not certify zero-shot performance.
Steps: fetch activities -> batch molecule SMILES -> 3-model predictions ->
metrics (r/RMSE/MAE) -> least-squares model weights + linear recalibration ->
apply to the delivered pairs (final_pairs.tsv) -> calibrated table.
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/pd_validation.json")
    ap.add_argument("--per-target-limit", type=int, default=200)
    args = ap.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())

    import numpy as np
    import pandas as pd
    from DeepPurpose import utils, DTI

    res = ROOT / cfg["results_dir"]
    inp = ROOT / cfg["inputs_dir"]
    inp.mkdir(parents=True, exist_ok=True)
    res.mkdir(parents=True, exist_ok=True)
    base = cfg["chembl_base"]

    # --- 1. holdout activities per branch target ---
    bp = pd.read_csv(
        ROOT / "product/chembl_branch/results/branch_pairs.tsv", sep="\t")
    tcs = sorted(bp["target_chembl_id"].unique())
    acc_of = dict(zip(bp["target_chembl_id"], bp["acc"]))
    rows = []
    for j, tc in enumerate(tcs, 1):
        a = get_json(f"{base}/activity.json?target_chembl_id={tc}"
                     f"&pchembl_value__isnull=false&limit="
                     f"{args.per_target_limit}")
        for act in (a or {}).get("activities", []):
            rows.append(dict(target_chembl_id=tc,
                             molecule=act.get("molecule_chembl_id"),
                             pchembl=float(act["pchembl_value"])))
        time.sleep(0.2)
        if j % 10 == 0:
            print(f"  activities {j}/{len(tcs)} targets", flush=True)
    ho = pd.DataFrame(rows).dropna()
    ho = (ho.sort_values("pchembl", ascending=False)
            .drop_duplicates(["target_chembl_id", "molecule"]))
    print(f"holdout pairs: {len(ho)}", flush=True)

    # --- 2. molecule SMILES in batches ---
    mols = sorted(ho["molecule"].unique())
    smi_of, ik_of = {}, {}
    for s in range(0, len(mols), 50):
        chunk = mols[s:s + 50]
        q = ",".join(chunk)
        rec = get_json(f"{base}/molecule.json?molecule_chembl_id__in={q}"
                       f"&limit=50")
        for m in (rec or {}).get("molecules", []):
            ms = m.get("molecule_structures") or {}
            if ms.get("canonical_smiles"):
                smi_of[m["molecule_chembl_id"]] = ms["canonical_smiles"]
                ik_of[m["molecule_chembl_id"]] = ms.get(
                    "standard_inchi_key", "")
        time.sleep(0.2)
        if (s // 50) % 10 == 9:
            print(f"  molecules {s + 50}/{len(mols)}", flush=True)
    ho = ho[ho["molecule"].isin(smi_of)].reset_index(drop=True)
    ho["smiles"] = ho["molecule"].map(smi_of)
    ho["acc"] = ho["target_chembl_id"].map(acc_of)
    ho.to_csv(inp / "calibration_holdout.tsv", sep="\t", index=False)
    print(f"holdout with SMILES: {len(ho)} | "
          f"pchembl range {ho.pchembl.min():.1f}-{ho.pchembl.max():.1f}",
          flush=True)

    seq_of = {}
    name = None
    for line in (ROOT / cfg["universe_fasta"]).read_text().splitlines():
        if line.startswith(">"):
            name = line.split("|")[1]
        elif name:
            seq_of.setdefault(name, "")
            seq_of[name] += line.strip()

    y = ho["pchembl"].to_numpy()
    P = {}
    for drug_enc, tgt_enc, name_m in MODELS:
        path = utils.download_pretrained_model(name_m)
        model = DTI.model_pretrained(path)
        X = utils.data_process(ho["smiles"].tolist(),
                               [seq_of[a] for a in ho["acc"]],
                               [0] * len(ho), drug_enc, tgt_enc,
                               split_method="no_split")
        P[name_m] = np.array([float(v) for v in model.predict(X)])
        r = float(np.corrcoef(P[name_m], y)[0, 1])
        rmse = float(np.sqrt(np.mean((P[name_m] - y) ** 2)))
        print(f"{name_m}: r={r:.3f} rmse={rmse:.2f}", flush=True)
    Xw = np.column_stack([P[m[2]] for m in MODELS])
    w, *_ = np.linalg.lstsq(Xw, y, rcond=None)
    ens = Xw @ w
    a1, b1 = np.polyfit(ens, y, 1)
    ens_cal = a1 * ens + b1
    metrics = dict(
        n_holdout=len(ho), pchembl_range=[float(y.min()), float(y.max())],
        per_model={m[2]: dict(r=round(float(np.corrcoef(P[m[2]], y)[0, 1]), 3),
                              rmse=round(float(np.sqrt(np.mean((P[m[2]] - y) ** 2))), 2))
                   for m in MODELS},
        ensemble_weighted=dict(r=round(float(np.corrcoef(ens, y)[0, 1]), 3),
                               rmse=round(float(np.sqrt(np.mean((ens - y) ** 2))), 2)),
        recalibration=dict(slope=round(float(a1), 3), intercept=round(float(b1), 3),
                           r=round(float(np.corrcoef(ens_cal, y)[0, 1]), 3),
                           rmse=round(float(np.sqrt(np.mean((ens_cal - y) ** 2))), 2)),
        weights={m[2]: round(float(wi), 3) for m, wi in zip(MODELS, w)},
        caveat="holdout shares source DB with DeepPurpose pretrained data; "
               "metrics may be mildly optimistic; used for ensemble "
               "calibration")

    # --- 3. re-score delivered pairs ---
    fp = pd.read_csv(res / "final_pairs.tsv", sep="\t")
    seqs = [seq_of.get(a, "") for a in fp["acc"]]
    em = pd.read_csv(
        ROOT / "product/execute_hiphop/results/exec_matrix.tsv", sep="\t")
    smi_of2 = dict(zip(em["inchikey"], em["smiles"]))
    smis = [smi_of2.get(i, "") for i in fp["inchikey"]]
    ok = [bool(s) and bool(q) for s, q in zip(smis, seqs)]
    Pf = {}
    for drug_enc, tgt_enc, name_m in MODELS:
        path = utils.download_pretrained_model(name_m)
        model = DTI.model_pretrained(path)
        X = utils.data_process(
            [s for s, o in zip(smis, ok) if o],
            [q for q, o in zip(seqs, ok) if o],
            [0] * sum(ok), drug_enc, tgt_enc, split_method="no_split")
        pred = np.array([float(v) for v in model.predict(X)])
        full = np.full(len(fp), np.nan)
        full[np.array(ok)] = pred
        Pf[name_m] = full
    Xf = np.column_stack([Pf[m[2]] for m in MODELS])
    fp["pd_pic50_calibrated"] = a1 * (Xf @ w) + b1
    m_ok = fp["pd_pic50_calibrated"].notna()
    if m_ok.sum() > 1:
        yy = fp.loc[m_ok, "pchembl_measured"].astype(float)
        pp = fp.loc[m_ok, "pd_pic50_calibrated"]
        metrics["delivered_pairs_calibrated"] = dict(
            n=int(m_ok.sum()),
            r=round(float(np.corrcoef(pp, yy)[0, 1]), 3),
            mae=round(float(np.mean(np.abs(pp - yy))), 2))
        raw = fp.loc[m_ok, "pd_pic50_mean"].astype(float)
        metrics["delivered_pairs_raw"] = dict(
            r=round(float(np.corrcoef(raw, yy)[0, 1]), 3),
            mae=round(float(np.mean(np.abs(raw - yy))), 2))
    fp.to_csv(res / "final_pairs_calibrated.tsv", sep="\t", index=False)
    (res / "calibration_report.json").write_text(
        json.dumps(metrics, indent=2))
    print(json.dumps(metrics, indent=2), flush=True)


if __name__ == "__main__":
    main()
