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
    ("Transformer", "CNN", "transformer_cnn_bindingdb"),
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
        try:
            path = utils.download_pretrained_model(name_m)
            model = DTI.model_pretrained(path)
            X = utils.data_process(smiles, seqs, [0] * len(smiles),
                                   drug_enc, tgt_enc,
                                   split_method="no_split")
            out[name_m] = np.array([float(v) for v in model.predict(X)])
            print(f"  {name_m}: done", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"  {name_m}: SKIPPED ({type(e).__name__}: "
                  f"{str(e)[:60]})", flush=True)
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
    holdout_path = inp / "calibration_holdout.tsv"
    if holdout_path.exists():
        ho = pd.read_csv(holdout_path, sep="\t")
        print(f"holdout cached: {len(ho)} pairs", flush=True)
    else:
        rows = []
        for j, tc in enumerate(sorted(acc_of), 1):
            a = get_json(f"{base}/activity.json?target_chembl_id={tc}"
                         f"&pchembl_value__isnull=false&limit="
                         f"{per_target_limit}")
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
            rec = get_json(f"{base}/molecule.json?"
                           f"molecule_chembl_id__in={q}&limit=50")
            for m in (rec or {}).get("molecules", []):
                ms = m.get("molecule_structures") or {}
                if ms.get("canonical_smiles"):
                    smi_of[m["molecule_chembl_id"]] = ms["canonical_smiles"]
            time.sleep(0.2)
        ho = ho[ho["molecule"].isin(smi_of)].reset_index(drop=True)
        ho["smiles"] = ho["molecule"].map(smi_of)
        ho["acc"] = ho["target_chembl_id"].map(acc_of)
        ho.to_csv(holdout_path, sep="\t", index=False)

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
    rs_of = {}
    ens = None
    for name_m, p in P.items():
        a1, b1 = np.polyfit(p, y, 1)
        calib[name_m] = [round(float(a1), 4), round(float(b1), 4)]
        rs = a1 * p + b1
        rs_of[name_m] = rs
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
    # non-negative least-squares weight optimization on the holdout
    from scipy.optimize import nnls

    names = sorted(rs_of)
    Xw = np.column_stack([rs_of[n] for n in names])
    w, _ = nnls(Xw, y)
    ens_w = Xw @ w
    report["weights"] = {n: round(float(wi), 4)
                         for n, wi in zip(names, w) if wi > 1e-6}
    report["ensemble_weighted"] = dict(
        r=round(float(np.corrcoef(ens_w, y)[0, 1]), 3),
        rmse=round(float(np.sqrt(np.mean((ens_w - y) ** 2))), 2),
        mae=round(float(np.mean(np.abs(ens_w - y))), 2))
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
    w_of = (cal.get("weights") or {}) if cfg.get("use_weights") else {}
    rescaled = {}
    for name_m, p in P.items():
        a1, b1 = cal["rescale"][name_m]
        full = np.full(len(df), np.nan)
        full[np.array(ok)] = a1 * p + b1
        df[f"pd_rescaled_{name_m}"] = full
        rescaled[name_m] = full
    wt = np.zeros(len(df))
    tot = 0.0
    for name_m, rs in rescaled.items():
        w = float(w_of.get(name_m, 0.0)) if w_of else 1.0
        wt = wt + w * np.nan_to_num(rs)
        tot += w
    df["pd_pic50_calibrated"] = wt / tot
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


def knn_stage(cfg):
    """Route 1 (literature-backed): per-target ECFP4 Tanimoto kNN potency.
    Leave-one-out evaluation on the calibration holdout, then applied to the
    delivered pairs. Delivered compounds are excluded from the reference set
    by ChEMBL id (leakage rule)."""
    import numpy as np
    import pandas as pd
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem

    res = ROOT / cfg["results_dir"]
    k = int(cfg.get("knn_k", 5))
    ho = pd.read_csv(
        ROOT / cfg["inputs_dir"] / "calibration_holdout.tsv", sep="\t")
    df = pd.read_csv(res / "final_pairs.tsv", sep="\t")

    ch = pd.read_csv(
        ROOT / "product/drug_annotation/results/chembl_targets.tsv",
        sep="\t").fillna("")
    cid_of = ch[ch["molecule_chembl_id"] != ""].groupby(
        "inchikey")["molecule_chembl_id"].first().to_dict()
    delivered_cids = {cid_of[i] for i in df["inchikey"] if i in cid_of}
    ho = ho[~ho["molecule"].isin(delivered_cids)].reset_index(drop=True)

    def fp_of(s):
        m = Chem.MolFromSmiles(s)
        return None if m is None else AllChem.GetMorganFingerprintAsBitVect(
            m, 2, nBits=2048)

    ho["fp"] = [fp_of(s) for s in ho["smiles"]]
    ho = ho[ho["fp"].notna()].reset_index(drop=True)
    refs = {}
    for tc, g in ho.groupby("target_chembl_id"):
        refs[tc] = (list(g["fp"]), g["pchembl"].to_numpy(),
                    list(g["molecule"]))

    def knn(qfp, ref, exclude_mol=None):
        fps, ys, mols = ref
        sims = np.array(DataStructs.BulkTanimotoSimilarity(qfp, fps))
        mask = np.array([m != exclude_mol for m in mols])
        sims = np.where(mask, sims, -1.0)
        order = np.argsort(sims)[::-1][:k]
        s = sims[order]
        if s.sum() <= 0:
            return np.nan, 0.0
        return float((s * ys[order]).sum() / s.sum()), float(s[0])

    preds, meas, nns = [], [], []
    for row in ho.itertuples():
        p, ns = knn(row.fp, refs[row.target_chembl_id],
                    exclude_mol=row.molecule)
        preds.append(p)
        meas.append(row.pchembl)
        nns.append(ns)
    preds, meas, nns = np.array(preds), np.array(meas), np.array(nns)
    loo = dict(
        n=int(len(preds)), k=k,
        r=round(float(np.corrcoef(preds, meas)[0, 1]), 3),
        rmse=round(float(np.sqrt(np.mean((preds - meas) ** 2))), 2),
        mae=round(float(np.mean(np.abs(preds - meas))), 2),
        median_nn_tanimoto=round(float(np.median(nns)), 3),
        excl_delivered_chembl_ids=len(delivered_cids))

    bp = pd.read_csv(
        ROOT / "product/chembl_branch/results/branch_pairs.tsv", sep="\t")
    tc_of = {}
    for r_ in bp.itertuples():
        tc_of.setdefault(r_.target_gene, r_.target_chembl_id)
    em = pd.read_csv(
        ROOT / "product/execute_hiphop/results/exec_matrix.tsv", sep="\t")
    smi_of = dict(zip(em["inchikey"], em["smiles"]))
    out_rows = []
    for r_ in df.itertuples():
        tc = tc_of.get(r_.target_gene)
        qfp = fp_of(smi_of.get(r_.inchikey, ""))
        if tc in refs and qfp is not None:
            p, ns = knn(qfp, refs[tc])
        else:
            p, ns = np.nan, 0.0
        out_rows.append(dict(inchikey=r_.inchikey, target_gene=r_.target_gene,
                             pd_pic50_knn=p, nn_tanimoto=ns,
                             n_ref=len(refs.get(tc, ([], [], []))[0])))
    ko = pd.DataFrame(out_rows)
    ko.to_csv(res / "final_pairs_knn.tsv", sep="\t", index=False)
    m = ko["pd_pic50_knn"].notna()
    if m.sum() > 1:
        yy = df.loc[m.values, "pchembl_measured"].astype(float)
        pp = ko.loc[m, "pd_pic50_knn"]
        loo["delivered_pairs"] = dict(
            n=int(m.sum()),
            mae=round(float(np.mean(np.abs(pp - yy))), 2),
            r_internal_only=round(float(np.corrcoef(pp, yy)[0, 1]), 3))
    (res / "knn_report.json").write_text(json.dumps(loo, indent=2))
    print(json.dumps(loo, indent=2), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/pd_validation.json")
    ap.add_argument("--stage", choices=["calibrate", "apply", "knn"],
                    default="apply")
    ap.add_argument("--per-target-limit", type=int, default=200)
    args = ap.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())
    if args.stage == "calibrate":
        calibrate(cfg, args.per_target_limit)
    elif args.stage == "knn":
        knn_stage(cfg)
    else:
        apply_stage(cfg)


if __name__ == "__main__":
    main()
