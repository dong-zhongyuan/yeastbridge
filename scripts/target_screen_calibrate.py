#!/usr/bin/env python
"""target_screen stage 5: calibration and hit calling.

(a) Retrospective gate: registered benchmark drugs must rank their known
universe targets in the top fraction of that drug's full-universe ranking
(>= min_evaluable evaluable and >= min_pass within top_frac to certify the
screen). (b) Per-ligand empirical null: the bulk of each ligand's
universe-wide affinities is treated as a non-binding background
(median/MAD Gaussian), BH-FDR controls the tail. (c) Secondary non-gating
decoy check: MW-matched random ChEMBL molecules docked on a
family-stratified target sample quantify how extreme the benchmark
known-target scores are (top pocket only). Registered in
product/target_screen/DESIGN.md; parameters in configs/target_screen.json.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import urllib.parse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from target_screen_prep_ligands import get_json, make_pdbqt  # noqa: E402


def read_gene_map(fasta: Path):
    acc_of_gene, gene_of_acc = {}, {}
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            parts = line[1:].split("|")
            if len(parts) >= 3:
                toks = parts[2].split()
                gene = next((t[3:] for t in toks if t.startswith("GN=")),
                            toks[0].split("_")[0])
                acc_of_gene[gene] = parts[1]
                gene_of_acc[parts[1]] = gene
    return acc_of_gene, gene_of_acc


def bh(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * n / rank)
        q[i] = min(1.0, prev)
    return q


def phi(z):
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def _dock_decoy_job(job, dec_paths, cfg):
    from vina import Vina

    acc, rec, pocket = job
    try:
        v = Vina(cpu=cfg["vina"]["cpu_per_job"])
        v.set_receptor(rec)
        v.compute_vina_maps(center=pocket["center"], box_size=pocket["size"])
        scores = {}
        for p in dec_paths:
            try:
                v.set_ligand_from_file(p)
                v.dock(exhaustiveness=cfg["vina"]["exhaustiveness"],
                       n_poses=cfg["vina"]["num_modes"])
                scores[p] = float(v.energies(n_poses=1)[0][0])
            except Exception:  # noqa: BLE001
                scores[p] = None
        return acc, scores
    except Exception:  # noqa: BLE001
        return acc, {}


def run_decoys(cfg, cal, scores, acc_of_gene):
    import numpy as np
    import pandas as pd

    base = ROOT / cfg["structures_dir"]
    dec_dir = ROOT / cfg["inputs_dir"] / "decoys"
    dec_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(cal["decoy_seed"])

    fam = pd.read_csv(ROOT / cfg["universe_targets"], sep="\t")
    fam_col = "target_family" if "target_family" in fam.columns else None
    inv = pd.read_csv(ROOT / cfg["results_dir"] / "inventory.tsv", sep="\t")
    _, gene_of_acc = read_gene_map(ROOT / cfg["universe_fasta"])
    fam_of_acc = {a: (dict(zip(fam["target_id"], fam[fam_col]))
                      .get(gene_of_acc.get(a, ""), "NA") if fam_col else "NA")
                  for a in inv["acc"]}

    sample = set()
    accs = sorted(inv["acc"])
    if fam_col:
        groups = {}
        for a in accs:
            groups.setdefault(fam_of_acc[a], []).append(a)
        total = len(accs)
        for _, lst in groups.items():
            k = max(1, round(cal["decoy_target_sample"] * len(lst) / total))
            sample.update(rng.choice(lst, size=min(k, len(lst)),
                                     replace=False).tolist())
    else:
        sample.update(rng.choice(
            accs, size=min(cal["decoy_target_sample"], len(accs)),
            replace=False).tolist())

    report = {}
    for b in cfg["benchmarks"]:
        lid = f"BENCH__{b['name']}"
        known_accs = [acc_of_gene[g] for g in b["known_targets"]
                      if g in acc_of_gene]
        sample_b = sorted(sample | set(known_accs))
        rec = get_json(f"{cfg['chembl_base']}/molecule.json?pref_name__"
                       f"iexact={urllib.parse.quote(b['name'])}&limit=1")
        mw = None
        if rec and rec.get("molecules"):
            mp = rec["molecules"][0].get("molecule_properties") or {}
            mw = mp.get("mw_free_base")
        if not mw:
            report[b["name"]] = {"skipped": "no_mw"}
            continue
        lo, hi = int(float(mw) * 0.75), int(float(mw) * 1.25)
        cand = get_json(f"{cfg['chembl_base']}/molecule.json?"
                        f"molecule_properties__mw_free_base__range={lo}-{hi}"
                        f"&limit=100")
        smis = [m["molecule_structures"]["canonical_smiles"]
                for m in (cand or {}).get("molecules", [])
                if m.get("molecule_structures")]
        if len(smis) < cal["decoys_per_benchmark"]:
            report[b["name"]] = {"skipped": f"only_{len(smis)}_candidates"}
            continue
        chosen = rng.choice(smis, size=cal["decoys_per_benchmark"],
                            replace=False).tolist()
        dec_paths = []
        for i, smi in enumerate(chosen):
            did = f"DECOY__{b['name']}__{i}"
            p = dec_dir / f"{did}.pdbqt"
            if not p.exists():
                try:
                    make_pdbqt(smi, cfg["ph"], cfg["conformer_seed"], p)
                except Exception:  # noqa: BLE001
                    continue
            if p.exists():
                dec_paths.append(str(p))
        jobs = []
        for a in sample_b:
            pj = base / "pockets" / f"{a}.json"
            rp = base / "receptors" / f"{a}.pdbqt"
            if pj.exists() and rp.exists():
                jobs.append((a, str(rp), json.loads(pj.read_text())[0]))
        with ProcessPoolExecutor(max_workers=cfg["max_workers"]) as ex:
            res = list(ex.map(partial(_dock_decoy_job, dec_paths=dec_paths,
                                      cfg=cfg), jobs))
        dec_vals = [v for _, sc in res for v in sc.values() if v is not None]
        sub = scores[(scores["lid"] == lid)
                     & (scores["acc"].isin(sample_b))]
        known_rows = sub[sub["acc"].isin(known_accs)]
        if not known_rows.empty and dec_vals:
            bench_aff = float(known_rows["affinity"].min())
            frac_worse = float(np.mean([d < bench_aff for d in dec_vals]))
            report[b["name"]] = {
                "n_decoys": len(dec_paths), "n_targets": len(jobs),
                "n_decoy_scores": len(dec_vals),
                "benchmark_best_known_affinity": bench_aff,
                "decoy_median": float(np.median(dec_vals)),
                "fraction_decoys_better": round(frac_worse, 4),
                "note": "fraction of decoy scores more negative than the "
                        "benchmark's best known-target affinity; low is good",
            }
        else:
            report[b["name"]] = {"skipped": "no_known_target_in_sample"}
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/target_screen.json")
    ap.add_argument("--skip-decoys", action="store_true")
    args = ap.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())
    cal = cfg["calibration"]

    import numpy as np
    import pandas as pd

    res = ROOT / cfg["results_dir"]
    scores = pd.read_csv(res / "target_scores.tsv", sep="\t")
    fam = pd.read_csv(ROOT / cfg["universe_targets"], sep="\t")
    fam_col = "target_family" if "target_family" in fam.columns else None
    acc_of_gene, _ = read_gene_map(ROOT / cfg["universe_fasta"])

    report = {}

    bench_rows = []
    for b in cfg["benchmarks"]:
        lid = f"BENCH__{b['name']}"
        sub = scores[scores["lid"] == lid].sort_values("affinity")
        if sub.empty:
            bench_rows.append({**b, "evaluable": False,
                               "reason": "no_scores"})
            continue
        n = len(sub)
        rank_of = {a: i + 1 for i, a in enumerate(sub["acc"])}
        fr = {g: rank_of[acc_of_gene[g]] / n
              for g in b["known_targets"] if acc_of_gene.get(g) in rank_of}
        if not fr:
            bench_rows.append({**b, "evaluable": False,
                               "reason": "known_targets_not_scored"})
            continue
        bench_rows.append({**b, "evaluable": True,
                           "best_rank_fraction": round(min(fr.values()), 5),
                           "per_target": {k: round(v, 5) for k, v in fr.items()}})
    n_eval = sum(1 for r in bench_rows if r.get("evaluable"))
    n_pass = sum(1 for r in bench_rows if r.get("evaluable")
                 and r["best_rank_fraction"] <= cal["top_frac"])
    report["benchmarks"] = bench_rows
    report["gate"] = {
        "rule": (f">={cal['min_evaluable']} evaluable and "
                 f">={cal['min_pass']} within top {cal['top_frac']:.0%}"),
        "n_evaluable": n_eval, "n_pass": n_pass,
        "passed": bool(n_eval >= cal["min_evaluable"]
                       and n_pass >= cal["min_pass"]),
    }

    hits = []
    for lid, sub in scores.groupby("lid"):
        s = sub["affinity"].to_numpy(float)
        med = float(np.median(s))
        mad = float(np.median(np.abs(s - med)))
        sd = 1.4826 * mad
        if sd < 1e-9:
            continue
        z = (s - med) / sd
        p = [phi(v) for v in z]
        q = bh(p)
        for (_, row), qi, pi, zi in zip(sub.iterrows(), q, p, z):
            if qi <= cal["fdr_alpha"]:
                hits.append({"lid": lid, "label": row["label"],
                             "source": row["source"], "acc": row["acc"],
                             "gene": row["gene"],
                             "affinity": float(row["affinity"]),
                             "z": round(float(zi), 2), "p": float(pi),
                             "fdr": float(qi)})
    hits_df = pd.DataFrame(hits)
    if fam_col and len(hits_df):
        fmap = dict(zip(fam["target_id"], fam[fam_col]))
        hits_df["family"] = hits_df["gene"].map(fmap)
    if len(hits_df):
        hits_df = hits_df.sort_values(["fdr", "affinity"])
    out_cols = [c for c in ["lid", "label", "source", "acc", "gene", "family",
                            "affinity", "z", "p", "fdr"] if c in hits_df]
    hits_df.to_csv(res / "screen_hits.tsv", sep="\t", index=False,
                   columns=out_cols)
    report["empirical_null"] = {
        "method": "per-ligand median/MAD Gaussian background, BH-FDR",
        "alpha": cal["fdr_alpha"], "n_hits": int(len(hits_df)),
        "n_ligands": int(scores["lid"].nunique()),
    }

    if not args.skip_decoys and cal.get("decoys_per_benchmark", 0):
        report["decoy_check"] = run_decoys(cfg, cal, scores, acc_of_gene)

    (res / "calibration_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps({"gate": report["gate"],
                      "n_hits": report["empirical_null"]["n_hits"]},
                     indent=2), flush=True)


if __name__ == "__main__":
    main()
