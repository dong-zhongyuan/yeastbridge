#!/usr/bin/env python
"""target_screen stage 4: inverse docking of all ligands vs the universe.

One job = one (target, pocket): build Vina maps for the pocket box once,
dock every prepared ligand, append per-ligand best affinities to a JSONL
file (resumable per ligand). results/target_scores.tsv aggregates the best
score over pockets for each ligand-target pair. Registered in
product/target_screen/DESIGN.md; parameters in configs/target_screen.json.
"""
from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read_gene_map(fasta: Path):
    gene_of = {}
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            parts = line[1:].split("|")
            if len(parts) >= 3:
                toks = parts[2].split()
                gene = next((t[3:] for t in toks if t.startswith("GN=")),
                            toks[0].split("_")[0])
                gene_of[parts[1]] = gene
    return gene_of


def load_ligands(cfg):
    import pandas as pd

    reg = pd.read_csv(ROOT / cfg["inputs_dir"] / "ligands" / "registry.tsv",
                      sep="\t").fillna("")
    reg = reg.drop_duplicates(subset="lid", keep="last")
    lig = []
    for r in reg.to_dict("records"):
        if r["status"] == "ok":
            p = ROOT / cfg["inputs_dir"] / "ligands" / f"{r['lid']}.pdbqt"
            if p.exists():
                lig.append({"lid": r["lid"], "path": str(p)})
    return lig


def pocket_jobs(cfg):
    base = ROOT / cfg["structures_dir"]
    n_screen = int(cfg.get("screen_pockets_per_target", 3))
    jobs = []
    for pj in sorted((base / "pockets").glob("*.json")):
        acc = pj.stem
        rec = base / "receptors" / f"{acc}.pdbqt"
        if not rec.exists():
            continue
        for pk in json.loads(pj.read_text())[:n_screen]:
            jobs.append({"acc": acc, "receptor": str(rec), "pocket": pk})
    return jobs


def run_job(job, lig, cfg, gene_of):
    from vina import Vina

    out = (ROOT / cfg["results_dir"] / "raw_pockets" /
           f"{job['acc']}__p{job['pocket']['pocket']}.jsonl")
    out.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if out.exists():
        for line in out.read_text().splitlines():
            try:
                done.add(json.loads(line)["lid"])
            except Exception:  # noqa: BLE001
                pass
    todo = [l for l in lig if l["lid"] not in done]
    if not todo:
        return 0
    v = Vina(cpu=cfg["vina"]["cpu_per_job"])
    v.set_receptor(job["receptor"])
    v.compute_vina_maps(center=job["pocket"]["center"],
                        box_size=job["pocket"]["size"])
    with out.open("a") as fh:
        for l in todo:
            row = {"acc": job["acc"], "gene": gene_of.get(job["acc"], ""),
                   "pocket": job["pocket"]["pocket"], "lid": l["lid"]}
            try:
                v.set_ligand_from_file(l["path"])
                v.dock(exhaustiveness=cfg["vina"]["exhaustiveness"],
                       n_poses=cfg["vina"]["num_modes"])
                row["affinity"] = float(v.energies(n_poses=1)[0][0])
            except Exception as e:  # noqa: BLE001
                row["error"] = f"{type(e).__name__}:{str(e)[:60]}"
            fh.write(json.dumps(row) + "\n")
        fh.flush()
    return len(todo)


def aggregate(cfg):
    import pandas as pd

    rows = []
    for f in (ROOT / cfg["results_dir"] / "raw_pockets").glob("*.jsonl"):
        for line in f.read_text().splitlines():
            try:
                rows.append(json.loads(line))
            except Exception:  # noqa: BLE001
                pass
    if not rows:
        print("no raw rows yet", flush=True)
        return
    df = pd.DataFrame(rows)
    if "affinity" not in df or df["affinity"].notna().sum() == 0:
        print("no scored rows yet", flush=True)
        return
    reg = pd.read_csv(ROOT / cfg["inputs_dir"] / "ligands" / "registry.tsv",
                      sep="\t").fillna("")
    reg = reg.drop_duplicates(subset="lid", keep="last")
    meta = reg.set_index("lid")[["label", "source"]].to_dict("index")
    df = df[df["affinity"].notna()]
    best = (df.sort_values("affinity")
              .groupby(["lid", "acc"], as_index=False).first())
    best["label"] = best["lid"].map(lambda x: meta.get(x, {}).get("label", ""))
    best["source"] = best["lid"].map(
        lambda x: meta.get(x, {}).get("source", ""))
    out = ROOT / cfg["results_dir"] / "target_scores.tsv"
    best[["lid", "label", "source", "acc", "gene", "pocket", "affinity"]] \
        .sort_values("affinity").to_csv(out, sep="\t", index=False)
    n_err = int(df["error"].notna().sum()) if "error" in df else 0
    (ROOT / cfg["results_dir"] / "docking_progress.json").write_text(
        json.dumps({"pairs": int(len(best)), "errors": n_err}))
    print(f"aggregated {len(best)} ligand-target pairs -> {out}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/target_screen.json")
    ap.add_argument("--aggregate-only", action="store_true")
    args = ap.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())

    gene_of = read_gene_map(ROOT / cfg["universe_fasta"])

    if not args.aggregate_only:
        lig = load_ligands(cfg)
        jobs = pocket_jobs(cfg)
        print(f"ligands: {len(lig)}, pocket jobs: {len(jobs)}", flush=True)
        if not jobs:
            print("no pocket jobs; run stage 3 first", flush=True)
            return
        with ProcessPoolExecutor(max_workers=cfg["max_workers"]) as ex:
            futs = [ex.submit(run_job, j, lig, cfg, gene_of) for j in jobs]
            for n, f in enumerate(futs, 1):
                f.result()
                if n % 20 == 0:
                    print(f"  {n}/{len(jobs)} pocket jobs finished", flush=True)
    aggregate(cfg)


if __name__ == "__main__":
    main()
