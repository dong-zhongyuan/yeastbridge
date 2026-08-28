#!/usr/bin/env python
"""target_screen stage 3: receptor preparation (clean -> fpocket -> PDBQT).

gemmi strips the raw structure to first-model polymer-only (waters,
heteroatoms, alternate conformations removed; all polymer chains kept);
fpocket (default parameters) defines pockets, top-N by Score retained;
OpenBabel writes the pH-protonated rigid receptor PDBQT. Pocket boxes span
the pocket-lining atoms plus padding (per-axis minimum enforced).
AlphaFold-sourced receptors additionally require mean pocket-residue
pLDDT >= threshold (pLDDT read from the cleaned PDB B-factor column).
Resumable per target via structures/statuses.tsv. Registered in
product/target_screen/DESIGN.md; parameters in configs/target_screen.json.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def bin_of(name):
    cand = Path(sys.executable).parent / name
    return str(cand) if cand.exists() else name


def clean_polymer(src: Path, dst: Path) -> None:
    import gemmi

    st = gemmi.read_structure(str(src))
    st.remove_alternative_conformations()
    if len(st) > 1:
        del st[1:]
    out = gemmi.Structure()
    out.name = st.name
    for model in st:
        m2 = gemmi.Model(model.name)
        for ch in model:
            c2 = gemmi.Chain(ch.name)
            for res in ch:
                info = gemmi.find_tabulated_residue(res.name)
                if info is not None and (info.is_amino_acid()
                                         or info.is_nucleic_acid()):
                    c2.add_residue(res)
            if len(c2):
                m2.add_chain(c2)
        out.add_model(m2)
    out.setup_entities()
    out.write_pdb(str(dst))


def mean_plddt_by_residue(pdb: Path) -> dict:
    out = {}
    for line in pdb.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and line[17:20] != "STP":
            key = (line[21], line[22:27])
            cur = out.setdefault(key, [0.0, 0])
            cur[0] += float(line[60:66])
            cur[1] += 1
    return {k: v[0] / v[1] for k, v in out.items()}


def parse_pqr_atoms(path: Path):
    atoms = []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and line[17:20] != "STP":
            atoms.append((line[21], line[22:27],
                          float(line[30:38]), float(line[38:46]),
                          float(line[46:54])))
    return atoms


def prep_target(row, cfg):
    acc = row["acc"]
    # structural ID includes pdb_id (state-aware: multiple structures per acc)
    raw_path = ROOT / row["path"]
    struct_id = raw_path.stem  # e.g. "P18089_6K41" or "P18089_af2"
    source = row["source"]
    base = ROOT / cfg["structures_dir"]
    work = base / "fpocket_out" / struct_id
    try:
        if not raw_path.exists():
            return struct_id, "failed", "raw_missing"
        work.mkdir(parents=True, exist_ok=True)
        clean = work / f"{struct_id}.pdb"
        if not clean.exists():
            clean_polymer(raw_path, clean)
        plddt = mean_plddt_by_residue(clean) if source == "af2" else None
        if not (work / f"{struct_id}_out" / f"{struct_id}_info.txt").exists():
            subprocess.run([bin_of("fpocket"), "-f", clean.name], cwd=work,
                           check=True, capture_output=True, timeout=3600)
        info = (work / f"{struct_id}_out" / f"{struct_id}_info.txt").read_text()
        scores = {}
        for m in re.finditer(
                r"Pocket\s+(\d+)\s*:\n(?:.*\n)*?\s*Score\s*:\s*([-\d.]+)",
                info):
            scores[int(m.group(1))] = float(m.group(2))
        order = sorted(scores, key=lambda k: -scores[k])[:cfg["fpocket_top_n"]]
        pad, mins = cfg["box_padding"], cfg["box_min_size"]
        kept = []
        for pk in order:
            pqr = work / f"{acc}_out" / "pockets" / f"pocket{pk}_atm.pdb"
            if not pqr.exists():
                continue
            atoms = parse_pqr_atoms(pqr)
            if len(atoms) < 10:
                continue
            site_plddt = None
            if plddt is not None:
                vals = [plddt.get((a[0], a[1])) for a in atoms]
                vals = [v for v in vals if v is not None]
                if not vals:
                    continue
                site_plddt = sum(vals) / len(vals)
                if site_plddt < cfg["plddt_site_min"]:
                    continue
            xs = [a[2] for a in atoms]
            ys = [a[3] for a in atoms]
            zs = [a[4] for a in atoms]
            kept.append({
                "pocket": pk, "score": scores[pk],
                "center": [(max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2,
                           (max(zs) + min(zs)) / 2],
                "size": [max(mins, max(xs) - min(xs) + 2 * pad),
                         max(mins, max(ys) - min(ys) + 2 * pad),
                         max(mins, max(zs) - min(zs) + 2 * pad)],
                "site_plddt": site_plddt,
            })
        if not kept:
            return struct_id, "failed", "no_valid_pocket"
        rec = base / "receptors" / f"{struct_id}.pdbqt"
        rec.parent.mkdir(parents=True, exist_ok=True)
        if not rec.exists():
            subprocess.run(
                [bin_of("obabel"), "-ipdb", str(clean), "-opdbqt", "-xr",
                 "-p", str(cfg["ph"]), "-O", str(rec)],
                check=True, capture_output=True, timeout=1800)
        pdir = base / "pockets"
        pdir.mkdir(parents=True, exist_ok=True)
        (pdir / f"{struct_id}.json").write_text(json.dumps(kept))
        return struct_id, "ok", f"{len(kept)}pockets"
    except Exception as e:  # noqa: BLE001
        return struct_id, "failed", f"{type(e).__name__}:{str(e)[:80]}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/target_screen.json")
    args = ap.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())

    import pandas as pd

    inv = pd.read_csv(ROOT / cfg["results_dir"] / "inventory.tsv", sep="\t")
    rows = inv.to_dict("records")
    st_path = ROOT / cfg["structures_dir"] / "statuses.tsv"
    done = set()
    if st_path.exists():
        for r in pd.read_csv(st_path, sep="\t").to_dict("records"):
            if r["status"] == "ok":
                done.add(r["acc"])  # now holds struct_id (acc_pdbId)
    # dedup by raw file stem (struct_id), not by acc
    seen = set()
    unique_rows = []
    for r in rows:
        sid = Path(r["path"]).stem
        if sid not in seen and sid not in done:
            seen.add(sid)
            unique_rows.append(r)
    todo = unique_rows
    print(f"targets: {len(rows)}, todo: {len(todo)}", flush=True)

    fresh = not st_path.exists()
    with st_path.open("a", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        if fresh:
            w.writerow(["acc", "status", "note"])
        with ProcessPoolExecutor(max_workers=cfg["max_workers"]) as ex:
            for n, (acc, status, note) in enumerate(
                    ex.map(partial(prep_target, cfg=cfg), todo), 1):
                w.writerow([acc, status, note])
                fh.flush()
                if n % 25 == 0:
                    print(f"  {n}/{len(todo)}", flush=True)
    print("FINISHED receptor prep", flush=True)


if __name__ == "__main__":
    main()
