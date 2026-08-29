#!/usr/bin/env python
"""Delivery-package docking self-consistency check (re side).

Registered: feasibility/docking_qualification/REGISTRATION.md (2026-08-29).
For each delivered package (YeastBridge_Mol001..004): re-dock the delivered
ligand (from ligand.smi, freshly prepared) into the receptor pocket extracted
from the delivered complex.pdb, and check the affinity against the delivered
scores.tsv value (tolerance 1.0 kcal/mol, pre-declared).

Reads product/delivery/v1 read-only; writes only under
feasibility/docking_qualification/.
Run under /public/home/mengxl/dzy/envs/yeastbridge/bin/python.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import time
from pathlib import Path

RE_ROOT = Path("/public/home/mengxl/dzy/yeastbridge_re")
DELIVERY = RE_ROOT / "product/delivery/v1"
OUT = RE_ROOT / "feasibility/docking_qualification/results"
DOCK_BIN = Path("/public/home/mengxl/dzy/envs/yeastbridge_vs_docking/bin")
OBABEL = DOCK_BIN / "obabel"
VINA = DOCK_BIN / "vina"
MEEKO = DOCK_BIN / "mk_prepare_ligand.py"
DOCK_PY = DOCK_BIN / "python"
TOLERANCE = 1.0
SEED = 20260829
EXHAUSTIVENESS = 8
NUM_MODES = 9
PACKAGES = ["YeastBridge_Mol001", "YeastBridge_Mol002", "YeastBridge_Mol003", "YeastBridge_Mol004"]


def run(cmd: list[str], timeout: int = 600) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


def pose_centroid_and_bbox(sdf: Path) -> tuple[list[float], list[float]]:
    xs, ys, zs = [], [], []
    lines = sdf.read_text(errors="replace").splitlines()
    if len(lines) < 4:
        raise ValueError(f"pose sdf too short: {sdf}")
    counts = lines[3]
    n_atoms = int(counts[0:3])
    for line in lines[4 : 4 + n_atoms]:
        if len(line) >= 30:
            xs.append(float(line[0:10])); ys.append(float(line[10:20])); zs.append(float(line[20:30]))
    center = [sum(v) / len(v) for v in (xs, ys, zs)]
    size = [max(b) - min(a) + 8.0 for a, b in ((xs, xs), (ys, ys), (zs, zs))]
    size = [max(s, 20.0) for s in size]
    return center, size


def delivered_affinity(scores_tsv: Path) -> float | None:
    lines = scores_tsv.read_text(encoding="utf-8-sig").splitlines()
    header = lines[0].split("\t")
    col = header.index("对接亲和能_kcal_mol")
    for row in lines[1:]:
        if row.strip():
            return float(row.split("\t")[col])
    return None


def check(pkg: str) -> dict:
    work = OUT / pkg
    work.mkdir(parents=True, exist_ok=True)
    prefix = f"{pkg}_"
    complex_pdb = DELIVERY / pkg / f"{prefix}complex.pdb"
    lig_smi = DELIVERY / pkg / f"{prefix}ligand.smi"
    pose_sdf = DELIVERY / pkg / f"{prefix}ligand_pose.sdf"
    scores_tsv = DELIVERY / pkg / f"{prefix}scores.tsv"
    smiles = lig_smi.read_text().splitlines()[0].split()[0]

    receptor_pdb = work / "receptor.pdb"
    with receptor_pdb.open("w") as fh:
        fh.write("\n".join(l for l in complex_pdb.read_text(errors="replace").splitlines() if l.startswith("ATOM")) + "\n")
    receptor_pdbqt = work / "receptor.pdbqt"
    p = run([str(OBABEL), str(receptor_pdb), "-O", str(receptor_pdbqt), "-xr", "-p"])
    if not receptor_pdbqt.exists():
        return {"package": pkg, "status": "receptor_prep_failed", "stderr": p.stderr[-300:]}

    lig3d = work / "ligand3d.sdf"
    p = run([str(OBABEL), f"-:{smiles}", "-O", str(lig3d), "--gen3d", "best"])
    if not lig3d.exists():
        return {"package": pkg, "status": "ligand_gen3d_failed", "stderr": p.stderr[-300:]}
    lig_pdbqt = work / "ligand.pdbqt"
    p = run([str(DOCK_PY), str(MEEKO), "-i", str(lig3d), "-o", str(lig_pdbqt), "--charge_model", "gasteiger"])
    if not lig_pdbqt.exists():
        return {"package": pkg, "status": "ligand_prep_failed", "stderr": p.stderr[-300:]}

    center, size = pose_centroid_and_bbox(pose_sdf)
    out_pdbqt = work / "redock_out.pdbqt"
    p = run(
        [str(VINA), "--receptor", str(receptor_pdbqt), "--ligand", str(lig_pdbqt),
         "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
         "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
         "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", str(NUM_MODES), "--seed", str(SEED),
         "--out", str(out_pdbqt)],
        timeout=1800,
    )
    m = re.search(r"^\s+1\s+(-?[\d.]+)", p.stdout, re.MULTILINE)
    if not m:
        return {"package": pkg, "status": "dock_failed", "stdout_tail": p.stdout[-300:], "stderr": p.stderr[-300:]}
    redocked = float(m.group(1))
    delivered = delivered_affinity(scores_tsv)

    # pose sanity: redocked best-pose centroid vs delivered pose centroid
    rx, ry, rz = [], [], []
    in_root = False
    for line in out_pdbqt.read_text(errors="replace").splitlines():
        if line.startswith("ROOT"):
            in_root = True
            continue
        if line.startswith("ENDROOT"):
            break
        if in_root and line.startswith(("ATOM", "HETATM")):
            try:
                rx.append(float(line[30:38])); ry.append(float(line[38:46])); rz.append(float(line[46:54]))
            except ValueError:
                continue
    centroid_dist = None
    if rx:
        cx = [sum(v) / len(v) for v in (rx, ry, rz)]
        centroid_dist = round(sum((a - b) ** 2 for a, b in zip(cx, center)) ** 0.5, 2)

    deviation = abs(redocked - delivered) if delivered is not None else None
    return {
        "package": pkg,
        "status": "ok",
        "redocked_affinity_kcal_mol": redocked,
        "delivered_affinity_kcal_mol": delivered,
        "abs_deviation": round(deviation, 3) if deviation is not None else None,
        "within_tolerance": deviation is not None and deviation <= TOLERANCE,
        "pose_centroid_distance_angstrom": centroid_dist,
        "box_center": [round(c, 2) for c in center],
        "box_size": [round(s, 2) for s in size],
        "seed": SEED,
        "exhaustiveness": EXHAUSTIVENESS,
    }


def main() -> int:
    results = [check(pkg) for pkg in PACKAGES]
    OUT.mkdir(parents=True, exist_ok=True)
    payload = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "tolerance_kcal_mol": TOLERANCE,
        "results": results,
        "all_within_tolerance": all(r.get("within_tolerance") for r in results),
    }
    (OUT / "delivery_redock_check.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
