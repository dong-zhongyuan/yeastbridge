#!/usr/bin/env python
"""EXP3 three-arm conformation screening reproduction (re side).

Registered: feasibility/conformation_selection/REGISTRATION.md (2026-08-29).
Upstream protocol (read-only, frozen at yeastbridge_vs@7bdaf4a):
  scripts/screen_v1.py + scripts/screen_aggregate.py
  library/receptors/ligand pdbqt under yeastbridge_vs reports+data.
Engine substitution declared in the registration: Vina-GPU 2.1 instead of
CPU Vina 1.2.7; inference is at the conformation-preference level only.

Stages (all outputs under feasibility/conformation_selection/):
  setup      enumerate the frozen upstream pair list, resolve boxes
             (incl. 7CR0 Kabsch transfer), symlink ligand dirs (vs read-only)
  dock       run Vina-GPU per receptor conformation into the feasibility sandbox
  aggregate  arms A/B/C, reference-panel preference, exact binomial p,
             enrichment; writes results/RESULTS.md + results.json

Usage:
  python scripts/conformation_repro.py --stage setup
  python scripts/conformation_repro.py --stage dock      (nohup; hours)
  python scripts/conformation_repro.py --stage aggregate
Run setup/aggregate under /public/home/mengxl/dzy/envs/yeastbridge/bin/python.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

RE_ROOT = Path("/public/home/mengxl/dzy/yeastbridge_re")
VS_ROOT = Path("/public/home/mengxl/dzy/yeastbridge_vs")
OUT = RE_ROOT / "feasibility/conformation_selection"
WORK = OUT / "work"
RESULTS = OUT / "results"

SCREEN = VS_ROOT / "reports/competition_v1/screen_v1"
MATRIX = SCREEN / "results/score_matrix.jsonl"
LIB_JSONL = SCREEN / "library/screen_library_v1.jsonl"
STEMS_JSON = SCREEN / "library/ligand_stems.json"
LIG_DIR = SCREEN / "ligands"

RECEPTORS = {  # verbatim from vs scripts/screen_v1.py
    "lpar1": {
        "4Z34": ("lpar1_4z34_on7_v1", "receptor_chain_A"),
        "7TD0": ("lpar1_7td0_nkp_v1", "receptor_chains_ABGR"),
        "7YU3": ("lpar1_7yu3_k6l_v1", "receptor_chains_ABGRS"),
    },
    "kcnq2": {
        "7CR0": ("kcnq2_7cr0_apo_v1", "receptor_chains_ABCD"),
        "7CR1": ("kcnq2_7cr1_gb9_a_v1", "receptor_chains_ABCD"),
        "7CR2": ("kcnq2_7cr2_fbx_a_v1", "receptor_chains_ABCD"),
    },
}
DEFAULT_CONF = {"lpar1": "4Z34", "kcnq2": "7CR0"}
CONFS = {"lpar1": ["4Z34", "7TD0", "7YU3"], "kcnq2": ["7CR0", "7CR1", "7CR2"]}
REFERENCE_CLASS = {  # verbatim from vs scripts/screen_aggregate.py
    "REF_LPA": ("lpar1", ["7TD0", "7YU3"]),
    "REF_ONO0740556": ("lpar1", ["7TD0", "7YU3"]),
    "REF_ONO9780307": ("lpar1", ["4Z34"]),
    "REF_ZTZ240": ("kcnq2", ["7CR1", "7CR2"]),
    "REF_RETIGABINE": ("kcnq2", ["7CR1", "7CR2"]),
}
SUPERPOSE_RMSD_GATE = 3.0


def _sha256(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _boxes() -> dict[str, dict]:
    boxes: dict[str, dict] = {}
    for target, conformations in RECEPTORS.items():
        for pdb_id, (run_dir, _stem) in conformations.items():
            manifest = json.loads(
                (VS_ROOT / "data/prepared/redocking" / run_dir / "protocol_manifest.json").read_text()
            )
            box = manifest["docking"].get("box")
            if box is None and pdb_id == "7CR0":
                frozen = SCREEN / "results/7cr0_transferred_box.json"
                box = json.loads(frozen.read_text()) if frozen.exists() else _transfer_box_7cr0()
            boxes[f"{target}/{pdb_id}"] = box
    return boxes


def _transfer_box_7cr0() -> dict | None:
    """Verbatim port of vs screen_v1._transfer_box_7cr0 (Kabsch on CA atoms)."""
    import numpy as np

    def ca_positions(path: Path) -> dict[int, tuple[float, float, float]]:
        positions: dict[int, tuple[float, float, float]] = {}
        for line in path.read_text(errors="replace").splitlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA" and len(line) >= 54:
                try:
                    num = int(line[22:26])
                    xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                except ValueError:
                    continue
                positions.setdefault(num, xyz)
        return positions

    source = VS_ROOT / "data/prepared/redocking/kcnq2_7cr1_gb9_a_v1/receptor_chains_ABCD_extracted.pdb"
    target = VS_ROOT / "data/prepared/redocking/kcnq2_7cr0_apo_v1/receptor_chains_ABCD_extracted.pdb"
    src_manifest = json.loads(
        (VS_ROOT / "data/prepared/redocking/kcnq2_7cr1_gb9_a_v1/protocol_manifest.json").read_text()
    )
    source_box = src_manifest["docking"]["box"]
    fixed_ca = ca_positions(target)
    mobile_ca = ca_positions(source)
    common = sorted(set(fixed_ca) & set(mobile_ca))
    if len(common) < 3:
        return None
    fixed = np.array([fixed_ca[n] for n in common], dtype=np.float64)
    mobile = np.array([mobile_ca[n] for n in common], dtype=np.float64)
    mobile_center = mobile.mean(axis=0)
    fixed_center = fixed.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (fixed - fixed_center)
    u, _s, vt = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[2, 2] = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ correction @ u.T
    translation = fixed_center - rotation @ mobile_center
    fitted = (rotation @ (mobile - mobile_center).T).T + fixed_center
    rmsd = float(np.sqrt(((fitted - fixed) ** 2).sum(axis=1).mean()))
    if rmsd > SUPERPOSE_RMSD_GATE:
        return None
    center = source_box["center_angstrom"]
    moved = rotation @ np.array([center["x"], center["y"], center["z"]]) + translation
    return {
        "center_angstrom": {
            "x": round(float(moved[0]), 3),
            "y": round(float(moved[1]), 3),
            "z": round(float(moved[2]), 3),
        },
        "size_angstrom": source_box["size_angstrom"],
        "transfer_note": f"kabsch_ca_superposition rmsd={rmsd:.3f}",
    }


def receptor_pdbqt(target: str, pdb_id: str) -> Path:
    run_dir, stem = RECEPTORS[target][pdb_id]
    return VS_ROOT / "data/prepared/redocking" / run_dir / f"{stem}.pdbqt"


def stage_setup() -> None:
    stems = json.loads(STEMS_JSON.read_text())
    roles = {}
    for line in LIB_JSONL.open(encoding="utf-8"):
        row = json.loads(line)
        roles[row["id"]] = row["role"]
    pairs = []
    missing_stem = []
    for line in MATRIX.open(encoding="utf-8"):
        row = json.loads(line)
        lid, receptor = row["id"], row["receptor"]
        stem = stems.get(lid)
        if stem is None or not (LIG_DIR / f"{stem}.pdbqt").exists():
            missing_stem.append(lid)
            continue
        pairs.append({"id": lid, "stem": stem, "receptor": receptor, "role": roles.get(lid, row.get("role", ""))})
    boxes = _boxes()
    for key, box in boxes.items():
        if box is None:
            sys.exit(f"box resolution failed for {key}")
    WORK.mkdir(parents=True, exist_ok=True)
    counts = {}
    for pair in pairs:
        conf_dir = WORK / "ligdirs" / pair["receptor"].replace("/", "_")
        conf_dir.mkdir(parents=True, exist_ok=True)
        dst = conf_dir / f"{pair['stem']}.pdbqt"
        src = LIG_DIR / f"{pair['stem']}.pdbqt"
        if not dst.exists():
            dst.symlink_to(src)
        counts[pair["receptor"]] = counts.get(pair["receptor"], 0) + 1
    receptor_hashes = {
        f"{t}/{p}": _sha256(receptor_pdbqt(t, p)) for t in RECEPTORS for p in RECEPTORS[t]
    }
    setup = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "upstream_matrix_sha256": _sha256(MATRIX),
        "library_sha256": _sha256(LIB_JSONL),
        "n_pairs_enumerated": len(pairs),
        "n_missing_ligand_pdbqt": len(missing_stem),
        "pairs_per_receptor": counts,
        "boxes": boxes,
        "receptor_pdbqt_sha256": receptor_hashes,
    }
    (OUT / "setup.json").write_text(json.dumps(setup, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with (OUT / "pairs.jsonl").open("w", encoding="utf-8") as handle:
        for pair in pairs:
            handle.write(json.dumps(pair) + "\n")
    print(json.dumps({"stage": "setup", **{k: setup[k] for k in ("n_pairs_enumerated", "n_missing_ligand_pdbqt", "pairs_per_receptor")}}, indent=2))


def _pick_gpu() -> int:
    proc = subprocess.run(
        ["nvidia-smi", "--query-gpu=index,memory.used", "--format=csv,noheader,nounits"],
        capture_output=True, text=True, check=True,
    )
    gpus = [(int(i), int(m)) for i, m in (ln.split(", ") for ln in proc.stdout.strip().splitlines())]
    return min(gpus, key=lambda g: g[1])[0]


def stage_dock() -> None:
    chembl_cfg = json.loads((RE_ROOT / "configs/chembl_branch.json").read_text())
    vg = chembl_cfg["vina_gpu"]
    boxes = json.loads((OUT / "setup.json").read_text())["boxes"]
    gpu = _pick_gpu()
    env = dict(os.environ)
    env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    env["LD_LIBRARY_PATH"] = "/public/home/mengxl/dzy/envs/gpubuild/lib:/usr/local/cuda-12.4/lib64"
    scores_path = RESULTS / "scores.tsv"
    RESULTS.mkdir(parents=True, exist_ok=True)
    done_receptors = set()
    if scores_path.exists():
        for line in scores_path.read_text().splitlines()[1:]:
            done_receptors.add(line.split("\t")[1])
    print(f"using GPU {gpu}", flush=True)
    for target in RECEPTORS:
        for pdb_id in RECEPTORS[target]:
            key = f"{target}/{pdb_id}"
            if key in done_receptors:
                print(f"skip {key} (already scored)", flush=True)
                continue
            box = boxes[key]
            c, s = box["center_angstrom"], box["size_angstrom"]
            ligdir = WORK / "ligdirs" / key.replace("/", "_")
            outdir = RESULTS / "gpu_out" / key.replace("/", "_")
            outdir.mkdir(parents=True, exist_ok=True)
            cfgfile = outdir / "vinagpu.cfg"
            cfgfile.write_text(
                f"receptor = {receptor_pdbqt(target, pdb_id)}\n"
                f"ligand_directory = {ligdir}\n"
                f"output_directory = {outdir}\n"
                f"opencl_binary_path = {vg['opencl_binary_path']}\n"
                f"center_x = {c['x']:.3f}\ncenter_y = {c['y']:.3f}\ncenter_z = {c['z']:.3f}\n"
                f"size_x = {s['x']:.3f}\nsize_y = {s['y']:.3f}\nsize_z = {s['z']:.3f}\n"
                f"thread = {vg['thread']}\n"
            )
            t0 = time.time()
            try:
                subprocess.run(
                    [vg["bin"], "--config", str(cfgfile)],
                    env=env, cwd=outdir, timeout=vg.get("run_timeout", 21600),
                    check=True, capture_output=True,
                )
            except Exception as exc:  # record and continue with other receptors
                (outdir / "error.log").write_text(f"{type(exc).__name__}: {exc}\n", encoding="utf-8")
                print(f"{key}: FAILED {type(exc).__name__}", flush=True)
                continue
            remark = re.compile(r"REMARK VINA RESULT:\s*(-?[\d.]+)")
            n = 0
            fresh = not scores_path.exists()
            with scores_path.open("a", encoding="utf-8") as fh:
                if fresh:
                    fh.write("ligand_id\treceptor\tscore\n")
                for f in sorted(outdir.glob("*_out.pdbqt")):
                    stem = f.name[: -len("_out.pdbqt")]
                    m = remark.search(f.read_text(errors="ignore"))
                    if m:
                        fh.write(f"{stem}\t{key}\t{m.group(1)}\n")
                        n += 1
            print(f"{key}: {n} ligands scored in {time.time()-t0:.0f}s", flush=True)
    print("DOCK_STAGE_DONE", flush=True)


def _rank_metrics(labels: list[int], scores: list[float]) -> dict:
    """EF1% + ROC-AUC (label 1 = active, lower Vina score = better)."""
    n = len(labels)
    order = sorted(range(n), key=lambda i: scores[i])
    n_act = sum(labels)
    k = max(1, math.ceil(0.01 * n))
    top = order[:k]
    hits_top = sum(labels[i] for i in top)
    ef1 = (hits_top / k) / (n_act / n) if n_act and n else 0.0
    # rank-based AUC with tie handling on scores
    auc_num, auc_den = 0.0, 0.0
    pos = [s for s, l in zip(scores, labels) if l == 1]
    neg = [s for s, l in zip(scores, labels) if l == 0]
    for p in pos:
        for q in neg:
            auc_den += 1
            auc_num += 1.0 if p < q else (0.5 if p == q else 0.0)
    auc = auc_num / auc_den if auc_den else 0.5
    return {"n": n, "actives": n_act, "ef1_percent": round(ef1, 3), "roc_auc": round(auc, 4)}


def _binomial_p_at_least(k: int, n: int, p: float = 0.5) -> float:
    return sum(math.comb(n, i) * (p ** i) * ((1 - p) ** (n - i)) for i in range(k, n + 1))


def stage_aggregate() -> None:
    import csv

    table: dict[str, dict[str, float]] = {}
    roles: dict[str, str] = {}
    with (RESULTS / "scores.tsv").open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            table.setdefault(row["ligand_id"], {})[row["receptor"]] = float(row["score"])

    reference_rows = []
    for ref_id, (target, correct_confs) in REFERENCE_CLASS.items():
        entry = table.get(ref_id, {})
        all_confs = CONFS[target]
        per_conf = {c: entry.get(f"{target}/{c}") for c in all_confs}
        correct_scores = [per_conf[c] for c in correct_confs if per_conf.get(c) is not None]
        incorrect_scores = [
            per_conf[c] for c in all_confs if c not in correct_confs and per_conf.get(c) is not None
        ]
        if not correct_scores or not incorrect_scores:
            continue
        best_correct, best_incorrect = min(correct_scores), min(incorrect_scores)
        argmin_conf = min(all_confs, key=lambda c: per_conf[c] if per_conf.get(c) is not None else math.inf)
        reference_rows.append(
            {
                "id": ref_id,
                "target": target,
                "correct_confs": correct_confs,
                "per_conf_scores": per_conf,
                "delta_best_correct_minus_incorrect": round(best_correct - best_incorrect, 3),
                "B_preference_correct": best_correct < best_incorrect,
                "C_top_conf": argmin_conf,
                "C_top_conf_correct": argmin_conf in correct_confs,
                "A_default_conf": DEFAULT_CONF[target],
                "A_default_correct": DEFAULT_CONF[target] in correct_confs,
            }
        )
    n_ref = len(reference_rows)
    acc = lambda key: sum(r[key] for r in reference_rows) / max(n_ref, 1)  # noqa: E731
    pref = {
        "rows": reference_rows,
        "armA_top_conf_accuracy": acc("A_default_correct"),
        "armC_top_conf_accuracy": acc("C_top_conf_correct"),
        "armB_delta_accuracy": acc("B_preference_correct"),
    }
    p_armC = _binomial_p_at_least(round(pref["armC_top_conf_accuracy"] * n_ref), n_ref)
    p_armB = _binomial_p_at_least(round(pref["armB_delta_accuracy"] * n_ref), n_ref)

    # enrichment per target per arm (A: default conf; C: best over all confs)
    lib_role = {}
    for line in LIB_JSONL.open(encoding="utf-8"):
        row = json.loads(line)
        lib_role[row["id"]] = row["role"]
    enrichment = {}
    for target in ("lpar1", "kcnq2"):
        ids = [i for i, r in lib_role.items() if r in ("active_" + target, "decoy") and i in table]
        for arm in ("A", "C"):
            labels, scores = [], []
            for i in ids:
                entry = table[i]
                if arm == "A":
                    v = entry.get(f"{target}/{DEFAULT_CONF[target]}")
                else:
                    vals = [entry.get(f"{target}/{c}") for c in CONFS[target]]
                    v = min([x for x in vals if x is not None], default=None)
                if v is None:
                    continue
                labels.append(1 if lib_role[i] == "active_" + target else 0)
                scores.append(v)
            enrichment[f"{target}_arm{arm}"] = _rank_metrics(labels, scores)

    upstream = {"armA": 0.2, "armC": 1.0, "armB": 1.0, "sign_p": 0.031}
    gate_pass = p_armB < 0.05 or p_armC < 0.05
    output = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "engine": "Vina-GPU 2.1 (substitution declared in REGISTRATION.md)",
        "n_ligands_scored": len(table),
        "reference_preference": pref,
        "binomial_p": {"armC": round(p_armC, 5), "armB": round(p_armB, 5)},
        "enrichment_secondary": enrichment,
        "upstream_reference": upstream,
        "gate": {
            "rule": "reference-panel conformation preference, exact binomial p < 0.05 (registered)",
            "pass": bool(gate_pass),
            "claim_boundary": "conformation-preference level only (same boundary as upstream EXP3)",
        },
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "results.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# 构象层复现结果（EXP3 re 版）",
        "",
        f"引擎：Vina-GPU 2.1（注册声明的替换；仅构象偏好层面推断）。打分配体数：{len(table)}。",
        "",
        "## 参照面板构象偏好（主判定）",
        "",
        "| 参照 | 靶点 | 正确构象 | 各构象得分 | delta | B 正确 | C 顶构象 | C 正确 | A 正确 |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in reference_rows:
        per = ", ".join(f"{c}={r['per_conf_scores'][c]}" for c in CONFS[r["target"]])
        lines.append(
            f"| {r['id']} | {r['target']} | {'/'.join(r['correct_confs'])} | {per} | "
            f"{r['delta_best_correct_minus_incorrect']} | {'✓' if r['B_preference_correct'] else '✗'} | "
            f"{r['C_top_conf']} | {'✓' if r['C_top_conf_correct'] else '✗'} | {'✓' if r['A_default_correct'] else '✗'} |"
        )
    lines += [
        "",
        f"- armA 默认构象正确率：{pref['armA_top_conf_accuracy']:.0%}（上游 20%）",
        f"- armC 顶构象正确率：{pref['armC_top_conf_accuracy']:.0%}（上游 100%），精确二项 p = {p_armC:.4f}",
        f"- armB delta 正确率：{pref['armB_delta_accuracy']:.0%}（上游 100%），精确二项 p = {p_armB:.4f}",
        f"- **Gate：{'PASS' if gate_pass else 'FAIL'}**（规则：参照面板偏好 p < 0.05；边界措辞与上游一致：仅构象偏好层面）",
        "",
        "## 次级富集指标（不作为判定）",
        "",
        "| 组 | n | actives | EF1% | ROC-AUC |",
        "|---|---|---|---|---|",
    ]
    for k2, v in enrichment.items():
        lines.append(f"| {k2} | {v['n']} | {v['actives']} | {v['ef1_percent']} | {v['roc_auc']} |")
    (RESULTS / "RESULTS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "aggregate", "gate_pass": bool(gate_pass), "armB": pref["armB_delta_accuracy"], "armC": pref["armC_top_conf_accuracy"], "armA": pref["armA_top_conf_accuracy"]}))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("setup", "dock", "aggregate"), required=True)
    args = parser.parse_args()
    {"setup": stage_setup, "dock": stage_dock, "aggregate": stage_aggregate}[args.stage]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
