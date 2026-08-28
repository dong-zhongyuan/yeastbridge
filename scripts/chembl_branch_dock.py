#!/usr/bin/env python
"""chembl_branch stage 4: GPU inverse-docking driver (Vina-GPU 2.1).

Runs on the idle GPU (auto-selected by free memory; never evicts existing
processes) so the CPU fleet stays with the main target_screen line. For
each branch target (structures from stage 2/3 reuse) and its top-N fpocket
pockets: symlink that target's branch compounds into a per-target ligand
directory, write a Vina-GPU config, run one GPU invocation per pocket, and
parse the best affinity from each {ligand}_out.pdbqt. Resumable per
(acc, pocket) via marker files. Registered in product/chembl_branch/
DESIGN.md; parameters in configs/chembl_branch.json.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import resource
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REMARK = re.compile(r"REMARK VINA RESULT:\s+([-\d.]+)")


def _stack_ok():
    # Vina-GPU needs >= 8 MB of stack (README requirement)
    resource.setrlimit(resource.RLIMIT_STACK,
                       (64 << 20, resource.RLIM_INFINITY))


def pick_gpu():
    out = subprocess.run(
        ["nvidia-smi", "--query-gpu=index,memory.used",
         "--format=csv,noheader,nounits"],
        capture_output=True, text=True, check=True).stdout.strip().splitlines()
    gpus = [(int(i.split(",")[0]), int(i.split(",")[1])) for i in out]
    return min(gpus, key=lambda g: g[1])[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/chembl_branch.json")
    args = ap.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())
    vg = cfg["vina_gpu"]

    import pandas as pd

    pairs = pd.read_csv(ROOT / cfg["results_dir"] / "branch_pairs.tsv",
                        sep="\t")
    pairs = pairs[pairs["ligand_ready"]]
    base = ROOT / cfg["structures_dir"]
    res = ROOT / cfg["results_dir"]
    res.mkdir(parents=True, exist_ok=True)
    out_tsv = res / "gpu_dock.tsv"
    done = set()
    fresh = not out_tsv.exists()
    if not fresh:
        for r in pd.read_csv(out_tsv, sep="\t").to_dict("records"):
            done.add((r["acc"], r["pocket"]))

    gpu = vg["gpu_id"] if str(vg["gpu_id"]).isdigit() else pick_gpu()
    print(f"using GPU {gpu}; pairs: {len(pairs)}", flush=True)

    env = dict(os.environ)
    env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    lig_of = {}
    for r in pairs.itertuples():
        lig_of.setdefault(r.acc, []).append((r.inchikey, r.target_gene))

    rows_written = 0
    for acc in sorted(lig_of):
        rec = base / "receptors" / f"{acc}.pdbqt"
        pj = base / "pockets" / f"{acc}.json"
        if not rec.exists() or not pj.exists():
            continue
        for pocket in json.loads(pj.read_text())[:cfg["fpocket_top_n"]]:
            key = (acc, pocket["pocket"])
            if key in done:
                continue
            center = pocket["center"]
            size = [min(vg["box_max"], s) for s in pocket["size"]]
            ligdir = ROOT / cfg["inputs_dir"] / "gpu_ligdirs" / f"{acc}_p{pocket['pocket']}"
            ligdir.mkdir(parents=True, exist_ok=True)
            for ik, _g in lig_of[acc]:
                src = ROOT / cfg["main_ligands_dir"] / f"{ik}.pdbqt"
                dst = ligdir / f"{ik}.pdbqt"
                if src.exists() and not dst.exists():
                    dst.symlink_to(src)
            outdir = res / "gpu_out" / f"{acc}__p{pocket['pocket']}"
            outdir.mkdir(parents=True, exist_ok=True)
            cfgfile = res / "gpu_out" / f"{acc}__p{pocket['pocket']}.cfg"
            cfgfile.write_text(
                f"receptor = {rec}\n"
                f"ligand_directory = {ligdir}\n"
                f"output_directory = {outdir}\n"
                f"opencl_binary_path = {vg['opencl_binary_path']}\n"
                f"center_x = {center[0]:.3f}\ncenter_y = {center[1]:.3f}\n"
                f"center_z = {center[2]:.3f}\n"
                f"size_x = {size[0]:.3f}\nsize_y = {size[1]:.3f}\n"
                f"size_z = {size[2]:.3f}\n"
                f"thread = {vg['thread']}\n")
            try:
                subprocess.run([vg["bin"], "--config", str(cfgfile)],
                               env=env, cwd=outdir, timeout=vg["run_timeout"],
                               check=True, capture_output=True,
                               preexec_fn=_stack_ok)
            except Exception as e:  # noqa: BLE001
                print(f"  {acc} p{pocket['pocket']}: {type(e).__name__}",
                      flush=True)
                continue
            aff = {}
            for f in outdir.glob("*_out.pdbqt"):
                ik = f.name[:-len("_out.pdbqt")]
                m = REMARK.search(f.read_text(errors="ignore"))
                if m:
                    aff[ik] = float(m.group(1))
            genes = {ik: g for ik, g in lig_of[acc]}
            with out_tsv.open("a", newline="") as fh:
                if fresh:
                    fh.write("acc\tpocket\tinchikey\ttarget_gene\taffinity\n")
                    fresh = False
                for ik, a in aff.items():
                    fh.write(f"{acc}\t{pocket['pocket']}\t{ik}\t"
                             f"{genes.get(ik, '')}\t{a}\n")
                    rows_written += 1
                fh.flush()
            print(f"  {acc} p{pocket['pocket']}: {len(aff)} ligands scored",
                  flush=True)
    print(f"FINISHED gpu dock, new rows: {rows_written}", flush=True)


if __name__ == "__main__":
    main()
