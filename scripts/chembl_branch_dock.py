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

    # state-aware structure selection: load inventory to map acc -> structures
    inv = pd.read_csv(ROOT / cfg["results_dir"] / "inventory.tsv",
                      sep="\t").fillna("")
    # for each acc, find best structure per state
    struct_of = {}  # (acc, required_state) -> struct_id
    for acc_id, grp in inv.groupby("acc"):
        states = {}
        for _, r in grp.iterrows():
            sid = Path(r["path"]).stem if r["path"] else ""
            if not sid:
                continue
            state = r.get("conformational_state", "unannotated")
            res_val = float(r["resolution"]) if r["resolution"] else 999.0
            if state not in states or res_val < states[state][1]:
                # prefer PDB over AF2 within same state category
                if r["source"] == "pdb" or state not in states:
                    states[state] = (sid, res_val)
        struct_of[acc_id] = states

    def select_structure(acc_id, action_types):
        """Select structure whose state matches compound's direction.
        Rule: INHIBITOR/ANTAGONIST/BLOCKER -> inactive;
              AGONIST/ACTIVATOR -> active;
              no_data -> best resolution any state."""
        states = struct_of.get(acc_id, {})
        if not states:
            return None
        at = str(action_types).upper() if isinstance(
            action_types, str) and action_types else ""
        is_inhibitor = any(k in at for k in
                           ["INHIBITOR", "ANTAGONIST", "BLOCKER",
                            "NEGATIVE MODULATOR"])
        is_agonist = any(k in at for k in
                         ["AGONIST", "ACTIVATOR", "POSITIVE MODULATOR"])
        if is_inhibitor and "inactive" in states:
            return states["inactive"][0]
        if is_inhibitor and "af2_inactive_like" in states:
            return states["af2_inactive_like"][0]
        if is_agonist and "active" in states:
            return states["active"][0]
        # fallback: any structure (best resolution)
        best = min(states.values(), key=lambda x: x[1])
        return best[0]

    gpu = vg["gpu_id"] if str(vg["gpu_id"]).isdigit() else pick_gpu()
    print(f"using GPU {gpu}; pairs: {len(pairs)}", flush=True)

    env = dict(os.environ)
    env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    env["LD_LIBRARY_PATH"] = (
        "/public/home/mengxl/dzy/envs/gpubuild/lib:"
        "/usr/local/cuda-12.4/lib64")
    # group by (acc, action_type) for state-aware receptor selection
    lig_of = {}
    for r in pairs.itertuples():
        acts = getattr(r, "action_types", "")
        lig_of.setdefault(r.acc, []).append(
            (r.inchikey, r.target_gene, acts))

    rows_written = 0
    for acc in sorted(lig_of):
        # select structure per compound direction (may differ per compound)
        # group compounds by their selected structure
        struct_groups = {}
        for ik, gene, acts in lig_of[acc]:
            sid = select_structure(acc, acts)
            if sid is None:
                continue
            struct_groups.setdefault(sid, []).append((ik, gene))
        for struct_id, compounds in struct_groups.items():
            rec = base / "receptors" / f"{struct_id}.pdbqt"
            pj = base / "pockets" / f"{struct_id}.json"
            if not rec.exists() or not pj.exists():
                continue
            for pocket in json.loads(pj.read_text())[
                    :cfg["fpocket_top_n"]]:
                key = (struct_id, pocket["pocket"])
                if key in done:
                    continue
                center = pocket["center"]
                size = [min(vg["box_max"], s) for s in pocket["size"]]
                ligdir = (ROOT / cfg["inputs_dir"] / "gpu_ligdirs" /
                          f"{struct_id}_p{pocket['pocket']}")
                ligdir.mkdir(parents=True, exist_ok=True)
                for ik, _g in compounds:
                    src = ROOT / cfg["main_ligands_dir"] / f"{ik}.pdbqt"
                    dst = ligdir / f"{ik}.pdbqt"
                    if src.exists() and not dst.exists():
                        dst.symlink_to(src)
                outdir = res / "gpu_out" / \
                    f"{struct_id}__p{pocket['pocket']}"
                outdir.mkdir(parents=True, exist_ok=True)
                cfgfile = res / "gpu_out" / \
                    f"{struct_id}__p{pocket['pocket']}.cfg"
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
                print(f"  {struct_id} p{pocket['pocket']}: "
                      f"{type(e).__name__}", flush=True)
                continue
            aff = {}
            for f in outdir.glob("*_out.pdbqt"):
                ik = f.name[:-len("_out.pdbqt")]
                m = REMARK.search(f.read_text(errors="ignore"))
                if m:
                    aff[ik] = float(m.group(1))
            genes = {ik: g for ik, g, _a in
                     [t for cs in lig_of[acc] for t in [(cs[0], cs[1], cs[2])]]}
            with out_tsv.open("a", newline="") as fh:
                if fresh:
                    fh.write("acc\tpocket\tinchikey\ttarget_gene\taffinity\n")
                    fresh = False
                for ik, a in aff.items():
                    fh.write(f"{struct_id}\t{pocket['pocket']}\t{ik}\t"
                             f"{genes.get(ik, '')}\t{a}\n")
                    rows_written += 1
                fh.flush()
            print(f"  {struct_id} p{pocket['pocket']}: "
                  f"{len(aff)} ligands scored", flush=True)
    print(f"FINISHED gpu dock, new rows: {rows_written}", flush=True)


if __name__ == "__main__":
    main()
