#!/usr/bin/env python
"""chembl_branch stage 5: membrane MD for top docked pairs.

Build (run with structscreen python; antechamber/tleap/packmol-memgen are
called from the mdenv binaries): per selected pair, rebuild the ligand from
the Vina-GPU pose with meeko (docking frame == protein frame), merge into
the cleaned protein PDB, parameterize with GAFF2/AM1-BCC, and build the
POPC bilayer system with packmol-memgen (ff14SB/lipid21/tip3p, 0.15 M KCl,
PPM3 orientation). Output: md_systems/{gene}/complex_v1.parm7/.rst7.

Run (mdenv python, --stage run): OpenMM CUDA, LangevinMiddle 310 K,
semiisotropic MonteCarlo barostat, PME, 2 fs; minimization, 1 ns NPT
equilibration (positional restraints first 0.5 ns), then N replicas x
production with registered seeds. Registered in product/chembl_branch/
DESIGN.md (MD 协议 section); parameters in configs/chembl_branch.json.
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def selected_pairs(cfg, n=None):
    import pandas as pd

    best = pd.read_csv(
        ROOT / cfg["results_dir"] / "gpu_dock_pairs.tsv", sep="\t")
    top = (best.sort_values("affinity")
              .drop_duplicates("target_gene"))
    return top.head(n or cfg["md"]["n_systems"])


def build_one(row, cfg, inv):
    md = cfg["md"]
    gene = row.target_gene
    sysdir = ROOT / cfg["results_dir"] / "md_systems" / gene
    sysdir.mkdir(parents=True, exist_ok=True)
    if (sysdir / "complex_v1.parm7").exists():
        print(f"{gene}: already built", flush=True)
        return
    prot = ROOT / "product/chembl_branch/structures/fpocket_out" / \
        f"{row.acc}" / f"{row.acc}.pdb"
    pose = (ROOT / cfg["results_dir"] / "gpu_out" /
            f"{row.acc}__p{int(row.pocket)}" / f"{row.inchikey}_out.pdbqt")
    if not prot.exists() or not pose.exists():
        print(f"{gene}: missing inputs", flush=True)
        return

    from rdkit import Chem
    from meeko import PDBQTMolecule, RDKitMolCreate

    pm = PDBQTMolecule.from_file(str(pose), skip_typing=True)
    mols = [m for m in RDKitMolCreate.from_pdbqt_mol(pm) if m is not None]
    mol = mols[0]
    charge = Chem.GetFormalCharge(mol)
    sdf = sysdir / "lig.sdf"
    w = Chem.SDWriter(str(sdf))
    w.write(mol)
    w.close()
    print(f"{gene}: ligand extracted, charge {charge}", flush=True)

    lig_pdb = sysdir / "lig.pdb"
    Chem.MolToPDBFile(mol, str(lig_pdb))
    het = []
    for line in lig_pdb.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            line = line.replace(" UNL ", " LIG ")
            het.append("HETATM" + line[6:])
    complex_pdb = sysdir / "complex.pdb"
    complex_pdb.write_text(prot.read_text() + "\n".join(het) + "\nTER\nEND\n")

    b = md["mdenv_bin"]
    tool_env = {"AMBERHOME": str(Path(b).parent),
                "PATH": b + ":/usr/bin:/bin"}
    lig_mol2 = sysdir / "lig.mol2"
    subprocess.run(
        [f"{b}/antechamber", "-i", str(sdf), "-fi", "sdf",
         "-o", str(lig_mol2), "-fo", "mol2", "-c", "bcc", "-at", "gaff2",
         "-nc", str(charge), "-pf", "y", "-dr", "no"],
        check=True, capture_output=True, timeout=3600, env=tool_env)
    frcmod = sysdir / "lig.frcmod"
    subprocess.run(
        [f"{b}/parmchk2", "-i", str(lig_mol2), "-f", "mol2",
         "-o", str(frcmod)], check=True, capture_output=True, timeout=600,
        env=tool_env)
    lib = sysdir / "lig.lib"
    (sysdir / "tleap_lig.in").write_text(
        "source leaprc.gaff2\n"
        f"LIG = loadmol2 {lig_mol2}\n"
        f"saveoff LIG {lib}\nquit\n")
    subprocess.run([f"{b}/tleap", "-f", str(sysdir / "tleap_lig.in")],
                   check=True, capture_output=True, cwd=sysdir, timeout=600,
                   env=tool_env)
    print(f"{gene}: ligand parameterized", flush=True)

    subprocess.run(
        [f"{b}/packmol-memgen", "--pdb", str(complex_pdb),
         "--ligand_param", f"{frcmod}:{lib}",
         "--lipids", md["lipids"], "--ffwat", md["ffwat"],
         "--ffprot", "ff14SB", "--fflip", "lipid21",
         "--dist_wat", str(md["dist_wat"]),
         "--salt", "--saltcon", str(md["saltcon"]),
         "--gaff2"],
        check=True, cwd=sysdir, timeout=7200, env=tool_env)
    print(f"{gene}: system built", flush=True)


def run_one(gene, cfg):
    import openmm
    from openmm import app, unit

    md = cfg["md"]
    sysdir = ROOT / cfg["results_dir"] / "md_systems" / gene
    parm7 = sysdir / "complex_v1.parm7"
    rst7 = sysdir / "complex_v1.rst7"
    if not parm7.exists():
        print(f"{gene}: not built", flush=True)
        return
    pdb = app.AmberPrmtopFile(str(parm7))
    crd = app.AmberInpcrdFile(str(rst7))
    system = pdb.createSystem(nonbondedMethod=app.PME,
                              constraints=app.HBonds,
                              rigidWater=True)
    # positional restraints on all heavy atoms (protein + ligand);
    # k switched per the registered equilibration schedule
    rest = openmm.CustomExternalForce(
        "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    rest.addGlobalParameter(
        "k", 0.0 * unit.kilojoule_per_mole / unit.nanometer ** 2)
    for p in ("x0", "y0", "z0"):
        rest.addPerParticleParameter(p)
    for atom in pdb.topology.atoms():
        if atom.element is not None and atom.element.symbol != "H":
            pos = crd.positions[atom.index]
            rest.addParticle(atom.index,
                             (pos[0], pos[1], pos[2]))
    system.addForce(rest)
    system.addForce(openmm.MonteCarloBarostat(
        1 * unit.atmosphere, md["temperature_k"] * unit.kelvin, 25))
    platform = openmm.Platform.getPlatformByName("CUDA")
    dt = md["timestep_fs"] * unit.femtoseconds
    prod_steps = int(md["production_ns"] * 500000)  # 1 ns = 5e5 steps @2fs

    for seed in md["seeds"][:md["replicas"]]:
        out = sysdir / f"rep{seed}"
        out.mkdir(exist_ok=True)
        if (out / "final.chk").exists():
            continue
        integ = openmm.LangevinMiddleIntegrator(
            md["temperature_k"] * unit.kelvin, 1 / unit.picosecond, dt)
        integ.setRandomNumberSeed(seed)
        sim = app.Simulation(pdb.topology, system, integ, platform)
        sim.context.setPositions(crd.positions)
        sim.context.setPeriodicBoxVectors(*crd.boxVectors)
        sim.minimizeEnergy()
        sim.context.setVelocitiesToTemperature(
            md["temperature_k"] * unit.kelvin, seed)
        sim.context.setParameter(
            "k", 42.0 * unit.kilojoule_per_mole / unit.nanometer ** 2)
        sim.step(250000)  # 0.5 ns restrained equilibration
        sim.context.setParameter(
            "k", 0.0 * unit.kilojoule_per_mole / unit.nanometer ** 2)
        sim.step(250000)  # 0.5 ns free equilibration
        print(f"{gene} rep{seed}: equil done", flush=True)
        sim.reporters.append(app.DCDReporter(
            str(out / "prod.dcd"), 100000))  # 0.2 ns
        sim.reporters.append(app.StateDataReporter(
            str(out / "log.txt"), 250000, step=True, time=True,
            speed=True, remainingTime=True))
        sim.step(prod_steps)
        sim.saveState(str(out / "final.chk"))
        print(f"{gene} rep{seed}: production done", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/chembl_branch.json")
    ap.add_argument("--stage", choices=["build", "run"], required=True)
    ap.add_argument("--gene")
    args = ap.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())

    import pandas as pd

    inv = pd.read_csv(
        ROOT / cfg["results_dir"] / "inventory.tsv", sep="\t")
    if args.stage == "build":
        for row in selected_pairs(cfg).itertuples():
            if args.gene and row.target_gene != args.gene:
                continue
            build_one(row, cfg, inv)
    else:
        for gene in ([args.gene] if args.gene else
                     selected_pairs(cfg)["target_gene"]):
            run_one(gene, cfg)


if __name__ == "__main__":
    main()
