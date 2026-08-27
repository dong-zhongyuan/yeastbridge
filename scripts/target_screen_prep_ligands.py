#!/usr/bin/env python
"""target_screen stage 2: ligand preparation (hits + registered benchmarks).

Resolve canonical SMILES from ChEMBL by InChIKey (hits) and pref_name
(benchmarks), protonate at target pH with dimorphite-dl, embed a single
ETKDG conformer (Vina searches torsions from this seed conformer), and
write Vina-ready PDBQT via meeko. Resumable per ligand via
inputs/ligands/registry.tsv (only status=ok rows count as done; pending/
failed rows are retried and re-appended, aggregation keeps the last row
per lid). Registered in product/target_screen/DESIGN.md; parameters in
configs/target_screen.json.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def get_json(url, tries=8):
    for t in range(tries):
        try:
            with urllib.request.urlopen(
                    urllib.request.Request(url), timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception as e:  # noqa: BLE001
            print(f"  retry {t + 1}/{tries}: {e}", flush=True)
            time.sleep(min(600, 5 * 2 ** t))
    return None


def get_text(url, tries=4):
    for t in range(tries):
        try:
            with urllib.request.urlopen(
                    urllib.request.Request(url), timeout=60) as r:
                return r.read().decode().strip()
        except Exception:  # noqa: BLE001
            time.sleep(min(120, 5 * 2 ** t))
    return None


def smiles_fallback(inchikey):
    """Resolve SMILES outside ChEMBL: CACTUS keeps stereochemistry from the
    InChIKey; PubChem is connectivity-only (flat) fallback. Returns
    (smiles, source) or ("", "pending_network")."""
    txt = get_text(
        f"https://cactus.nci.nih.gov/chemical/structure/"
        f"{urllib.parse.quote(inchikey)}/smiles")
    if txt and not txt.startswith("<") and "\n" not in txt and len(txt) < 2000:
        return txt, "cactus"
    rec = get_json(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
        f"{inchikey}/property/CanonicalSMILES/JSON")
    if rec and "PropertyTable" in rec:
        props = rec["PropertyTable"]["Properties"][0]
        smi = props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")
        if smi:
            return smi, "pubchem_flat"
    if (txt is None) and (rec is None):
        return "", "pending_network"
    return "", "not_resolvable"


def protonated(smiles, ph):
    try:
        from dimorphite_dl import protonate_smiles
        out = protonate_smiles(smiles, ph_min=ph, ph_max=ph)
        if out:
            return out[0], "dimorphite"
    except Exception as e:  # noqa: BLE001
        print(f"  dimorphite unavailable/failed ({e}); neutral fallback",
              flush=True)
    return smiles, "neutral_fallback"


def make_pdbqt(smiles, ph, seed, out_path: Path):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    prot, note = protonated(smiles, ph)
    mol = Chem.MolFromSmiles(prot)
    if mol is None:
        prot, note = smiles, "neutral_fallback_sanitize"
        mol = Chem.MolFromSmiles(prot)
    if mol is None:
        raise ValueError("rdkit_sanitize_failed")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=seed,
                             useRandomCoords=True) != 0:
        raise ValueError("etkdg_embedding_failed")
    try:
        AllChem.MMFFOptimizeMolecules([mol], maxIters=500)
    except Exception:  # noqa: BLE001
        pass
    setup = MoleculePreparation()(mol)[0]
    string, ok, err = PDBQTWriterLegacy.write_string(setup)
    if not ok:
        raise ValueError(f"meeko_write_failed:{err}")
    out_path.write_text(string)
    return prot, note


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/target_screen.json")
    args = ap.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())

    import pandas as pd

    lig_dir = ROOT / cfg["inputs_dir"] / "ligands"
    lig_dir.mkdir(parents=True, exist_ok=True)
    reg_path = lig_dir / "registry.tsv"

    sig = pd.read_csv(ROOT / cfg["exec_matrix"], sep="\t")
    hits = sorted(set(sig.loc[sig["q"] < cfg["fdr_alpha"], "inchikey"]))
    print(f"hit compounds: {len(hits)}", flush=True)

    jobs = [("hit", ik, ik) for ik in hits]
    for b in cfg["benchmarks"]:
        jobs.append(("benchmark", f"BENCH__{b['name']}", b["name"]))

    done = set()
    if reg_path.exists():
        for r in pd.read_csv(reg_path, sep="\t").fillna("").to_dict("records"):
            if r["status"] == "ok":
                done.add(r["lid"])
    print(f"already prepared: {len(done)}", flush=True)

    fresh = not reg_path.exists()
    base = cfg["chembl_base"]
    with reg_path.open("a", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        if fresh:
            w.writerow(["lid", "label", "source", "inchikey", "smiles",
                        "status", "note"])
        for n, (source, lid, label) in enumerate(jobs, 1):
            if lid in done:
                continue
            if source == "hit":
                rec = get_json(
                    f"{base}/molecule.json?molecule_structures__"
                    f"standard_inchi_key={label}&limit=1")
                if rec is None:
                    w.writerow([lid, label, source, label, "",
                                "pending_network", ""])
                    continue
                src_note = "chembl"
                if not rec.get("molecules"):
                    smiles, src_note = smiles_fallback(label)
                    if not smiles:
                        w.writerow([lid, label, source, label, "",
                                    "pending_network"
                                    if src_note == "pending_network"
                                    else "failed", src_note])
                        continue
                    ik = label
                else:
                    ms = rec["molecules"][0].get("molecule_structures") or {}
                    smiles = ms.get("canonical_smiles", "")
                    ik = ms.get("standard_inchi_key", label)
            else:
                rec = get_json(
                    f"{base}/molecule.json?pref_name__iexact="
                    f"{urllib.parse.quote(label)}&limit=10")
                if rec is None:
                    w.writerow([lid, label, source, "", "",
                                "pending_network", ""])
                    continue
                cands = [m for m in rec.get("molecules", [])
                         if m.get("molecule_structures")]
                if not cands:
                    w.writerow([lid, label, source, "", "",
                                "failed", "name_not_resolved"])
                    continue
                ms = cands[0]["molecule_structures"]
                smiles = ms.get("canonical_smiles", "")
                ik = ms.get("standard_inchi_key", "")
            if not smiles:
                w.writerow([lid, label, source, ik, "", "failed", "no_smiles"])
                continue
            try:
                prot, note = make_pdbqt(smiles, cfg["ph"],
                                        cfg["conformer_seed"],
                                        lig_dir / f"{lid}.pdbqt")
                w.writerow([lid, label, source, ik, prot, "ok",
                            f"{src_note}:{note}"])
            except Exception as e:  # noqa: BLE001
                w.writerow([lid, label, source, ik, smiles, "failed",
                            f"prep:{type(e).__name__}"])
            fh.flush()
            if n % 20 == 0:
                print(f"  {n}/{len(jobs)} processed", flush=True)
    print("FINISHED ligand prep", flush=True)


if __name__ == "__main__":
    main()
