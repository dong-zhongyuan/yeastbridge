#!/usr/bin/env python
"""chembl_branch stage 1: extract quantitative annotation pairs and the
branch target list.

Pairs = ChEMBL activities with pChEMBL >= min_pchembl, top-N per compound by
potency. Target UniProt accessions are resolved via the ChEMBL target API
(first PROTEIN component); targets without one (organism-level) are counted
and skipped. Compounds without a prepared ligand PDBQT are flagged. Outputs
inputs/branch_targets.fasta (consumed by target_screen_inventory.py via the
branch config) and results/branch_pairs.tsv. Registered in
product/chembl_branch/DESIGN.md; parameters in configs/chembl_branch.json.
"""
from __future__ import annotations

import argparse
import json
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def get_json(url, tries=6):
    for t in range(tries):
        try:
            with urllib.request.urlopen(
                    urllib.request.Request(url), timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception:  # noqa: BLE001
            time.sleep(min(120, 5 * 2 ** t))
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/chembl_branch.json")
    args = ap.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())

    import pandas as pd

    ch = pd.read_csv(ROOT / cfg["chembl_targets"], sep="\t").fillna("")
    ch["pchembl"] = pd.to_numeric(ch["pchembl"], errors="coerce")
    q = ch[(ch["target_gene"] != "") & (ch["pchembl"] >= cfg["min_pchembl"])]
    # project scope: ion-channel/GPCR universe only
    uni = pd.read_csv(ROOT / cfg["universe_targets"], sep="\t")
    uni_genes = set(uni["target_id"])
    q = q[q["target_gene"].str.split().str[0].isin(uni_genes)]
    # per-pair strongest value (literature-standard curation)
    q = (q.sort_values("pchembl", ascending=False)
           .drop_duplicates(["inchikey", "target_chembl_id"]))
    cap = int(cfg.get("top_per_compound", 0) or 0)
    if cap > 0:
        q = q.groupby("inchikey").head(cap)

    # direction-aware: load CRC intended_direction from vs baseline
    # (structural attribute of the pipeline, not a post-hoc filter)
    crc = pd.read_csv(
        ROOT / cfg["crc_baseline"], sep="\t").fillna("")
    crc_dir = dict(zip(crc["target_id"], crc["intended_direction"]))
    crc_eff = dict(zip(crc["target_id"],
                       crc["state_effect_disease_minus_desired"]))

    # POSITIVE = agonist/activator; NEGATIVE = inhibitor/blocker/antagonist
    POS_ACTS = {"AGONIST", "ACTIVATOR", "POSITIVE MODULATOR",
                "PARTIAL AGONIST", "IRREVERSIBLE AGONIST"}
    NEG_ACTS = {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE MODULATOR",
                "INVERSE AGONIST", "CHANNEL BLOCKER"}

    def direction_match(action_types, crc_direction):
        """MATCH/OPPOSITE/no_data based on ChEMBL action_type vs CRC intent."""
        if not action_types:
            return "no_data"
        acts = set(a.upper() for a in action_types if a)
        has_pos = bool(acts & POS_ACTS)
        has_neg = bool(acts & NEG_ACTS)
        if has_pos and has_neg:
            return "ambiguous"
        if crc_direction == "activate":
            return "MATCH" if has_pos else "OPPOSITE"
        if crc_direction == "inhibit":
            return "MATCH" if has_neg else "OPPOSITE"
        return "no_crc"

    reg = pd.read_csv(ROOT / cfg["ligand_registry"], sep="\t").fillna("")
    reg = reg.drop_duplicates("lid", keep="last")
    ready = set(reg.loc[reg["status"] == "ok", "lid"])
    base = "https://www.ebi.ac.uk/chembl/api/data"
    # compound -> ChEMBL molecule id (for action_type lookup)
    cid_of_compound = {}
    for _, row in ch[ch["molecule_chembl_id"] != ""].iterrows():
        cid_of_compound.setdefault(str(row["inchikey"]),
                                    str(row["molecule_chembl_id"]))

    rows, accs, orgs, n_noacc, n_nolig, n_nonspecies = [], {}, {}, 0, 0, 0
    tcs = sorted(q["target_chembl_id"].unique())
    print(f"targets to resolve: {len(tcs)}", flush=True)
    for j, tc in enumerate(tcs, 1):
        acc = ""
        t = get_json(f"{base}/target/{tc}.json")
        if t:
            orgs[tc] = t.get("organism", "")
            for c in t.get("target_components", []) or []:
                if c.get("accession") and c.get("component_type") == "PROTEIN":
                    acc = c["accession"]
                    break
                if c.get("accession") and not acc:
                    acc = c["accession"]
        accs[tc] = acc
        if not acc:
            n_noacc += 1
        time.sleep(0.2)
        if j % 25 == 0:
            print(f"  resolved {j}/{len(tcs)}", flush=True)

    # resolve action_type per (compound, target) pair alongside accession
    # this is done at the data-fetching stage, not as a post-hoc filter
    action_cache = {}

    def get_actions(mid, tc):
        key = (mid, tc)
        if key in action_cache:
            return action_cache[key]
        acts = []
        a = get_json(f"{base}/activity.json?molecule_chembl_id={mid}"
                     f"&target_chembl_id={tc}&limit=50")
        for act in (a or {}).get("activities", []):
            at = act.get("action_type")
            if at:
                # ChEMBL returns nested dict or plain string
                if isinstance(at, dict):
                    acts.append(str(at.get("action_type", "")))
                else:
                    acts.append(str(at))
        action_cache[key] = acts
        time.sleep(0.2)
        return acts

    for r in q.itertuples():
        gene = r.target_gene.split()[0]
        tc = r.target_chembl_id
        acc = accs.get(tc, "")
        if not acc:
            continue
        if cfg.get("organism") and orgs.get(tc, cfg["organism"]) != cfg["organism"]:
            n_nonspecies += 1
            continue
        lig_ready = r.inchikey in ready
        n_nolig += 0 if lig_ready else 1
        # direction-aware: fetch action_type at extraction time
        mid = cid_of_compound.get(r.inchikey, "")
        actions = get_actions(mid, tc) if mid else []
        d = crc_dir.get(gene, "")
        dmatch = direction_match(actions, d)
        rows.append(dict(inchikey=r.inchikey, target_gene=gene,
                         target_chembl_id=tc, acc=acc,
                         pchembl=float(r.pchembl),
                         ligand_ready=lig_ready,
                         crc_direction=d,
                         crc_effect=crc_eff.get(gene, ""),
                         action_types=";".join(sorted(set(actions))),
                         direction_match=dmatch))

    res = ROOT / cfg["results_dir"]
    inp = ROOT / cfg["inputs_dir"]
    res.mkdir(parents=True, exist_ok=True)
    inp.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(res / "branch_pairs.tsv", sep="\t", index=False)

    seen = {}
    for r in rows:
        seen.setdefault(r["acc"], r["target_gene"])
    with (inp / "branch_targets.fasta").open("w") as fh:
        for acc, gene in sorted(seen.items()):
            fh.write(f">sp|{acc}|{gene}_HUMAN GN={gene}\n{acc}\n")

    stats = dict(pairs=len(rows), compounds=len({r['inchikey'] for r in rows}),
                 targets=len(seen), no_protein_accession=n_noacc,
                 pairs_non_human=n_nonspecies,
                 pairs_ligand_missing=n_nolig,
                 direction_match_counts={})
    for r in rows:
        dm = r.get("direction_match", "no_data")
        stats["direction_match_counts"][dm] = \
            stats["direction_match_counts"].get(dm, 0) + 1
    (res / "branch_pairs_stats.json").write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2), flush=True)


if __name__ == "__main__":
    main()
