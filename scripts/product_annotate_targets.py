#!/usr/bin/env python
"""Product step 4: compound x target annotation + convergence tiers
(registered: product/drug_annotation/DESIGN.md; wet-lab interface v1).

Joins the executed yeast-task pairs (exec_matrix, q < fdr_alpha) with the
fetched ChEMBL pharmacology (results/chembl_targets.tsv, genes backfilled).
Three registered cases per (compound, executed-task target) pair:
tier 1 convergent (ChEMBL knows this compound acting on this very target),
tier 2 divergent (compound annotated on other targets only),
tier 3 novel (compound has no ChEMBL annotation; task target stands as
hypothesis, carried to the structure channel).

Outputs: results/annotation_pairs.tsv (all significant pairs, tiered),
results/compound_summary.tsv (per-compound overview),
results/wetlab_interface_v1.tsv (actionable shortlist for the wet lab:
all tier-1 pairs by potency, top-50 tier-3 novel pairs by yeast evidence),
results/annotation_summary.json.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/drug_annotation.json")
    args = ap.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())
    res = ROOT / cfg["results_dir"]

    em = pd.read_csv(ROOT / cfg["exec_matrix"], sep="\t")
    sig = em[em["q"] < cfg["fdr_alpha"]][
        ["target_id", "inchikey", "dose", "spearman_rho", "q", "smiles"]].copy()

    ch = pd.read_csv(res / "chembl_targets.tsv", sep="\t").fillna("")
    ch = ch[ch["target_gene"] != ""]
    ch["pchembl"] = pd.to_numeric(ch["pchembl"], errors="coerce")
    ch_exp = ch.assign(target_gene=ch["target_gene"].str.split()).explode(
        "target_gene")
    per_cg = ch_exp.groupby(["inchikey", "target_gene"]).agg(
        pchembl_max=("pchembl", "max"), n_act=("pchembl", "size")
    ).reset_index()
    genes_of = ch_exp.groupby("inchikey")["target_gene"].apply(set).to_dict()
    n_genes = {k: len(v) for k, v in genes_of.items()}
    chembl_id_of = ch.groupby("inchikey")["molecule_chembl_id"].first().to_dict()
    pchembl_of = per_cg.set_index(["inchikey", "target_gene"])[
        "pchembl_max"].to_dict()

    fam = pd.read_csv(
        ROOT / "product/transfer_route_b/inputs/universe_targets.tsv", sep="\t")
    fam_of = dict(zip(fam["target_id"], fam["target_family"]))

    rows = []
    for r in sig.itertuples():
        ik, tgt = r.inchikey, r.target_id
        g = genes_of.get(ik, set())
        tier = 1 if tgt in g else (2 if g else 3)
        rows.append(dict(
            inchikey=ik, chembl_id=chembl_id_of.get(ik, ""), task_target=tgt,
            family=fam_of.get(tgt, ""), dose=r.dose, rho=r.spearman_rho,
            q=r.q, tier=tier, pchembl_task=pchembl_of.get((ik, tgt)),
            n_chembl_genes=n_genes.get(ik, 0), smiles=r.smiles))
    df = pd.DataFrame(rows)
    df.to_csv(res / "annotation_pairs.tsv", sep="\t", index=False)

    comp = df.groupby("inchikey").agg(
        n_sig_tasks=("task_target", "size"),
        n_tier1=("tier", lambda s: int((s == 1).sum())),
        n_tier2=("tier", lambda s: int((s == 2).sum())),
        n_tier3=("tier", lambda s: int((s == 3).sum())),
        max_rho=("rho", "max"), min_q=("q", "min")).reset_index()
    comp["chembl_id"] = comp["inchikey"].map(chembl_id_of)
    comp["n_chembl_genes"] = comp["inchikey"].map(lambda k: n_genes.get(k, 0))
    topg = per_cg.sort_values("pchembl_max", ascending=False).groupby(
        "inchikey")["target_gene"].apply(
        lambda s: " ".join(s.head(5))).to_dict()
    comp["top_chembl_genes"] = comp["inchikey"].map(topg)
    comp["smiles"] = df.drop_duplicates("inchikey").set_index("inchikey")[
        "smiles"]
    comp.to_csv(res / "compound_summary.tsv", sep="\t", index=False)

    t1 = df[df.tier == 1].sort_values(
        ["pchembl_task", "q"], ascending=[False, True], na_position="last")
    t1 = t1.assign(priority=range(1, len(t1) + 1))
    t3 = df[df.tier == 3].sort_values(
        ["rho", "q"], ascending=[False, True]).head(50)
    t3 = t3.assign(priority=range(1, len(t3) + 1))
    cols = ["tier", "priority", "inchikey", "chembl_id", "task_target",
            "family", "rho", "q", "pchembl_task", "n_chembl_genes", "dose",
            "smiles"]
    wl = pd.concat([t1[cols], t3[cols]])
    assay = {"ion_channel": "patch clamp / ion flux assay",
             "gpcr": "second messenger (cAMP/IP1/Ca2+) or beta-arrestin"}
    wl["suggested_assay"] = wl["family"].map(assay).fillna("per target class")
    wl.to_csv(res / "wetlab_interface_v1.tsv", sep="\t", index=False)

    stats = dict(
        n_sig_pairs=int(len(df)),
        n_compounds=int(df.inchikey.nunique()),
        n_targets=int(df.task_target.nunique()),
        tier1_pairs=int((df.tier == 1).sum()),
        tier1_compounds=int(df[df.tier == 1].inchikey.nunique()),
        tier2_compounds=int(df[df.tier == 2].inchikey.nunique()),
        tier3_compounds=int(df[df.tier == 3].inchikey.nunique()),
        n_annotated_compounds=int(len(genes_of)),
    )
    (res / "annotation_summary.json").write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2), flush=True)


if __name__ == "__main__":
    main()
