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
results/wetlab_interface_v1.tsv (PRIMARY wet-lab deliverable, yeast system:
per prioritized pair the task's top-15 strains plus bottom-5 specificity
controls with the compound's strain z-profile at the task dose and
direction predictions — the validation object is the human-goal-to-yeast
transfer and its execution, not a human-target assay),
results/target_context_v1.tsv (human-target context for the same pairs),
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
    # quantitative-only gene sets: convergence must rest on activities that
    # carry a pChEMBL value (value-less records are often counterscreen or
    # inactive deposits, e.g. gemcitabin "activities" on adrenergic
    # receptors)
    quant_pairs = per_cg[per_cg["pchembl_max"].notna()]
    quant_genes_of = quant_pairs.groupby("inchikey")["target_gene"].apply(
        set).to_dict()
    n_genes = {k: len(v) for k, v in quant_genes_of.items()}
    chembl_id_of = ch.groupby("inchikey")["molecule_chembl_id"].first().to_dict()
    pchembl_of = per_cg.set_index(["inchikey", "target_gene"])[
        "pchembl_max"].to_dict()

    fam = pd.read_csv(
        ROOT / "product/transfer_route_b/inputs/universe_targets.tsv", sep="\t")
    fam_of = dict(zip(fam["target_id"], fam["target_family"]))

    rows = []
    for r in sig.itertuples():
        ik, tgt = r.inchikey, r.target_id
        g_all = genes_of.get(ik, set())
        g_quant = quant_genes_of.get(ik, set())
        if tgt in g_quant:
            tier = 1
        elif g_all:
            tier = 2
        else:
            tier = 3
        rows.append(dict(
            inchikey=ik, chembl_id=chembl_id_of.get(ik, ""), task_target=tgt,
            family=fam_of.get(tgt, ""), dose=r.dose, rho=r.spearman_rho,
            q=r.q, tier=tier, pchembl_task=pchembl_of.get((ik, tgt)),
            n_chembl_genes=n_genes.get(ik, 0), smiles=r.smiles,
            chembl_task_mention_only=bool(
                tgt in g_all and tgt not in g_quant)))
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
    pd.concat([t1[cols], t3[cols]]).to_csv(
        res / "target_context_v1.tsv", sep="\t", index=False)

    # --- 酵母湿实验接口(主交付):验证对象是"人功能目标 -> 酵母任务"迁移
    # 与化合物任务执行;每对给任务 top15 菌株 + bottom5 特异性对照,
    # 附化合物在任务最优剂量下的菌株 z 谱与方向预测 ---
    import numpy as np

    cfg_exec = json.loads((ROOT / "configs/product_execute.json").read_text())
    npz = np.load(cfg_exec["response_npz"], allow_pickle=True)
    orfs = npz["strain_orfs"]
    iks_arr = npz["compound_inchikeys"]
    doses_arr = npz["doses"]
    units_arr = npz["dose_units"]
    veh_arr = npz["is_vehicle"]
    Z = npz["z_score"]
    orf_ix = {o: i for i, o in enumerate(orfs)}
    task_dir = ROOT / "product/transfer_route_b/results"

    def pair_zcolumn(ik, dose):
        """Library array for this compound: the dose recorded in exec_matrix
        when the pair was called; fallback = its strongest-effect array."""
        best, best_score = None, -1.0
        for j in range(len(iks_arr)):
            if iks_arr[j] != ik or veh_arr[j] or units_arr[j] != "micromolar":
                continue
            if dose is not None and \
                    abs(float(doses_arr[j]) - float(dose)) < 0.01:
                return j, doses_arr[j]
            sc = float(np.abs(Z[:, j]).mean())
            if sc > best_score:
                best, best_score = j, sc
        return (best, doses_arr[best]) if best is not None else (None, None)

    rows = []
    for r in pd.concat([t1, t3]).itertuples():
        tf = task_dir / f"yeast_task_{r.task_target}.tsv"
        if not tf.exists():
            continue
        task = pd.read_csv(tf, sep="\t")
        sel = pd.concat([task.head(15).assign(control_role="top_task"),
                         task.tail(5).assign(
                             control_role="bottom_task_control")])
        j, dose_used = pair_zcolumn(r.inchikey, r.dose)
        for t in sel.itertuples():
            z = None
            if j is not None and t.yeast_gene in orf_ix:
                z = float(Z[orf_ix[t.yeast_gene], j])
            rows.append(dict(
                tier=r.tier, priority=r.priority, inchikey=r.inchikey,
                chembl_id=r.chembl_id, task_target=r.task_target,
                family=r.family, rho=r.rho, q=r.q, dose_uM=dose_used,
                smiles=r.smiles, orf=t.yeast_gene, task_rank=t.rank,
                task_cosine=t.cosine, compound_z=z,
                predicted_direction=(
                    None if z is None or z == 0 else
                    ("hypersensitive" if z < 0 else "resistant")),
                control_role=t.control_role))
    yl = pd.DataFrame(rows)
    yl.to_csv(res / "wetlab_interface_v1.tsv", sep="\t", index=False)

    stats = dict(
        n_sig_pairs=int(len(df)),
        n_compounds=int(df.inchikey.nunique()),
        n_targets=int(df.task_target.nunique()),
        tier1_pairs=int((df.tier == 1).sum()),
        tier1_compounds=int(df[df.tier == 1].inchikey.nunique()),
        tier2_compounds=int(df[df.tier == 2].inchikey.nunique()),
        tier3_compounds=int(df[df.tier == 3].inchikey.nunique()),
        n_annotated_compounds=int(len(genes_of)),
        n_valueless_mention_pairs=int(df["chembl_task_mention_only"].sum()),
        n_interface_pairs=int(yl["inchikey"].nunique()),
        n_interface_strain_rows=int(len(yl)),
    )
    (res / "annotation_summary.json").write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2), flush=True)


if __name__ == "__main__":
    main()
