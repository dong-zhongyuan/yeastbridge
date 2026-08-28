#!/usr/bin/env python
"""Delivery package builder: candidate small-molecule target pairs that
passed every step of the fixed route (yeast-screen hit -> ChEMBL
quantitative annotation -> docking gate -> PD gate).

Outputs (configs/delivery.json):
- one summary CSV (one row per pair, all pipeline scores, utf-8-sig for Excel)
- one folder per pair: complex PDB (clean receptor + docked pose ligand),
  ligand pose SDF, ligand SMILES, per-pair README (submission-format fields:
  basic info / computational design / algorithm statement / design logic /
  honest 'not performed' experiment section), per-pair scores TSV.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/delivery.json")
    args = ap.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())

    import pandas as pd
    from rdkit import Chem
    from meeko import PDBQTMolecule, RDKitMolCreate

    df = pd.read_csv(ROOT / cfg["pairs_calibrated_tsv"], sep="\t")
    df = df[df["final_pass"]].sort_values("docking_affinity").reset_index(
        drop=True)
    knn = pd.read_csv(
        ROOT / "product/pd_validation/results/final_pairs_knn.tsv",
        sep="\t")
    knn_of = {(r.inchikey, r.target_gene): r for r in knn.itertuples()}

    em = pd.read_csv(
        ROOT / "product/execute_hiphop/results/exec_matrix.tsv", sep="\t")
    em_sig = em[em["q"] < 0.1]
    yeast_ev = em_sig.groupby("inchikey").agg(
        yeast_max_rho=("spearman_rho", "max"),
        yeast_min_q=("q", "min"),
        yeast_n_tasks=("target_id", "nunique"))
    smi_of = dict(zip(em["inchikey"], em["smiles"]))

    fam = pd.read_csv(
        ROOT / "product/transfer_route_b/inputs/universe_targets.tsv",
        sep="\t")
    fam_of = dict(zip(fam["target_id"], fam["target_family"]))
    bp = pd.read_csv(
        ROOT / "product/chembl_branch/results/branch_pairs.tsv", sep="\t")
    tchembl_of = {(r.inchikey, r.target_gene): r.target_chembl_id
                  for r in bp.itertuples()}
    gd = pd.read_csv(
        ROOT / "product/chembl_branch/results/gpu_dock.tsv", sep="\t")
    pocket_of = {(r.inchikey, r.acc): int(r.pocket)
                 for r in gd.sort_values("affinity")
                 .groupby(["inchikey", "acc"]).first().reset_index()
                 .itertuples()}
    ch = pd.read_csv(
        ROOT / "product/drug_annotation/results/chembl_targets.tsv",
        sep="\t").fillna("")
    cid_of = ch[ch["molecule_chembl_id"] != ""].groupby(
        "inchikey")["molecule_chembl_id"].first().to_dict()

    out = ROOT / cfg["out_dir"]
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, r in enumerate(df.itertuples(), 1):
        cand = f"{cfg['candidate_prefix']}_Mol{i:03d}"
        folder = out / cand
        folder.mkdir(exist_ok=True)
        smi = smi_of.get(r.inchikey, "")
        ev = yeast_ev.loc[r.inchikey] if r.inchikey in yeast_ev.index else None
        kn = knn_of.get((r.inchikey, r.target_gene))

        # complex PDB: clean polymer receptor + docked pose as HETATM LIG
        pk = pocket_of.get((r.inchikey, r.acc), 1)
        pose = (ROOT / "product/chembl_branch/results/gpu_out" /
                f"{r.acc}__p{pk}" / f"{r.inchikey}_out.pdbqt")
        prot = (ROOT / "product/chembl_branch/structures/fpocket_out" /
                f"{r.acc}" / f"{r.acc}.pdb")
        lig_block, sdf_text = "", ""
        if pose.exists():
            pm = PDBQTMolecule.from_file(str(pose), skip_typing=True)
            mols = [m for m in RDKitMolCreate.from_pdbqt_mol(pm)
                    if m is not None]
            if mols:
                lig_block = Chem.MolToPDBBlock(mols[0])
                lig_block = "\n".join(
                    "HETATM" + ln[6:] for ln in lig_block.splitlines()
                    if ln.startswith(("ATOM", "HETATM")))
                s = Chem.SDWriter(str(folder / f"{cand}_ligand_pose.sdf"))
                s.write(mols[0])
                s.close()
        if prot.exists() and lig_block:
            (folder / f"{cand}_complex.pdb").write_text(
                prot.read_text() + lig_block + "TER\nEND\n")
        (folder / f"{cand}_ligand.smi").write_text(
            f"{smi}\t{cand}\t{r.target_gene}\n")

        row = dict(
            候选ID=cand, 靶点名称=r.target_gene, 靶点Uniprot=r.acc,
            靶点家族=fam_of.get(r.target_gene, ""),
            主要功能="结合调节（激动/拮抗未实验判定）",
            化合物InChIKey=r.inchikey,
            化合物ChEMBL_ID=cid_of.get(r.inchikey, ""),
            靶点ChEMBL_ID=tchembl_of.get((r.inchikey, r.target_gene), ""),
            SMILES=smi,
            酵母筛选_q=(None if ev is None else round(float(ev["yeast_min_q"]), 4)),
            酵母筛选_最高rho=(None if ev is None else round(float(ev["yeast_max_rho"]), 4)),
            酵母筛选_显著任务数=(None if ev is None else int(ev["yeast_n_tasks"])),
            ChEMBL实测_pChEMBL=r.pchembl_measured,
            对接亲和能_kcal_mol=r.docking_affinity,
            PD预测_pIC50_DTI集成=round(r.pd_pic50_calibrated, 2),
            PD预测_pIC50_kNN近邻=(None if kn is None else round(
                float(getattr(kn, "pd_pic50_knn")), 2)),
            PD适用域_最近邻Tanimoto=(None if kn is None else round(
                float(getattr(kn, "nn_tanimoto")), 2)),
            对接门_通过="是", PD门_通过="是")
        rows.append(row)

        readme = f"""# {cand}

## 1. 基础信息
- 候选ID：{cand}
- 靶点名称：{r.target_gene}（UniProt {r.acc}，{fam_of.get(r.target_gene, '')}；ChEMBL {tchembl_of.get((r.inchikey, r.target_gene), '')}）
- 主要功能：结合调节（激动/拮抗未实验判定）
- 化合物：InChIKey {r.inchikey}；ChEMBL {cid_of.get(r.inchikey, '')}

## 2. 计算设计方案
- 分子结构（SMILES）：`{smi}`（另见 {cand}_ligand.smi / _ligand_pose.sdf）
- 小分子-靶点复合物结构预测文件：{cand}_complex.pdb（蛋白取自 {r.acc} 实验结构的清洁聚合物链；配体为 Vina-GPU 最优口袋 pose）
- 算法说明：化合物源自酵母 HIP/HOP 化学基因组筛选（Lee et al. 2014, E-MTAB-2391；车辆对照 z 谱与迁移任务排名全谱 Spearman，菌株标签置换检验 BH-FDR）。靶点身份来自 ChEMBL 定量注释（pChEMBL≥6、人源、离子通道/GPCR）。分子对接：Vina-GPU 2.1（thread 8000），口袋由 fpocket 默认参数预测（top-3，Site Score 排序），受体经 gemmi 清洗 + OpenBabel pH7.4 质子化，配体经 dimorphite-dl pH7.4 + meeko 制备。药效学评估：DeepPurpose 预训练模型集成（MPNN/CNN/Morgan × BindingDB IC50），SMILES+蛋白序列输入。
- 设计逻辑：人功能靶点的任务经 scFoundation 酵母表示迁移（B 路线，ESM2 注入）定义酵母可执行任务；化合物在酵母中显著执行任务后，其人靶点身份由定量药理注释确立，再以对接（结合合理性，{r.docking_affinity:.1f} kcal/mol）与双引擎药效学评估确认——DTI 预训练集成（校准后 pIC50 {r.pd_pic50_calibrated:.2f}）与同靶点分子近邻法（kNN pIC50 {(float(getattr(kn, 'pd_pic50_knn')) if kn is not None else float('nan')):.2f}，最近邻 Tanimoto {(float(getattr(kn, 'nn_tanimoto')) if kn is not None else float('nan')):.2f}，即与该靶点已知配体的结构相似度）。实测 pChEMBL {r.pchembl_measured}。

## 3-4. 实验验证结果 / 关键实验记录
未开展湿实验。计划验证：酵母删除菌株 ± 化合物生长曲线（判定标准见项目湿实验接口文档）。

方法学细节见项目技术文档（含各计算层的适用范围说明）。
"""
        (folder / f"{cand}_README.md").write_text(readme)
        pd.DataFrame([row]).to_csv(folder / f"{cand}_scores.tsv",
                                   sep="\t", index=False)

    pd.DataFrame(rows).to_csv(out / "summary_candidates.csv",
                              index=False, encoding="utf-8-sig")
    print(json.dumps(dict(pairs=len(rows), out=str(out)), indent=2),
          flush=True)


if __name__ == "__main__":
    main()
