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


def _get(url, tries=4):
    import time as _t
    import urllib.request as _u
    for a in range(tries):
        try:
            with _u.urlopen(_u.Request(url), timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception:  # noqa: BLE001
            _t.sleep(min(60, 5 * 2 ** a))
    return None


FUNC_TYPES = {"EC50", "AC50"}
BIND_TYPES = {"IC50", "Ki", "Kd"}


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
    # direction-aware: only MATCH pairs enter delivery (structural rule)
    bp = pd.read_csv(
        ROOT / "product/chembl_branch/results/branch_pairs.tsv",
        sep="	", dtype=str).fillna("")
    # delivery csv has candidate ID + inchikey; join via target + direction
    # rebuild inchikey from previous delivery (stable mapping)
    # merge direction info from branch_pairs (proper DataFrame join)
    bp_dir = bp[["inchikey", "target_gene", "direction_match",
                 "crc_direction", "crc_effect", "action_types"]].copy()
    bp_dir = bp_dir.drop_duplicates(subset=["inchikey", "target_gene"])
    n_before = len(df)
    df = df.merge(bp_dir, on=["inchikey", "target_gene"], how="left")
    df = df[df["direction_match"] == "MATCH"].reset_index(drop=True)
    print(f"direction filter: {len(df)}/{n_before} pairs (MATCH only)",
          flush=True)

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

    # assay-type split (functional EC50/AC50 vs binding IC50/Ki/Kd),
    # cached; functional potency leads the metric columns per user order
    at_path = ROOT / "product/pd_validation/results/assay_types.tsv"
    if at_path.exists():
        at = pd.read_csv(at_path, sep="\t").fillna("")
        at_of = {(r.inchikey, r.target_gene):
                 dict(func=float(r.func_pchembl) if r.func_pchembl else None,
                      bind=float(r.bind_pchembl) if r.bind_pchembl else None,
                      types=r.types,
                      func_src=getattr(r, "func_src", ""),
                      bind_src=getattr(r, "bind_src", ""))
                 for r in at.itertuples()}
    else:
        at_of, rows_at = {}, []
        doc_doi = {}
        base = "https://www.ebi.ac.uk/chembl/api/data"

        def doi_of(doc):
            if not doc:
                return ""
            if doc not in doc_doi:
                r2 = _get(f"{base}/document/{doc}.json")
                doc_doi[doc] = (r2 or {}).get("doi", "") or ""
            return doc_doi[doc]

        for r in df.itertuples():
            cid = cid_of.get(r.inchikey, "")
            tc = tchembl_of.get((r.inchikey, r.target_gene), "")
            func = bind = None
            func_src = bind_src = ""
            types = []
            if cid and tc:
                rec = _get(f"{base}/activity.json?molecule_chembl_id={cid}"
                           f"&target_chembl_id={tc}"
                           f"&pchembl_value__isnull=false&limit=100")
                for act in (rec or {}).get("activities", []):
                    st = act.get("standard_type") or ""
                    pv = act.get("pchembl_value")
                    if pv is None:
                        continue
                    pv = float(pv)
                    types.append(st)
                    desc = (act.get("assay_description") or "")[:90]
                    src = f"{desc} [{doi_of(act.get('document_chembl_id'))}]"
                    if st in FUNC_TYPES and (func is None or pv > func):
                        func, func_src = pv, src
                    elif st in BIND_TYPES and (bind is None or pv > bind):
                        bind, bind_src = pv, src
            at_of[(r.inchikey, r.target_gene)] = dict(
                func=func, bind=bind, types=" ".join(sorted(set(types))),
                func_src=func_src, bind_src=bind_src)
            rows_at.append(dict(inchikey=r.inchikey, target_gene=r.target_gene,
                                func_pchembl=func, bind_pchembl=bind,
                                types=" ".join(sorted(set(types))),
                                func_src=func_src, bind_src=bind_src))
        pd.DataFrame(rows_at).to_csv(at_path, sep="\t", index=False)

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
        at = at_of.get((r.inchikey, r.target_gene),
                       dict(func=None, bind=None, types="",
                            func_src="", bind_src=""))
        func_desc = (f"功能性调节（EC50/AC50 测定，p={at['func']:.2f}；"
                     f"激动/拮抗方向未实验判定）" if at["func"] is not None
                     else "结合调节（激动/拮抗未实验判定）")

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
            CRC方向=getattr(r, "crc_direction", ""),
            CRC状态效应=getattr(r, "crc_effect", ""),
            化合物药理方向=getattr(r, "action_types", ""),
            主要功能=func_desc,
            化合物InChIKey=r.inchikey,
            化合物ChEMBL_ID=cid_of.get(r.inchikey, ""),
            靶点ChEMBL_ID=tchembl_of.get((r.inchikey, r.target_gene), ""),
            SMILES=smi,
            功能效力_pEC50_AC50=at["func"],
            结合效力_pIC50_Ki_Kd=at["bind"],
            功能效力_测定来源=at["func_src"],
            结合效力_测定来源=at["bind_src"],
            测定类型=at["types"],
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
            构象态=("非活性态" if any(k in str(r.action_types).upper() for k in ["INHIBIT","ANTAGON","BLOCK","NEGATIVE"]) else "活性态" if any(k in str(r.action_types).upper() for k in ["AGONIST","ACTIV","POSITIVE"]) else "未判定") if hasattr(r, "action_types") else "未判定",
            药理作用=("抑制剂" if any(k in str(r.action_types).upper() for k in ["INHIBIT","ANTAGON","BLOCK","NEGATIVE"]) else "激动剂" if any(k in str(r.action_types).upper() for k in ["AGONIST","ACTIV","POSITIVE"]) else "未判定") if hasattr(r, "action_types") else "未判定",
            对接门_通过="是", PD门_通过="是")
        rows.append(row)

        readme = f"""# {cand}

## 1. 基础信息
- 候选ID：{cand}
- 靶点名称：{r.target_gene}（UniProt {r.acc}，{fam_of.get(r.target_gene, '')}；ChEMBL {tchembl_of.get((r.inchikey, r.target_gene), '')}）
- CRC 方向：{getattr(r, 'crc_direction', '未知')}（状态效应 {getattr(r, 'crc_effect', '未知')}）
- 化合物药理方向：{getattr(r, 'action_types', '未记录')}
- 主要功能：{func_desc}
- 化合物：InChIKey {r.inchikey}；ChEMBL {cid_of.get(r.inchikey, '')}

## 2. 计算设计方案
- 分子结构（SMILES）：`{smi}`（另见 {cand}_ligand.smi / _ligand_pose.sdf）
- 小分子-靶点复合物结构预测文件：{cand}_complex.pdb（蛋白取自 {r.acc} 实验结构的清洁聚合物链；配体为 Vina-GPU 最优口袋 pose）
- 算法说明：化合物源自酵母 HIP/HOP 化学基因组筛选（Lee et al. 2014, E-MTAB-2391；车辆对照 z 谱与迁移任务排名全谱 Spearman，菌株标签置换检验 BH-FDR）。靶点身份来自 ChEMBL 定量注释（pChEMBL≥6、人源、离子通道/GPCR）。分子对接：Vina-GPU 2.1（thread 8000），口袋由 fpocket 默认参数预测（top-3，Site Score 排序），受体经 gemmi 清洗 + OpenBabel pH7.4 质子化，配体经 dimorphite-dl pH7.4 + meeko 制备。药效学评估：DeepPurpose 预训练模型集成（MPNN/CNN/Morgan × BindingDB IC50），SMILES+蛋白序列输入。
- 设计逻辑：人功能靶点的任务经 scFoundation 酵母表示迁移（B 路线，ESM2 注入）定义酵母可执行任务；化合物在酵母中显著执行任务后，其人靶点身份由定量药理注释确立，再以对接（结合合理性，{r.docking_affinity:.1f} kcal/mol）与双引擎药效学评估确认——DTI 预训练集成（校准后 pIC50 {r.pd_pic50_calibrated:.2f}）与同靶点分子近邻法（kNN pIC50 {(float(getattr(kn, 'pd_pic50_knn')) if kn is not None else float('nan')):.2f}，最近邻 Tanimoto {(float(getattr(kn, 'nn_tanimoto')) if kn is not None else float('nan')):.2f}，即与该靶点已知配体的结构相似度）。实测 pChEMBL {r.pchembl_measured}。

## 3-4. 实验验证结果 / 关键实验记录
未开展本项目湿实验。已引用的实测效力均来自 ChEMBL 策展的发表实验，可溯源：
- 功能效力（EC50/AC50）来源：{at["func_src"] or "无该类测定"}
- 结合效力（IC50/Ki/Kd）来源：{at["bind_src"] or "无该类测定"}
计划验证：酵母删除菌株 ± 化合物生长曲线（判定标准见项目湿实验接口文档）。

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
