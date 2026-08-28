# YeastBridge_Mol019

## 1. 基础信息
- 候选ID：YeastBridge_Mol019
- 靶点名称：DRD3（UniProt P35462，gpcr；ChEMBL CHEMBL234）
- 主要功能：结合调节（激动/拮抗未实验判定）
- 化合物：InChIKey LNEPOXFFQSENCJ-UHFFFAOYSA-N；ChEMBL CHEMBL54

## 2. 计算设计方案
- 分子结构（SMILES）：`C1CN(CCC1(C2=CC=C(C=C2)Cl)O)CCCC(=O)C3=CC=C(C=C3)F`（另见 YeastBridge_Mol019_ligand.smi / _ligand_pose.sdf）
- 小分子-靶点复合物结构预测文件：YeastBridge_Mol019_complex.pdb（蛋白取自 P35462 实验结构的清洁聚合物链；配体为 Vina-GPU 最优口袋 pose）
- 算法说明：化合物源自酵母 HIP/HOP 化学基因组筛选（Lee et al. 2014, E-MTAB-2391；车辆对照 z 谱与迁移任务排名全谱 Spearman，菌株标签置换检验 BH-FDR）。靶点身份来自 ChEMBL 定量注释（pChEMBL≥6、人源、离子通道/GPCR）。分子对接：Vina-GPU 2.1（thread 8000），口袋由 fpocket 默认参数预测（top-3，Site Score 排序），受体经 gemmi 清洗 + OpenBabel pH7.4 质子化，配体经 dimorphite-dl pH7.4 + meeko 制备。药效学评估：DeepPurpose 预训练模型集成（MPNN/CNN/Morgan × BindingDB IC50），SMILES+蛋白序列输入。
- 设计逻辑：人功能靶点的任务经 scFoundation 酵母表示迁移（B 路线，ESM2 注入）定义酵母可执行任务；化合物在酵母中显著执行任务后，其人靶点身份由定量药理注释确立，再以对接（结合合理性，-8.1 kcal/mol）与预训练效价模型（集成 pIC50 5.74，实测 pChEMBL 9.16）双重确认。

## 3-4. 实验验证结果 / 关键实验记录
未开展湿实验（如实标注）。计划验证：酵母删除菌株 ± 化合物生长曲线（任务 top15 + 对照 bottom5，判定标准见湿实验接口 v1 文档）。

## 局限（如实）
对接打分为真空近似；DeepPurpose 绝对效价精度有限（本集合实测值处于 6-8 窄区间，预测-实测相关弱，PD 结果作合理性确认而非效价定量）；AF2 来源结构无寡聚/脂环境。
