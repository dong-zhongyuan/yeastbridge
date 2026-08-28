# YeastBridge_Mol032

## 1. 基础信息
- 候选ID：YeastBridge_Mol032
- 靶点名称：KCNH2（UniProt Q12809，ion_channel；ChEMBL CHEMBL240）
- 主要功能：功能性调节（EC50/AC50 测定，p=7.40；激动/拮抗方向未实验判定）
- 化合物：InChIKey QOYHHIBFXOOADH-UHFFFAOYSA-N；ChEMBL CHEMBL46516

## 2. 计算设计方案
- 分子结构（SMILES）：`C1CN(CCC12C(=O)NCN2C3=CC=CC=C3)CCCC(C4=CC=C(C=C4)F)C5=CC=C(C=C5)F`（另见 YeastBridge_Mol032_ligand.smi / _ligand_pose.sdf）
- 小分子-靶点复合物结构预测文件：YeastBridge_Mol032_complex.pdb（蛋白取自 Q12809 实验结构的清洁聚合物链；配体为 Vina-GPU 最优口袋 pose）
- 算法说明：化合物源自酵母 HIP/HOP 化学基因组筛选（Lee et al. 2014, E-MTAB-2391；车辆对照 z 谱与迁移任务排名全谱 Spearman，菌株标签置换检验 BH-FDR）。靶点身份来自 ChEMBL 定量注释（pChEMBL≥6、人源、离子通道/GPCR）。分子对接：Vina-GPU 2.1（thread 8000），口袋由 fpocket 默认参数预测（top-3，Site Score 排序），受体经 gemmi 清洗 + OpenBabel pH7.4 质子化，配体经 dimorphite-dl pH7.4 + meeko 制备。药效学评估：DeepPurpose 预训练模型集成（MPNN/CNN/Morgan × BindingDB IC50），SMILES+蛋白序列输入。
- 设计逻辑：人功能靶点的任务经 scFoundation 酵母表示迁移（B 路线，ESM2 注入）定义酵母可执行任务；化合物在酵母中显著执行任务后，其人靶点身份由定量药理注释确立，再以对接（结合合理性，-7.2 kcal/mol）与双引擎药效学评估确认——DTI 预训练集成（校准后 pIC50 6.64）与同靶点分子近邻法（kNN pIC50 7.01，最近邻 Tanimoto 0.47，即与该靶点已知配体的结构相似度）。实测 pChEMBL 7.4。

## 3-4. 实验验证结果 / 关键实验记录
未开展湿实验。计划验证：酵母删除菌株 ± 化合物生长曲线（判定标准见项目湿实验接口文档）。

方法学细节见项目技术文档（含各计算层的适用范围说明）。
