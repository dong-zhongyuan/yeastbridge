# YeastBridge_Mol021

## 1. 基础信息
- 候选ID：YeastBridge_Mol021
- 靶点名称：HRH1（UniProt P35367，gpcr；ChEMBL CHEMBL231）
- 主要功能：功能性调节（EC50/AC50 测定，p=5.20；激动/拮抗方向未实验判定）
- 化合物：InChIKey LNEPOXFFQSENCJ-UHFFFAOYSA-N；ChEMBL CHEMBL54

## 2. 计算设计方案
- 分子结构（SMILES）：`C1CN(CCC1(C2=CC=C(C=C2)Cl)O)CCCC(=O)C3=CC=C(C=C3)F`（另见 YeastBridge_Mol021_ligand.smi / _ligand_pose.sdf）
- 小分子-靶点复合物结构预测文件：YeastBridge_Mol021_complex.pdb（蛋白取自 P35367 实验结构的清洁聚合物链；配体为 Vina-GPU 最优口袋 pose）
- 算法说明：化合物源自酵母 HIP/HOP 化学基因组筛选（Lee et al. 2014, E-MTAB-2391；车辆对照 z 谱与迁移任务排名全谱 Spearman，菌株标签置换检验 BH-FDR）。靶点身份来自 ChEMBL 定量注释（pChEMBL≥6、人源、离子通道/GPCR）。分子对接：Vina-GPU 2.1（thread 8000），口袋由 fpocket 默认参数预测（top-3，Site Score 排序），受体经 gemmi 清洗 + OpenBabel pH7.4 质子化，配体经 dimorphite-dl pH7.4 + meeko 制备。药效学评估：DeepPurpose 预训练模型集成（MPNN/CNN/Morgan × BindingDB IC50），SMILES+蛋白序列输入。
- 设计逻辑：人功能靶点的任务经 scFoundation 酵母表示迁移（B 路线，ESM2 注入）定义酵母可执行任务；化合物在酵母中显著执行任务后，其人靶点身份由定量药理注释确立，再以对接（结合合理性，-7.9 kcal/mol）与双引擎药效学评估确认——DTI 预训练集成（校准后 pIC50 6.88）与同靶点分子近邻法（kNN pIC50 6.16，最近邻 Tanimoto 0.37，即与该靶点已知配体的结构相似度）。实测 pChEMBL 6.49。

## 3-4. 实验验证结果 / 关键实验记录
未开展本项目湿实验。已引用的实测效力均来自 ChEMBL 策展的发表实验，可溯源：
- 功能效力（EC50/AC50）来源：Binding affinity towards human HRH1 in an in vitro assay measured by filtration [10.1038/s41467-023-40064-9]
- 结合效力（IC50/Ki/Kd）来源：DRUGMATRIX: Histamine H1, Central radioligand binding (ligand: [3H] Pyrilamine) []
计划验证：酵母删除菌株 ± 化合物生长曲线（判定标准见项目湿实验接口文档）。

方法学细节见项目技术文档（含各计算层的适用范围说明）。
