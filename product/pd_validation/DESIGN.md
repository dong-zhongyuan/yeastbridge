# product/pd_validation — 计算药效学验证包（零湿实验，用户指令 2026-08-28）

## 定位
化合物中心的药效学（PD）验证层：对管线分配的人靶点-化合物对给出计算药效证据。约束：不新增任何湿实验；预训练模型直接使用，不自训（用户指令）；权重通道 hf-mirror.com（服务器实测可达）与本地 GitHub 中转。

## 引擎
1. **效价引擎（主力）DeepPurpose**：SMILES + 蛋白序列 → 结合亲和力（BindingDB/DAVIS/KIBA 预训练）。架构性零样本。多个编码器组合的预训练模型用于集成。
2. **签名引擎（补充）PerturbNet / chemCPA**：权重取自 cyclopeta/PerturbNet_reproduce（HF，经 hf-mirror），chemicalVAE 对任意 SMILES 出嵌入。注册限制：训练语境为癌症细胞系，通道/GPCR 相关语境覆盖弱，仅作方向性签名旁证。

## 防泄露协议（用户质询后注册）
- 预训练模型的训练集含 BindingDB/ChEMBL，48 个注释化合物大概率在训练集内。
- **校准主判据 = 时间切分**：从 ChEMBL 拉本管线靶点家族（GPCR/离子通道）`document_year ≥ 2023` 的定量活性（模型训练截止后入库），评估 DeepPurpose 预测的 RMSE/相关/分类 AUC——这是无泄露的泛化测量。
- 48 注释化合物组作为"可能已在训练集"参考组单列，不作为闸门。
- 应用对象（tier-3：ChEMBL/BindingDB 均无记录）天然不在任何训练集。
- kNN-SAR（ECFP4 + 适用域）保留为可解释对照引擎，训练侧同样遵守时间切分。

## 应用范围
1. 支线注释对（已知靶点，效价预测 vs 实测 pChEMBL 一致性）。
2. 主线 tier-3 前导（湿实验接口 16 化合物优先）在其分配靶点上的预测效价 + 适用域。
3. 产出列：预测 pIC50（引擎集成）、最近邻 Tanimoto（适用域）、签名方向旁证（PerturbNet 可用时）。

## 主张边界
计算药效学证据链（效价预测 + 适用域 + 签名旁证），报告如实标注不替代实测药效；预训练模型的训练分布限制（癌症细胞系语境）如实入档。

## 治理
本文件 + configs/pd_validation.json 先行提交；校准数据（时间切分集）拉取脚本与结果入 results/；所有模型文件不入 git（.gitignore）。
