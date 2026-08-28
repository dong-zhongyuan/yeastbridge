# 模型选择证据链（feasibility）

## Norman 多基座扰动预测基准（vs 项目 2026-08-17，导入 re 2026-08-29）

来源：`feasibility/norman_benchmark/summary.json`（原 vs `reports/human_foundation_norman/full/`）。

**结果**：Norman 扰动数据上，linear_additive（Pearson 0.887, RMSE 0.045）胜出全部基座模型——scfoundation 0.824、scgpt 0.823、geneformer 0.800。Gate：`foundation_mean_over_linear = NO_GO`。

**对 re 的含义（如实）**：
- scF 在扰动响应预测上不优于线性基线——这是已注册的负结果，不因 re 使用 scF 而消失。
- re 选 scF 的依据不是 Norman 基准，而是 re 自己的五路线五任务评估（feasibility/transfer_routes/RESULTS_v3.md）：B2（scF + ESM2 注入）在迁移任务上胜出（T2 AUC 0.820 / T3 Spearman 0.168）。
- 两层证据串联：Norman 定义了基座模型在扰动预测上的天花板；五路线评估证明了 scF 基因表作为迁移表示空间的功能——后者是 re 管线的基础，前者不否定它。

## Legacy 六路线冻结表征审计（vs 项目完成，导入 re 2026-08-29）

来源：`yeastbridge_vs/docs/LEGACY_EFFECT_AUDIT.md`。

**结果**：六种冻结表征（esm2_mean / routea_init / routea_scgpt_ft / routeb_protein / routec_scgpt / routed_scyeast）在 T2/T3 上的公平比较。esm2_mean 单模型即接近最强，routeb_protein（旧蛋白桥接路线）在多任务上稳定。融合增益通过防作弊检验（OOF/val 锁定→测试集一次评分）。

**对 re 的含义**：B 路线（蛋白桥接）的可行性独立于 scF 选择而有 legacy 审计支撑。

## 方向感知管线设计（2026-08-29 注册）

vs 的 CRC target scan 为每个靶点输出了 `intended_direction`（activate/inhibit，来自 Becker GSE201348 疾病态 vs 期望态对比）和 `state_effect_disease_minus_desired` 数值。此方向信息在本管线中作为**结构属性贯穿全程**：

1. **任务定义层**：每个任务携带其来源靶点的 CRC 方向和状态效应值。
2. **ChEMBL 靶点注释层**：提取 `action_type`（AGONIST/INHIBITOR/BLOCKER 等），与 CRC 期望方向比对——只有方向兼容的注释对进入下游。
3. **对接层**：对接不区分激动/拮抗（如实标注为结构限制），仅作物理合理性门。
4. **药效学层**：pIC50 不编码方向（结构限制），仅确认结合亲和力。
5. **交付层**：只有化合物药理方向与 CRC 期望方向**匹配**的对进入交付包。方向不可判定的对如实排除。

**方向判定规则**：
- `MATCH`：action_type 为 AGONIST/ACTIVATOR/POSITIVE MODULATOR 且 CRC 方向为 activate；或 action_type 为 INHIBITOR/ANTAGONIST/BLOCKER/NEGATIVE MODULATOR 且 CRC 方向为 inhibit。
- `OPPOSITE`：方向相反。
- `no_data`：ChEMBL 无 action_type 记录——**不进入交付包**（方向未确认）。
