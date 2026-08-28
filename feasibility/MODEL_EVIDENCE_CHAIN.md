# 模型选择证据链（feasibility）

## 第一层：基座选择 — Norman retrieval 基准（re 项目自有结果）

来源：`feasibility/norman_foundation/retrieval/retrieval_result.json`。

**任务**：哪个基座的预训练基因表示能正确检索未见扰动的程序方向（6 程序，80 基因，随机基线 0.167）。

| 基座 | accuracy | permutation p |
|---|---|---|
| **scfoundation** | **0.400** | 0.000 |
| scgpt | 0.350 | 0.0004 |
| geneformer | 0.325 | 0.0018 |

**Gate：PASS。scF 胜出——这是 re 选择 scFoundation 的直接依据。**

## 第二层：微调策略选择 — 五路线评估（re 项目自有结果）

来源：`feasibility/transfer_routes/RESULTS_v3.md`。

在同一 scF 基座上比较三种微调路线（A2 正交初始化 / B2 ESM2 注入 / C2 随机初始化），在原框架五任务上公平竞争。

**结果**：B2（scF + ESM2 注入）胜出（T2 AUC 0.820 / T3 Spearman 0.168）；A2/C2 处于随机水平（scF 的 MaeAutobin 值重构不组织基因表，只有 ESM2 注入有效）。

## 第三层：扰动响应预测 — Norman full regression（re 项目复现了 vs 结果）

来源：`feasibility/norman_foundation/full/summary.json`。

linear_additive（Pearson 0.887）胜出全部基座模型（scF 0.824、scgpt 0.823、geneformer 0.800）。Gate：`foundation_mean_over_linear = NO_GO`。

**含义（如实）**：扰动响应预测任务上线性加法模型足够强，基座模型在该任务上无显著增益。这与 scF 在 retrieval 任务上的胜出不矛盾——前者测的是"预测响应值"，后者测的是"基因表示携带程序方向信息"。re 用 scF 做的是后者（蛋白→酵母基因表迁移），不是前者。

## Legacy 六路线冻结表征审计（vs 项目完成，导入 re 2026-08-29）

来源：`yeastbridge_vs/docs/LEGACY_EFFECT_AUDIT.md`。

六种冻结表征公平比较（T2/T3，防作弊检验）。esm2_mean 单模型即接近最强，routeb_protein 在多任务上稳定。

## 方向感知管线设计（2026-08-29 注册）

vs 的 CRC target scan 为每个靶点输出了 `intended_direction`（activate/inhibit）和 `state_effect_disease_minus_desired`。此方向在本管线中作为**结构属性贯穿全程**：

1. **任务定义层**：每个任务携带其来源靶点的 CRC 方向和状态效应值。
2. **ChEMBL 靶点注释层**：提取 `action_type`，与 CRC 期望方向比对——只有方向兼容的注释对进入下游。
3. **对接层**：不区分激动/拮抗（结构限制），仅作物理合理性门。
4. **药效学层**：不编码方向（结构限制），仅确认结合亲和力。
5. **交付层**：只有化合物药理方向与 CRC 期望方向**匹配**的对进入交付包。方向不可判定的对如实排除。

**方向判定规则**：
- `MATCH`：action_type 为 AGONIST/ACTIVATOR/POSITIVE MODULATOR 且 CRC 方向为 activate；或 INHIBITOR/ANTAGONIST/BLOCKER/NEGATIVE MODULATOR 且 CRC 方向为 inhibit。
- `OPPOSITE`：方向相反。
- `no_data`：ChEMBL 无 action_type——**不进入交付包**。
