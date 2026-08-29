# 模型选择证据链（feasibility）

## 管线顺序索引

| 管线阶段 | 证据层 | 状态 |
|---|---|---|
| 0. 数据来源与无变化证明 | 数据层（L0/L1 vendoring + repro_attestation） | 已执行，登记中 |
| 1. 基座选择 | 第一层 Norman retrieval（scF 胜出） | 已执行（re 自有） |
| 2. 迁移策略选择 | 第二层 五路线评估（B2 ESM2 注入） | 已执行（re 自有） |
| 3. 靶点全集与疾病签名 | 数据层：crc_scan L1 重跑 | 已执行，登记中 |
| 4. HIP/HOP 执行 | HIP/HOP 执行层（置换 + BH-FDR） | 已执行（re 自有），登记中 |
| 5. 构象状态选择 | 构象层（EXP3 受控重对接复现） | 已注册，执行中 |
| 6. 对接可信度 | 对接认证层（工具链 + 交付包自一致性） | 已注册，执行中 |
| 7. 方向贯穿 | 方向感知管线设计 | 已注册（2026-08-29） |

Legacy 六路线冻结表征审计作为历史导入层保留于文末。

## 第一层：基座选择 — Norman retrieval 基准（re 项目自有结果）

来源：`feasibility/norman_foundation/retrieval/retrieval_result.json`。

**任务**：哪个基座的预训练基因表示能正确检索未见扰动的程序方向（6 程序，80 基因，随机基线 0.167）。

| 基座 | accuracy | permutation p |
|---|---|---|
| **scfoundation** | **0.400** | 0.000 |
| scgpt | 0.350 | 0.0004 |
| geneformer | 0.325 | 0.0018 |

**Gate：PASS。scF 胜出——这是 re 选择 scFoundation 的直接依据。**

输入依赖：`raw_oof_full.npz`（`product/norman_assets/`，L0 manifest 入册）与 `features/*.npz`。驱动脚本 `scripts/norman_retrieval.py`。

## 第二层：微调策略选择 — 五路线评估（re 项目自有结果）

来源：`feasibility/transfer_routes/RESULTS_v3.md`。

在同一 scF 基座上比较微调路线（A2 正交初始化 / B2 ESM2 注入 / C2 随机初始化 / D' scYeast / E' SGA 图原生），在原框架五任务上公平竞争。

**结果**：B2（scF + ESM2 注入）胜出（T2 AUC 0.820 / T3 Spearman 0.168）；A2/C2 处于随机水平（scF 的 MaeAutobin 值重构不组织基因表，只有 ESM2 注入有效）；D' scYeast 单独作转移路线落败但保留酵母状态表征角色；E' 保留图原生用途。

## 数据层：vendoring 与无变化证明（登记中）

来源：`product/repro/crc_scan_v1/`（内容寻址重跑，seed 42，2026-08-29）与 `feasibility/repro_attestation/`（待产出）。

L0 原始输入（Becker h5ad、UniProt 快照、X-Atlas 批次、Norman h5ad、orthodb）以 SHA-256 入册；crc_scan L1 重跑产出 state_signature / candidate_baseline / target_universe（1,177 靶点）。无变化证明协议见 `feasibility/repro_attestation/REGISTRATION.md`。

## HIP/HOP 执行层（re 项目自有结果，登记中）

来源：`product/execute_hiphop/RESULTS.md` 与 `execute_summary.json`。

全谱 Spearman + 菌株标签置换检验 + BH-FDR：7,100 对 q<0.1，覆盖 1,038 靶点、484 化合物。此层为 re 自有执行，作为筛选统计证据登记。

## 构象层：EXP3 受控重对接复现（已注册，执行中）

上游证据：yeastbridge_vs@7bdaf4a EXP3 三臂筛选（参照面板构象偏好 arms B/C 5/5 vs arm A 1/5，sign test p=0.031，解锁「direction-aware conformation selection」，边界为构象偏好层面）。

re 侧复现协议：`feasibility/conformation_selection/REGISTRATION.md`。产出待 `RESULTS.md`。

## 对接认证层（已注册，执行中）

协议：`feasibility/docking_qualification/REGISTRATION.md`。工具链探针 + 交付包自一致性重对接，产出待 `RESULTS.md`。

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
