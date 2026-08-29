# 模型选择证据链（feasibility）

## 管线顺序索引

| 管线阶段 | 证据层 | 状态 | 裁决 |
|---|---|---|---|
| 0. 数据来源与无变化证明 | 数据层（L0/L1 vendoring + repro_attestation） | 已执行 | **通过**（117/117 manifest；交付链输入 bit 稳定） |
| 1. 基座选择 | 第一层 Norman retrieval（scF 胜出） | 已执行（re 自有） | PASS |
| 2. 迁移策略选择 | 第二层 五路线评估（B2 ESM2 注入） | 已执行（re 自有） | PASS |
| 3. 靶点全集与疾病签名 | 数据层：crc_scan L1 重跑 | 已执行 | 两份 bit 一致；candidate_baseline 数值等价（偏差已归档） |
| 4. HIP/HOP 执行 | HIP/HOP 执行层（置换 + BH-FDR） | 已执行（re 自有） | 7,100 对 q<0.1 / 484 化合物 |
| 5. 构象状态选择 | 构象层（EXP3 受控重对接复现） | 已执行（re 自有） | **PASS**（armB/C 5/5，p=0.031） |
| 6. 对接可信度 | 对接认证层（工具链 + 交付包自一致性） | 已执行 | **PASS**（4/4） |
| 7. 方向贯穿 | 方向感知管线设计 | 已注册（2026-08-29） | — |

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

溯源闭环（2026-08-29）：输入 `raw_oof_full.npz`（`product/norman_assets/`，L0 manifest 入册，内容为真值矩阵与条件/基因清单）；三基座特征由原始 checkpoint 在 re 内重新提取（`features_reextract/`，scgpt/scF/geneformer **3/3 与 vendored 副本 bit 级一致**，scF 并与 vs 原始 provider 记录三方吻合）；检索计算重跑结果 bit 级一致（sha256 `ceb025eb…`）。驱动：`scripts/norman_retrieval.py`。

## 第二层：微调策略选择 — 五路线评估（re 项目自有结果）

来源：`feasibility/transfer_routes/RESULTS_v3.md`。

在同一 scF 基座上比较微调路线（A2 正交初始化 / B2 ESM2 注入 / C2 随机初始化 / D' scYeast / E' SGA 图原生），在原框架五任务上公平竞争。

**结果**：B2（scF + ESM2 注入）胜出（T2 AUC 0.820 / T3 Spearman 0.168）；A2/C2 处于随机水平（scF 的 MaeAutobin 值重构不组织基因表，只有 ESM2 注入有效）；D' scYeast 单独作转移路线落败（T2 0.738 / T3 0.149）但保留酵母状态表征角色（second_round 双基座设计）；E' 保留图原生用途。

## 数据层：vendoring 与无变化证明（2026-08-29 执行）

来源：`product/repro/crc_scan_v1/`（内容寻址重跑，seed 42）与 `feasibility/repro_attestation/`（协议 + attestation.json + RESULTS.md）。

- **L0 manifest**：117/117 项通过（含 44GB X-Atlas 全批次、NORMAN h5ad、Becker h5ad、UniProt 快照、orthodb、RAW_OOF）。
- **交付链输入 bit 稳定**：state_signature / candidate_baseline / target_universe 三文件与 vs@7bdaf4a production_v1 输出 bit 级一致——交付所依赖的输入自上游冻结以来零变化。
- **L1 重跑保真度**：state_signature、target_universe 与交付输入 bit 一致；candidate_baseline 数值等价（浮点末位噪声 max 2.2×10⁻¹⁶ + 16 行酵母 ortholog 注释差异，根因为 vendored orthodb 版本不同，与交付 4 靶点零交集）。
- 运行时登记：`feasibility/runtime_probes/runtime.json`。

## HIP/HOP 执行层（re 项目自有结果）

来源：`product/execute_hiphop/RESULTS.md` 与 `execute_summary.json`。

全谱 Spearman + 菌株标签置换检验 + BH-FDR：7,100 对 q<0.1，覆盖 1,038 靶点、484 化合物。此层为 re 自有执行，作为筛选统计证据登记。

## 构象层：EXP3 受控重对接复现（re 项目自有结果，2026-08-29 执行）

上游证据：yeastbridge_vs@7bdaf4a EXP3 三臂筛选（armA 1/5、arms B/C 5/5，sign test p=0.031）。协议、库（4,196 配体）、受体（6 构象）、box（含 7CR0 Kabsch 转移冻结产物）全部只读复用上游冻结资产；引擎按注册声明替换为 Vina-GPU 2.1。

结果（`feasibility/conformation_selection/results/RESULTS.md`）：

| 臂 | 上游 | re 复现 |
|---|---|---|
| A 默认构象正确率 | 1/5 | **1/5** |
| B 方向匹配 delta 正确率 | 5/5 | **5/5** |
| C 顶构象正确率 | 5/5 | **5/5** |

精确二项 p = 0.0312。**Gate：PASS。**这是 re 管线构象状态选择（抑制→非活性构象、激活→活性构象）的直接受控实验依据，声明边界与上游一致（仅构象偏好层面）。次级富集不增益（与上游定性一致）。引擎损耗（大配体被拒，保留 58–61%，MW 偏向，构象间一致）已量化归档。

## 对接认证层（re 项目自有结果，2026-08-29 执行）

来源：`feasibility/docking_qualification/RESULTS.md`。工具链版本与 checkpoint 哈希登记 + 交付包自一致性重对接：4 个交付对（ADRA2C −9.822/−9.8、ADRA2B −8.892/−8.8、KCNK2 −8.282/−8.2、OPRM1 −7.871/−7.4，重对接/交付），最大偏差 0.471 kcal/mol ≤ 1.0 容差。**Gate：PASS（4/4）**——交付对接结果可被独立重算复现。

## Legacy 六路线冻结表征审计（vs 项目完成，导入 re 2026-08-29）

来源：`yeastbridge_vs/docs/LEGACY_EFFECT_AUDIT.md`。

六种冻结表征公平比较（T2/T3，防作弊检验）。esm2_mean 单模型即接近最强，routeb_protein 在多任务上稳定。

## 方向感知管线设计（2026-08-29 注册）

vs 的 CRC target scan 为每个靶点输出了 `intended_direction`（activate/inhibit）和 `state_effect_disease_minus_desired`。此方向在本管线中作为**结构属性贯穿全程**：

1. **任务定义层**：每个任务携带其来源靶点的 CRC 方向和状态效应值。
2. **ChEMBL 靶点注释层**：提取 `action_type`，与 CRC 期望方向比对——只有方向兼容的注释对进入下游。
3. **对接层**：构象状态选择依据方向执行（见构象层证据）。
4. **药效学层**：不编码方向（结构限制），仅确认结合亲和力。
5. **交付层**：只有化合物药理方向与 CRC 期望方向**匹配**的对进入交付包。方向不可判定的对如实排除。

**方向判定规则**：
- `MATCH`：action_type 为 AGONIST/ACTIVATOR/POSITIVE MODULATOR 且 CRC 方向为 activate；或 INHIBITOR/ANTAGONIST/BLOCKER/NEGATIVE MODULATOR 且 CRC 方向为 inhibit。
- `OPPOSITE`：方向相反。
- `no_data`：ChEMBL 无 action_type——**不进入交付包**。
