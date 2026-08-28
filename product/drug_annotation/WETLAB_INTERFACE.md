# 湿实验接口 v1（2026-08-28）

第一批可执行交付物，来自药理注释通道（ChEMBL 直查，已完成）。结构通道（反向对接 + MD）进行中，完成后发 v2 补充物理证据列。

## 交付文件

- `results/wetlab_interface_v1.tsv` — 主表（见列说明）
- `results/annotation_pairs.tsv` — 全部 7,100 个显著化合物×任务靶点对及分档（不截断）
- `results/compound_summary.tsv` — 484 个化合物概览（档位构成、ChEMBL 已知靶点 top5、SMILES）
- `results/annotation_summary.json` — 汇总统计

## 主表列

tier（1=收敛/2=偏移/3=新假设）、priority（档内优先级）、inchikey、chembl_id（可订购索引）、task_target（人靶点基因）、family（gpcr / ion_channel）、rho / q（酵母任务执行强度与显著性）、pchembl_task（ChEMBL 上该对的最大 pChEMBL，仅 tier1）、n_chembl_genes（化合物已知靶点数，涉多靶提示）、dose（HIP/HOP 测试剂量 µM）、smiles（来自原始筛库）、suggested_assay（按家族的检测建议）。

## 三档定义与数量

| 档 | 定义 | 数量 |
|---|---|---|
| 1 收敛 | ChEMBL 记录该化合物作用于其执行任务的同一靶点 | **6 对 / 2 化合物** |
| 2 偏移 | 化合物有 ChEMBL 注释，但与执行任务靶点无交集 | 135 化合物 |
| 3 新假设 | 化合物无任何 ChEMBL 蛋白靶点注释，任务靶点即为新假设 | 349 化合物 |

## tier-1 全部六对（最硬验证候选）

| 化合物 | 靶点 | 家族 | rho | q | pchembl |
|---|---|---|---|---|---|
| CHEMBL888 | ADRA1D | gpcr | 0.049 | 0 | 已注释 |
| CHEMBL888 | ADRA2C | gpcr | 0.058 | 0 | 已注释 |
| CHEMBL888 | ADRB1 | gpcr | 0.053 | 0 | 已注释 |
| CHEMBL888 | ADRB3 | gpcr | 0.055 | 0 | 已注释 |
| CHEMBL53463 | ADRB3 | gpcr | 0.045 | 0 | 已注释 |
| CHEMBL888 | DRD4 | gpcr | 0.048 | 0 | 已注释 |

如实说明：tier-1 的酵母任务执行强度（rho 0.045–0.058）偏弱；任务执行最强的对都在 tier-3（rho 最高 0.273，STING1/GPR149/RYR3 等）。目标身份证据与作用强度在此数据中解耦——tier-1 适合直接做受体药理验证，tier-3 的强对等待结构通道交叉后再进入实验。

## 优先级与检测建议

- tier-1 按 pchembl 降序、q 升序；tier-3 按 rho 降序、q 升序取 top50（全量在 annotation_pairs.tsv）。
- ion_channel：膜片钳（手动/自动）或离子流检测；gpcr：第二信使（cAMP / IP1 / Ca²⁺，按受体偶联）或 β-arrestin 募集。

## 主张边界（注册）

tier-1 为药理学与任务执行收敛的验证候选；tier-3 为新假设，结构通道对接/FDR 与 MD 物理合理性评估在 v2 补充；化合物是 transferred task 的候选效应因子，靶点身份以湿实验为最终仲裁。

## 复现

`/public/home/mengxl/dzy/envs/yeastbridge/bin/python scripts/product_annotate_targets.py`
（config: configs/drug_annotation.json；上游：execute_hiphop 显著对 + chembl_targets.tsv 回填完成版）
