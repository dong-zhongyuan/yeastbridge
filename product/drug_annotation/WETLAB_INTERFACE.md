# 湿实验接口 v1（2026-08-28，酵母体系）

**验证对象：高等生物功能目标迁移到酵母这一步**——即化合物在酵母菌株上是否产生其任务（由人靶点迁移而来）所预测的敏感谱。人靶点注释只作溯源背景，不是检测对象。

主交付 `results/wetlab_interface_v1.tsv`；来源为已完成的 ChEMBL 药理注释通道，结构通道（反向对接）完成后发 v2 补充物理证据列。

## 实验设计

每个接口对（化合物 × 任务）：

1. **菌株组**：任务的 top15 酵母菌株（ORF 删除菌株，按迁移排名）+ bottom5 特异性对照（任务排名末位菌株，预测无/弱效应）。
2. **处理**：化合物 3–5 个浓度，覆盖库内测试剂量（表内 `dose_uM` 为证据剂量）；每株 ± 药、含 vehicle 对照。
3. **读出**：生长曲线（96/384 孔板读数）或点滴稀释；报告相对生长率（处理/vehicle）。
4. **判定**：
   - 方向：top 菌株的效应方向应与 `predicted_direction` 一致（z<0 = 超敏感 = 生长受抑制更强；z>0 = 抗性方向）；
   - 强度排序：跨 top15 菌株的实测效应排名与 `task_rank` 的一致性（Spearman）；
   - 特异性：bottom5 对照株应显著弱于 top15（对照成立时迁移步骤获得验证）。

## 主表列（wetlab_interface_v1.tsv，长表：每行 = 对 × 菌株）

tier（1=ChEMBL 注释与任务靶点收敛 / 3=无注释新假设）、priority、inchikey / chembl_id / smiles（订购与身份）、task_target（任务来源的人靶点基因，仅溯源）、family、rho / q（任务执行强度与显著性）、dose_uM（证据剂量）、orf、task_rank / task_cosine（该菌株在任务中的排名与迁移得分）、compound_z（该菌株对化合物的实测 z，来自 HIP/HOP 库数据）、predicted_direction、control_role（top_task / bottom_task_control）。

## 优先级与档位（收敛要求定量活性）

| 档 | 定义 | 数量 |
|---|---|---|
| 1 收敛 | ChEMBL **定量活性**（有 pChEMBL 值）记录该化合物作用于其任务靶点 | **0 对** |
| 2 偏移 | 有 ChEMBL 注释但与任务靶点无定量交集 | 135 化合物 |
| 3 新假设 | 无任何 ChEMBL 蛋白靶点注释 | 349 化合物 |

**诚实披露**：初版曾把 6 对计为收敛（吉西他滨 CHEMBL888 对 ADRA1D/ADRA2C/ADRB1/ADRB3/DRD4；多柔比星 CHEMBL53463 对 ADRB3），核查发现这些活性记录**全部无数值**（典型为反筛/无活性沉积），已降为 tier-2，并在 `annotation_pairs.tsv` 的 `chembl_task_mention_only` 列标记。两个药物的真实定量药理（如吉西他滨→SLC29A1 pChEMBL 6.72）均不在任务靶点上。

**如实结论**：没有任何化合物的已知定量药理与其任务靶点收敛——这正是"迁移筛选面向无同源/无注释化合物"设定的直接体现。接口实验重心 = **tier-3 top50**（按 rho 降序，领头 STING1 / GPR149 / RYR3 来源任务，rho 至 0.273）；全部 7,100 对与 484 化合物明细在 `annotation_pairs.tsv` / `compound_summary.tsv`。

## 溯源背景（非检测对象）

`target_context_v1.tsv` 给同批对的人靶点信息（基因、家族、pchembl、涉多靶计数）。人靶点上的直接药理实验（膜片钳/第二信使等）属于后续可选验证，不在本接口范围。

## 主张边界（注册）

化合物是迁移任务的候选效应因子；本接口检验的是迁移与执行两步在酵母里的可重复性与特异性；人靶点身份以注释/结构通道为旁证；湿实验为最终仲裁。

## 复现

`/public/home/mengxl/dzy/envs/yeastbridge/bin/python scripts/product_annotate_targets.py`
（config: configs/drug_annotation.json；上游：execute_hiphop 显著对、chembl_targets.tsv 回填完成版、transfer_route_b 任务排名、hiphop strain_response.npz 菌株 z 谱）
