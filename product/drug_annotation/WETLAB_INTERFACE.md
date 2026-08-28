# 湿实验接口 v1（2026-08-28，酵母体系）

**验证对象：高等生物功能目标迁移到酵母这一步**——即化合物在酵母菌株上是否产生其任务（由人靶点迁移而来）所预测的敏感谱。人靶点注释只作溯源背景，不是检测对象。

主交付 `results/wetlab_interface_v1.tsv`；来源为已完成的 ChEMBL 药理注释通道，结构通道（反向对接）完成后发 v2 补充物理证据列。

## 实验设计（完整操作参数）

每个接口对（`pair_id` 标识，共 50 对 / 16 个化合物）：

1. **菌株组**：任务的 top15 酵母删除菌株（ORF，按迁移排名）+ bottom5 特异性对照。`strain_collection` 列已标注每个 ORF 应使用的删除库类型：essential → 杂合二倍体删除株（HIP 型）；nonessential → 同合二倍体删除株（HOP 型）——与 Lee et al. 2014 筛选平台同源。`gene_common_name` 给通用基因名。
2. **化合物**：按 `compounds_order_v1.tsv` 订购（16 个，名称/ChEMBL ID/SMILES/InChIKey；无名称的按 SMILES 订购）。
3. **处理**：每个化合物 4 个浓度 = 证据剂量（`dose_uM`）的 0.5×、1×、2×、4×；每个浓度 × 每个菌株 ≥3 个生物学重复；含 vehicle（DMSO 对应体积）与无药对照。
4. **条件与读出**：YPD 或相应缺陷培养基（如删除标记需要），30 °C，96/384 孔板生长曲线（OD600，每 20–30 min，24–48 h）或系列稀释点滴。报告**相对生长** = 药物孔曲线下面积 / 同批 vehicle 孔曲线下面积。
5. **判定标准**：
   - 主判据（与管线 v2 端点一致）：每对内跨 20 个菌株，|1 − 相对生长| 的实测排名与 `task_rank` 的 Spearman 相关；置换检验 p<0.05 为通过。
   - 次判据：top 菌株效应方向与 `predicted_direction` 一致性。
   - 特异性：top15 与 bottom5 的 |效应| 差异（如 Wilcoxon）。

## 数据回收

结果按 `results_template.csv` 的列填写（pair_id, inchikey, task_target, orf, dose_uM, replicate, relative_growth），返回后由管线侧完成上述全部统计判定并出验证报告。

## 库参考值缺失的说明

`compound_z` 是 HIP/HOP 库数据中该菌株×化合物的参考 z；barcode 数据存在结构性缺失，平均每对约 9/20 菌株无库量化值（`library_reference = not_quantified_in_library`）。**这些菌株不是错误，恰恰是湿实验要新鲜测量的部分**——主判据（实测效应排名 vs task_rank）在实验室实际测得数据的菌株上计算；有库参考的菌株可同时做方向一致性检查。

## 主表列（wetlab_interface_v1.tsv，长表：每行 = 对 × 菌株）

pair_id、tier、priority、inchikey / chembl_id / compound_name / smiles（订购与身份）、task_target（任务来源的人靶点基因，仅溯源）、family、rho / q（任务执行强度与显著性）、dose_uM（证据剂量）、orf / gene_common_name、essentiality / strain_collection（删除库选型）、task_rank / task_cosine（该菌株在任务中的排名与迁移得分）、compound_z（库数据中该菌株对化合物的实测 z）、predicted_direction、control_role（top_task / bottom_task_control）。

## 优先级与档位（收敛要求定量活性）

| 档 | 定义 | 数量 |
|---|---|---|
| 1 收敛 | ChEMBL **定量活性**（有 pChEMBL 值）记录该化合物作用于其任务靶点 | **0 对** |
| 2 偏移 | 有 ChEMBL 注释但与任务靶点无定量交集 | 135 化合物 |
| 3 新假设 | 无任何 ChEMBL 蛋白靶点注释 | 349 化合物 |

**诚实披露**：初版曾把 6 对计为收敛（吉西他滨 CHEMBL888 对 ADRA1D/ADRA2C/ADRB1/ADRB3/DRD4；多柔比星 CHEMBL53463 对 ADRB3），核查发现这些活性记录**全部无数值**（典型为反筛/无活性沉积），已降为 tier-2，并在 `annotation_pairs.tsv` 的 `chembl_task_mention_only` 列标记。两个药物的真实定量药理（如吉西他滨→SLC29A1 pChEMBL 6.72）均不在任务靶点上。

**如实结论**：没有任何化合物的已知定量药理与其任务靶点收敛——这正是"迁移筛选面向无同源/无注释化合物"设定的直接体现。接口实验重心 = **tier-3 top50**（按 rho 降序，领头 STING1 / GPR149 / RYR3 来源任务，rho 至 0.273）；全部 7,100 对与 484 化合物明细在 `annotation_pairs.tsv` / `compound_summary.tsv`。

**队列冻结**：本 top50 即交付实验队列。结构通道（对接/校准）完成后的 v2 只补充证据列，不重排此队列——湿验证依据是酵母侧证据，不与尚在运行的计算线耦合。

## 溯源背景（非检测对象）

`target_context_v1.tsv` 给同批对的人靶点信息（基因、家族、pchembl、涉多靶计数）。人靶点上的直接药理实验（膜片钳/第二信使等）属于后续可选验证，不在本接口范围。

## 主张边界（注册）

化合物是迁移任务的候选效应因子；本接口检验的是迁移与执行两步在酵母里的可重复性与特异性；人靶点身份以注释/结构通道为旁证；湿实验为最终仲裁。

## 复现

`/public/home/mengxl/dzy/envs/yeastbridge/bin/python scripts/product_annotate_targets.py`
（config: configs/drug_annotation.json；上游：execute_hiphop 显著对、chembl_targets.tsv 回填完成版、transfer_route_b 任务排名、hiphop strain_response.npz 菌株 z 谱）
