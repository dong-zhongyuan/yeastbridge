# product/target_screen — 任务定向结构验证通道（docking）

## 定位与架构（2026-08-28 用户决策重定）
化合物中心管线的结构验证通道。**验证对象是 execute 步已指认的显著（化合物 × 任务靶点）对**：对接给出物理合理性，与 pd_validation（DeepPurpose 效价引擎）并联评估——**docking 门与 PD 门双过的对进入最终结果**。

架构变更记录：本模块原设计为"全库反向筛选"（465 配体 × 1177 universe 泛搜索，复刻海报方法并修复其缺陷）。2026-08-28 经用户决策废除该模式：任务靶点分配已由 execute 步给出，全库搜索重复回答同一问题且耗时/CPU 占用不成比例。部分对接产物已删除（git 历史留档）；基础设施（结构盘点、配体制备、受体制备、对接脚本）全部复用，改为 `pairs_file` 驱动的任务定向模式。

## 流程（参数全部在 configs/target_screen.json）
1. inventory（已完成，复用）：靶点结构盘点（沉积优先，AF2 兜底，口袋 pLDDT 过滤）。477 沉积 + 698 AF2。
2. 配体制备（已完成，复用）：SMILES 解析链（exec_matrix/ChEMBL/CACTUS/PubChem）→ 质子化 → PDBQT。465/489。
3. 受体制备（已完成，复用）：gemmi 清洗 → fpocket top-3 口袋 → obabel PDBQT。1169 可对接。
4. **任务定向对接**：仅对 `pairs_file`（annotation_pairs.tsv，7,100 显著对）中配体已制备的 (化合物, 任务靶点) 对接；每靶点 top-3 口袋，exhaustiveness=8，逐口袋 JSONL 断点续跑；配体×靶点得分 = 各口袋最优。判定门：`dock_gate_kcal`（默认 −7.0 kcal/mol，Vina 文献常规命中阈值）。
5. 基准药参考：注册基准药（nifedipine 等 5 个，见 config benchmarks）同流程对接其已知通道靶点，作为打分参考系（不再适用全库经验零分布——该层随全库模式一并废除）。
6. 与 pd_validation 汇合：双门判定产出最终 (化合物, 靶点) 结果表。

## 注册的结构性限制（如实报告）
- AF2 单体模型无寡聚组装与脂环境；通道界面口袋（Kv/Nav/LGIC 中央孔）在单体上不成立——受体制备的 pLDDT 过滤缓解但不消除；超大蛋白（RYR/PIEZO/ADGRV1 等 6 个）受体制备失败，相关对无对接分，如实留空。
- Vina 真空打分无膜静电；门限为经验阈值并标注基准药参考系。
- 对接门为物理合理性过滤，不是效价预测——效价归 pd_validation。

## 计算引擎
CPU Vina 1.2.7（本通道）；GPU Vina-GPU 2.1 保留用于支线（chembl_branch）与精修轮，不与主线任务争用。

## 治理
注册制；各阶段断点续跑；结果入 results/；湿实验接口 v1 队列冻结不受本通道任何结果影响（2026-08-28 注册）。
