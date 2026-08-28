# product/chembl_branch — ChEMBL 注释对的物理验证支线（GPU docking + MD）

## 定位与来源（用户指令 2026-08-28）
- 主线（target_screen 全库反向对接）占用 CPU 期间，本支线用 **GPU**（Vina-GPU 2.1）运行，避免资源竞争。
- 验证对象：有 ChEMBL **定量注释**的化合物×已知靶点对（pChEMBL≥5，每化合物取 top5）。
- 主线的 MD 在主线对接完成后另行补充（注册待执行，见主线 DESIGN 4b）。

## 流程（参数全部在 configs/chembl_branch.json）
1. pairs：从 `product/drug_annotation/results/chembl_targets.tsv` 取定量对，按文献标准筛选（文献依据 2026-08-28：pChEMBL≥5 为公认 active 线、≥6 为高置信集惯例；无文献采用每化合物 top-N 筛实验注释对——原 top5 规则已废除）：**pChEMBL≥6 + 物种限定 Homo sapiens + 每对取最强值 + 靶点限定 ion-channel/GPCR universe**（2026-08-28 修正：初版未限定 universe，导致 CYP 代谢酶注释涌入打分榜——范围修正后废除）；解析靶点 UniProt accession（ChEMBL target API），生成 `inputs/branch_targets.fasta` + `results/branch_pairs.tsv`；无蛋白组分、非人源、无制备配体的对如实跳过并计数。
2. 结构盘点：复用 `scripts/target_screen_inventory.py --config configs/chembl_branch.json`（同一脚本，独立目录）。
3. 受体制备：复用 `scripts/target_screen_prep_receptors.py`（max_workers 压低，避免与主线抢 CPU）。
4. GPU 对接：`scripts/chembl_branch_dock.py` —— Vina-GPU 2.1（OpenCL，GPU 选卡规则：自动选空闲显存最大的卡，不挤既有进程），逐靶点×top-3 口袋，批量对接该靶点的全部支线化合物；结果 `results/gpu_dock.tsv`。
5. MD：对接后按亲和能选取 top 对，走 EXP2b 多副本协议（后续执行）。

## 注册边界与限制
- 支线引擎为 Vina-GPU 2.1：与主线 CPU Vina 同打分函数系但搜索范式不同（thread/search_depth，无 exhaustiveness）；支线内部分析自洽，不与主线分数混用。
- 注释靶点不限于离子通道/GPCR universe（ChEMBL 注释如实），家族构成在结果中报告；与 universe 交集的对单列。
- docking/MD 只评估物理合理性，不构成活性结论；湿实验为最终仲裁。

## 治理
注册制：本文件 + config 先行提交；各阶段断点续跑；产物入 results/。
