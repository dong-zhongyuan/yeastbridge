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

## MD 协议（2026-08-28 注册，阶段 5）
- **体系选取**：`gpu_dock_pairs.tsv` 按"每靶点取亲和能最优对"取前 5 个不同靶点（ADRA2C / DRD3 / CHRM1 / HTR2B / OPRL1，全部为沉积实验结构 2.17–2.80 Å）。
- **构建**：配体 pose 从 Vina-GPU 输出经 meeko 重建分子（对接坐标系即蛋白坐标系，直接并入蛋白 PDB）；antechamber AM1-BCC + GAFF2 电荷/参数，parmchk2 补缺；packmol-memgen 建膜体系（POPC 双层、ff14SB/lipid21/tip3p、0.15 M KCl、水层 15 Å、PPM3 跨膜定向、二硫键自动检测），产出 Amber parm7/rst7。
- **运行**：OpenMM 8.6（CUDA，空闲卡规则同支线对接），LangevinMiddle 310 K、MonteCarlo 半各向同性 1 atm、PME、2 fs 步长；能量最小化 → NPT 平衡 1 ns（蛋白/配体重原子位置限制 42 kJ/mol/nm² 前 0.5 ns，后 0.5 ns 释放）→ 3 副本 × **20 ns** 生产（种子注册于 config；v1 长度按单卡聚合吞吐 ~40-60 ns/天校准，2026-08-28 由 50 ns 下调，可从 checkpoint 扩展）。
- **分析**：蛋白/配体 RMSD、配体接触持久度（4 Å 逐帧）、姿势漂移；逐体系报告。
- **边界**：MD 评估物理合理性（姿势稳定性/接触持续性），不产出结合自由能排名；湿实验为最终仲裁。
- **执行环境**：`/public/home/mengxl/dzy/envs/mdenv`（conda-forge ambertools+openmm，pip 无 ambertools），构建用 structscreen env（meeko）+ mdenv 二进制混合调用。

## Scheduling (2026-08-28)
MD builds/runs PAUSED: CPU contention with the main-line docking and
future GPU contention with Tahoe inference outweigh MD urgency (wet-lab
interface already delivered; v2 waits on the main line only). Resume
builds when the main docking finishes; run MD on GPU0 after Tahoe
calibration or share with explicit accounting. Built artifacts under
results/md_systems/ remain valid; packmol phase restarts on resume.
