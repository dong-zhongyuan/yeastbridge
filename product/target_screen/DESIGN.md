# product/target_screen — 结构驱动的化合物→人靶点注释通道（反向筛选）

## 定位与范围
- 化合物中心管线的第二注释通道，与 `product/drug_annotation`（ChEMBL 药理注释）并行；对 ChEMBL 无记录的化合物（含后续接入的 RL 生成分子）是唯一可用通道。
- 靶点范围 = 迁移步的 1177 膜蛋白 universe（config 指定 fasta/tsv）。
- 配体输入 = `product/execute_hiphop` q<0.1 命中化合物 + 注册基准药；接口为 SMILES 列表，后续接入 RL 分子不改代码。
- 对海报版反向筛选方法的缺陷修复：删掉无校准的 AI 粗筛层（1177 规模直接全库对接）；阴性对照与阈值全部经验校准；AF2 结构加口袋级 pLDDT 过滤；结果经回顾性基准门禁。

## 流程（参数全部在 configs/target_screen.json）
1. inventory：按 UniProt accession 查 RCSB 沉积结构（entity 精确匹配，仅实验结构，取分辨率最高者）；无沉积 → AlphaFold DB 预测模型（prediction API 解析文件 URL，优先 PDB 格式直链；B 因子=pLDDT）。产物 `results/inventory.tsv`，结构文件入 `structures/raw/`（gitignore，脚本可再生）。
2. 配体制备：ChEMBL canonical SMILES（命中按 InChIKey、基准药按 pref_name 解析）→ dimorphite-dl pH 7.4 质子化 → RDKit ETKDG 单构象（注册种子，Vina 柔性搜索只需输入种子构象）→ meeko PDBQT。产物 `inputs/ligands/`。
3. 受体制备：gemmi 只保留聚合物（去水/杂原子/替代构象，首模型）→ fpocket 默认参数 Score top-N 口袋；AF2 来源口袋要求口袋残基平均 pLDDT ≥ 阈值；口袋盒 = 口袋原子包围盒 + padding（各维下限）；OpenBabel `-p 7.4 -xr` 转刚性受体 PDBQT。产物 `structures/receptors|pockets/`（gitignore）。
4. 反向对接：Vina python API，逐口袋建图后批对接全部配体（exhaustiveness=2、9 poses、cpu_per_job）；配体×靶点得分 = 各口袋最优；多进程按口袋并行，逐口袋 JSONL 断点续跑。产物 `results/target_scores.tsv`（全矩阵）+ `results/raw_pockets/`。
5. 校准（门禁，海报版缺失的核心层）：
   - 回顾性基准：5 个注册基准药，其已知 universe 靶点须进入该药全 universe 排名 top 5%；≥3 个可评估且 ≥2 达标才放行正式结果。
   - 每配体经验零分布：同配体 1177 靶点得分主体（median/MAD 高斯）为非结合背景，BH-FDR 定命中。
   - 次级 decoy 校验（非门禁）：基准药已知靶点得分 vs MW 匹配随机 ChEMBL 分子在家族分层靶点子集上的分位报告。
   - 产物 `results/screen_hits.tsv`、`results/calibration_report.json`。
6. 短名单 → CONF-01 制备 / EXP2b 多副本 MD（复用既有协议）。

## 注册的结构性限制（如实报告，不构成减分）
- AF2 单体模型无寡聚组装与脂环境；通道孔腔外侧开口等功能位点几何可能失真。口袋级 pLDDT 过滤缓解但不消除。
- Vina 真空打分无膜静电项；跨膜深口袋得分存在系统性偏差，以校准层经验阈值兜底。
- 沉积结构覆盖率在盘点前未知，盘点后如实写入 inventory.tsv。
- 经验零分布假设每配体的 1177 个得分中绝大多数为非结合背景。

## 治理
本文件与 configs/target_screen.json 先行提交后执行；五阶段各自断点续跑；网络失败目标不落盘、重跑补齐；执行环境 `/public/home/mengxl/dzy/envs/structscreen`（micromamba 建，fpocket 来自 bioconda，其余 pip）。
