# 对接认证层预注册（docking_qualification）

Registered 2026-08-29。执行前冻结。上游对应：vs CONF-01 面板认证（11/11 cognate 重对接通过）与工具链 PARTIAL_GO。re 版按交付实际收缩为两项。

## 项 1：工具链探针

登记 re 管线实际使用的工具与版本：Vina-GPU 2.1 二进制、CPU Vina 1.2.7、Meeko、RDKit、OpenBabel、gemmi；scF 与 ESM2 checkpoint 的 SHA-256。产出并入 `feasibility/runtime_probes/runtime.json`。判定：版本可执行且哈希与 configs 登记一致 → PASS。

## 项 2：交付包自一致性重对接

对 4 个交付对（ADRA2C、ADRA2B、KCNK2/TREK-1、OPRM1）逐一执行：

- 输入：交付包内 `ligand.smi`（重新制备 pdbqt，制备参数与原管线登记一致）+ 交付包 `complex.pdb` 提取的受体口袋（只读）。
- 引擎：CPU Vina 1.2.7，固定种子，exhaustiveness 与原管线登记一致。
- 输出：只写 `feasibility/docking_qualification/results/`。
- 判定（预声明）：重对接亲和能与交付 `scores.tsv` 记录值的偏差 ≤ 1.0 kcal/mol，且重对接 pose 与交付 `ligand_pose.sdf` 的重原子 RMSD 有可比性（同口袋内）→ 自一致性 PASS。偏差超过阈值 → 如实记录并报告，不修改交付包。

## 边界

本层检验的是「交付对接结果可被独立复算」，不是新的生物学裁决；`product/` 零写入。
