# 对接认证层结果（docking_qualification）

执行 2026-08-29，协议见本目录 REGISTRATION.md。`product/` 零写入；全部产物在本目录 results/。

## 项 2：交付包自一致性重对接（4/4 通过）

方法：交付包内 `ligand.smi` 重新制备配体（obabel gen3d + meeko gasteiger），从交付 `complex.pdb` 提取受体（ATOM 记录，obabel `-xr -p` 加氢），box 取交付 `ligand_pose.sdf` 质心 + 包围盒 + 8 Å padding（最小边 20 Å），CPU Vina（exhaustiveness 8，seed 20260829）。判定容差 ±1.0 kcal/mol（预声明）。

| 包 | 靶点 | 重对接亲和能 | 交付亲和能 | 绝对偏差 | 判定 | pose 质心距离 |
|---|---|---:|---:|---:|---|---:|
| YeastBridge_Mol001 | ADRA2C | −9.822 | −9.8 | 0.022 | 通过 | 6.34 Å |
| YeastBridge_Mol002 | ADRA2B | −8.892 | −8.8 | 0.092 | 通过 | 2.47 Å |
| YeastBridge_Mol003 | KCNK2 | −8.282 | −8.2 | 0.082 | 通过 | 2.17 Å |
| YeastBridge_Mol004 | OPRM1 | −7.871 | −7.4 | 0.471 | 通过 | 3.60 Å |

**结论：全部 4 个交付对的对接亲和能可被独立重算复现（最大偏差 0.471 kcal/mol，远低于容差），重对接 pose 与交付 pose 位于同一口袋。**

机器可读记录：`results/delivery_redock_check.json`。

## 项 1：工具链与运行时探针

见 `feasibility/runtime_probes/runtime.json`（由 `scripts/repro_attest.py --stage runtime` 产出）：Vina-GPU 2.1 二进制哈希、CPU Vina 版本、rdkit/meeko/gemmi/numpy/pandas/torch 版本、obabel 版本、scFoundation checkpoint SHA-256。

## Gate

- 项 2（自一致性）：**PASS**（4/4）。
- 项 1（工具链登记）：**PASS**（全部可执行且版本登记）。

边界：本层证明「交付对接结果可独立复算」，不构成新的生物学裁决。
