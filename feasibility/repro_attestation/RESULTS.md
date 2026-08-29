# 无变化证明结果（repro_attestation）

执行 2026-08-29，协议见本目录 REGISTRATION.md。`product/` 零写入。机器可读记录：`attestation.json`。

## 裁决总表

| 项 | 结果 |
|---|---|
| L0 manifest 全项校验 | **117/117 通过，0 失败**（含 44GB X-Atlas 全批次、NORMAN h5ad、Becker h5ad、UniProt 快照、orthodb） |
| 交付链输入 bit 稳定性 | **通过**：state_signature / candidate_baseline / target_universe 三文件与 vs@7bdaf4a production_v1 输出 bit 级一致 |
| L1 重跑保真度 | state_signature、target_universe：**bit 级一致**；candidate_baseline：数值等价 + 已归档偏差（见下） |
| 检索层（Layer 1） | 重跑结果 bit 级一致（`retrieval_result.json` sha256 `ceb025eb…` 前后不变） |
| 三基座特征重提取 | scgpt / scfoundation / geneformer **3/3 bit 级一致**（scF 与 vs 原始 provider 记录三方吻合） |

## candidate_baseline 重跑偏差（如实归档）

重跑（`product/repro/crc_scan_v1/`）与交付链输入（`product/target_scan/inputs/`）在 1,177 行 × 45 列上：

1. **浮点末位噪声**：17 个数值列 + 4 个数值字符串列共 569 个格子，最大绝对偏差 **2.22×10⁻¹⁶**（机器精度量级；多线程浮点求和顺序差异）。
2. **16 行酵母 ortholog 注释差异**：重跑映射到了 ortholog（如 CACNA1C→YHL017W;YKL039W），交付版未映射。根因：**vendored orthodb 版本不同**（re 入册 `9daadfbd…` vs vs production 使用 `d7f7680f…`）。涉及 16 靶点：TMEM87A、SLC26A6、CACNA1C/D/F/S、GPHRA/B、KCNAB3、NALCN、SCN11A、SCN4A、TMEM63A/B/C、UCP1。
3. **与交付 4 靶点（ADRA2C/ADRA2B/KCNK2/OPRM1）交集为空**——偏差不影响任何已交付结果。

## 结论

交付链消费的全部输入与其上游来源 bit 级一致且持续可校验；re 的 L1 重跑完整复现两份产物、第三份数值等价（偏差已定位、量化、归档，且与交付无关）。**无变化证明：通过。**

## 运行时登记

`feasibility/runtime_probes/runtime.json`：obabel 3.2.1、meeko 0.7.1、rdkit 2026.03.1、gemmi 0.6.5、numpy 1.26.4、pandas 2.3.2、torch 2.6.0+cu124、CPU Vina f458505-mod（1.2.x）、Vina-GPU 2.1 二进制哈希、scFoundation checkpoint SHA-256、ESM2 checkpoint 清单。
