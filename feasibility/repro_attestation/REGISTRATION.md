# 无变化证明预注册（repro_attestation）

Registered 2026-08-29。执行前冻结。

## 目标

证明 re 的数据层与已交付产品结果之间不存在任何变化：原始输入未被改动，L1 重跑与交付链使用的输入 bit 级一致，上游（vs@7bdaf4a）原输出与 re 重跑输出一致。

## 协议

1. **L0 manifest 全项校验**：`data/external/crc_scan_raw/MANIFEST.sha256`（含补齐后的 X-Atlas 全批次、NORMAN_PERTURB_PINNED.h5ad、orthodb、Becker h5ad、UniProt 快照）逐文件重算 SHA-256 并比对。任何一项不匹配 → FAIL 并停止，不自动修复。
2. **L1 位级比对（三方）**：
   - A = `product/repro/crc_scan_v1/{state_signature,candidate_baseline,target_universe}.tsv`（re 重跑输出）
   - B = `product/target_scan/inputs/` 同名三文件（交付链实际消费的输入）
   - C = vs 上游 crc_scan 原输出（vs@7bdaf4a，只读）
   判定：A ≡ B ≡ C（SHA-256 相等）。B 与 C 不一致而 A 与 C 一致 → 记录差异并单独报告（说明交付链输入与重跑输出的关系），不修改任何交付产物。
3. **删除层一致性核验**：确认阴性 Norman full 层删除后，Layer 1 检索结果可从存活输入（`product/norman_assets/raw_oof_full.npz` + `features/*.npz`）完整重算且与 `retrieval_result.json` bit 级一致。
4. **运行时探针**（附带产出）：Vina-GPU/CPU、Meeko、RDKit、OpenBabel、gemmi 版本；scF/ESM2 checkpoint SHA-256。落 `feasibility/runtime_probes/runtime.json`。

## 输出

`feasibility/repro_attestation/attestation.json`（机器可读，逐项 pass/fail + 哈希）+ 本目录 `RESULTS.md`（人读裁决）。失败项如实记录，不做任何自动修改。
