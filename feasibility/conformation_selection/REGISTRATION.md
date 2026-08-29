# 构象层复现预注册 — EXP3 三臂受控重对接（re 版）

Registered 2026-08-29。上游协议：`yeastbridge_vs@7bdaf4a` `experiments/competition_v1/EXP3_three_arm_screening.md`（只读引用）。本注册声明 re 侧复现的协议、替换与判定规则，执行前冻结。

## 目标

复现上游 EXP3 的核心裁决：**方向匹配的构象选择在受控比对中优于默认单构象**（上游结果：参照面板构象偏好 arms B/C 5/5 vs arm A 1/5，sign test p=0.031；声明边界为「构象偏好层面」）。该裁决是 re 管线构象状态选择步骤（抑制方向对接非活性构象、激活方向对接活性构象）的直接证据依据。

## 库与受体（冻结来源）

- 配体库：vs 冻结库 + 1:10 性质匹配 decoys（MW/cLogP 分箱匹配、对 actives 最大 Tanimoto < 0.4），来源 `yeastbridge_vs/reports/competition_v1/screen_v1/{library,ligands}/`，逐文件 SHA-256 登记。
- 受体：vs 认证构象（LPAR1: 4Z34/7YU3/7TD0；KCNQ2: 7CR0/7CR1/7CR2，pdbqt 来自 vs 制备产物），只读引用 + SHA-256 登记。不足处在本目录下补制备。

## 臂定义（与上游一致）

- **Arm A（行业默认）**：单一默认构象（LPAR1→4Z34；KCNQ2→7CR0 apo）。
- **Arm B（方向匹配）**：方向匹配构象集成，delta-score = score(正确构象) − score(错误构象)。
- **Arm C（朴素集成）**：全部认证构象对接取最优，无方向逻辑。

## 真值（与上游一致）

文献确立的配体药理（上游 reference_panel_evidence.md）：LPA agonist→active；ONO 系列 antagonist→inactive；retigabine/ZTZ240/flupirtine opener→open；linopirdine/XE991 blocker→apo/closed。

## 引擎替换声明

上游：Vina 1.2.7 CPU，exhaustiveness 32，固定种子与 CPU 预算。re 版：**Vina-GPU 2.1**（同一 Vina 打分函数的 GPU 实现，`configs/chembl_branch.json` 登记的二进制路径），固定种子、空卡执行。分数与 CPU 版不保证 bit 级一致；判定只看偏好方向性，不看分数绝对值。

## 指标与 Gate（预声明）

- **主判定（对应 re 实际使用的声明）**：参照面板配体的构象偏好——方向匹配构象得分为优的配体比例，binomial/sign test，p < 0.05 判 PASS；边界措辞与上游相同（conformation-preference level only）。
- 次级报告：EF1%、BEDROC、ROC-AUC（正确方向富集）、delta-score AUC、配对 bootstrap CI。次级指标不通过只如实报告，不推翻主判定。

## 输出边界

全部产物只写 `feasibility/conformation_selection/{results,logs}/`；`product/` 与 vs 目录零写入。
