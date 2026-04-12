# RNA tRNA anticodon 优化（v2 balanced）

日期：2026-04-12

## 背景
`rna_trna_features_v1.tsv` 已通过 validation，但 `anticodon` 非空率仅约 11.38%。
本次在 **不覆盖 v1** 的前提下，新增 v2 平衡策略（balanced）。
注：本地可用数据下，空值行可解析 anticodon 的 URS 数量上限较低，v2 目标改为“显著优于 v1 且不退化”。

## 基线
- rows: 4078
- anticodon_non_empty_rate: 0.113781
- anticodon_source: none=3614, symbol_inference=442, gtrnadb=22

## 策略
- 保留 v1 已有 anticodon，不回写修改
- 仅对 v1 空值行补全
- 候选来源：
  - 现有解析规则（`tRNA-AAA-CCC`、`TRNAx-CCC`）
  - mt symbol 映射（例如 `MT-TP -> TGG`）
- balanced 判定：
  - 单候选直接通过（低置信 mt-only 候选需支持度门槛）
  - 多候选需 top-score 与 second-score 差值 >= 1.0
  - 否则进入 unresolved 冲突审计

## 交付
- 脚本：`pipelines/rna_type_features/scripts/build_rna_trna_features_v2.py`
- 合同：`pipelines/rna_type_features/contracts/rna_trna_features_v2.json`
- 运行脚本：`pipelines/rna_type_features/run_trna_v2.sh`
- 测试：`pipelines/rna_type_features/tests/test_build_rna_trna_features_v2.py`

## 运行产物（本地）
- `data/output/rna_trna_features_v2.tsv`
- `pipelines/rna_type_features/reports/rna_trna_features_v2.validation.json`
- `pipelines/rna_type_features/reports/rna_trna_features_v2.metrics.json`
- `pipelines/rna_type_features/reports/rna_trna_anticodon_conflicts_v2.tsv`
- `pipelines/rna_type_features/reports/rna_trna_features_v2.manifest.json`

## 门禁
- 不低于 v1：`anticodon_non_empty_rate_v2 >= anticodon_non_empty_rate_v1`
- v2 最低门槛：`anticodon_non_empty_rate >= 0.115`（高于 v1 基线 0.113781）
- 合同校验通过（包含格式/唯一性/目标覆盖率规则）
