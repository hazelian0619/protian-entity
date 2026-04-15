# INTERACTION L2 Local Readiness (2026-04-12)

- 生成时间(UTC): 2026-04-12T09:46:42.379734+00:00
- 扫描模式: sample
- 总体结论: **FAIL**

## 1) 三类互作同口径状态表

| 互作 | 组件 | 属性字段 | 数据类型 | 获取方式 | 存储形式 | 示例 | 用途/指标 | 指标 | 层次结构 | 语义关系 | 数据来源 | 维护机构 | 更新频率 | 证据路径 | 当前状态 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PPI | 已完成 | 已完成 | 已完成 | 缺失 | 已完成 | 已完成 | 缺失 | 缺失 | 已完成 | 已完成 | 已完成 | 已完成 | 部分 | 已完成 | 部分 |
| PSI | 已完成 | 已完成 | 已完成 | 缺失 | 已完成 | 已完成 | 缺失 | 缺失 | 已完成 | 已完成 | 已完成 | 已完成 | 部分 | 已完成 | 部分 |
| RPI | 已完成 | 已完成 | 已完成 | 已完成 | 已完成 | 已完成 | 已完成 | 已完成 | 部分 | 已完成 | 已完成 | 已完成 | 部分 | 已完成 | 部分 |

## 2) 关键指标摘要

### PPI
- 组件完成度: 2/2 (100.00%)
- 字段完成度: 2/2 (100.00%)
- evidence-edge join_rate: 1.0000
- 来源预览: string(400)
- 证据路径: {"edges": "data/output/edges/edges_ppi_v1.tsv", "evidence": "data/output/evidence/ppi_evidence_v1.tsv", "contract_edges": "pipelines/edges_ppi/contracts/edges_ppi_v1.json", "contract_evidence": "pipelines/edges_ppi/contracts/ppi_evidence_v1.json", "report_validation_edges": "pipelines/edges_ppi/reports/edges_ppi_v1.validation.json", "report_validation_evidence": "pipelines/edges_ppi/reports/ppi_evidence_v1.validation.json"}

### PSI
- 组件完成度: 3/3 (100.00%)
- 字段完成度: 2/2 (100.00%)
- evidence-edge join_rate: 1.0000
- 来源预览: drugbank(400)
- 证据路径: {"drug_target_edges": "data/output/edges/drug_target_edges_v1.tsv", "drug_target_evidence": "data/output/evidence/drug_target_evidence_v1.tsv", "molecules_m3_psi": "data/output/molecules/molecules_m3_psi_edges_v1.tsv", "contract_edges": "pipelines/drugbank/contracts/drug_target_edges_v1.json", "contract_evidence": "pipelines/drugbank/contracts/drug_target_evidence_v1.json", "report_validation_edges": "pipelines/drugbank/reports/drug_target_edges_v1.validation.json", "report_validation_evidence": "pipelines/drugbank/reports/drug_target_evidence_v1.validation.json"}

### RPI
- 组件完成度: 2/2 (100.00%)
- 字段完成度: 2/2 (100.00%)
- evidence-edge join_rate: 0.0100
- 来源预览: starbase(400)
- 证据路径: {"edges": "data/output/edges/rna_protein_edges_v1.tsv", "evidence": "data/output/evidence/rna_protein_evidence_v1.tsv", "contract_edges": "pipelines/rna_rpi/contracts/rna_protein_edges_v1.json", "contract_evidence": "pipelines/rna_rpi/contracts/rna_protein_evidence_v1.json", "report_metrics": "pipelines/rna_rpi/reports/rna_rpi_v1.metrics.json", "report_gates": "pipelines/rna_rpi/reports/rna_rpi_v1.gates.json", "report_manifest": "pipelines/rna_rpi/reports/rna_rpi_v1.manifest.json", "report_validation_edges": "pipelines/rna_rpi/reports/rna_protein_edges_v1.validation.json", "report_validation_evidence": "pipelines/rna_rpi/reports/rna_protein_evidence_v1.validation.json"}

## 3) Top10 缺口（阻塞级/可延期 + 工作量）

| # | 互作 | 缺口 | 分级 | 工作量 | 建议 |
|---|---|---|---|---|---|
| 1 | RPI | 层次结构为部分 | 阻塞级 | M | 按统一 L2 字段规范补齐并在本地复跑 QA。 |
| 2 | PPI | 更新频率为部分 | 可延期 | S | 按统一 L2 字段规范补齐并在本地复跑 QA。 |
| 3 | PSI | 更新频率为部分 | 可延期 | S | 按统一 L2 字段规范补齐并在本地复跑 QA。 |
| 4 | RPI | 更新频率为部分 | 可延期 | S | 按统一 L2 字段规范补齐并在本地复跑 QA。 |

## 4) 本地通过后再上传：前置门槛

| Gate | 说明 | 结果 |
|---|---|---|
| ppi_required_artifacts_ready | PPI 必需产物（edges/evidence/核心衍生）齐全 | PASS |
| psi_required_artifacts_ready | PSI 必需产物（edges/evidence/核心衍生）齐全 | PASS |
| rpi_required_artifacts_ready | RPI 必需产物（edges/evidence/核心衍生）齐全 | PASS |
| rpi_edge_evidence_hierarchy_ge_0_99 | RPI evidence.edge_id 与 edges.edge_id 连接率 >= 0.99 | FAIL |
| interaction_schema_uniformity | 三类互作均满足统一 L2 字段口径 | PASS |
| zero_blocking_gaps | Top gaps 中阻塞级缺口为 0 | FAIL |

## 5) 收口建议

1. 先补齐 PPI/PSI 的本地 edges+evidence+contract+validation 产物，再做统一上传。
2. 保持三类互作统一 L2 字段（method/reference/score/source_version/fetch_date）。
3. 所有缺口关闭后，重新运行 interaction_readiness 全量核查并锁定 manifest。

