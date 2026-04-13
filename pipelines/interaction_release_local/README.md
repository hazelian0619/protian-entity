# 互作本地收口发布包（interaction_release_local）

目标：基于已实现并通过 QA 的互作数据集（A/B/C/D/E），生成一份可提审的本地发布收口包：

- release summary（JSON）
- gates（JSON）
- summary contract validation（JSON）
- checklist（Markdown）
- manifest（JSON）

> 本 pipeline 聚焦“已有成果可提审”，不新增研发缺口填补。

## 输入（必须存在）

### 数据表
- PPI: `edges_ppi_v1.tsv`, `ppi_evidence_v1.tsv`, `ppi_method_context_v2.tsv`, `ppi_function_context_v2.tsv`
- PSI: `drug_target_edges_v1.tsv`, `drug_target_evidence_v1.tsv`, `psi_activity_context_v2.tsv`, `psi_structure_evidence_v2.tsv`
- RPI: `rna_protein_edges_v1.tsv`, `rna_protein_evidence_v1.tsv`, `rpi_site_context_v2.tsv`, `rpi_domain_context_v2.tsv`, `rpi_function_context_v2.tsv`
- 融合层: `interaction_cross_validation_v2.tsv`, `interaction_aggregate_score_v2.tsv`, `interaction_ontology_mapping_v2.tsv`

### QA/Gates 报告
- A: `pipelines/ppi_semantic_enrichment/reports/ppi_semantic_enrichment_v2.qa.json`
- B: `pipelines/psi_activity_structure_enrichment/reports/psi_activity_structure_enrichment_v2.qa.json`
- C: `pipelines/rpi_site_domain_enrichment/reports/rpi_site_domain_enrichment_v2.gates.json`
- D: `pipelines/interaction_cross_validation/reports/interaction_cross_validation_v2.qa.json`
- E: `pipelines/interaction_ontology_mapping/reports/interaction_ontology_mapping_v2.gates.json`

## 运行

```bash
bash pipelines/interaction_release_local/run.sh
```

可选：

```bash
DATE_TAG=2026-04-13 bash pipelines/interaction_release_local/run.sh
```

## 输出

- `pipelines/interaction_release_local/reports/interaction_release_local_<DATE>.summary.json`
- `pipelines/interaction_release_local/reports/interaction_release_local_<DATE>.gates.json`
- `pipelines/interaction_release_local/reports/interaction_release_local_<DATE>.validation.json`
- `pipelines/interaction_release_local/reports/interaction_release_local_<DATE>.manifest.json`
- `../_coord/互作补洗_本地收口Checklist_<DATE>.md`

## 收口判定（PASS 条件）

- required_tables_present = true
- all_qa_reports_present = true
- all_qa_pass = true
