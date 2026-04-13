# 互作本体标准化层管线（interaction_ontology_mapping_v2）

目标：输出统一本体映射表：

- PSI-MI 方法映射
- GO 术语映射
- 统一谓词标准化（`detected_by / supported_by / participates_in / standardized_as`）

产物：

- `data/output/evidence/interaction_ontology_mapping_v2.tsv`

## 输入

必需：

- `data/output/edges/edges_ppi_v1.tsv`
- `data/output/edges/drug_target_edges_v1.tsv`
- `data/output/edges/rna_protein_edges_v1.tsv`

可选：

- `data/output/evidence/rpi_function_context_v2.tsv`（用于 GO 术语映射）

## 输出字段

`mapping_id,interaction_type,edge_id,mapping_scope,raw_term,standardized_predicate,ontology_namespace,ontology_id,ontology_label,ontology_uri,mapping_confidence,source,source_version,fetch_date`

## 运行

```bash
bash pipelines/interaction_ontology_mapping/run.sh
```

流程：
1. preflight 检查
2. 样本运行（`--limit-per-type 2000`）
3. 全量运行
4. contract 校验
5. manifest 生成
6. gates 汇总

## 报告

- `interaction_ontology_mapping_v2.blocked_or_ready.json`
- `interaction_ontology_mapping_v2.sample.metrics.json`
- `interaction_ontology_mapping_v2.metrics.json`
- `interaction_ontology_mapping_v2.validation.json`
- `interaction_ontology_mapping_v2.manifest.json`
- `interaction_ontology_mapping_v2.gates.json`

## 验收门槛（任务 E）

- 可映射记录占比 >= 0.85
- 本体 URI 规范率 >= 0.99
- 合同 + gates PASS
