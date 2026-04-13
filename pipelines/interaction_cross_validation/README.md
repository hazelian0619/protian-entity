# interaction_cross_validation（v2）

跨库一致性与统一评分层（任务 D）：

- `data/output/evidence/interaction_cross_validation_v2.tsv`
- `data/output/evidence/interaction_aggregate_score_v2.tsv`

目标：对 **PPI / PSI / RPI** 三类互作统一计算：

1. `consistent_across_n`（多来源一致性计数）
2. `conflict_flag` + `conflict_reason`（跨库冲突）
3. `aggregate_score`（0-1 聚合评分）

---

## 输入

### 必需（核心）

- PPI：`edges_ppi_v1.tsv` + `ppi_evidence_v1.tsv`
- PSI：`drug_target_edges_v1.tsv` + `drug_target_evidence_v1.tsv`
- **PSI 任务 B 更新表**（强依赖）：
  - `psi_activity_context_v2.tsv`
  - `psi_structure_evidence_v2.tsv`
- RPI：`rna_protein_edges_v1.tsv` + `rna_protein_evidence_v1.tsv`
  - 若本仓库无 RPI v1，`run.sh` 会自动尝试 `../protian-entity/...`
  - 再无则回退到 `gene_to_rna_*_v1`（仅兜底）
- 小分子 xref：`molecule_xref_core_v2.tsv`（无则回退 v1）
  - 用于 DrugBank ID ↔ InChIKey 对齐，支撑 PSI 跨库一致性（DrugBank 与 ChEMBL/PDB）

### 可选（上下文）

- PPI：`ppi_method_context_v2.tsv`, `ppi_function_context_v2.tsv`
- RPI：`rpi_site_context_v2.tsv`, `rpi_domain_context_v2.tsv`, `rpi_function_context_v2.tsv`

---

## 输出字段

### 1) cross-validation

`record_id, interaction_type, edge_id, entity_pair_key, consistent_across_n, source_list, evidence_count, distinct_methods_n, reference_count, predicate_set, direction_set, conflict_flag, conflict_reason, source_version_list, fetch_date`

### 2) aggregate score

`aggregate_id, interaction_type, edge_id, consistent_across_n, evidence_count, distinct_methods_n, numeric_score_max, numeric_score_norm, reference_coverage, context_coverage, conflict_flag, aggregate_score, score_bucket, fetch_date`

---

## 运行

```bash
cd /Users/pluviophile/graph/1218
bash pipelines/interaction_cross_validation/run.sh
```

流程：

1. preflight（阻塞检查 + blocked_or_ready 报告）
2. smoke（`SMOKE_LIMIT_PER_TYPE=5000`）
3. full
4. 合同校验（两张表）
5. QA gates
6. manifest

可选参数：

```bash
SMOKE_LIMIT_PER_TYPE=8000 DATA_VERSION=kg-data-local FETCH_DATE=2026-04-12 \
  bash pipelines/interaction_cross_validation/run.sh
```

---

## 验收口径（任务 D）

- 三类互作都能计算 cross-validation（PPI/PSI/RPI）
- 聚合评分分布合理（非全 0 / 全 1，且非常数）
- 合同 + QA gates 全 PASS

---

## 缺外部数据时处理

本 pipeline 会先跑 preflight；若必需输入缺失，输出：

- `pipelines/interaction_cross_validation/reports/interaction_cross_validation_v2.blocked_or_ready.json`

并给出下载/放置/sha256 校验提示，避免静默失败。
