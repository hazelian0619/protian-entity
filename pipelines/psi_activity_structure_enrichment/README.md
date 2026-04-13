# PSI 活性 + 结构证据增强（v2, P0）

本 pipeline 针对 PSI（protein-small molecule interaction）补齐两类增强表：

- `data/output/evidence/psi_activity_context_v2.tsv`
- `data/output/evidence/psi_structure_evidence_v2.tsv`

目标对齐任务 B 验收口径：

1. 活性关系符规范化：`=, <, >, <=, >=, ~`
2. assay 类型与语境：`B/F/A/T/P` + 描述语境
3. 条件对象抽取：`pH / 温度 / 体系关键词`
4. 结构证据：PDB 复合体语境 + 结构关联亲和力分数

---

## 输入

必需：

- `data/output/molecules/chembl_m3.sqlite`
  - 表：`psi_edges_v1`, `psi_evidence_v1`
- `data/output/protein/pdb_structures_v1.tsv`

说明：

- 本版优先使用本地可复现输入（ChEMBL PSI + RCSB PDB 快照）。
- BindingDB / PDBbind / BioLiP 可在后续轮次叠加为更细粒度“化合物-位点级”结构证据；当前版本先完成 P0 验收所需覆盖率与语义字段。

---

## 输出 schema（核心字段）

### 1) `psi_activity_context_v2.tsv`

核心字段：

- 锚点：`edge_id`
- 活性：`standard_type, standard_relation, standard_value, standard_units, standard_value_nM, pchembl_value_eff`
- assay：`assay_type, assay_type_desc, assay_description, assay_context, assay_confidence_score, bao_format`
- 条件：`condition_pH, condition_temperature_c, condition_system, condition_context`
- 追溯：`source, source_version, fetch_date`

### 2) `psi_structure_evidence_v2.tsv`

核心字段：

- 锚点：`edge_id`
- 结构证据：`pdb_id, pdb_experimental_method, pdb_resolution, pdb_release_date, pdb_ligand_count`
- 语义：`structure_evidence_type=pdb_target_complex_context, complex_match_level=target_level`
- 结构关联亲和力：`standard_type, standard_relation, standard_value_nM, pchembl_value_eff, structure_affinity_score`
- 追溯：`source, source_version, fetch_date`

---

## 执行

在仓库根目录（`1218/`）运行：

```bash
bash pipelines/psi_activity_structure_enrichment/run.sh
```

默认流程：

1. smoke（`SMOKE_ROWS=50000`）
2. full（全量）
3. contract validation
4. QA gates（验收指标）
5. manifest

可选环境变量：

- `FETCH_DATE`（默认 UTC 当天）
- `SMOKE_ROWS`（默认 50000）
- `MIN_LIGAND_COUNT`（默认 1）
- `DATA_VERSION`

---

## 报告

- `pipelines/psi_activity_structure_enrichment/reports/psi_activity_structure_enrichment_v2.build.json`
- `pipelines/psi_activity_structure_enrichment/reports/psi_activity_context_v2.validation.json`
- `pipelines/psi_activity_structure_enrichment/reports/psi_structure_evidence_v2.validation.json`
- `pipelines/psi_activity_structure_enrichment/reports/psi_activity_structure_enrichment_v2.qa.json`
- `pipelines/psi_activity_structure_enrichment/reports/psi_activity_structure_enrichment_v2.manifest.json`

---

## 验收阈值（QA gates）

- `edge_id_join_rate_activity >= 0.99`
- `activity_type_coverage >= 0.95`
- `assay_field_coverage >= 0.90`
- `structure_edge_coverage >= 0.20`
- 两张表 contract 全 PASS
