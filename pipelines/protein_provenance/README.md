# Protein provenance / source_version sidecar（v1）

本 pipeline 只做“版本可追溯层”，**不改动** `data/processed/protein_master_v6_clean.tsv` 的语义与列结构。

## 输入

- `data/processed/protein_master_v6_clean.tsv`
- `data/processed/protein_edges.tsv`
- `data/processed/ptm_sites.tsv`
- `data/processed/pathway_members.tsv`
- `data/processed/alphafold_quality.tsv`
- `data/processed/hgnc_core.tsv`

## 输出

- `data/output/protein/protein_source_versions_v1.tsv`
  - 字段：`dataset,primary_source,source_version,fetch_date,evidence_field,notes`
- `pipelines/protein_provenance/reports/protein_source_versions_v1.validation.json`
- `pipelines/protein_provenance/reports/protein_source_versions_v1.qa.json`
- `pipelines/protein_provenance/reports/protein_source_versions_v1.manifest.json`

## 运行

```bash
cd /Users/pluviophile/graph/1218
bash pipelines/protein_provenance/run.sh
```

### 运行策略（符合协作约束）

`run.sh` 会先执行 smoke（默认每个输入扫描前 2000 行），再执行 full：

1. smoke build + validation + QA
2. full build + validation + QA + manifest

可调参数：

```bash
SMOKE_ROWS=5000 DATA_VERSION=kg-data-local bash pipelines/protein_provenance/run.sh
```

## 关键口径

### 为什么主表暂不强行新增 `source_version` 列？

1. **避免破坏现有 contract 与下游 join 假设**：`protein_master_v6_clean.tsv` 已被 gene/spine/ppi/drugbank 等 pipeline 直接消费。
2. **主表是多源融合产物**：单列 `source_version` 难表达 UniProt / AlphaFold / HGNC / STRING 等异构版本。
3. **sidecar 更稳妥**：将版本锚点独立到 `protein_source_versions_v1.tsv`，可追溯且不改变主表语义。

## 验收映射

- 版本字段非空率 100%：由 contract `source_version_non_empty` + QA gate `source_version_non_empty_100pct` 保证。
- 每个数据集恰好一行版本锚点：由 `dataset_unique` + 各 `has_*_row` 规则 + QA gate `each_dataset_exactly_one_row` 保证。
- 主表多版本规则明确：`protein_master_v6_clean.tsv` 行采用 composite anchor，并在 `notes` 与政策文档中说明。
