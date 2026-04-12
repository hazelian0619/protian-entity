# rna_pdb：RNA 3D 结构映射补齐（v1）

把 RNA 与 PDB 结构建立可复跑映射，产出主线表：

- `data/output/rna_pdb_structures_v1.tsv`

本管线默认执行：**sample → full → contract validation → QA → manifest**。

---

## 输入

必需：

- `data/output/rna_master_v1.tsv`
- `data/raw/rna/rnacentral/id_mapping.tsv(.gz)`

可选增强：

- `data/output/rna_xref_mrna_enst_urs_v2.tsv`
  - 若存在，优先用于 ENST→URS 映射；否则回退到 `id_mapping.tsv` 的 ENSEMBL* 记录投票映射。

RCSB API（在线可用时）：

- GraphQL：`https://data.rcsb.org/graphql`
- fallback：`https://data.rcsb.org/rest/v1/core/entry/{pdb_id}`

---

## 输出

- 主表：`data/output/rna_pdb_structures_v1.tsv`
- 合约：`pipelines/rna_pdb/contracts/rna_pdb_structures_v1.json`
- 报告：
  - `pipelines/rna_pdb/reports/rna_pdb_structures_v1.build.json`
  - `pipelines/rna_pdb/reports/rna_pdb_structures_v1.api_audit.json`
  - `pipelines/rna_pdb/reports/rna_pdb_structures_v1.validation.json`
  - `pipelines/rna_pdb/reports/rna_pdb_structures_v1.qa.json`
  - `pipelines/rna_pdb/reports/rna_pdb_structures_v1.manifest.json`
- sample 报告：
  - `pipelines/rna_pdb/reports/rna_pdb_structures_v1.sample.*.json`

---

## 主表字段

- `rna_id`
- `rna_type`
- `urs_id`
- `pdb_id`
- `pdb_entity_id`
- `mapping_strategy`
- `experimental_method`
- `resolution`
- `release_date`
- `source`
- `source_version`
- `fetch_date`

---

## 运行

在仓库根目录执行：

```bash
bash pipelines/rna_pdb/run.sh
```

可选环境变量：

- `INPUT_MASTER`（默认 `data/output/rna_master_v1.tsv`）
- `INPUT_ID_MAPPING`（默认 `data/raw/rna/rnacentral/id_mapping.tsv.gz`）
- `INPUT_XREF`（默认 `data/output/rna_xref_mrna_enst_urs_v2.tsv`）
- `SOURCE_VERSION`（默认 `RNAcentral:25;RCSB:live`）
- `DATA_VERSION`（默认 `kg-data-local`）

sample 参数：

- `SAMPLE_MASTER_ROWS`（默认 `0`，即不截断；>0 时仅扫描前 N 行 master）
- `SAMPLE_IDMAP_LINES`（默认 `0`，即不截断；>0 时仅扫描前 N 行 id_mapping）
- `SMOKE_MAX_UNIQUE_PDB`（默认 `200`）

RCSB API 参数：

- `RCSB_API_FAIL_THRESHOLD`（默认 `0.05`）
- `RCSB_BATCH_SIZE`（默认 `500`）
- `RCSB_TIMEOUT`（默认 `40`）
- `RCSB_RETRIES`（默认 `3`）
- `RCSB_SLEEP_BETWEEN_BATCHES`（默认 `0.05`）

---

## QA 验收口径

`02_qa_rna_pdb_structures_v1.py` 主要检查：

- `rna_id` 回链到 `rna_master_v1.tsv`
- `pdb_id` 格式：`^[0-9A-Z]{4}$`
- `release_date`：日期格式与范围
- `resolution`：数值格式与范围（EM 与非 EM 分档）
- `api_failure_rate <= fail_threshold`

---

## API 失败阈值与手动下载

若 `api_error_unique_pdb / selected_unique_pdb > RCSB_API_FAIL_THRESHOLD`：

1. 构建中断（退出码 `3`）
2. 报告写入手动下载清单（不假成功）：
   - `pipelines/rna_pdb/reports/rna_pdb_structures_v1.api_audit.json`
   - `pipelines/rna_pdb/reports/rna_pdb_structures_v1.build.json`

清单建议放置目录：`data/raw/rna/pdb_bulk/`。
