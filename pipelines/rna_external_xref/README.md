# RNA External XREF v1

目标：生成 `data/output/rna_external_xref_v1.tsv`，提供 RNA 主表到外部库（RNAcentral / Ensembl / RefSeq / miRBase）的统一映射，并输出 coverage、conflicts、validation、manifest。

## 输入

- `data/output/rna_master_v1.tsv`
- `data/raw/rna/rnacentral/id_mapping.tsv.gz`
- `data/raw/rna/aux_xref/*`（至少一个包含 ENST + RefSeq 列的 TSV/TXT）

> 若缺输入，`run.sh` 会中断并生成：
> `pipelines/rna_external_xref/reports/rna_external_xref_v1.manual_download.json`

## 输出

- `data/output/rna_external_xref_v1.tsv`
- `pipelines/rna_external_xref/reports/rna_external_xref_v1.build.json`
- `pipelines/rna_external_xref/reports/rna_external_xref_v1.coverage.json`
- `pipelines/rna_external_xref/reports/rna_external_xref_v1.conflicts.json`
- `pipelines/rna_external_xref/reports/rna_external_xref_v1.validation.json`
- `pipelines/rna_external_xref/reports/rna_external_xref_v1.manifest.json`

## 表结构

`rna_external_xref_v1.tsv`

- `rna_id`
- `xref_db`
- `xref_id`
- `xref_type`
- `source`
- `source_version`
- `fetch_date`

## 运行

```bash
bash pipelines/rna_external_xref/run.sh
```

可选环境变量：

- `INPUT_MASTER`
- `INPUT_ID_MAPPING`
- `INPUT_AUX_DIR`
- `SAMPLE_MASTER_ROWS`（默认 `20000`）
- `SAMPLE_IDMAP_LINES`（默认 `3000000`）
- `MIN_BACKLINK_RATE`（默认 `0.99`）
- `SOURCE_VERSION`（默认 `RNAcentral:25`）
- `DATA_VERSION`（默认 `kg-rna-external-xref-v1`）

## 质量门

- `rna_id` 回链主表比例：`coverage.summary.backlink_rate >= 0.99`
- 合同校验通过：`validation.passed == true`
- 每个 `xref_db` 输出独立覆盖率和冲突率：
  - `coverage.per_xref_db`
  - `conflicts.per_xref_db`
