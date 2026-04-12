# rna_xref_enst_urs pipeline（mRNA ENST↔URS 映射修复 v2）

本 pipeline 产出独立 xref 表，不改写 `rna_master` 主表。

目标：将 mRNA 的 `ENST_*_9606` 映射到 RNAcentral `URS_*_9606`，同时输出覆盖率与冲突审计报告，并在覆盖不足时自动给出手动下载清单。

## Inputs

- `data/output/rna_master_v1.tsv`（默认）
- `data/raw/rna/rnacentral/id_mapping.tsv.gz`（默认，可改为 `.tsv`）

## Outputs

- 主输出（本地保留）  
  `data/output/rna_xref_mrna_enst_urs_v2.tsv`
- 合约  
  `pipelines/rna_xref_enst_urs/contracts/rna_xref_mrna_enst_urs_v2.json`
- 报告  
  `pipelines/rna_xref_enst_urs/reports/rna_xref_mrna_enst_urs_v2.validation.json`  
  `pipelines/rna_xref_enst_urs/reports/mrna_enst_urs_coverage_v2.json`  
  `pipelines/rna_xref_enst_urs/reports/mrna_enst_urs_conflicts_v2.json`
- 覆盖不足时自动生成  
  `pipelines/rna_xref_enst_urs/reports/mrna_enst_urs_manual_download_v2.json`

## Run

在仓库根目录执行：

```bash
bash pipelines/rna_xref_enst_urs/run.sh
```

可选环境变量：

- `SAMPLE_MASTER_ROWS`（默认 `5000`）
- `SAMPLE_IDMAP_LINES`（默认 `2000000`）
- `MIN_COVERAGE`（默认 `0.70`）
- `TARGET_COVERAGE`（默认 `0.90`）
- `STRICT_COVERAGE_GATE`（默认 `0`，设为 `1` 时低覆盖率会阻塞退出）
- `MASTER_PATH`（默认 `data/output/rna_master_v1.tsv`）
- `ID_MAPPING_PATH`（默认 `data/raw/rna/rnacentral/id_mapping.tsv.gz`）

## Match strategy

- `exact_rna_id_base`：`rna_id` 直接去 `_9606` 后匹配 `id_mapping` 中 ENST
- `fallback_ensembl_transcript_version`：回退到 `ensembl_transcript_id`（带版本）
- `fallback_ensembl_transcript_base`：回退到 `ensembl_transcript_id`（去版本）
- `fallback_rna_name_base`：回退到 `rna_name`（去版本）

## Gates

- Contract validation：字段与格式校验
- Coverage KPI：默认最低 `70%`（目标 `90%`）
- **默认非阻塞**：覆盖率低时仅告警并输出手动下载清单，pipeline 继续完成
- **严格模式可选**：`STRICT_COVERAGE_GATE=1` 时，低覆盖率会以非 0 退出

## 统一核查表（RNA）

| 组件 | 属性字段 | 数据类型 | 获取方式 | 存储形式 | 示例 | 用途/指标 | 指标 | 层次结构 | 语义关系 | 数据来源 | 维护机构 | 更新频率 | 当前状态 | 证据路径 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mRNA xref v2 | `rna_id`,`xref_id`,`match_strategy` 等 10 列 | TSV（字符串列） | 本地脚本流式 join（master × RNAcentral id_mapping） | `data/output/rna_xref_mrna_enst_urs_v2.tsv` | `ENST00000514884_9606 -> URS00002EBCD8_9606` | ENST↔URS 桥接；覆盖率与冲突审计 | 覆盖率、冲突数、策略分布 | `rna_master_mrna_v1` 子层 xref | mRNA transcript 到 RNAcentral sequence ID | RNAcentral + Ensembl 导出 | RNAcentral / Ensembl | 跟随数据快照重跑 | 已落地（v2） | `pipelines/rna_xref_enst_urs/reports/*.json` |
