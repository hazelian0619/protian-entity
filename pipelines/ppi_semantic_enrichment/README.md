# ppi_semantic_enrichment pipeline（PPI 语义增强，v2）

将 `edges_ppi_v1` 的 PPI 边补充为两张语义上下文证据表：

- `data/output/evidence/ppi_method_context_v2.tsv`
  - 检测方法标准化（Y2H/Co-IP/AP-MS/PCA 等）
  - HT/LT 通量标签
  - PMID / DOI 拆分
  - STRING 的实验分与文本挖掘分拆分
- `data/output/evidence/ppi_function_context_v2.tsv`
  - GO / Reactome / KEGG 功能重叠上下文

## 输入

- `pipelines/edges_ppi/data/output/edges/edges_ppi_v1.tsv`
- `data/processed/protein_master_v6_clean.tsv`
- `data/processed/pathway_members.tsv`
- `data/output/protein/protein_kegg_pathway_v1.tsv`

## 外部数据（自动下载到 cache）

缓存目录：`pipelines/ppi_semantic_enrichment/data/raw/`

- STRING detailed links（human）
- IntAct MITAB（优先 human）
- BioGRID TAB3（优先 latest organism/human）

> 说明：若 IntAct/BioGRID 下载失败，pipeline 会继续执行（保留 STRING fallback）；但 QA 仍会严格检查验收阈值。
> 若外部源缺失，会产出 `blocked_missing_inputs` 报告（下载 URL / 版本提示 / sha256）。

## 运行

在 `1218/` 根目录：

```bash
bash pipelines/ppi_semantic_enrichment/run.sh
```

可选环境变量：

- `SAMPLE_ROWS`：样本预跑行数（默认 `5000`）
- `DATA_VERSION`：manifest 的 data_version（默认 `kg-data-local`）

## 产物与报告

- 合同：
  - `pipelines/ppi_semantic_enrichment/contracts/ppi_method_context_v2.json`
  - `pipelines/ppi_semantic_enrichment/contracts/ppi_function_context_v2.json`
- 报告：
  - `pipelines/ppi_semantic_enrichment/reports/ppi_semantic_enrichment_v2.sample.json`
  - `pipelines/ppi_semantic_enrichment/reports/ppi_semantic_enrichment_v2.build.json`
  - `pipelines/ppi_semantic_enrichment/reports/ppi_method_context_v2.validation.json`
  - `pipelines/ppi_semantic_enrichment/reports/ppi_function_context_v2.validation.json`
  - `pipelines/ppi_semantic_enrichment/reports/ppi_semantic_enrichment_v2.qa.json`
  - `pipelines/ppi_semantic_enrichment/reports/ppi_semantic_enrichment_v2.manifest.json`
  - `pipelines/ppi_semantic_enrichment/reports/ppi_semantic_enrichment_v2.result_evidence.json`（行数/覆盖率/source 分布/抽样 20 行路径）
  - `pipelines/ppi_semantic_enrichment/reports/ppi_semantic_enrichment_v2.blocked_missing_inputs.json`
  - `pipelines/ppi_semantic_enrichment/reports/ppi_method_context_v2.sample20.tsv`
  - `pipelines/ppi_semantic_enrichment/reports/ppi_function_context_v2.sample20.tsv`

## 验收门禁（QA）

- `edge_id` 映射覆盖率 >= 0.98（method / function 两张表）
- `method` 非空率 >= 0.95
- `pmid/doi` 至少一项覆盖率 >= 0.80
- 合同校验 PASS
