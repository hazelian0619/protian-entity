# RNA 表达证据 + RBP 锚点管线（rna_expression_rbp_v1）

先做最小可用证据层，输出两张表：

- `data/output/rna_expression_evidence_v1.tsv`
- `data/output/rna_rbp_sites_v1.tsv`

## 输入

- `data/output/rna_master_v1.tsv`
- `data/raw/rna/encode/expression_evidence.tsv(.gz)`（由预处理脚本生成）
- `data/raw/rna/encode/rbp_sites.tsv(.gz)`（由预处理脚本生成）
- 预处理依赖：
  - `.tmp/ENCFF937GQE.tsv`（ENCODE RNA-seq gene quant）
  - `.tmp/ENCFF593RED.bed.gz`（ENCODE eCLIP peaks）
  - `data/raw/rna/ensembl/Homo_sapiens.GRCh38.115.chr.gtf.gz`（转录本坐标）

> 输入支持 `TSV/CSV`，脚本自动识别分隔符；支持 `.gz`。

## 输出字段

### 1) 表达证据：`rna_expression_evidence_v1.tsv`

`rna_id,biosample_tissue,biosample_cell_type,expression_value,expression_unit,assay,sample_id,experiment_id,source_dataset,source_file,evidence_level,source,source_version,fetch_date`

### 2) RBP 锚点：`rna_rbp_sites_v1.tsv`

`rna_id,rbp_symbol,biosample_cell_type,assay,peak_count,score,evalue,sample_id,experiment_id,source_dataset,source_file,evidence_level,source,source_version,fetch_date`

## 目录结构

- `pipelines/rna_expression_rbp/run.sh`
- `pipelines/rna_expression_rbp/scripts/prepare_encode_minimal_inputs.py`
- `pipelines/rna_expression_rbp/scripts/build_rna_expression_rbp.py`
- `pipelines/rna_expression_rbp/contracts/rna_expression_evidence_v1.json`
- `pipelines/rna_expression_rbp/contracts/rna_rbp_sites_v1.json`
- `pipelines/rna_expression_rbp/reports/*.json`

## 运行

```bash
bash pipelines/rna_expression_rbp/run.sh
```

运行流程：
1. 预处理 ENCODE 输入（生成 `expression_evidence.tsv.gz` / `rbp_sites.tsv.gz`）
2. preflight 检查（缺输入即中断并输出下载清单）
3. 最小样本 (`--limit 200`)
4. 全量构建
5. 两张表各自 contract 校验
6. 产出 manifest

## 报告

- `rna_expression_rbp_v1.blocked_or_ready.json`
- `rna_expression_rbp_v1.prep.metrics.json`
- `rna_expression_rbp_v1.sample.metrics.json`
- `rna_expression_rbp_v1.metrics.json`
- `rna_expression_evidence_v1.validation.json`
- `rna_rbp_sites_v1.validation.json`
- `rna_expression_rbp_v1.manifest.json`

## 手动下载中断条件

若 ENCODE 数据缺失或接口不稳定，立即停止并按下列目录放置离线文件：

- `data/raw/rna/encode/expression_evidence.tsv.gz`
- `data/raw/rna/encode/rbp_sites.tsv.gz`
- `data/raw/rna/encode/metadata.tsv`

校验：

```bash
sha256sum data/raw/rna/encode/expression_evidence.tsv.gz
sha256sum data/raw/rna/encode/rbp_sites.tsv.gz
sha256sum data/raw/rna/encode/metadata.tsv
```

## 统一核查模板（本管线）

| 组件 | 属性字段 | 数据类型 | 获取方式 | 存储形式 | 示例 | 用途/指标 | 指标 | 层次结构 | 语义关系 | 数据来源 | 维护机构 | 更新频率 | 当前状态 | 证据路径 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 证据层 | `expression_value`,`rbp_symbol`,`evidence_level` | `Numeric/String/Enum` | ENCODE 表达量与 eCLIP peak 汇总（离线文件） | `rna_expression_evidence_v1.tsv` + `rna_rbp_sites_v1.tsv` | `12.3`, `ELAVL1`, `high` | 回链 RNA 主表并补充证据 | `rna_id`回链率、元数据齐全率、evidence_level 非空率 | L2 证据层 | RNA `supported_by` 表达证据；RNA `bound_by` RBP | ENCODE | ENCODE Consortium | 滚动更新 | 缺输入时阻塞（按约束中断） | `pipelines/rna_expression_rbp/reports/rna_expression_rbp_v1.blocked_or_ready.json` |
