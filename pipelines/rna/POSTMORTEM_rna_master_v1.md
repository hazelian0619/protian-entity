# RNA L1（v1）复盘 + 决策说明（工业化闭环 v1）

> 目标：把 RNA（miRNA + mRNA/transcript）做成可发布、可复跑、可回归的 L1 实体表交付。
> 核心约束：产物体积较大，因此数据通过 Release 分发；仓库内保留 **pipeline + contract + QA 报告 + manifest**。

## 1) 背景与主线对齐（为什么要做 RNA L1）

在 KG 执行流中（见 `docs/KG_EXECUTION_FLOW.md`），RNA 与 Protein 同属于 L1 实体层：

- 为后续 edge（例如 RPI、PSI）提供节点与属性
- 为 cross-entity 对齐（gene_master）提供基因映射字段

RNA 的难点不是“脚本能不能跑”，而是：

- 数据源多且大（xref/FASTA），本地环境差异会导致复现成本高
- 产物大，直接 commit 会造成仓库膨胀与 GitHub 限制问题

因此 v1 的主线是：

- 把复现路径写清楚（能跑/能定位）
- 把交付形态标准化（contract/QA/manifest）
- 把大文件分发从 git 中剥离（Release）

## 2) 输入与产物（做什么）

### 输入（本地重跑所需）

详见 `docs/rna/RNA_SOURCES_AND_VERSIONS.md`，v1 主要依赖：

- RNAcentral Release 25（`id_mapping.tsv.gz` 等）
- Ensembl cDNA FASTA（`Homo_sapiens.GRCh38.cdna.all.fa.gz`）
- Protein 主表（作为种子 gene 映射来源）：`data/processed/protein_master_v6_clean.tsv`

约定路径：`data/raw/rna/`。

### 输出（本地产物）

- `data/output/rna_master_v1.tsv`
- `data/output/rna_master_mirna_v1.tsv`
- `data/output/rna_master_mrna_v1.tsv`

### 仓库内保留（小文件，可审计）

- pipeline：`pipelines/rna/run.sh`
- contract：`pipelines/rna/contracts/rna_master_v1.json`
- QA 报告：`pipelines/rna/reports/rna_master_v1.validation.json`
- manifest：`pipelines/rna/reports/rna_v1.manifest.json`

## 3) 范围界定（v1 不做什么）

- v1 不追求导入 Rfam 的家族注释（文档中保留升级触发条件）
- v1 先聚焦 Human 9606 的 miRNA + mRNA/transcript 口径

这样做的原因是“抓大放小”：先把 L1 主表交付闭环跑通，并保证字段稳定。

## 4) 复现闭环（如何重跑 & 如何不重跑）

### 4.1 不重跑（推荐给大多数同事）

直接下载 Release（包含 `.tsv.gz` + QA + manifest）：

- https://github.com/hazelian0619/protian-entity/releases/tag/rna-l1-v1

适用场景：只需要用数据做 join/分析/建图，不需要改 ETL。

### 4.2 本地重跑（开发者）

前置条件：把 v1 所需 raw 文件放到 `data/raw/rna/`（见 `docs/rna/RNA_SOURCES_AND_VERSIONS.md`）。

运行：

```bash
cd /Users/pluviophile/graph/1218
bash pipelines/rna/run.sh
```

重跑完成后可执行 QA + manifest（与 release 同范式）：

```bash
python3 tools/kg_validate_table.py \
  --contract pipelines/rna/contracts/rna_master_v1.json \
  --table data/output/rna_master_v1.tsv \
  --out pipelines/rna/reports/rna_master_v1.validation.json

python3 tools/kg_make_manifest.py \
  --data-version kg-data-local \
  --out pipelines/rna/reports/rna_v1.manifest.json \
  data/output/rna_master_v1.tsv \
  data/output/rna_master_mirna_v1.tsv \
  data/output/rna_master_mrna_v1.tsv
```

## 5) 常见卡点 / 隐藏点（最省时间的部分）

1) **缺 raw 文件**：这是最常见原因。
   - 现象：脚本报找不到 `id_mapping.tsv.gz` 或 Ensembl FASTA。
   - 解决：按 `docs/rna/RNA_SOURCES_AND_VERSIONS.md` 放到 `data/raw/rna/`。

2) **产物太大不该进 git**：
   - `data/output/*.tsv` 属于可重建产物，默认不提交。
   - 发布给同事使用走 Release（带 manifest/QA）。

3) **版本口径漂移**：
   - v1 锁定 RNAcentral Release 25（见 `docs/rna/RNA_SOURCES_AND_VERSIONS.md`）。
   - 如果 raw 版本不一致，可能造成行数/覆盖率变化，进而影响 QA。

4) **下游对齐问题**：
   - 先用 manifest 确认“是不是同一份文件”。
   - 再看 validation 的 coverage/格式规则是否退化。

## 6) 已知边界与 v2 方向

- v2 触发条件已在 `docs/rna/RNA_SOURCES_AND_VERSIONS.md` 明确（RNAcentral 升级、Rfam 真正导入等）。
- v2 建议引入 CHANGELOG，固化版本间差异与回归解释。

