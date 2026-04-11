# protein_kegg pipeline（KEGG pathway 映射，和 Reactome 并存）

本 pipeline 将 `protein_master_v6_clean.tsv` 中的 `ncbi_gene_id` 映射到 KEGG 人类基因（`hsa:*`），再展开到 KEGG pathway（`hsa\d{5}`），产出独立表：

- `data/output/protein/protein_kegg_pathway_v1.tsv`

该表不会覆盖或改写现有 `data/processed/pathway_members.tsv`（Reactome），只做并存补充。

## Inputs

- `data/processed/protein_master_v6_clean.tsv`
  - required: `uniprot_id`, `ncbi_gene_id`
- `data/processed/pathway_members.tsv`（用于 Reactome 互补性统计）

## KEGG resources（REST）

- `https://rest.kegg.jp/info/pathway`
- `https://rest.kegg.jp/conv/hsa/ncbi-geneid`
- `https://rest.kegg.jp/link/pathway/hsa`
- `https://rest.kegg.jp/list/pathway/hsa`

缓存位置：`pipelines/protein_kegg/.cache/`

## Outputs

- **主输出（大文件，本地保留）**
  - `data/output/protein/protein_kegg_pathway_v1.tsv`
- **小文件（QA / 审计）**
  - `pipelines/protein_kegg/reports/protein_kegg_pathway_v1.sample.json`
  - `pipelines/protein_kegg/reports/protein_kegg_pathway_v1.build.json`
  - `pipelines/protein_kegg/reports/protein_kegg_pathway_v1.validation.json`
  - `pipelines/protein_kegg/reports/protein_kegg_pathway_v1.qa.json`
  - `pipelines/protein_kegg/reports/protein_kegg_pathway_v1.manifest.json`

## How to run

在仓库根目录（`1218/`）执行：

```bash
bash pipelines/protein_kegg/run.sh
```

可选环境变量：

- `SAMPLE_ROWS`：最小样本行数（默认 `200`）
- `DATA_VERSION`：写入 manifest（默认 `kg-data-local`）

## Validation / QA

- Contract: `pipelines/protein_kegg/contracts/protein_kegg_pathway_v1.json`
- 关键 gate：
  - join 到主表（`uniprot_id` 外键）
  - `kegg_pathway_id` 格式校验：`^hsa\d{5}$`
  - 与 Reactome 的互补统计（overlap / kegg_only / reactome_only / jaccard）

## 手动下载中断条件（按要求）

若 KEGG API 无法稳定访问（重试后仍失败），脚本会直接报错并中断，错误信息会包含：

1. 失败 URL
2. 建议的本地放置路径（`pipelines/protein_kegg/.cache/...`）
3. 校验命令（`sha256sum <file>`）

即：先输出“手动下载/授权需求说明”，不继续生成结果，避免产出不完整数据。
