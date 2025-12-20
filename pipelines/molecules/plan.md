# M3（PSI：蛋白-小分子活性证据边）最小可用设计 + 实现准备（v1）

> 目标：从 ChEMBL 抽取 **PSI 活性证据边**（protein–small molecule interaction evidence），形成可追溯（assay/doc）、可 QC、可聚合的关系层数据，为后续图查询/训练集提供稳定入口。
>
> 重要约束：M3 **不对** `data/output/molecules/molecules_m1.sqlite` 做任何写操作；M3 输出写到新的库 `data/output/molecules/chembl_m3.sqlite`。

## 1) 现状对齐（M1/M2 与 check 文档）

- `pipelines/molecules/POSTMORTEM_molecules_m1_v1.md`：M1 已把小分子结构锚点（InChIKey）与 ChEMBL ID 固化，且提供 strict 子集口径。
- `check/项目大纲.pdf`：L2 关系层必须“边+证据+置信度+溯源”。PSI 是核心关系之一（与 PPI/RPI 并列）。
- `check/组件细节.docx`（PSI 部分要点）：活性数据需要带 `standard_type/value/units/relation`、`assay` 上下文、置信度（confidence_score）、引用（doi/pubmed/source）。
- `check/数据库进度.pdf`：PSI 优先数据源是 ChEMBL（其次 BindingDB / PDBbind / DrugBank…），强调 assay 上下文。

## 2) M3 最小可用 Schema（v1）

建议采用“两层表”：
- **`psi_evidence_v1`**：证据层（粒度 = ChEMBL activity）。保留原始测定上下文与引用，方便回钻。
- **`psi_edges_v1`**：边层（聚合后的 PSI 关系）。用于图查询/训练数据入口。

### 2.1 `psi_evidence_v1`（证据层，必须字段）

满足你提出的“必须包含”项：

- target 标识：
  - `target_uniprot_accession`（优先；为空时可用 `target_chembl_id` 兜底）
  - `target_chembl_id`, `target_tid`
- compound 锚点：
  - `compound_inchikey`（优先）
  - `compound_chembl_id`, `molregno`
- assay 上下文：
  - `assay_id`, `assay_type`, `assay_confidence_score`, `assay_description`, `bao_format`
- 活性核心：
  - `standard_type`, `standard_value`, `standard_units`, `standard_relation`, `pchembl_value`
- reference：
  - `doc_id`, `doi`, `pubmed_id`, `source`

v1 额外推荐字段（便于 QA/筛选/聚合）：
- `standard_value_nM`：常见浓度单位换算到 nM（不可换算则 NULL）
- `pchembl_value_eff`：`COALESCE(pchembl_value, computed_from_standard_value_nM)`
- `evidence_score_v1`：证据评分（见 2.3）

### 2.2 `psi_edges_v1`（聚合边层，最小字段）

v1 聚合键（edge key）建议：
- 优先：`(compound_inchikey, target_uniprot_accession, standard_type)`
- 兜底：若 `target_uniprot_accession` 为空，则使用 `target_chembl_id` 作为 target key。

边表建议输出：
- `compound_inchikey`, `compound_chembl_id`
- `target_uniprot_accession`（可空）, `target_chembl_id`
- `standard_type`
- `n_evidence`, `n_docs`
- `pchembl_max`, `pchembl_mean`
- `assay_confidence_score_max`, `evidence_score_max`
- 代表证据溯源：`best_activity_id`, `best_assay_id`, `best_doc_id`, `best_doi`, `best_pubmed_id`

### 2.3 证据评分（evidence_score_v1，v1 可解释方案）

目的：把 assay 上下文与测定可比性压缩为 [0,1] 分数，支持“筛边/选代表证据”。

v1 评分建议（离线、可解释）：
- `conf_term = assay_confidence_score/9`（无则 0，截断到 [0,1]）
- `assay_type_term`（经验权重）：B=1.0, F=0.9, A=0.6, T=0.4, P=0.2, 其他=0.5
- `relation_term`：`=`=1.0，`<=`/`>=`=0.9，`<`/`>`=0.7，`~`=0.5，其他/NULL=0.6
- `units_term`：nM/pM=1.0，uM=0.8，mM=0.6，M=0.5，其他/NULL=0.5

组合：

`evidence_score_v1 = 0.45*conf_term + 0.25*assay_type_term + 0.2*relation_term + 0.1*units_term`

### 2.4 去重与聚合策略（v1）

- 证据层（`psi_evidence_v1`）尽量 1:1 保留 activity（便于追溯与审计）。
- 边层聚合：同一 `compound-target-type` 多测定聚合。
- 代表证据选择（tie-break）：
  1) `pchembl_value_eff` 最大
  2) `evidence_score_v1` 最大
  3) `assay_confidence_score` 最大
  4) `activity_id` 最小（稳定）

## 3) 推荐数据源形态 + 手动下载清单（你需要做的事）

### 3.1 首选：ChEMBL SQLite dump（推荐）

- 推荐文件：`chembl_36_sqlite.tar.gz`
- 下载目录（EBI 官方 FTP）：
  - `https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_36/`

放置与解压（建议固定路径，脚本默认读这里）：

```bash
mkdir -p data/raw/molecules/chembl_36
# 把 chembl_36_sqlite.tar.gz 放到 data/raw/molecules/chembl_36/
cd data/raw/molecules/chembl_36
tar -xzf chembl_36_sqlite.tar.gz

# 最终保证存在：
# data/raw/molecules/chembl_36/chembl_36.db
```

### 3.2 次选：TSV 分表（不建议作为 v1 主路径）

- 联表口径与外键一致性更容易出错；体量更难管理。

## 4) ETL 步骤（v1）

1) 打开输入 ChEMBL SQLite（只读）
2) 创建输出 `data/output/molecules/chembl_m3.sqlite`，写入 `meta_run`
3) 抽取 `psi_evidence_v1`
   - 过滤：`standard_value` 非空，`standard_type` 在白名单（默认 IC50/Ki/Kd/EC50）
   - join：`activities` × `assays` × `target_dictionary` × `molecule_dictionary` × `compound_structures` × `docs` ×（`target_components` × `component_sequences`）
   - 计算：`standard_value_nM`、`pchembl_value_eff`、`evidence_score_v1`
4) 构建 `psi_edges_v1` 聚合边

## 5) QC（验收清单）

建议脚本输出：
- `pipelines/molecules/reports/m3/qc_m3_chembl.md`
- `pipelines/molecules/reports/m3/qc_m3_chembl.json`

最小 QC：
- 规模：evidence 行数、edge 行数
- 覆盖：distinct compound(inchikey) / target(uniprot、chembl)
- 分布：`standard_type`、`standard_units`、`assay_type`、`assay_confidence_score`
- 缺失率：inchikey/uniprot/doi/pubmed/pchembl 缺失

## 6) “抓大放小”路线建议

- 推荐：**先做 PSI 子集 → 再算属性**
  - 先用 M3 得到“有活性证据的化合物子集”，再对该子集做 RDKit/3D/更重特征；计算量与存储更可控，且更贴近 L2/L4。
- 备选：strict 全量属性 → 再做 PSI
  - 得到通用化学空间底座，但短期成本更高。

## 7) 对应脚本

- `pipelines/molecules/scripts/m3_extract_chembl_psi.py`：ChEMBL SQLite → `data/output/molecules/chembl_m3.sqlite`
