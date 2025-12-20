# M3（ChEMBL PSI 活性证据边）交付复盘（闭环 v1）

> 目标：从 ChEMBL SQLite（release 36）抽取 **PSI（protein–small molecule interaction）活性证据边**，形成可复跑、可追溯、可 QC 的关系层数据。
>
> 关键约束：M3 全程 **不写入** `data/output/molecules/molecules_m1.sqlite`；M3 输出写到独立库 `data/output/molecules/chembl_m3.sqlite`。

## 1) 输入与产物

### 输入

- ChEMBL SQLite：`data/raw/molecules/chembl_36/chembl_36.db`
  - 预检脚本：`pipelines/molecules/scripts/m3_check_chembl_db.py`

### 输出

- M3 SQLite：`data/output/molecules/chembl_m3.sqlite`
  - `psi_evidence_v1`：证据层（粒度=ChEMBL activity）
  - `psi_edges_v1`：聚合边（粒度=compound–target–type）
  - `target_uniprot_map_v1`：target→component→Uniprot 映射（用于审计/对齐）

- QC（smoke）：`pipelines/molecules/reports/m3/qc_m3_smoke.md`、`pipelines/molecules/reports/m3/qc_m3_smoke.json`
- QC（full）：`pipelines/molecules/reports/m3/qc_m3_full.md`、`pipelines/molecules/reports/m3/qc_m3_full.json`

## 2) v1 口径（过滤规则）

v1 默认以“抓大放小”为目标，强约束对齐与可比性：

- `standard_type IN (IC50, Ki, Kd, EC50)`
- `assay.confidence_score >= 6`
- `compound_structures.standard_inchi_key IS NOT NULL`（`--require-inchikey`）
- `target_uniprot_accession IS NOT NULL`（`--require-uniprot`）
  - 实现上选取每个 target 的一个代表 accession：优先 Swiss-Prot（`component_sequences.db_source='SWISS-PROT'`），否则取最小 accession（稳定且可复现）。

单位与强度字段：

- `pchembl_value_eff = COALESCE(pchembl_value, computed_from_standard_value_nM)`
- `standard_value_nM`：仅对常见浓度单位做换算（pM/nM/uM/mM/M），其余单位保留原值但 `standard_value_nM=NULL`。

## 3) 运行闭环（如何复现）

### 3.1 预检（只读输入库）

```bash
python3 pipelines/molecules/scripts/m3_check_chembl_db.py \
  --chembl-db data/raw/molecules/chembl_36/chembl_36.db
```

### 3.2 Smoke-run（1 万条）

```bash
python3 pipelines/molecules/scripts/m3_extract_chembl_psi.py \
  --chembl-db data/raw/molecules/chembl_36/chembl_36.db \
  --out-db data/output/molecules/chembl_m3.sqlite \
  --outdir pipelines/molecules/reports/m3 \
  --qc-prefix qc_m3_smoke \
  --types IC50,Ki,Kd,EC50 \
  --min-confidence 6 \
  --require-inchikey \
  --limit 10000 \
  --progress-every 2000
```

### 3.3 Full-run（全量）

```bash
python3 pipelines/molecules/scripts/m3_extract_chembl_psi.py \
  --chembl-db data/raw/molecules/chembl_36/chembl_36.db \
  --out-db data/output/molecules/chembl_m3.sqlite \
  --outdir pipelines/molecules/reports/m3 \
  --qc-prefix qc_m3_full \
  --types IC50,Ki,Kd,EC50 \
  --min-confidence 6 \
  --require-inchikey \
  --require-uniprot \
  --limit 0 \
  --progress-every 200000
```

### 3.4 断点/恢复（只做聚合 + QC，不重抽 evidence）

全量最耗时步骤通常是聚合边（`psi_edges_v1`）。如果 evidence 已写入但聚合被中断，可用：

```bash
python3 pipelines/molecules/scripts/m3_extract_chembl_psi.py \
  --mode aggregate \
  --chembl-db data/raw/molecules/chembl_36/chembl_36.db \
  --out-db data/output/molecules/chembl_m3.sqlite \
  --outdir pipelines/molecules/reports/m3 \
  --qc-prefix qc_m3_full \
  --types IC50,Ki,Kd,EC50 \
  --min-confidence 6 \
  --require-inchikey \
  --require-uniprot
```

## 4) 结果摘要（full）

来自 `pipelines/molecules/reports/m3/qc_m3_full.md` 与 `pipelines/molecules/reports/m3/qc_m3_full.json`：

- `psi_evidence_v1`：3,113,487
- `psi_edges_v1`：2,360,023
- `distinct_compounds_inchikey`：1,231,861
- `distinct_targets_uniprot`：8,795
- `missing_inchikey`：0
- `missing_uniprot`：0

与 M1/M2 的只读 join 覆盖率：

- `join_props_rate`（`psi_evidence_v1` → `m1.molecule_props_rdkit_v1`）：0.958236
  - 说明：这依赖 M2 全量进度，后续可能继续上升。

## 5) 性能与工程经验

- 证据抽取阶段（写入 `psi_evidence_v1`）是“顺序扫描 + 批量插入”，总体可控。
- 全量聚合阶段（生成 `psi_edges_v1`）如果用 Python 按 edge 回表选代表证据，会在百万级边上非常慢。
  - v1 最终实现采用 SQLite 的窗口函数（`ROW_NUMBER() OVER (...)`）在 SQL 内完成“每组挑 best evidence + 聚合统计”，把全量聚合耗时降到分钟级。

## 6) 已知边界与后续 v2 方向

- 单位长尾：`standard_units` 存在少量非摩尔浓度或异常单位（例如 `ug.mL-1` 等），导致 `standard_value_nM`/`pchembl_value_eff` 无法计算。
  - v2 可增加可选过滤：仅保留可换算到 nM 的单位，或引入分子量后把质量浓度转为摩尔浓度。
- Target 代表 Uniprot：复杂 target（复合物/家族）可能对应多个 components；v1 选一个代表 accession（Swiss-Prot 优先）。
  - v2 可把多 component 展开为多条 evidence（或保存 component-level 边），用于更精细的对齐。
