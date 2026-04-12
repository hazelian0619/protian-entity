# molecule_activity_fusion（ChEMBL + BindingDB + PDBbind 活性证据融合）

目标：将小分子-靶点活性证据从 ChEMBL 单源扩展为多源融合，统一活性类型与单位，并保留可审计的原始来源字段。

---

## 输入

必需：

- `chembl_m3.sqlite`（默认优先：`data/output/molecules/chembl_m3.sqlite`，fallback：`../12182/out/m3/chembl_m3.sqlite`）
  - 表：`psi_evidence_v1`、`psi_edges_v1`

可选（自动下载）：

- BindingDB PubChem TSV ZIP（默认落地：`pipelines/molecule_activity_fusion/.cache/BindingDB_PubChem_YYYYMM_tsv.zip`）
- PDBbind v2020 plain-text index（默认落地：`pipelines/molecule_activity_fusion/.cache/PDBbind_v2020_plain_text_index.tar.gz`）

> 若外部下载失败，流程继续（ChEMBL-only + 已可用源），并在报告里输出手动下载清单与校验命令。

---

## 输出

- `data/output/evidence/molecule_activity_evidence_v2.tsv`
- `data/output/edges/molecule_target_edges_v2.tsv`
- `pipelines/molecule_activity_fusion/contracts/*.json`
- `pipelines/molecule_activity_fusion/reports/*.json`

核心报告：

- `molecule_activity_fusion_v2.build.json`
- `molecule_activity_fusion_v2.qa.json`
- `molecule_activity_evidence_v2.validation.json`
- `molecule_target_edges_v2.validation.json`
- `molecule_activity_conflict_audit_v2.tsv`
- `molecule_activity_conflict_audit_v2.json`
- `molecule_activity_fusion_v2.manifest.json`

---

## 证据表字段（关键）

证据表包含并门禁以下字段：

- `source_db`
- `assay_type`
- `standard_type`（`IC50/Ki/Kd/EC50`）
- `standard_value`
- `standard_unit`
- `normalized_nM`
- `confidence`（`high/medium/low`）

并保留来源与参考字段：

- `source_record_id`
- `reference_doi` / `reference_pubmed_id` / `reference_pdb_id`
- `source_raw`
- `source_version`

---

## 冲突审计

在 `(compound_inchikey, target_uniprot_accession, standard_type)` 粒度上聚合多源记录：

- 输出 `source_count > 1` 的多源组到 `molecule_activity_conflict_audit_v2.tsv`
- 给出 `min_nM / max_nM / fold_change` 与 `conflict_flag`（默认阈值 `fold_change >= 10`）

---

## 运行

```bash
bash pipelines/molecule_activity_fusion/run.sh
```

可选：

```bash
SMOKE_ROWS=50000 DATA_VERSION=kg-data-local bash pipelines/molecule_activity_fusion/run.sh
```

---

## smoke -> full

`run.sh` 固定执行：

1. smoke build（限制 ChEMBL 行数）
2. smoke contract + smoke QA
3. full build
4. full contract + full QA（严格门禁）
5. manifest

---

## 质量门禁（full）

- 证据表包含必需字段并通过 contract
- 冲突审计表存在且可解析
- ChEMBL 基线不退化（与 `psi_evidence_v1` 同口径比较）
  - `v2` 中 `source_db=chembl` 的证据条数不低于基线
  - `v2` 中包含 chembl 来源的 edge 条数不低于基线
