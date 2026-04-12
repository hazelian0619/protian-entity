# molecule_xref_enrichment_v2（DrugBank / ChEMBL / PubChem 统一回填）

目标：把 `molecule_xref_core_v1` 升级为可审计的 `v2`，并满足：

- `inchikey` 唯一率 100%
- `chembl_id` 覆盖率相对 v1 提升（报告给绝对增量）
- `pubchem_cid` 覆盖率相对 v1 提升（报告给绝对增量）
- DrugBank 缺 InChIKey 条目执行二阶段回填，输出 `match_strategy` + `confidence`

---

## 输入

默认输入（可用环境变量覆盖）：

- `data/output/molecules/molecule_xref_core_v1.tsv` (`V1_CORE_PATH`)
- `data/output/drugbank/drug_master_v1.tsv` (`DRUG_MASTER_PATH`)
- `data/output/drugbank/drug_xref_molecules_v1.tsv` (`DRUG_XREF_PATH`)
- `data/output/molecules/molecules_m1.sqlite` (`M1_DB_PATH`)
- `data/raw/molecules/chembl_36/chembl_36.db` (`CHEMBL_DB_PATH`)

`run.sh` 内置了当前协作工作区的 fallback 路径（`../1218/...`, `../12182/...`），缺默认输入时会自动尝试。

---

## 输出

- 主表：`data/output/molecules/molecule_xref_core_v2.tsv`
- Contract：`pipelines/molecule_xref_enrichment_v2/contracts/molecule_xref_core_v2.json`
- 报告：
  - `pipelines/molecule_xref_enrichment_v2/reports/molecule_xref_core_v2.build.json`
  - `pipelines/molecule_xref_enrichment_v2/reports/molecule_xref_core_v2.backfill_audit.json`
  - `pipelines/molecule_xref_enrichment_v2/reports/molecule_xref_core_v2.missing_audit.json`
  - `pipelines/molecule_xref_enrichment_v2/reports/molecule_xref_core_v2.validation.json`
  - `pipelines/molecule_xref_enrichment_v2/reports/molecule_xref_core_v2.qa.json`
  - `pipelines/molecule_xref_enrichment_v2/reports/molecule_xref_core_v2.manifest.json`

---

## 字段（v2）

- `inchikey`
- `chembl_id`（`;` 分隔）
- `drugbank_id`（`;` 分隔）
- `pubchem_cid`（`;` 分隔）
- `match_strategy`（`;` 分隔，含 backfill 策略）
- `confidence`（`high/medium/low`）
- `xref_source`
- `fetch_date`
- `source_version`

---

## 二阶段回填策略

1. `backfill_stage1_name_exact_unique`
   - DrugBank `name` 与 ChEMBL `pref_name` 做大小写不敏感精确匹配
   - 且候选唯一（`chembl_id + inchikey` 唯一）
   - `confidence=high`

2. `backfill_stage2_synonym_exact_unique`
   - 对 stage1 未命中的条目，使用 DrugBank `synonyms` 与 ChEMBL `molecule_synonyms` 匹配
   - 同样要求候选唯一
   - `confidence=medium`

未命中或歧义条目进入 `missing_audit`。

---

## 运行

在仓库根目录：

```bash
bash pipelines/molecule_xref_enrichment_v2/run.sh
```

可选：

```bash
SMOKE_MAX_ROWS=8000 DATA_VERSION=kg-data-local bash pipelines/molecule_xref_enrichment_v2/run.sh
```

---

## smoke -> full

`run.sh` 固定执行：

1. smoke build（抽样）
2. smoke contract validation + smoke QA
3. full build
4. full contract validation + full QA（强制 coverage 正增量 gate）
5. manifest

---

## 缺数据处理

如 5 个输入中任一缺失且 fallback 仍找不到：

- 直接失败退出
- 输出 `pipelines/molecule_xref_enrichment_v2/reports/molecule_xref_core_v2.manual_download.json`
- 含下载清单、目标放置路径、`sha256sum` 校验命令
