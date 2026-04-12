# molecule_semantic_layer（v1）

构建小分子语义层标签表：

- **ChEBI 本体映射**（基于 `chebi.obo` 的 InChIKey 对齐）
- **ChEBI 父子层级展开**（`is_a` 路径）
- **ATC 分类映射**（ChEMBL ATC 表）
- **DrugBank 分类字段**（`groups`）
- **ZINC 标签**（若有 `molecule_3d_registry_v1.tsv` 则读取）

---

## 输入

必需：

- `data/output/molecules/molecule_xref_core_v2.tsv`（若无自动降级到 v1）
- `data/raw/molecules/chebi/chebi.obo`
- `data/output/drugbank/drug_master_v1.tsv`
- `data/raw/molecules/chembl_36/chembl_36.db`

可选：

- `data/output/molecules/molecule_3d_registry_v1.tsv`（用于 ZINC 标签增强）

---

## 输出

- 主表：`data/output/molecules/molecule_semantic_tags_v1.tsv`
- Contract：`pipelines/molecule_semantic_layer/contracts/molecule_semantic_tags_v1.json`
- 报告：
  - `molecule_semantic_tags_v1.build.json`
  - `molecule_semantic_tags_v1.hierarchy.json`
  - `molecule_semantic_tags_v1.coverage.json`
  - `molecule_semantic_tags_v1.validation.json`
  - `molecule_semantic_tags_v1.qa.json`
  - `molecule_semantic_tags_v1.manifest.json`

---

## 关键字段

- `inchikey`（主键）
- `chebi_id / chebi_name / chebi_ancestor_ids / chebi_paths`
- `atc_code / atc_level1_code / atc_level1_description / ...`
- `drugbank_group`
- `zinc_id / zinc_label`
- `semantic_tags`
- `chebi_source / atc_source / drugbank_source / zinc_source / source_version`

---

## 运行

```bash
bash pipelines/molecule_semantic_layer/run.sh
```

可选参数：

```bash
SMOKE_MAX_ROWS=5000 bash pipelines/molecule_semantic_layer/run.sh
```

---

## smoke → full

`run.sh` 固定流程：

1. smoke build
2. smoke contract validation + smoke QA
3. full build
4. full contract validation + full QA
5. manifest

---

## 验收映射

- `inchikey` 回链率（QA gate）`>= 0.99`
- ChEBI 层级可展开（`hierarchy` + QA gate）
- 分类字段来源与版本锚点（source columns + `source_version`）

---

## 缺输入处理

若缺必需输入，pipeline 会生成：

- `pipelines/molecule_semantic_layer/reports/molecule_semantic_tags_v1.manual_download.json`

包含下载清单、放置路径和 `sha256sum` 校验命令。
