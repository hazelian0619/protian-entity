# molecule_3d_registry（v1）

构建小分子 3D 结构**引用层（registry）**，不做结构重建模：

- PubChem 3D：基于 CID 批量查询 `ConformerCount3D`
- ZINC 元数据：基于 ChEMBL→ZINC 同义词映射（本地 `chembl_36.db`）

并统一输出来源、格式、可用性、类型（`computational/experimental`）等字段。

---

## 输入

- `data/output/molecules/molecule_xref_core_v2.tsv`（优先）
  - 若无，自动回退 `molecule_xref_core_v1.tsv`
- `data/raw/molecules/chembl_36/chembl_36.db`

`run.sh` 在当前协作工作区支持 fallback 路径（`../1218/...`、`../12182/...`）。

---

## 输出

- 表：`data/output/molecules/molecule_3d_registry_v1.tsv`
- Contract：`pipelines/molecule_3d_registry/contracts/molecule_3d_registry_v1.json`
- 报告：
  - `molecule_3d_registry_v1.build.json`
  - `molecule_3d_registry_v1.coverage.json`
  - `molecule_3d_registry_v1.validation.json`
  - `molecule_3d_registry_v1.qa.json`
  - `molecule_3d_registry_v1.manifest.json`

---

## 字段说明（核心）

- `inchikey`：主键，可回链 xref core
- `pubchem_3d_available` / `pubchem_conformer_count` / `pubchem_download_url`
- `zinc_3d_available` / `zinc_id` / `zinc_download_url`
- `source_providers`：来源聚合（`PubChem`/`ZINC`/`none`）
- `structure_formats`：格式聚合（如 `SDF;JSON`）
- `availability_status`：`available` / `unavailable`
- `source_types`：`computational` / `experimental` / `mixed` / `unknown`

> 当前 v1 主要是计算结构来源（PubChem/ZINC），因此 `source_types` 以 `computational` 为主。

---

## 运行

```bash
bash pipelines/molecule_3d_registry/run.sh
```

可选参数：

```bash
SMOKE_MAX_ROWS=5000 PUBCHEM_CHUNK_SIZE=100 PUBCHEM_TIMEOUT=30 bash pipelines/molecule_3d_registry/run.sh
```

---

## 流程

`run.sh` 固定执行：

1. smoke build（小样本）
2. smoke contract validation + smoke QA
3. full build
4. full contract validation + full QA
5. manifest

---

## 缺数据处理

若缺必需输入，pipeline 退出并产出：

- `pipelines/molecule_3d_registry/reports/molecule_3d_registry_v1.manual_download.json`

含下载清单、放置路径和 `sha256sum` 校验命令。
