# molecule_3d_experimental_linker（v1）

新增小分子 **实验 3D 证据层（registry-only）**，不做结构重建。

核心链路：

1. 从 RCSB GraphQL 读取 PDB entry 的 nonpolymer ligand 信息（含 `InChIKey`）
2. 通过 `pdb_structures_v1.tsv` 的 PDB→UniProt 映射拿到 target
3. 与 `molecule_xref_core_v2.tsv` 用 InChIKey 回链
4. 输出实验结构证据表（`source_type=experimental`）及版本锚点

---

## 输入

- `data/output/molecules/molecule_xref_core_v2.tsv`（优先，缺失时回退 v1）
- `data/output/protein/pdb_structures_v1.tsv`
- RCSB Data API（GraphQL）

---

## 输出

- 表：`data/output/molecules/molecule_3d_experimental_v1.tsv`
- Contract：`pipelines/molecule_3d_experimental_linker/contracts/molecule_3d_experimental_v1.json`
- 报告：
  - `molecule_3d_experimental_v1.build.json`
  - `molecule_3d_experimental_v1.coverage.json`
  - `molecule_3d_experimental_v1.validation.json`
  - `molecule_3d_experimental_v1.qa.json`
  - `molecule_3d_experimental_v1.manifest.json`

---

## 运行

```bash
bash pipelines/molecule_3d_experimental_linker/run.sh
```

可选参数：

```bash
SMOKE_MAX_PDB=800 FULL_MAX_PDB=5000 RCSB_BATCH_SIZE=100 RCSB_HTTP_TIMEOUT=35 bash pipelines/molecule_3d_experimental_linker/run.sh
```

---

## 验收对应

- `source_type=experimental` 行数 > 0：由 QA gate `experimental_rows_present` 保证
- InChIKey 回链率 >= 0.99（可映射子集）：QA gate `inchikey_backlink_rate_mappable_subset`
- `pdb_id/experimental_method/resolution` 非空率：QA gate `pdb_method_resolution_non_empty_rate`
- 明确区分 experimental vs computational：`source_type` 严格为 `experimental`

