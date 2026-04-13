# molecule_zinc_direct_xref（v1）

把 ZINC 从“桥接补丁”升级为“直连映射层”：

- 使用 ZINC 官方 catalogs 可下载文件（`catalogs/<tier>/<vendor>/<vendor>.info.txt(.gz)`）
- 基于 `InChIKey -> ZINC ID` 直接映射
- 补充 `SMILES`（`catalogs/source/<vendor>.src.txt(.gz)`，若可得）
- 输出可购买/层级标签与冲突审计

## 输入

- `data/output/molecules/molecule_xref_core_v2.tsv`（优先，缺失回退 v1）
- ZINC 官方公开目录（默认 `http://files.docking.org/catalogs/`）

## 输出

- `data/output/molecules/molecule_zinc_xref_v1.tsv`
- `pipelines/molecule_zinc_direct_xref/contracts/molecule_zinc_xref_v1.json`
- `pipelines/molecule_zinc_direct_xref/reports/*`
  - `molecule_zinc_xref_v1.build.json`
  - `molecule_zinc_xref_v1.coverage.json`
  - `molecule_zinc_xref_v1.conflict_audit.json`
  - `molecule_zinc_xref_v1.validation.json`
  - `molecule_zinc_xref_v1.qa.json`
  - `molecule_zinc_xref_v1.manifest.json`

## 运行

```bash
bash pipelines/molecule_zinc_direct_xref/run.sh
```

可选参数：

```bash
SMOKE_MAX_VENDORS=10 FULL_MAX_VENDORS=80 ZINC_HTTP_TIMEOUT=45 bash pipelines/molecule_zinc_direct_xref/run.sh
```

## 验收门禁（对应需求）

- 覆盖率显著高于基线（默认基线取 `molecule_3d_registry_v1.coverage.json` 的 `zinc_3d_rate`）
- 每条映射带 `source/source_version/fetch_date`
- 多映射冲突有审计报告（`conflict_audit.json`）

