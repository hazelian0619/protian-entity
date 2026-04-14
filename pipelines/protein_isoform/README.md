# Protein Isoform Layer (`protein_isoform_v1`)

该管线在 **不改动 `protein_master_v6_clean.tsv`** 的前提下，新增 Protein 异构体层（isoform layer）：

- `data/output/protein/protein_isoform_v1.tsv`
- `data/output/protein/protein_isoform_map_v1.tsv`（canonical ↔ isoform）

## 口径定义

### canonical 与 isoform 边界

- **canonical**：`protein_master_v6_clean.tsv` 的 `uniprot_id` 行（主表仍保持 canonical 粒度）。
- **isoform**：从 `protein_master_v6_clean.tsv.isoforms` 里解析出的 `IsoId=...` 标记。
- **边界规则**：`Sequence=Displayed` 视为 canonical 边界；非 Displayed（如 `VSP_*`）视为 isoform 边界。

### 映射规则（canonical ↔ isoform）

- `isoform_uniprot_id` 以 `canonical_uniprot_id-` 开头 ⇒ `mapping_scope=self_accession`
- 其他前缀（跨 accession）⇒ `mapping_scope=cross_accession`
- relation_type：
  - `canonical_displayed`
  - `canonical_isoform`
  - `cross_accession_displayed`
  - `cross_accession_isoform`

### 缺失处理策略

- `isoforms` 为空或仅 `ALTERNATIVE PRODUCTS:` 占位、无 `IsoId=` 时：
  - **不生成行**（不伪造 isoform）
  - 在 QA 报告中记录 `total_without_isoid / placeholder_rows`，保证可审计

## 运行

从仓库根目录执行：

```bash
bash pipelines/protein_isoform/run.sh
```

## 产物

- Contracts
  - `pipelines/protein_isoform/contracts/protein_isoform_v1.json`
  - `pipelines/protein_isoform/contracts/protein_isoform_map_v1.json`
- QA
  - `pipelines/protein_isoform/reports/protein_isoform_v1.qa.json`
  - `pipelines/protein_isoform/reports/protein_isoform_map_v1.qa.json`
- Validation
  - `pipelines/protein_isoform/reports/protein_isoform_v1.validation.json`
  - `pipelines/protein_isoform/reports/protein_isoform_map_v1.validation.json`
- Manifest
  - `pipelines/protein_isoform/reports/protein_isoform_v1.manifest.json`
  - `pipelines/protein_isoform/reports/protein_isoform_map_v1.manifest.json`

## 审计字段

两个输出表均包含：

- `source`
- `source_version`
- `fetch_date`

用于追溯到主表来源与版本锚点。
