# Protein Data Dictionary (canonical / isoform / physchem)

> 适用范围：**Protein 线**。本文件不定义 RNA / Molecule / Interaction 合同。

## 1) 分层总览

| Layer | 表 | 粒度 | 主键 | 兼容性 |
|---|---|---|---|---|
| Canonical Core | `data/processed/protein_master_v6_clean.tsv` | canonical protein（每行一个 `uniprot_id`） | `uniprot_id` | 既有主表，保持不变 |
| Isoform Extended | `data/output/protein/protein_isoform_v1.tsv` | canonical × isoform | `isoform_record_id` | 新增层，不替换主表 |
| Isoform Map | `data/output/protein/protein_isoform_map_v1.tsv` | canonical ↔ isoform 映射 | `map_id` | 新增层，用于 join |
| Physchem Core | `data/output/protein/protein_physchem_v1.tsv` | canonical protein | `uniprot_id` | 既有 core，保持不变 |
| Physchem Extended | `data/output/protein/protein_physchem_extended_v1.tsv` | canonical protein | `uniprot_id` | 新增层，向后兼容 |

---

## 2) Canonical 主表：`protein_master_v6_clean.tsv`

### 主键与定位
- **主键**：`uniprot_id`
- **定位**：Protein 线 canonical 主实体，不展开 isoform。

### 关键字段
- 标识：`uniprot_id`, `entry_name`, `protein_name`, `gene_names`, `symbol`, `hgnc_id`
- 序列：`sequence`, `sequence_len`, `mass`
- 功能：`function`, `go_biological_process`, `go_molecular_function`, `go_cellular_component`
- 结构：`pdb_ids`, `has_alphafold`, `alphafold_*`
- 异构体原始信息：`isoforms`（文本，不是标准化表）
- 血缘：`source`, `date_modified`, `fetch_date`

### 质量门槛（以合同为准）
- 合同：`pipelines/protein/contracts/protein_master_v6.json`
- 校验报告：`pipelines/protein/reports/protein_master_v6.validation.json`

---

## 3) Isoform 标准化层：`protein_isoform_v1.tsv`

### 粒度与边界
- **canonical 定义**：主表中的 `uniprot_id`
- **isoform 定义**：从主表 `isoforms` 字段解析的 `IsoId=...`
- **边界规则**：
  - `Sequence=Displayed` → `boundary_class=canonical`
  - 其他（如 `VSP_*`）→ `boundary_class=isoform`

### 主键
- `isoform_record_id = canonical_uniprot_id|isoform_uniprot_id`

### 关键字段
- 键与关系：`canonical_uniprot_id`, `isoform_uniprot_id`, `mapping_scope`, `relation_type`
- 解析上下文：`event_types`, `named_isoforms_declared`, `isoform_name_token`, `isoform_synonyms`, `sequence_status_raw`
- 审计字段：`source`, `source_version`, `fetch_date`

### 缺失处理策略
- `isoforms` 为空或仅占位文本、无 `IsoId=`：**不产出 isoform 行**。
- 缺失规模在 QA 报告披露，不伪造补全。

### 合同与报告
- 合同：`pipelines/protein_isoform/contracts/protein_isoform_v1.json`
- QA：`pipelines/protein_isoform/reports/protein_isoform_v1.qa.json`
- Validation：`pipelines/protein_isoform/reports/protein_isoform_v1.validation.json`

---

## 4) Canonical ↔ Isoform 映射层：`protein_isoform_map_v1.tsv`

### 作用
提供轻量映射表，供下游 join 与图关系建模。

### 主键
- `map_id = canonical_uniprot_id|isoform_uniprot_id`

### 关键字段
- 映射：`canonical_uniprot_id`, `isoform_uniprot_id`
- 关系：`boundary_class`, `relation_type`, `mapping_scope`, `sequence_status_class`
- 审计字段：`source`, `source_version`, `fetch_date`

### 合同与报告
- 合同：`pipelines/protein_isoform/contracts/protein_isoform_map_v1.json`
- QA：`pipelines/protein_isoform/reports/protein_isoform_map_v1.qa.json`
- Validation：`pipelines/protein_isoform/reports/protein_isoform_map_v1.validation.json`

---

## 5) Physchem Core 层：`protein_physchem_v1.tsv`

### 定位
稳定 core 指标层（不破坏历史下游）。

### 主键
- `uniprot_id`

### 指标
- `mass_recalc`
- `isoelectric_point`
- `gravy`
- `aromaticity`
- `instability_index`

---

## 6) Physchem Extended 层：`protein_physchem_extended_v1.tsv`

### 定位
在 core 基础上新增扩展指标；core 字段保留并对齐。

### 主键
- `uniprot_id`

### 字段（新增/对齐）
- 对齐 core：`mass_recalc`, `isoelectric_point`, `gravy`, `aromaticity`, `instability_index`
- 新增：
  - `extinction_coefficient`（M^-1 cm^-1，280nm，oxidized）
  - `extinction_coefficient_reduced`（M^-1 cm^-1，280nm，reduced）
  - `aliphatic_index`（unitless）
- 审计字段：`source`, `source_version`, `fetch_date`

### 计算口径
- 引擎：Biopython `ProteinAnalysis`
- `instability_index/gravy/isoelectric_point` 与 core 同口径
- `aliphatic_index`：Ikai 公式
- 特殊残基：`U→C`, `O→K`

### 合同与报告
- 合同：`pipelines/protein_physchem_extended/contracts/protein_physchem_extended_v1.json`
- QA：`pipelines/protein_physchem_extended/reports/protein_physchem_extended_v1.qa.json`
- Validation：`pipelines/protein_physchem_extended/reports/protein_physchem_extended_v1.validation.json`

---

## 7) 审计与质量门禁（Protein）

统一检查维度：
- 主键唯一率（100%）
- 必填字段非空率（100%）
- `fetch_date` 日期格式（`YYYY-MM-DD`）
- `source/source_version/fetch_date` 可审计

执行入口（仓库根目录）：

```bash
bash pipelines/protein/run.sh
bash pipelines/protein_isoform/run.sh
bash pipelines/protein_physchem/run.sh
bash pipelines/protein_physchem_extended/run.sh
```

---

## 8) 向后兼容声明

- 不修改 `protein_master_v6_clean.tsv` 结构与合同。
- 不修改 `protein_physchem_v1.tsv` 结构与合同。
- 新增层均为 additive，可被选择性消费。
- 本文档不引入 RNA / Molecule / Interaction 的数据逻辑或合同变更。
