# PSI 条件字段增强 v3（P0）

目标：在 `psi_activity_context_v2.tsv` 基础上增强条件字段覆盖率，输出：

- `data/output/evidence/psi_activity_context_v3.tsv`
- `data/output/evidence/psi_condition_parse_audit_v3.tsv`

并保证 v2 主干字段不退化（`edge_id / standard_type / standard_relation / assay_type`）。

---

## 输入

必需：

- `data/output/evidence/psi_activity_context_v2.tsv`

可选（仅记录可用性，不作为阻塞）：

- `data/output/molecules/chembl_m3.sqlite`
- `CHEMBL_M3_DB` 指向的外部 ChEMBL SQLite
- `data/input/bindingdb/bindingdb_conditions.tsv`

---

## 抽取增强（v3）

### 1) 多字段抽取（至少四个源字段）

- `assay_description`
- `assay_context`
- `activity_comment`
- `data_validity_comment`

### 2) 规则抽取 + 标准化

- **pH**：支持 `pH 7.4` / `pH 7-8` / `between pH 7 and 8`
- **温度**：支持 `°C/°F/K`，统一到 `condition_temperature_c`
- **结构化上下文**：抽取 system / cell line / buffer

并带有推断规则（低置信度）：

- buffer 默认 pH（如 PBS/HEPES/Tris）
- cell culture + incubation 推断 pH=7.4、temperature=37

### 3) 新增字段

- `condition_context_json`
- `condition_extract_confidence`（0-1）
- `condition_extract_source_field`

### 4) 冲突审计

- 同一 `activity_id` 多来源条件值冲突时输出审计行：
  - `conflict_flag=true`
  - `conflict_fields`
  - `raw_condition_candidates_json`
  - `resolution_strategy`

---

## 执行

```bash
bash pipelines/psi_condition_enrichment/run.sh
```

流程：

1. preflight + `blocked_missing_inputs.json`
2. extractor 单测
3. smoke（默认 50k）
4. full
5. contract validation
6. QA gates
7. result evidence + sample20
8. manifest + summary

可选参数：

- `SMOKE_ROWS`（默认 50000）
- `FETCH_DATE`
- `DATA_VERSION`
- `CHEMBL_M3_DB`

---

## 验收门禁（硬门禁）

- edge_id 覆盖率保持 1.0
- standard_type / standard_relation / assay_type 覆盖率不低于 v2
- condition_context >= 0.55
- condition_pH >= 0.22
- condition_temperature_c >= 0.15
- contract validation PASS

目标门禁（尽量）：

- condition_context >= 0.65
- condition_pH >= 0.30
- condition_temperature_c >= 0.22
