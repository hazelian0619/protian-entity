# Protein 主表基线（`protein_master_v6_clean.tsv`）

本目录只负责对既有主表做 **Contract Validation + Manifest**，不重建主表、不改动主表内容。

- 输入：`data/processed/protein_master_v6_clean.tsv`
- 输出：
  - `pipelines/protein/reports/protein_master_v6.validation.json`
  - `pipelines/protein/reports/protein_master_v6.manifest.json`

## 向后兼容声明

- `protein_master_v6_clean.tsv` 仍是 Protein 线 canonical 主表（兼容既有下游）。
- isoform 与 physchem 扩展均以**新增层**提供，不替换主表。

## 分层关系（Protein only）

- **Canonical 主层**：`data/processed/protein_master_v6_clean.tsv`
- **Isoform 层（新增）**：
  - `data/output/protein/protein_isoform_v1.tsv`
  - `data/output/protein/protein_isoform_map_v1.tsv`
- **Physchem Core 层**：`data/output/protein/protein_physchem_v1.tsv`
- **Physchem Extended 层（新增）**：`data/output/protein/protein_physchem_extended_v1.tsv`

## 运行

从仓库根目录执行：

```bash
bash pipelines/protein/run.sh
```

成功标志：`[PASS] protein_master_v6 ...`

## Contract 与门禁

- Contract：`pipelines/protein/contracts/protein_master_v6.json`
- 校验工具：`python3 tools/kg_validate_table.py`
- 失败定位：查看 `pipelines/protein/reports/protein_master_v6.validation.json` 中 `passed=false` 的规则

## 常见误解

1. **“会不会生成新的 protein 主表？”**
   - 不会。本 pipeline 只做 QA + manifest。

2. **“source_version 在主表里吗？”**
   - 主表合同保持不变；版本锚点维护在 `data/output/protein/protein_source_versions_v1.tsv`。

3. **“如何跑完整 Protein 层？”**
   - 主表基线：`bash pipelines/protein/run.sh`
   - isoform 层：`bash pipelines/protein_isoform/run.sh`
   - physchem core：`bash pipelines/protein_physchem/run.sh`
   - physchem extended：`bash pipelines/protein_physchem_extended/run.sh`
