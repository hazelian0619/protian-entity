# Protein Physchem Extended Layer (`protein_physchem_extended_v1`)

该管线在保留 `protein_physchem_v1`（core）不变的前提下，新增扩展理化层：

- `data/output/protein/protein_physchem_extended_v1.tsv`

## Core vs Extended 分层

- **Core (`protein_physchem_v1`)**：
  - `mass_recalc`, `isoelectric_point`, `gravy`, `aromaticity`, `instability_index`
- **Extended (`protein_physchem_extended_v1`)**：
  - 完整保留 core 指标（用于对齐比对）
  - 额外新增：
    - `extinction_coefficient`（氧化态，M^-1 cm^-1）
    - `extinction_coefficient_reduced`（还原态，M^-1 cm^-1）
    - `aliphatic_index`（无量纲）

## 计算口径

- `instability_index`：`Bio.SeqUtils.ProtParam.ProteinAnalysis.instability_index`
- `gravy`：`ProteinAnalysis.gravy`（Kyte-Doolittle）
- `isoelectric_point`：`ProteinAnalysis.isoelectric_point`
- `extinction_coefficient`：`ProteinAnalysis.molar_extinction_coefficient()` 的 **oxidized** 返回值
- `aliphatic_index`：Ikai 公式 `X(Ala)+2.9*X(Val)+3.9*(X(Ile)+X(Leu))`，X 为摩尔百分比

特殊残基处理：
- `U -> C`
- `O -> K`

用于保持与 core 一致的可解析率策略。

## 运行

从仓库根目录执行：

```bash
bash pipelines/protein_physchem_extended/run.sh
```

## 产物

- Contract：`pipelines/protein_physchem_extended/contracts/protein_physchem_extended_v1.json`
- QA：`pipelines/protein_physchem_extended/reports/protein_physchem_extended_v1.qa.json`
- Validation：`pipelines/protein_physchem_extended/reports/protein_physchem_extended_v1.validation.json`
- Manifest：`pipelines/protein_physchem_extended/reports/protein_physchem_extended_v1.manifest.json`

## 审计字段

输出包含：

- `source`
- `source_version`
- `fetch_date`

满足 lineage 审计与可追溯要求。
