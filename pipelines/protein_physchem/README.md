# Protein Physchem Core Layer (`protein_physchem_v1`)

从 `protein_master_v6_clean.tsv` 的序列计算核心理化指标，生成 **core 层**：

- `data/output/protein/protein_physchem_v1.tsv`

## Core 指标

- `mass_recalc`
- `isoelectric_point`
- `gravy`
- `aromaticity`
- `instability_index`

> `protein_physchem_v1` 保持稳定，不被扩展层替换。

## 与 Extended 层关系

- Extended 层：`data/output/protein/protein_physchem_extended_v1.tsv`
- Extended 会复算并对齐 core 指标，同时新增：
  - `extinction_coefficient`
  - `extinction_coefficient_reduced`
  - `aliphatic_index`

## 依赖

```bash
python3 -m pip install biopython pandas
```

## 运行

从仓库根目录执行：

```bash
bash pipelines/protein_physchem/run.sh
```

## 输出与报告

- Contract：`pipelines/protein_physchem/contracts/protein_physchem_v1.json`
- Metrics：`pipelines/protein_physchem/reports/protein_physchem_v1.metrics.json`
- Validation：`pipelines/protein_physchem/reports/protein_physchem_v1.validation.json`
- Manifest：`pipelines/protein_physchem/reports/protein_physchem_v1.manifest.json`

## 特殊残基策略

为保证可解析率：
- `U -> C`
- `O -> K`

替换计数记录在 metrics 报告中。
