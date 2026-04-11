# Protein 序列理化属性管线（protein_physchem_v1）

从 `data/processed/protein_master_v6_clean.tsv` 的 `uniprot_id, sequence` 计算标准理化属性，输出辅助表：

- `data/output/protein/protein_physchem_v1.tsv`

## 输出字段

- `uniprot_id`
- `sequence_len`（由 `len(sequence)` 重算）
- `mass_recalc`
- `isoelectric_point`
- `gravy`
- `aromaticity`
- `instability_index`

## 目录结构

- `pipelines/protein_physchem/run.sh`
- `pipelines/protein_physchem/scripts/build_protein_physchem.py`
- `pipelines/protein_physchem/contracts/protein_physchem_v1.json`
- `pipelines/protein_physchem/reports/*.json`

## 依赖

```bash
python3 -m pip install biopython pandas
```

> 若环境缺包会直接报错并退出，请先安装后重试。

## 运行

从仓库根目录执行：

```bash
bash pipelines/protein_physchem/run.sh
```

该脚本会：
1. 先跑最小样本（`--limit 20`）
2. 再跑全量
3. 按 contract 做验证并输出 JSON 报告

## 报告说明

- `protein_physchem_v1.sample.metrics.json`：最小样本验收
- `protein_physchem_v1.metrics.json`：全量覆盖率/数值可解析率/`sequence_len` 差异率
- `protein_physchem_v1.validation.json`：contract 校验结果

## 特殊残基处理

Biopython 的 ProtParam 不支持 `U/O` 直接计算。为保证 100% 可解析率，本管线在计算前做替换：

- `U -> C`
- `O -> K`

替换计数会记录在 `metrics.json` 的 `special_residue_substitution.counts`。
