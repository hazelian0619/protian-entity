# Unified Protein Schema v2（历史草案，已弃用）

> 状态：**Deprecated / 历史设计稿**。本文件不再作为当前数据合同来源。

## 为什么弃用

旧版内容与当前仓库实现存在漂移（例如：行数规模、字段命名、是否“每行一异构体”等），若继续作为执行依据会导致误导。

## 当前生效的 Protein 合同与字典

- Canonical 主表合同：`pipelines/protein/contracts/protein_master_v6.json`
- Isoform 层合同：
  - `pipelines/protein_isoform/contracts/protein_isoform_v1.json`
  - `pipelines/protein_isoform/contracts/protein_isoform_map_v1.json`
- Physchem core 合同：`pipelines/protein_physchem/contracts/protein_physchem_v1.json`
- Physchem extended 合同：`pipelines/protein_physchem_extended/contracts/protein_physchem_extended_v1.json`
- 当前数据字典（Protein）：`docs/DATA_DICTIONARY.md`

## 迁移说明

若你在下游文档中引用本文件，请迁移到 `docs/DATA_DICTIONARY.md` 与对应 pipeline README，避免口径漂移。
