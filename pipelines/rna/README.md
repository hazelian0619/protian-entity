# RNA L1（v1）：miRNA + mRNA/transcript（Human 9606）

这条管线产出知识图谱的 **RNA 实体表（L1）**。

## 同事怎么用（最短路径）

1) 直接下载数据（不需要跑代码）：
- Release: https://github.com/hazelian0619/protian-entity/releases/tag/rna-l1-v1

2) 你会拿到这些文件（都在 Release assets 里）：
- `rna_master_v1.tsv.gz`：全量 RNA 主表（miRNA + mRNA/transcript）
- `rna_master_mirna_v1.tsv.gz`：miRNA 子表
- `rna_master_mrna_v1.tsv.gz`：mRNA/transcript 子表
- `rna_master_v1.validation.json`：质量门禁（QA）报告（应为 PASS）
- `manifest.json`：每个产物的 sha256/大小/生成时间/commit（用于校验与追溯）

3) 读取数据（示例）：

```python
import pandas as pd

# 直接读取压缩 TSV
rna = pd.read_csv('rna_master_v1.tsv.gz', sep='\t', compression='gzip')
print(rna.head())
```

## 表结构与规则（必读）

- 数据字典：`docs/rna/DATA_DICTIONARY_RNA.md`
- 数据源与版本锚点：`docs/rna/RNA_SOURCES_AND_VERSIONS.md`

## 需要本地重跑（开发者）

从仓库根目录运行：

```bash
bash pipelines/rna/run.sh
```

说明：

- 本地输出默认写到 `data/output/`（大文件，默认不提交到 git）
- 输入数据放在 `data/raw/rna/`（具体需要哪些文件见 `docs/rna/RNA_SOURCES_AND_VERSIONS.md`）

## QA 与 manifest（本地跑完后）

```bash
python3 tools/kg_validate_table.py \
  --contract pipelines/rna/contracts/rna_master_v1.json \
  --table data/output/rna_master_v1.tsv \
  --out pipelines/rna/reports/rna_master_v1.validation.json

python3 tools/kg_make_manifest.py \
  --data-version kg-data-local \
  --out pipelines/rna/reports/rna_v1.manifest.json \
  data/output/rna_master_v1.tsv \
  data/output/rna_master_mirna_v1.tsv \
  data/output/rna_master_mrna_v1.tsv
```
