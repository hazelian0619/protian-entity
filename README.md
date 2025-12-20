# 人类知识图谱数据集（Protein + RNA）

这个仓库按“工业级数据产品”的方式组织：

- **代码 / 规范 / 质量报告**进入仓库（可审计、可复现）
- **体积大的数据产物**通过 **GitHub Releases** 发布（可下载、可校验、可回滚）

## 快速入口（给同事看这一段就够）

- **仓库怎么用**：`docs/REPO_USAGE.md`
- **Protein（L1）数据集**：`data/processed/protein_master_v6_clean.tsv`（仓库内可直接下载）
- **Protein（L1）工业化入口**：`pipelines/protein/README.md`（QA + manifest，一键验证与追溯）
- **RNA（L1, v1）数据集**：Release `rna-l1-v1`（包含 `.tsv.gz` + `manifest.json` + QA 报告）
  - Release: https://github.com/hazelian0619/protian-entity/releases/tag/rna-l1-v1
  - RNA 使用说明：`pipelines/rna/README.md`
  - RNA 规范：`docs/rna/README.md`

---

## 🧬 Protein 实体（L1）

构建以蛋白质为中心的高质量数据集，整合 UniProt、AlphaFold、HGNC、STRING 等多源数据。

### 📊 数据概览

| 项目 | 数量/覆盖率 | 说明 |
|------|------------|------|
| **蛋白质总数** | 19,135 | 去重后的人类蛋白质 |
| **字段数** | 33 | 完整信息字段 |
| **基因ID映射** | 99.6% | NCBI + Ensembl |
| **AlphaFold结构** | 99.7% | 含质量评分 |
| **功能注释** | 86% | 含证据代码+文献 |
| **GO注释** | 82-94% | 三个维度 |
| **PDB实验结构** | 44.3% | 实验解析结构 |

**主数据文件**：`data/processed/protein_master_v6_clean.tsv` (60MB, 19,135 行 × 33 列)

---

### ✅ 核心字段（33列）

#### 基础信息
`uniprot_id` | `protein_name` | `gene_names` | `sequence` | `mass`

#### 功能注释
`function` | `subcellular_location` | `diseases` | `ptms`

#### GO注释
`go_biological_process` | `go_molecular_function` | `go_cellular_component`

#### 基因ID
`ncbi_gene_id` | `ensembl_gene_id` | `hgnc_id` | `symbol` | `gene_synonyms`

#### 结构信息
`alphafold_id` | `alphafold_mean_plddt` | `pdb_ids` | `domains`

#### 交互数据
`string_ids` | `keywords`

---

### 🚀 快速使用

```python
import pandas as pd

df = pd.read_csv('data/processed/protein_master_v6_clean.tsv', sep='\t')

tp53 = df[df['gene_names'].str.contains('TP53', na=False)]
print(tp53[['uniprot_id', 'ncbi_gene_id', 'alphafold_mean_plddt']])
```

---

### 📁 辅助数据

```
data/processed/
├── alphafold_quality.tsv       # AlphaFold 每残基质量
├── protein_edges.tsv           # STRING 交互网络（约 88 万条）
├── ptm_sites.tsv               # 翻译后修饰（约 23 万条）
└── pathway_members.tsv         # 通路成员（约 12 万条）
```

---

### 🎯 设计原则

- ✅ **一级信息为主**：提取原始数据，不做推断
- ✅ **以蛋白质为主体**：每个 UniProt ID 一行
- ✅ **保留完整原文**：功能描述含证据代码和 PubMed 引用
- ✅ **多源整合**：7 个主要生物数据库

---

### 📋 项目状态

**阶段**：✅ 完成  \
**版本**：v6_clean  \
**更新**：2025-10-27

**数据源**：UniProt | AlphaFold | HGNC | STRING | GO | PDB  \
**时效性**：截止 2025-10-26

---

## 🧬 RNA 实体（L1, v1）

RNA（miRNA + mRNA/transcript）属于 L1 实体表；输出体积较大（单文件 >100MB），所以：

- 数据产物不直接 commit 进仓库（避免触发 GitHub 单文件限制、避免仓库膨胀）
- 统一通过 **GitHub Releases** 发布，并附带可核验的 `manifest.json` 与 QA 报告

- Release: https://github.com/hazelian0619/protian-entity/releases/tag/rna-l1-v1
- 使用说明：`pipelines/rna/README.md`
- 规范：`docs/rna/README.md`

---

## 💊 小分子（Molecules）实体（L1）与 PSI（L2, v1）

小分子按同样的“工业级数据产品”方式交付：

- 代码/规范/QA/manifest 进仓库
- 大体积数据（SQLite）通过 GitHub Releases 发布

入口：

- 使用说明：`pipelines/molecules/README.md`
- 规范：`docs/molecules/README.md`

Releases：

- Molecules L1（M1+M2）：https://github.com/hazelian0619/protian-entity/releases/tag/molecules-l1-v1
- Molecules PSI L2（M3）：https://github.com/hazelian0619/protian-entity/releases/tag/molecules-psi-l2-v1
