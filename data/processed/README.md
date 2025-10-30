# 1025 项目 - 处理后数据说明

## 概述

本目录包含人类蛋白质知识图谱的核心数据集，整合了多个权威生物信息学数据库的信息。

**数据更新时间**：2025-10-26
**数据版本**：v6
**物种**：Homo sapiens (人类，Taxonomy ID: 9606)

---

## 主表文件

### protein_master_v6_clean.tsv
**人类蛋白质核心信息表（推荐使用）**

- **行数**：19,135条蛋白质记录
- **列数**：33个字段
- **文件大小**：约58 MB
- **编码**：UTF-8
- **分隔符**：制表符（\t）

**数据来源**：
- UniProt (Swiss-Prot + TrEMBL) - 核心序列、功能、定位、疾病
- HGNC - 官方基因命名与ID映射
- AlphaFold Database v6 - 结构预测与置信度
- Gene Ontology (GO) - 功能注释（BP/MF/CC）
- STRING v12.0 - 交互网络ID
- PhosphoSite - 修饰位点摘要
- Reactome - 通路信息

**核心字段**：
- 基础标识：`uniprot_id`, `entry_name`, `protein_name`, `symbol`, `hgnc_id`
- 序列信息：`sequence`, `sequence_len`, `mass`
- 基因关联：`ncbi_gene_id`, `ensembl_gene_id`, `ensembl_transcript_id`, `gene_synonyms`
- 功能注释：`function`, `go_biological_process`, `go_molecular_function`, `go_cellular_component`
- 结构信息：`pdb_ids`, `alphafold_pdb_url`, `alphafold_mean_plddt`
- 其他信息：`subcellular_location`, `ptms`, `diseases`, `domains`, `isoforms`, `keywords`
- 血缘追踪：`source`, `date_modified`, `fetch_date`

**数据完整性**：
- 必填字段（uniprot_id/sequence/symbol）：100%
- function：约86%
- subcellular_location：约86%
- pdb_ids：约44%（其余由AlphaFold覆盖）
- AlphaFold覆盖率：约99.7%

**使用示例**：
```python
import pandas as pd

# 读取主表
df = pd.read_csv('protein_master_v6_clean.tsv', sep='\t')
print(f"Total proteins: {len(df)}")

# 查询特定蛋白
tp53 = df[df['symbol'] == 'TP53']
print(tp53[['uniprot_id', 'protein_name', 'function']])

# 筛选有PDB结构的蛋白
pdb_proteins = df[df['pdb_ids'].notna()]
print(f"Proteins with PDB structure: {len(pdb_proteins)}")
```

---

## 辅助表文件

### hgnc_core.tsv
**HGNC基因命名核心数据**

- **行数**：19,278个基因
- **数据源**：HGNC (2025-10版)
- **用途**：提供权威基因命名、别名、ID映射（NCBI/Ensembl/UniProt）

**关键字段**：HGNC ID、symbol、基因名称、locus类型、Ensembl/NCBI/UniProt交叉引用、别名

---

### protein_edges.tsv
**蛋白质交互网络（STRING）**

- **行数**：884,555条交互边
- **数据源**：STRING v12.0
- **更新日期**：2025-10-26

**字段说明**：
- `source_uniprot`：源蛋白UniProt ID
- `target_uniprot`：目标蛋白UniProt ID
- `combined_score`：交互置信分数（0-1000，通常>400为可信）
- `source`：STRING_v12.0
- `fetch_date`：数据获取日期

**网络规模**：覆盖约19,000个蛋白节点，平均每蛋白约46条边

**使用示例**：
```python
import pandas as pd
import networkx as nx

edges = pd.read_csv('protein_edges.tsv', sep='\t')
high_conf = edges[edges['combined_score'] > 700]

G = nx.from_pandas_edgelist(high_conf,
                             source='source_uniprot',
                             target='target_uniprot',
                             edge_attr='combined_score')
print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
```

---

### ptm_sites.tsv
**翻译后修饰位点（PhosphoSite）**

- **行数**：235,502个修饰位点
- **数据源**：PhosphoSitePlus
- **更新日期**：2025-10-26

**字段说明**：`uniprot_id`, `gene_name`, `ptm_type`, `residue`, `position`, `site_sequence`, `domain`, `lit_references`, `ms_references`

**覆盖情况**：磷酸化位点约80%，其他修饰包括乙酰化、泛素化、甲基化等

---

### pathway_members.tsv
**通路成员关系（Reactome）**

- **行数**：127,096条蛋白-通路关系
- **数据源**：Reactome
- **更新日期**：2025-10-26

**字段说明**：`uniprot_id`, `pathway_id`, `pathway_name`, `species`, `source`, `fetch_date`

**用途**：蛋白质通路注释、功能富集分析、免疫/代谢通路查询

---

### alphafold_quality.tsv
**AlphaFold结构质量详情**

- **行数**：23,586条记录
- **数据源**：AlphaFold Database v6
- **更新日期**：2025-10-26

**字段说明**：pLDDT统计（mean/median/min/max）、置信度分布、质量等级、残基数、版本

**质量分级标准**：
- Very high (>90 pLDDT)：主链和侧链均高度准确
- High (70-90 pLDDT)：主链准确，侧链可能偏差
- Low (50-70 pLDDT)：低置信度，需谨慎使用
- Very low (<50 pLDDT)：无序区域，不可靠

---

### uniprot_seed.json
**HGNC→UniProt种子映射**

- **记录数**：约19,278个映射
- **格式**：JSON数组
- **用途**：ETL起点，确保覆盖官方命名的蛋白质编码基因

---

## 数据质量说明

### 整体评估
- **当前评分**：⭐⭐⭐⭐½ (4.5/5)
- **覆盖范围**：主要人类蛋白质（约占Swiss-Prot的95%）
- **数据新鲜度**：2025年最新数据
- **血缘可追溯**：所有数据标注来源和获取时间

### 已知限制
1. **异构体未完全展开**：当前每个基因1条记录，未展开为所有异构体变体（目标约100k条）
2. **部分字段空值**：function约14%空，subcellular_location约14%空（可通过GO补充）
3. **PDB结构覆盖率**：约44%（其余可通过AlphaFold补充）

### 改进计划
详见 `/docs/Data_Quality_Assessment.md` 和 `/docs/Protein_Entity_Completeness_Audit.md`

---

## 数据使用注意事项

1. **字符编码**：所有TSV文件使用UTF-8编码
2. **分隔符**：制表符（\t），部分字段内容可能包含逗号
3. **空值表示**：空字符串或NaN
4. **ID唯一性**：uniprot_id为主键，保证唯一
5. **列表字段**：部分字段（如gene_synonyms）以JSON列表形式存储

---

## 引用与致谢

本数据集整合了以下公开数据库：

- **UniProt**: The UniProt Consortium (2025). UniProt: the Universal Protein Knowledgebase in 2025.
- **HGNC**: HUGO Gene Nomenclature Committee at the European Bioinformatics Institute
- **AlphaFold**: Jumper, J. et al. (2021). Highly accurate protein structure prediction with AlphaFold. Nature.
- **Gene Ontology**: The Gene Ontology Consortium (2025). The Gene Ontology resource.
- **STRING**: Szklarczyk, D. et al. (2025). STRING v12: protein-protein association networks.
- **PhosphoSitePlus**: Hornbeck, P.V. et al. (2025). PhosphoSitePlus: mutations, PTMs and recalibrations.
- **Reactome**: Gillespie, M. et al. (2025). The Reactome Pathway Knowledgebase.

---

## 更新历史

- **v6 (2025-10-26)**：新增ncbi_gene_id、ensembl_gene_id、ensembl_transcript_id、gene_synonyms字段
- **v5 (2025-10-26)**：整合v3-v4所有字段
- **v4 (2025-10-26)**：增强疾病信息
- **v3 (2025-10-26)**：基础版（UniProt + HGNC + AlphaFold + GO）

---

**文档版本**：v1.0
**最后更新**：2025-10-26
