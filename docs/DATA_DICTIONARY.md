# 数据字典（Data Dictionary）

## 概述

本文档详细描述了protein_master_v6_clean.tsv主表的所有字段，包括数据类型、来源、说明和空值情况。

**主表**：protein_master_v6_clean.tsv
**行数**：19,135条
**列数**：33列
**更新日期**：2025-10-26

---

## 字段详细说明

### 一、基础标识字段

| 字段名 | 数据类型 | 说明 | 数据来源 | 空值率 | 示例值 |
|--------|---------|------|---------|--------|--------|
| **uniprot_id** | VARCHAR(20) | UniProt Accession（主键），蛋白质唯一标识符 | UniProt | 0% | P04637 |
| **entry_name** | VARCHAR(50) | UniProt Entry Name，格式为{PROTEIN}_{ORGANISM} | UniProt | 0% | TP53_HUMAN |
| **protein_name** | VARCHAR(500) | 推荐蛋白质名称（全称） | UniProt | 0% | Cellular tumor antigen p53 |
| **gene_names** | TEXT | 基因名称（可能包含多个别名，空格分隔） | UniProt | 0% | TP53 P53 |
| **symbol** | VARCHAR(50) | HGNC官方基因符号（权威命名） | HGNC | 0% | TP53 |
| **hgnc_id** | VARCHAR(20) | HGNC ID（格式：HGNC:数字） | HGNC | 0% | HGNC:11998 |

**说明**：
- `uniprot_id`是主键，保证全局唯一
- `symbol`使用HGNC官方命名，优先级高于`gene_names`
- 查询时建议使用`symbol`字段

---

### 二、序列信息字段

| 字段名 | 数据类型 | 说明 | 数据来源 | 空值率 | 示例值 |
|--------|---------|------|---------|--------|--------|
| **sequence** | TEXT | 完整氨基酸序列（单字母代码，无FASTA头） | UniProt | 0% | MEEPQSDPSVEPPLSQ... |
| **sequence_len** | INTEGER | 序列长度（氨基酸残基数） | UniProt | 0% | 393 |
| **mass** | INTEGER | 分子量（单位：道尔顿 Da） | UniProt | 0% | 43653 |

**说明**：
- `sequence`仅包含20种标准氨基酸字母（ACDEFGHIKLMNPQRSTVWY）
- `sequence_len` = len(sequence)，可用于校验
- `mass`为理论分子量，不含翻译后修饰

---

### 三、基因关联字段

| 字段名 | 数据类型 | 说明 | 数据来源 | 空值率 | 示例值 |
|--------|---------|------|---------|--------|--------|
| **ncbi_gene_id** | VARCHAR(20) | NCBI Gene ID（Entrez Gene ID） | HGNC | ~5% | 7157 |
| **ensembl_gene_id** | VARCHAR(20) | Ensembl基因ID（格式：ENSG开头） | HGNC | ~3% | ENSG00000141510 |
| **ensembl_transcript_id** | VARCHAR(20) | Ensembl转录本ID（格式：ENST开头） | HGNC | ~10% | ENST00000269305 |
| **gene_synonyms** | TEXT | 基因同义词列表（JSON数组格式） | HGNC | ~20% | ["P53", "TRP53", "BCC7"] |

**说明**：
- 这些字段在v6版本中新增，用于多数据库交叉引用
- `gene_synonyms`可用于模糊查询和别名匹配
- 部分空值是由于HGNC映射不完整

---

### 四、功能注释字段

| 字段名 | 数据类型 | 说明 | 数据来源 | 空值率 | 示例值 |
|--------|---------|------|---------|--------|--------|
| **function** | TEXT | 功能描述原文（完整保留，可能包含文献引用） | UniProt | ~14% | Acts as a tumor suppressor... {ECO:0000269\|PubMed:123456} |
| **go_biological_process** | TEXT | GO生物过程注释（多个term用分号分隔） | Gene Ontology | ~13% | apoptotic process [GO:0006915]; cell cycle arrest [GO:0007050] |
| **go_molecular_function** | TEXT | GO分子功能注释 | Gene Ontology | ~18% | DNA binding [GO:0003677]; protein binding [GO:0005515] |
| **go_cellular_component** | TEXT | GO细胞组分注释 | Gene Ontology | ~6% | nucleus [GO:0005634]; cytoplasm [GO:0005737] |

**说明**：
- `function`字段保留UniProt原文，包括证据代码（ECO）和PubMed引用
- GO注释字段格式为：term名称 [GO:ID]，多个term用分号分隔
- 空值可通过补充GO数据库或文本挖掘填补

---

### 五、定位与修饰字段

| 字段名 | 数据类型 | 说明 | 数据来源 | 空值率 | 示例值 |
|--------|---------|------|---------|--------|--------|
| **subcellular_location** | TEXT | 亚细胞定位描述原文 | UniProt | ~14% | Cytoplasm. Nucleus. Endoplasmic reticulum. |
| **ptms** | TEXT | 翻译后修饰摘要（概览） | UniProt | ~72% | Phosphorylated; Ubiquitinated; Acetylated |
| **diseases** | TEXT | 相关疾病描述原文 | UniProt | ~73% | Li-Fraumeni syndrome (LFS) [MIM:151623]... |

**说明**：
- `subcellular_location`描述蛋白质在细胞中的定位
- `ptms`为修饰概览，详细信息见`ptm_sites.tsv`辅助表
- `diseases`空值率高是正常现象（多数蛋白无已知疾病关联）

---

### 六、结构域与异构体字段

| 字段名 | 数据类型 | 说明 | 数据来源 | 空值率 | 示例值 |
|--------|---------|------|---------|--------|--------|
| **isoforms** | TEXT | 异构体信息（未展开，描述性文本） | UniProt | ~1% | Named isoforms=3; Name=1; IsoId=P04637-1... |
| **domains** | TEXT | 结构域信息（多个结构域用分号分隔） | UniProt | ~57% | DNA-binding domain 102-292; Tetramerization domain 320-355 |

**说明**：
- `isoforms`字段包含异构体描述，但未展开为独立行（待改进）
- `domains`空值率高是因为部分蛋白缺少明确的结构域注释
- 详细结构域信息可补充InterPro数据库

---

### 七、结构信息字段

| 字段名 | 数据类型 | 说明 | 数据来源 | 空值率 | 示例值 |
|--------|---------|------|---------|--------|--------|
| **pdb_ids** | TEXT | PDB结构ID列表（多个ID用逗号或分号分隔） | UniProt | ~56% | 1TUP,1TSR,2OCJ |
| **has_alphafold** | BOOLEAN | 是否有AlphaFold预测结构 | AlphaFold DB | 0% | TRUE / FALSE |
| **alphafold_id** | VARCHAR(50) | AlphaFold条目ID（格式：AF-{uniprot_id}-F1） | AlphaFold DB | ~0.3% | AF-P04637-F1 |
| **alphafold_entry_url** | VARCHAR(500) | AlphaFold网页链接 | AlphaFold DB | ~0.3% | https://alphafold.ebi.ac.uk/entry/P04637 |
| **alphafold_pdb_url** | VARCHAR(500) | AlphaFold结构文件下载链接 | AlphaFold DB | ~0.3% | https://alphafold.ebi.ac.uk/.../AF-P04637-F1-model_v6.pdb |
| **alphafold_mean_plddt** | FLOAT | AlphaFold平均pLDDT置信度（0-100） | AlphaFold DB | ~0.3% | 85.3 |

**说明**：
- `pdb_ids`为实验解析的三维结构（X-ray/NMR/Cryo-EM）
- AlphaFold覆盖率>99%，可作为PDB的补充
- `alphafold_mean_plddt` > 70为高置信度，> 90为非常高置信度
- 详细AlphaFold质量分析见`alphafold_quality.tsv`

---

### 八、交互与关键词字段

| 字段名 | 数据类型 | 说明 | 数据来源 | 空值率 | 示例值 |
|--------|---------|------|---------|--------|--------|
| **string_ids** | TEXT | STRING数据库ID（可能包含多个，逗号分隔） | UniProt交叉引用 | ~3% | 9606.ENSP00000269305 |
| **keywords** | TEXT | UniProt关键词（功能/定位/疾病等标签） | UniProt | ~2% | Acetylation; Activator; Apoptosis; DNA-binding |

**说明**：
- `string_ids`用于关联`protein_edges.tsv`交互网络数据
- `keywords`为分号分隔的关键词列表，便于快速分类和检索

---

### 九、血缘追踪字段

| 字段名 | 数据类型 | 说明 | 数据来源 | 空值率 | 示例值 |
|--------|---------|------|---------|--------|--------|
| **source** | VARCHAR(100) | 数据主要来源 | - | 0% | UniProt/Swiss-Prot |
| **date_modified** | DATE | UniProt数据最后修改日期 | UniProt | 0% | 2025-06-18 |
| **fetch_date** | DATE | 数据抓取日期 | - | 0% | 2025-10-26 |

**说明**：
- `source`标注数据来自Swiss-Prot（人工审核）或TrEMBL（自动注释）
- `date_modified`来自UniProt官方，反映源数据更新时间
- `fetch_date`为本地ETL处理时间

---

## 数据类型说明

### 文本字段处理
- **TEXT类型**：无长度限制，保留完整原文（包括Unicode字符、换行符）
- **分隔符**：多个值通常使用分号（;）分隔，部分字段使用逗号（,）
- **引用格式**：文献引用格式为`{ECO:证据代码|PubMed:ID}`或`[MIM:疾病ID]`

### 空值表示
- 空字符串：`""`
- pandas读取后显示为：`NaN`
- 建议使用`pd.isna()`或`.notna()`判断

### JSON格式字段
- `gene_synonyms`：`["P53", "TRP53"]`（JSON数组）
- 使用`json.loads()`或`ast.literal_eval()`解析

---

## 字段优先级建议

### 蛋白质查询
1. **首选**：`symbol`（HGNC官方）
2. **备选**：`uniprot_id`, `gene_names`, `gene_synonyms`

### 功能注释
1. **首选**：`function`（最详细）
2. **补充**：`go_biological_process`, `go_molecular_function`
3. **关键词筛选**：`keywords`

### 结构信息
1. **实验结构**：`pdb_ids`（优先）
2. **预测结构**：`alphafold_pdb_url`（覆盖率高）
3. **质量判断**：`alphafold_mean_plddt` > 70

---

## 字段关联关系

### 外键关系
- `uniprot_id` → `protein_edges.tsv` (source_uniprot / target_uniprot)
- `uniprot_id` → `ptm_sites.tsv` (uniprot_id)
- `uniprot_id` → `pathway_members.tsv` (uniprot_id)
- `uniprot_id` → `alphafold_quality.tsv` (uniprot_id)

### ID映射
- `symbol` ↔ `hgnc_core.tsv` (symbol)
- `ncbi_gene_id` ↔ NCBI Gene数据库
- `ensembl_gene_id` ↔ Ensembl数据库
- `string_ids` ↔ STRING数据库

---

## 数据质量指标

| 指标 | 值 | 说明 |
|------|-----|------|
| **必填字段完整率** | 100% | uniprot_id, sequence, symbol全部非空 |
| **核心功能字段** | 86% | function字段覆盖率 |
| **结构覆盖率** | 99%+ | PDB + AlphaFold综合覆盖 |
| **GO注释覆盖率** | 82-87% | 三个维度平均 |
| **ID唯一性** | 100% | uniprot_id无重复 |

---

## 使用示例

### 示例1：查询特定蛋白
```python
import pandas as pd

df = pd.read_csv('protein_master_v6_clean.tsv', sep='\t')

# 按symbol查询
tp53 = df[df['symbol'] == 'TP53'].iloc[0]
print(f"UniProt ID: {tp53['uniprot_id']}")
print(f"Function: {tp53['function'][:200]}...")
print(f"PDB: {tp53['pdb_ids']}")
print(f"AlphaFold pLDDT: {tp53['alphafold_mean_plddt']}")
```

### 示例2：筛选高质量结构
```python
# 有实验结构的蛋白
exp_structure = df[df['pdb_ids'].notna()]

# 高置信度AlphaFold结构
high_conf_af = df[(df['alphafold_mean_plddt'] > 80) & (df['pdb_ids'].isna())]

print(f"Experimental structures: {len(exp_structure)}")
print(f"High-confidence AlphaFold: {len(high_conf_af)}")
```

### 示例3：功能关键词筛选
```python
# 查找转录因子
tf = df[df['keywords'].str.contains('Transcription', case=False, na=False)]

# 查找DNA结合蛋白
dna_binding = df[df['go_molecular_function'].str.contains('DNA binding', case=False, na=False)]

print(f"Transcription factors: {len(tf)}")
print(f"DNA-binding proteins: {len(dna_binding)}")
```

### 示例4：空值处理
```python
# 填充空的function字段（使用GO注释）
df['function_filled'] = df['function'].fillna(df['go_biological_process'])

# 统计空值
print(df.isnull().sum())
```

---

## 数据更新说明

**当前版本**：v6 (2025-10-26)

**主要改进**：
- ✅ 新增4个基因关联字段（ncbi_gene_id, ensembl_gene_id等）
- ✅ 完整保留功能描述原文（TEXT类型）
- ✅ AlphaFold v6数据更新

**待改进项**：
- ⏳ 异构体展开（将isoforms字段拆分为独立行）
- ⏳ 补充空值功能字段（通过GO/InterPro）
- ⏳ 添加audit_status质量审计字段

详见项目文档：`/docs/Protein_Entity_Completeness_Audit.md`

---

**文档版本**：v1.0
**最后更新**：2025-10-26
**维护者**：项目团队
