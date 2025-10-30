# 数据质量评估与完善建议（2025-10-26）

## 📊 当前数据概览

### 已有数据表
| 文件名 | 行数 | 大小 | 用途 |
|--------|------|------|------|
| **protein_master_v2.tsv** | 19,253 | 58 MB | 主表（蛋白质实体） |
| **protein_edges.tsv** | 884,555 | 37 MB | 蛋白质交互边（STRING） |
| **ptm_sites.tsv** | 235,502 | 20 MB | 翻译后修饰位点 |
| **pathway_members.tsv** | 127,096 | 11 MB | 通路成员 |
| **hgnc_core.tsv** | 19,278 | 9.5 MB | HGNC核心数据 |
| **uniprot_seed.json** | - | 1.7 MB | UniProt种子数据 |

---

## ✅ 当前优点

### 1. 数据量合理
- **主表19,253条**：覆盖了主要的人类蛋白质（Swiss-Prot ~20k）
- **边数据88万条**：STRING交互网络规模适中
- **PTM位点23万条**：修饰数据相当完整

### 2. 字段质量较好
抽样检查发现：
- ✅ **function字段**：保留完整原文，包含文献引用（如 ECO:0000269|PubMed:xxx）
- ✅ **subcellular_location字段**：详细的定位信息，包含注释
- ✅ **序列信息**：完整的氨基酸序列
- ✅ **GO注释**：三个维度齐全（BP/MF/CC）

### 3. 多源整合
- ✅ UniProt主干数据
- ✅ STRING交互网络
- ✅ PTM修饰数据
- ✅ HGNC基因命名
- ✅ 通路信息

### 4. 时间戳管理
- ✅ `date_modified`: 2025-06-18（UniProt数据日期）
- ✅ `fetch_date`: 2025-10-26（抓取日期）

---

## ⚠️ 关键问题与差距

### 问题1：数据量远低于目标（严重）

**现状**：
- 当前：19,253条蛋白质记录
- 目标：~100,000条蛋白质变体

**原因分析**：
- ❌ **异构体未完全展开**：
  - 从抽样看，`isoforms`字段有内容（如"Named isoforms=3"），但这些异构体并未展开为独立行
  - 例如：A6NFL2有2个异构体，但只有1行记录
  - 导致数据丢失：~80%的蛋白质变体未录入

**影响**：
- 违背核心原则："以蛋白质为主体，每行一变体"
- 与1017项目犯了同样的错误（虽然程度较轻）

**预期改进后数据量**：
- Swiss-Prot: ~20k主蛋白 × 平均3个异构体 = **~60k条**
- TrEMBL: ~20k补充 × 平均2个异构体 = **~40k条**
- 合计：**~100k条**（符合目标）

---

### 问题2：字段结构与Schema不一致（中等）

**缺失关键字段**（对比 Unified_Protein_Schema_v2.md）：

#### 基本标识字段缺失
| 缺失字段 | 重要性 | 用途 |
|---------|--------|------|
| `entry_type` | ⭐⭐⭐⭐⭐ | 区分Swiss-Prot（审核）/TrEMBL（自动） |
| `is_canonical` | ⭐⭐⭐⭐⭐ | 标记代表性异构体（关键！） |
| `isoform_name` | ⭐⭐⭐⭐ | 异构体具体名称（如"Isoform 2"） |
| `protein_name_alternative` | ⭐⭐⭐ | 其他蛋白质名称（JSON数组） |

#### 基因关联字段缺失
| 缺失字段 | 重要性 | 用途 |
|---------|--------|------|
| `ncbi_gene_id` | ⭐⭐⭐⭐⭐ | NCBI Gene ID（外键） |
| `ensembl_gene_id` | ⭐⭐⭐⭐ | Ensembl基因ID |
| `ensembl_transcript_id` | ⭐⭐⭐⭐ | 对应转录本ID（异构体来源） |

#### 结构信息字段缺失
| 缺失字段 | 重要性 | 用途 |
|---------|--------|------|
| `structure_available` | ⭐⭐⭐⭐⭐ | 结构类型标识（PDB/AlphaFold/None） |
| `pdb_best_resolution` | ⭐⭐⭐⭐ | 最佳分辨率（Å） |
| `alphafold_model_url` | ⭐⭐⭐⭐ | AlphaFold预测结构下载链接 |
| `alphafold_mean_plddt` | ⭐⭐⭐⭐⭐ | 平均置信度（判断结构可靠性） |
| `alphafold_version` | ⭐⭐⭐ | AlphaFold版本 |

#### 功能注释字段缺失
| 缺失字段 | 重要性 | 用途 |
|---------|--------|------|
| `catalytic_activity` | ⭐⭐⭐⭐ | 催化活性（独立字段） |
| `cofactor` | ⭐⭐⭐ | 辅因子信息 |
| `activity_regulation` | ⭐⭐⭐ | 活性调控机制 |
| `subunit_structure` | ⭐⭐⭐ | 亚基组成 |
| `pathway` | ⭐⭐⭐⭐ | 代谢通路（独立字段） |

#### 数据血缘字段缺失
| 缺失字段 | 重要性 | 用途 |
|---------|--------|------|
| `primary_source` | ⭐⭐⭐ | 主要数据源 |
| `uniprot_version` | ⭐⭐⭐⭐ | UniProt版本号 |
| `audit_status` | ⭐⭐⭐⭐⭐ | 质量状态（OK/MULTI/MISSING） |
| `audit_comments` | ⭐⭐⭐⭐ | 补漏日志 |

**字段命名不一致**：
| 当前字段名 | Schema设计字段名 | 建议 |
|-----------|----------------|------|
| `protein_name` | `protein_name_recommended` | 保持一致 |
| `gene_names` | `gene_synonyms` | 拆分为`hgnc_symbol`和`gene_synonyms` |
| `mass` | `sequence_mass` | 改名 |

---

### 问题3：辅助表不完整（中等）

**当前有**：
- ✅ `protein_edges.tsv`（交互）
- ✅ `ptm_sites.tsv`（PTM修饰）

**缺失**（对比Schema设计）：
- ❌ `protein_aliases.tsv`：别名映射表
- ❌ `protein_features.tsv`：结构特征表（domain/region/binding site）
- ❌ `pdb_structures.tsv`：PDB结构详情表
- ❌ `disease_variants.tsv`：疾病变异表

**影响**：
- 别名检索不便（需要解析主表的`gene_names`字段）
- 结构特征信息未拆分（domain/region等混在一起）
- PDB详情缺失（分辨率、配体、实验方法等）
- 疾病变异未结构化存储

---

### 问题4：数据空值问题（轻微-中等）

**抽样检查发现**：
| 字段 | 空值情况 | 影响 |
|------|---------|------|
| `diseases` | ~70%为空 | 可接受（多数蛋白无疾病关联） |
| `isoforms` | ~30%为空 | 正常（单一异构体） |
| `pdb_ids` | ~70%为空 | **需补AlphaFold结构** |
| `domains` | ~20%为空 | 可用InterPro补充 |

**缺少结构信息的蛋白质**：
- 当前`pdb_ids`空值~70%
- 但AlphaFold DB覆盖率>95%
- **建议**：补充`alphafold_model_url`和`alphafold_mean_plddt`字段

---

### 问题5：时效性问题（轻微）

**当前数据日期**：
- `date_modified`: 2025-06-18（UniProt数据）
- `fetch_date`: 2025-10-26（抓取日期）

**目标**：≤2025-10-31

**建议**：
- ✅ UniProt数据较新（2025-06-18），可接受
- ⚠️ 但缺少`uniprot_version`字段（如"2025_04"）
- ⚠️ 建议补充其他数据源的版本信息：
  - PDB版本：2025-10周
  - AlphaFold版本：v4
  - GO版本：2025-10-01
  - STRING版本：v12.5

---

## 📋 完善建议清单（按优先级排序）

### 🔴 高优先级（必须）

#### 1. 异构体展开（最关键！）
**任务**：将`isoforms`字段中的异构体信息展开为独立行

**实施步骤**：
```python
# 伪代码
for protein in protein_master:
    if protein['isoforms']:
        # 解析isoforms字段：'Named isoforms=3; Name=1; IsoId=...'
        isoform_list = parse_isoforms(protein['isoforms'])
        for isoform in isoform_list:
            new_row = protein.copy()
            new_row['uniprot_id'] = isoform['IsoId']  # 如 P04637-2
            new_row['isoform_name'] = isoform['Name']  # 如 "Isoform 2"
            new_row['is_canonical'] = (isoform['Sequence'] == 'Displayed')
            # 从UniProt API获取该异构体的具体序列
            new_row['sequence'] = fetch_isoform_sequence(isoform['IsoId'])
            output_row(new_row)
```

**预期结果**：
- 行数从19,253 → ~100,000
- 每个异构体独立成行
- `uniprot_id`包含异构体后缀（如P04637-2）

**数据源**：
- UniProt API: `/uniprotkb/{accession}/isoforms`

---

#### 2. 补充结构信息（AlphaFold）
**任务**：为缺少PDB结构的蛋白质补充AlphaFold预测结构

**新增字段**：
- `structure_available`: ENUM('PDB', 'AlphaFold', 'None')
- `alphafold_model_url`: AlphaFold模型下载链接
- `alphafold_mean_plddt`: 平均置信度（0-100）
- `alphafold_version`: v4
- `pdb_best_resolution`: 最佳分辨率（从现有PDB提取）

**实施步骤**：
```python
for protein in protein_master:
    if protein['pdb_ids']:
        protein['structure_available'] = 'PDB'
        # 调用PDB API获取分辨率
        protein['pdb_best_resolution'] = get_best_resolution(protein['pdb_ids'])
    else:
        # 调用AlphaFold API
        af_data = fetch_alphafold(protein['uniprot_id'])
        if af_data:
            protein['structure_available'] = 'AlphaFold'
            protein['alphafold_model_url'] = af_data['pdbUrl']
            protein['alphafold_mean_plddt'] = calculate_mean_plddt(af_data)
        else:
            protein['structure_available'] = 'None'
```

**数据源**：
- AlphaFold DB API: `https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}`

**预期结果**：
- 结构覆盖率：30% (PDB) + 65% (AlphaFold) = **95%**

---

#### 3. 补充基因关联字段
**任务**：添加NCBI Gene ID和Ensembl ID

**新增字段**：
- `ncbi_gene_id`
- `ensembl_gene_id`
- `ensembl_transcript_id`

**数据源**：
- UniProt交叉引用：`uniProtKBCrossReferences.GeneID`
- UniProt交叉引用：`uniProtKBCrossReferences.Ensembl`
- 或使用HGNC映射表

**实施步骤**：
```python
# 方案1：从UniProt API获取
uniprot_data = fetch_uniprot(uniprot_id)
ncbi_gene_id = uniprot_data['uniProtKBCrossReferences']['GeneID']
ensembl_id = uniprot_data['uniProtKBCrossReferences']['Ensembl']

# 方案2：从HGNC映射
hgnc_data = hgnc_core[hgnc_core['symbol'] == gene_symbol]
ncbi_gene_id = hgnc_data['entrez_id']
ensembl_id = hgnc_data['ensembl_gene_id']
```

---

#### 4. 添加质量审计字段
**任务**：记录数据质量状态和补漏日志

**新增字段**：
- `audit_status`: ENUM('OK', 'MULTI', 'MISSING', 'PARTIAL')
- `audit_comments`: TEXT（记录补漏操作）

**审计规则**：
```python
audit_status = 'OK'
audit_comments = []

# 检查必填字段
if not sequence or not function:
    audit_status = 'MISSING'
    audit_comments.append('缺少核心字段')

# 检查异构体
if '|' in uniprot_id:  # 多个ID未拆分
    audit_status = 'MULTI'
    audit_comments.append('异构体需展开')

# 检查结构信息
if not pdb_ids and not alphafold_model_url:
    audit_status = 'PARTIAL'
    audit_comments.append('缺少结构信息')

# 记录补漏操作
if function_补充了GO:
    audit_comments.append('补充GO功能注释')
```

---

### 🟡 中优先级（重要）

#### 5. 创建辅助表

**protein_aliases.tsv**（别名表）：
```
alias_id | uniprot_id | alias_type | alias_value | source | fetch_date
1        | P04637     | gene_alt   | P53         | HGNC   | 2025-10-26
2        | P04637     | gene_alt   | TRP53       | HGNC   | 2025-10-26
```

**protein_features.tsv**（特征表）：
```
feature_id | uniprot_id | feature_type | description | start_pos | end_pos | evidence | source
1          | P04637     | Domain       | DNA-binding | 102       | 292     | IDA      | UniProt
```

**pdb_structures.tsv**（PDB详情表）：
```
pdb_id | uniprot_id | title | method | resolution | r_factor | ligands | release_date
1TUP   | P04637     | p53.. | X-RAY  | 2.10       | 0.195    | [ZN]    | 1994-10-31
```

**disease_variants.tsv**（变异表）：
```
variant_id | uniprot_id | position | original_aa | variant_aa | disease | dbsnp_id | significance
1          | P04637     | 175      | R           | H          | LFS     | rs28934576 | Pathogenic
```

---

#### 6. 字段拆分与规范化

**拆分`gene_names`字段**：
```python
# 当前：gene_names = "CYP2D7"
# 改为：
hgnc_symbol = "CYP2D7"  # 主字段
gene_synonyms = []  # 从HGNC获取同义词
```

**拆分`isoforms`字段**：
- 当前混合在一个TEXT字段中
- 改为结构化存储在独立行中

**规范化PTM字段**：
- 当前`ptms`字段为TEXT
- 已有独立的`ptm_sites.tsv`表（很好）
- 主表保留概览，详情在辅助表

---

#### 7. 补充详细功能注释字段

**从UniProt API提取以下独立字段**：
- `catalytic_activity`: 催化活性
- `cofactor`: 辅因子
- `activity_regulation`: 活性调控
- `subunit_structure`: 亚基组成
- `pathway`: 代谢通路

**当前问题**：
- 这些信息可能混在`function`字段中
- 需要从UniProt API的`comments`中提取对应部分

---

### 🟢 低优先级（可选）

#### 8. 补充其他数据源

**PhosphoSite**：
- 补充PTM修饰细节（已有ptm_sites.tsv，可增强）

**Reactome**：
- 补充通路信息（已有pathway_members.tsv）

**Human Protein Atlas**：
- 补充组织表达数据

---

#### 9. 数据版本管理

**添加`source_version`字段**：
```python
source_version = "UniProt:2025_04, HGNC:2025-10, STRING:v12.5"
```

**记录各数据源的具体版本**：
- UniProt: 2025_04
- PDB: 2025-10
- AlphaFold: v4
- HGNC: 2025-10
- STRING: v12.5
- GO: 2025-10-01

---

#### 10. 索引优化

**为主表添加索引**（如果导入数据库）：
```sql
CREATE INDEX idx_uniprot_id ON protein_master(uniprot_id);
CREATE INDEX idx_hgnc_symbol ON protein_master(hgnc_symbol);
CREATE INDEX idx_entry_type ON protein_master(entry_type);
CREATE INDEX idx_is_canonical ON protein_master(is_canonical);
CREATE INDEX idx_structure_available ON protein_master(structure_available);
```

---

## 📈 完善后的预期结果

| 指标 | 当前 | 目标 | 改善 |
|------|------|------|------|
| **记录数** | 19,253 | ~100,000 | ✅ 5倍提升 |
| **字段数** | 24 | ~40 | ✅ 信息更全面 |
| **结构覆盖率** | ~30% (PDB) | ~95% (PDB+AlphaFold) | ✅ 3倍提升 |
| **异构体覆盖** | 0%（未展开） | 100% | ✅ 完全覆盖 |
| **辅助表数量** | 3 | 7 | ✅ 更结构化 |
| **数据完整率** | ~85% | >99% | ✅ 达标 |
| **audit追踪** | 无 | 有 | ✅ 可追溯 |

---

## 🎯 实施路线图

### 第一阶段（1-2天）：核心修复
1. ✅ 异构体展开脚本
2. ✅ AlphaFold结构补充
3. ✅ 基因ID映射

### 第二阶段（2-3天）：字段补全
4. ✅ 添加audit字段
5. ✅ 拆分gene_names
6. ✅ 补充详细功能字段

### 第三阶段（2-3天）：辅助表创建
7. ✅ protein_aliases表
8. ✅ protein_features表
9. ✅ pdb_structures表
10. ✅ disease_variants表

### 第四阶段（1天）：质量检查
11. ✅ 运行QA脚本
12. ✅ 修复空值
13. ✅ 生成质量报告

**总计**：约6-9天

---

## 📝 总结

**当前数据质量评分**：⭐⭐⭐⭐ (4/5)

**优点**：
- ✅ 数据量合理，覆盖主要蛋白质
- ✅ 字段内容完整，保留原文
- ✅ 多源整合较好
- ✅ 时间戳管理到位

**主要问题**：
- ❌ 异构体未展开（最严重）
- ❌ 字段结构与Schema设计不完全一致
- ❌ 缺少AlphaFold结构信息
- ❌ 辅助表不完整

**完善后评分预期**：⭐⭐⭐⭐⭐ (5/5)

---

**文档版本**：v1.0
**创建时间**：2025-10-26
**作者**：项目团队
**下次审核**：2025-10-27（完成第一阶段后）
