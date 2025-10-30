# 蛋白质实体信息完整性排查报告

## 📋 当前数据盘点

### 主表：protein_master_v3_final_fixed.tsv
**行数**：19,253条
**列数**：29列
**文件大小**：58 MB

### 辅助表
| 表名 | 行数 | 用途 | 状态 |
|------|------|------|------|
| alphafold_quality.tsv | ~19k | AlphaFold详细质量数据 | ✅ 已完成 |
| protein_edges.tsv | 884,555 | STRING交互网络 | ✅ 已完成 |
| ptm_sites.tsv | 235,502 | PTM修饰位点 | ✅ 已完成 |
| pathway_members.tsv | 127,096 | Reactome通路成员 | ✅ 已完成 |
| hgnc_core.tsv | 19,278 | HGNC基因命名 | ✅ 已完成 |

---

## 🔍 字段空值率分析（主表）

### ✅ 完整字段（0%空值）
核心标识、序列、时间戳、AlphaFold基础字段全部完整

### 🟢 低空值率字段（<5%）
| 字段 | 空值率 | 空值数 |
|------|--------|--------|
| alphafold_pdb_url | 0.34% | 66 |
| alphafold_mean_plddt | 0.34% | 66 |
| isoforms | 1.08% | 208 |
| string_ids | 2.61% | 502 |
| go_cellular_component | 5.80% | 1,117 |

### 🟡 中等空值率字段（10-20%）
| 字段 | 空值率 | 空值数 | ⚠️ 需补充 |
|------|--------|--------|----------|
| function | 13.98% | 2,691 | GO功能注释 |
| subcellular_location | 13.95% | 2,685 | GO CC |
| go_biological_process | 13.22% | 2,546 | GO API |
| go_molecular_function | 18.39% | 3,541 | GO API |

### 🔴 高空值率字段（>50%）
| 字段 | 空值率 | 空值数 | 正常性 |
|------|--------|--------|--------|
| pdb_ids | 55.73% | 10,729 | ✅ 正常 |
| domains | 57.09% | 10,992 | 可补充InterPro |
| ptms | 71.89% | 13,841 | ✅ 正常（有独立表） |
| diseases | 72.89% | 14,033 | ✅ 正常 |

---

## ⚠️ 缺失的关键字段（对比Schema v3.0）

### 🔴 高优先级缺失（必须补充）

#### A. 基因ID映射（4个字段）⭐⭐⭐⭐⭐
```
❌ ncbi_gene_id          - NCBI Gene ID
❌ ensembl_gene_id       - Ensembl基因ID
❌ ensembl_transcript_id - Ensembl转录本ID
❌ gene_synonyms         - 基因同义词（JSON）
```
**数据源**：hgnc_core.tsv（已有）或UniProt交叉引用
**工作量**：⭐ (低，1天)

#### B. 条目质量标识（2个字段）⭐⭐⭐⭐⭐
```
❌ entry_type     - Swiss-Prot vs TrEMBL
❌ is_canonical   - 是否代表性异构体
```
**数据源**：UniProt API或从entry_name推断
**工作量**：⭐ (低)

#### C. 结构信息补充（2个字段）⭐⭐⭐⭐⭐
```
❌ structure_available  - ENUM('PDB', 'AlphaFold', 'None')
❌ pdb_best_resolution  - 最佳PDB分辨率（Å）
```
**实施逻辑**：
```python
if pdb_ids:
    structure_available = "PDB"
    pdb_best_resolution = min(get_resolutions(pdb_ids))
elif has_alphafold and alphafold_mean_plddt > 70:
    structure_available = "AlphaFold"
else:
    structure_available = "None"
```
**工作量**：⭐⭐ (中，需调PDB API)

#### D. 数据审计（2个字段）⭐⭐⭐⭐⭐
```
❌ audit_status    - OK/MULTI/MISSING/PARTIAL
❌ audit_comments  - 补漏日志
```
**工作量**：⭐ (低，脚本生成)

---

### 🟡 中优先级缺失（推荐补充）

#### E. 详细功能注释（6个字段）⭐⭐⭐⭐
```
❌ catalytic_activity   - 催化活性
❌ cofactor             - 辅因子
❌ activity_regulation  - 活性调控
❌ subunit_structure    - 亚基组成
❌ pathway              - 代谢通路
❌ interaction_partners - 相互作用蛋白
```
**当前状态**：这些信息可能混在`function`字段中
**数据源**：UniProt API comments部分
**工作量**：⭐⭐⭐ (中高，2-3天)
**必要性**：如需结构化查询 → 补充；否则暂缓

#### F. 序列校验（3个字段）⭐⭐⭐
```
❌ sequence_crc64   - 序列校验和
❌ sequence_md5     - 序列MD5
❌ sequence_version - 序列版本号
```
**用途**：序列去重、版本追踪
**工作量**：⭐⭐ (中)

#### G. 定位补充（3个字段）⭐⭐⭐
```
❌ tissue_specificity  - 组织特异性
❌ developmental_stage - 发育阶段
❌ induction           - 诱导条件
```
**工作量**：⭐⭐⭐ (中高)

---

### 🟢 低优先级缺失（可选）

#### H. 蛋白质名称细化（2个字段）⭐⭐
```
❌ protein_name_alternative - 其他名称（JSON）
❌ protein_name_short       - 简称
```

#### I. 版本追踪（2个字段）⭐⭐⭐
```
❌ source_version   - 多源版本号
❌ alphafold_version - AlphaFold版本
```

---

## 🎯 补全计划（分阶段实施）

### Phase 1：基础映射补全（1天）⭐⭐⭐⭐⭐
**新增8个字段**：
1. `ncbi_gene_id` ← 从hgnc_core.tsv
2. `ensembl_gene_id` ← 从hgnc_core.tsv
3. `ensembl_transcript_id` ← 从UniProt/HGNC
4. `gene_synonyms` ← 从hgnc_core.tsv
5. `entry_type` ← 从UniProt推断
6. `is_canonical` ← 设为TRUE
7. `structure_available` ← 综合pdb_ids和has_alphafold
8. `pdb_best_resolution` ← 从PDB API

**实施步骤**：
```python
# 1. 读取HGNC映射
hgnc = pd.read_csv('hgnc_core.tsv', sep='\t')
mapping = hgnc.set_index('symbol')[['entrez_id', 'ensembl_gene_id', 'alias_symbol']].to_dict('index')

# 2. 读取主表
df = pd.read_csv('protein_master_v3_final_fixed.tsv', sep='\t')

# 3. 映射基因ID
df['ncbi_gene_id'] = df['symbol'].map(lambda x: mapping.get(x, {}).get('entrez_id'))
df['ensembl_gene_id'] = df['symbol'].map(lambda x: mapping.get(x, {}).get('ensembl_gene_id'))
df['gene_synonyms'] = df['symbol'].map(lambda x: mapping.get(x, {}).get('alias_symbol'))

# 4. 推断entry_type
df['entry_type'] = df['source'].apply(lambda x: 'Swiss-Prot' if 'reviewed' in x.lower() else 'TrEMBL')

# 5. 设置is_canonical
df['is_canonical'] = True  # 当前每行都是主要蛋白

# 6. 生成structure_available
def get_structure_flag(row):
    if pd.notna(row['pdb_ids']) and row['pdb_ids']:
        return 'PDB'
    elif row['has_alphafold'] and row['alphafold_mean_plddt'] > 70:
        return 'AlphaFold'
    else:
        return 'None'

df['structure_available'] = df.apply(get_structure_flag, axis=1)

# 7. 获取PDB分辨率（仅针对有pdb_ids的）
# 需调用PDB API...
```

**输出**：`protein_master_v4.tsv`（37列）

---

### Phase 2：质量审计（0.5天）⭐⭐⭐⭐⭐
**新增2个字段**：
9. `audit_status`
10. `audit_comments`

**逻辑**：
```python
def audit_record(row):
    status = "OK"
    comments = []

    if pd.isna(row['function']) or not row['function']:
        status = "PARTIAL"
        comments.append("缺少function")

    if pd.isna(row['go_biological_process']) and pd.isna(row['go_molecular_function']):
        comments.append("缺少GO注释")

    if row['structure_available'] == 'None':
        comments.append("无结构信息")

    return status, "; ".join(comments)

df[['audit_status', 'audit_comments']] = df.apply(audit_record, axis=1, result_type='expand')
```

---

### Phase 3：功能字段拆分（2-3天，按需）⭐⭐⭐⭐
**新增6个字段**：详细功能注释

**决策点**：
- ✅ 当前`function`字段已包含完整原文
- ❓ 是否需要结构化查询？
- 建议：暂缓，Phase 1/2完成后再评估

---

### Phase 4：序列与定位（1-2天，可选）⭐⭐⭐
**新增6个字段**：序列校验、组织表达等

---

### Phase 5：版本追踪（0.5天）⭐⭐⭐
**新增2个字段**：source_version、alphafold_version

---

## 📊 补全效果预测

| 指标 | 当前 | Phase 1 | Phase 2 | 全部完成 |
|------|------|---------|---------|----------|
| **字段数** | 29 | 37 | 39 | 51 |
| **基因ID映射** | 50% | **100%** | 100% | 100% |
| **结构标识** | 部分 | **完整** | 完整 | 完整 |
| **质量追溯** | 无 | 无 | **有** | 有 |
| **功能细粒度** | 混合 | 混合 | 混合 | 拆分 |

---

## 🚀 立即行动建议

### 最小可行方案（推荐）：
**完成 Phase 1 + Phase 2（1.5天）**

**收益**：
- 基因ID映射完整率：50% → 100%
- 新增结构类型标识字段
- 新增数据质量审计
- 字段数：29 → 39

**成本**：
- 工作量：1.5天
- 技术难度：低-中

---

## ❓ 需要你的决策

### 问题1：功能字段是否拆分？
**现状**：`function`字段包含完整原文（可能混合催化活性、辅因子等）

**选项A**：保持现状 ✅ 推荐
- 优点：简单，原文完整
- 缺点：无法结构化查询

**选项B**：拆分为6个字段
- 优点：结构化查询
- 缺点：+2-3天工作量

### 问题2：组织表达数据是否需要？
**Schema有**：tissue_specificity、developmental_stage
**当前**：缺失
**建议**：如研究不关注组织特异性 → 暂缓

### 问题3：异构体是否展开？
**当前**：19k条（1蛋白1行）
**Schema**：100k条（每异构体1行）
**建议**：如不区分异构体 → **保持现状**

---

## ✅ 核心结论

**当前数据质量**：⭐⭐⭐⭐½ (4.5/5)

**主要缺失**：
1. 🔴 基因ID映射（ncbi_gene_id、ensembl_gene_id等）
2. 🔴 条目质量标识（entry_type、is_canonical）
3. 🔴 结构类型字段（structure_available）
4. 🔴 数据审计字段（audit_status）

**补全建议**：
- **立即实施**：Phase 1 + Phase 2（1.5天）
- **按需实施**：Phase 3-5（取决于下游需求）

**补全后评分**：⭐⭐⭐⭐⭐ (5/5)

---

**文档版本**：v1.0
**创建时间**：2025-10-27
**作者**：项目团队
**下次审核**：Phase 1完成后
