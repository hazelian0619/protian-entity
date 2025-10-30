# 1025 蛋白质知识图谱项目

## 🎯 项目目标

收集截止到2025年11月前更新的**全部人类蛋白质实体信息**，构建以蛋白质为中心的高质量数据集。

### 核心原则
- ✅ **以蛋白质为主体**：主键为`uniprot_id`（包含异构体）
- ✅ **保留完整原文**：TEXT字段无截断，保留文献引用
- ✅ **多源整合**：UniProt + PDB + AlphaFold + HGNC + Ensembl + GO + STRING
- ✅ **时效性**：数据≤2025-10-31
- ✅ **完整性**：目标覆盖率>99%，结构覆盖>95%

---

## 📊 当前数据概览

### 数据文件
```
data/processed/
├── protein_master_v2.tsv      # 19,253条  - 主表（蛋白质实体）
├── protein_edges.tsv           # 884,555条 - STRING交互网络
├── ptm_sites.tsv              # 235,502条 - 翻译后修饰位点
├── pathway_members.tsv         # 127,096条 - 通路成员
├── hgnc_core.tsv              # 19,278条  - HGNC基因命名
└── uniprot_seed.json          # 1.7MB     - UniProt种子数据
```

### 文档
```
docs/
├── Project_Goal_Revision.md          # 项目目标矫正（vs 1017问题诊断）
├── Data_Source_Field_Comparison.md   # 8大数据源字段对比
├── Unified_Protein_Schema_v2.md      # 统一Schema设计
└── Data_Quality_Assessment.md        # 数据质量评估与完善建议 ⭐
```

---

## ⚠️ 核心问题

### 🔴 问题1：异构体未展开（最严重）
- **现状**：19,253条记录
- **目标**：~100,000条（包含所有异构体）
- **原因**：`isoforms`字段有数据，但未展开为独立行
- **影响**：违背"以蛋白质为主体"原则，丢失80%变体数据

**示例**：
```
当前：P04637 (TP53, 1行, 含3个异构体信息但未拆分)
应为：
  - P04637    (TP53, Isoform 1, is_canonical=true)
  - P04637-2  (TP53, Isoform 2, is_canonical=false)
  - P04637-3  (TP53, Isoform 3, is_canonical=false)
```

### 🟡 问题2：字段不完整
**缺失关键字段**：
- `entry_type` (Swiss-Prot/TrEMBL)
- `is_canonical` (是否代表性异构体)
- `ncbi_gene_id`, `ensembl_gene_id`
- `structure_available`, `alphafold_mean_plddt`
- `audit_status`, `audit_comments`

### 🟡 问题3：结构覆盖率低
- PDB覆盖：~30%
- AlphaFold字段：**缺失**
- 目标：95%（PDB 30% + AlphaFold 65%）

---

## ✅ 当前优点

1. **字段质量好**：
   - ✅ `function`字段保留完整原文（>1000字）
   - ✅ 包含文献引用（ECO:0000269|PubMed:xxx）
   - ✅ `subcellular_location`详细准确

2. **多源整合**：
   - ✅ UniProt主干
   - ✅ STRING交互网络（88万边）
   - ✅ PTM修饰数据（23万位点）
   - ✅ 通路信息（12万成员）

3. **时间戳管理**：
   - ✅ `date_modified`: 2025-06-18
   - ✅ `fetch_date`: 2025-10-26

---

## 📋 完善清单（优先级排序）

### 🔴 高优先级（必须，1-2天）

#### 1. 异构体展开
- [ ] 解析`isoforms`字段
- [ ] 调用UniProt API获取每个异构体的序列
- [ ] 生成独立行，添加`is_canonical`字段
- [ ] 目标：19k → 100k行

#### 2. 补充AlphaFold结构
- [ ] 调用AlphaFold API
- [ ] 添加字段：
  - `structure_available` (PDB/AlphaFold/None)
  - `alphafold_model_url`
  - `alphafold_mean_plddt`
  - `pdb_best_resolution`
- [ ] 目标：结构覆盖率>95%

#### 3. 补充基因ID映射
- [ ] 添加字段：
  - `ncbi_gene_id`
  - `ensembl_gene_id`
  - `ensembl_transcript_id`
- [ ] 数据源：UniProt交叉引用 或 HGNC映射

#### 4. 添加审计字段
- [ ] `audit_status` (OK/MULTI/MISSING/PARTIAL)
- [ ] `audit_comments` (补漏日志)

---

### 🟡 中优先级（重要，2-3天）

#### 5. 创建辅助表
- [ ] `protein_aliases.tsv` (别名映射)
- [ ] `protein_features.tsv` (domain/region/binding site)
- [ ] `pdb_structures.tsv` (PDB详情：分辨率/配体/方法)
- [ ] `disease_variants.tsv` (疾病变异)

#### 6. 字段拆分
- [ ] `gene_names` → `hgnc_symbol` + `gene_synonyms`
- [ ] `mass` → `sequence_mass`
- [ ] 补充详细功能字段：
  - `catalytic_activity`
  - `cofactor`
  - `activity_regulation`
  - `subunit_structure`
  - `pathway`

---

### 🟢 低优先级（可选，1-2天）

#### 7. 数据版本管理
- [ ] 添加`source_version`字段
- [ ] 记录各数据源版本：
  ```
  UniProt:2025_04, PDB:2025-10, AlphaFold:v4,
  HGNC:2025-10, STRING:v12.5, GO:2025-10-01
  ```

#### 8. 质量检查
- [ ] 空值率统计
- [ ] ID唯一性检查
- [ ] 序列长度一致性验证
- [ ] 生成QA报告

---

## 🎯 预期改善

| 指标 | 当前 | 目标 | 改善 |
|------|------|------|------|
| **记录数** | 19,253 | ~100,000 | 5倍 ✅ |
| **字段数** | 24 | ~40 | 全面 ✅ |
| **结构覆盖** | 30% | 95% | 3倍 ✅ |
| **异构体** | 0% | 100% | 完全 ✅ |
| **辅助表** | 3 | 7 | 结构化 ✅ |
| **完整率** | 85% | >99% | 达标 ✅ |

---

## 🚀 快速开始

### 查看详细评估报告
```bash
cat docs/Data_Quality_Assessment.md
```

### 查看Schema设计
```bash
cat docs/Unified_Protein_Schema_v2.md
```

### 查看数据源对比
```bash
cat docs/Data_Source_Field_Comparison.md
```

---

## 📚 关键文档索引

| 文档 | 用途 | 重要性 |
|------|------|--------|
| **Data_Quality_Assessment.md** | 数据质量评估与完善建议 | ⭐⭐⭐⭐⭐ |
| **Unified_Protein_Schema_v2.md** | 目标Schema设计 | ⭐⭐⭐⭐⭐ |
| **Project_Goal_Revision.md** | 项目目标矫正 | ⭐⭐⭐⭐ |
| **Data_Source_Field_Comparison.md** | 数据源字段对比 | ⭐⭐⭐⭐ |

---

## 🔗 下一步行动

1. **阅读** `docs/Data_Quality_Assessment.md`（最重要）
2. **确认** 完善清单中的优先级
3. **开始** 异构体展开脚本开发
4. **补充** AlphaFold结构数据
5. **创建** 辅助表

---

**项目状态**：数据收集完成 → 正在完善阶段
**当前评分**：⭐⭐⭐⭐ (4/5)
**目标评分**：⭐⭐⭐⭐⭐ (5/5)
**预计完成**：6-9天

---

**最后更新**：2025-10-26
**责任人**：项目团队
