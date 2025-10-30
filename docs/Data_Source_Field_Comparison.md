## 📄 精简版报告（v3.0 - 实战版）


# 蛋白质知识图谱数据源整合方案（精简版 v3.0）

## 1. 整合原则

### 目标
从蛋白质实体角度，采集完整一级信息。1个UniProt ID = 1行蛋白。

### 核心策略
```
HGNC（权威基因符号）→ UniProt（主干数据）→ 去重 → 多源补漏
                          ↓
            [基因ID] [结构] [功能] [定位] [疾病] [PTM]
```

---

## 2. 5层数据源模型

| 层级 | 主源 | 用途 | 覆盖 | 状态 |
|------|------|------|------|------|
| **L1** | **UniProt** | 序列/功能/定位/疾病/PTM | 100% | ✅完成 |
| **L2** | HGNC+Ensembl | 基因命名/ID映射 | 99.6% | ✅完成 |
| **L3** | PDB+AlphaFold | 实验/预测结构 | 95% | ✅完成 |
| **L4** | GO | 功能补(function空时) | 82-94% | ✅完成 |
| **L5** | STRING+Reactome | 蛋白交互/通路(KG二期) | - | ⏳规划 |

---

## 3. 关键字段清单（33列现状→46列目标）

### L1: UniProt核心字段（29列，v3基准）
```
基础: uniprot_id, entry_name, protein_name, gene_names
序列: sequence(TEXT完整), sequence_len, mass
注释: function(TEXT), subcellular_location(TEXT), ptms, diseases(TEXT)
功能: go_biological_process, go_molecular_function, go_cellular_component
交叉: pdb_ids, string_ids, keywords, domains, isoforms(TEXT)
元数据: date_modified, source, fetch_date
```

### L2: 基因映射字段（4列新增，v6现有）
```
✅ ncbi_gene_id (99.6%)
✅ ensembl_gene_id (99.4%)
✅ gene_synonyms (80.9%)
⚠️ ensembl_transcript_id (0% - HGNC无, 可跳)
```

### L3: 结构质量字段（2列已有）
```
✅ alphafold_id, alphafold_mean_plddt (99.7%)
✅ has_alphafold, alphafold_pdb_url, alphafold_entry_url
```

### L4: GO补充字段（需补，Phase 2+4）
```
⏳ go_id_bp / go_id_mf / go_id_cc (GO term ID列表)
⏳ go_evidence_bp / go_evidence_mf / go_evidence_cc (证据代码: IDA/IEA/etc)
```

### L5: 互作字段（需补，KG二期）
```
⏳ string_id, string_score (交互蛋白列表 → 单独edges表)
⏳ reactome_pathway_id (通路 → 单独pathways表)
```

---

## 4. 整合流程（Pandas伪代码）

### Step A: 基准 (v3_clean)
```
df = pd.read_csv('protein_master_v3_clean.tsv')  # 19,135 unique uniprot_id
```

### Step B: 补L2 (基因ID, ✅完成)
```
df_hgnc = pd.read_csv('hgnc_core.tsv')
mapping = df_hgnc.drop_duplicates('symbol').set_index('symbol')[['entrez_id', 'ensembl_gene_id', 'alias_symbol']]
df = df.merge(mapping.loc[df['primary_symbol']], how='left', on='symbol')
# Result: ncbi_gene_id 99.6%, ensembl_gene_id 99.4%
```

### Step C: 补L4 (GO)
```
# UniProt GO已有 (go_biological_process + term ID)
# 补缺: df.loc[df['go_biological_process'].isnull(), 'go_biological_process'] = GO_API_query(ncbi_gene_id)
# → 83% → 95%
```

### Step D: 补L5 (STRING/Reactome, 下一步)
```
# STRING: create edges_merged.tsv (uniprot_id_A, uniprot_id_B, score, evidence)
# Reactome: create pathways_merged.tsv (uniprot_id, reactome_id, role)
```

---

## 5. 缺失补救规则

| 字段 | 空值率 | 触发条件 | 补源 | 预期覆盖 |
|------|--------|---------|------|----------|
| **function** | 14% | cc_function空 | GO definition | 95% |
| **go_ids** | 18% | 新增计算 | UniProt GO term extract | 95% |
| **pdb_ids** | 56% | 可空 | 仅补if structure_quality>3Å | - |
| **diseases** | 73% | 可空 | UniProt已完整 | 27% |
| **ptms** | 5% | 补细节 | PhosphoSite (下一步) | >95% |

**触发逻辑**：
```
if df['function'].isnull():
    df['function'] += GO_API(ncbi_gene_id)['definition']
    df['function_source'] = 'GO'
else:
    df['function_source'] = 'UniProt'
```

---

## 6. 版本与质控

### 数据版本记录
```
protein_master_v6_clean.tsv
├─ source_versions:
│  ├─ UniProt: 2025_04 (10/2025)
│  ├─ HGNC: 2025-10
│  ├─ Ensembl: GRCh38.p14
│  ├─ AlphaFold: v4 (2025版)
│  └─ PDB: 2025-10 (weekly)
├─ row_count: 19,135 (unique uniprot_id)
└─ column_count: 33
```

### 质控指标（Audit)
```
✅ 唯一性: uniprot_id = 100%
✅ 空值: function<1%, go_ids<18% (补漏目标<5%)
✅ ID一致: ncbi_gene_id vs ensembl_gene_id 无矛盾 >99%
✅ 数据类型: TEXT/Array/Float 各字段正确
```

---

## 7. 下一步行动（优先级）

| 优先级 | 任务 | 时间 | 产出 |
|--------|------|------|------|
| **P0** | ⏳ Phase 2+4: GO_ids + go_evidence 补充 | 1h | v7 (37列) |
| **P1** | ⏳ 异构体附表: protein_isoforms.tsv | 2h | 50k行 |
| **P2** | ⏳ STRING互作: edges_merged.tsv | 2h | 50k边 |
| **P3** | ⏳ Reactome通路: pathways_merged.tsv | 1h | 30k通路 |
| **P4** | 📖 知识图谱构建 (KG nodes/edges) | - | 超图 |

---

## 总体现状

✅ **蛋白质基础实体完成** (v6_clean: 19,135行, 33列)
- UniProt序列/功能/定位
- 基因ID映射（99%+）
- 结构补充（95%）

⏳ **功能/交互补充中**
- Phase 2+4: GO IDs + 证据 → v7
- 然后：STRING/Reactome (KG二期)

---

**文档版本**: v3.0  
**最后更新**: 2025-10-30  
**对标**: 实战进度（vs理想化方案）
```
