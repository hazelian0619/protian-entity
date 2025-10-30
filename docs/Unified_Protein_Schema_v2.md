## 文档分析
主表protein_master覆盖静态实体核心（uniprot_id PK、序列/功能/结构/PTM全TEXT+JSON），辅助表（aliases/features/pdb_structures/disease_variants）支持扩展/查询优化，原则对齐蛋白中心（异构体展开、血缘追溯），校验/索引/实施实用，示例TP53详实，差异对比突出1017缺陷修复。这匹配项目点位：完整原文无截断、多源整合（UniProt主+PDB/AlphaFold/HGNC/Ensembl）、时效（版本/日期追踪），总字段40，主表行估100k，支持KG节点导入（Neo4j/Pandas）。 


# 统一蛋白质信息Schema设计（v3.0）

## 1. 设计原则

### 1.1 核心原则
- ✅ **以蛋白质为主体**：主键为 `uniprot_id`，展开所有异构体（isoforms），避免基因聚合丢失变体。
- ✅ **保留完整信息**：文本字段统一 TEXT 类型，无长度截断，支持原文引用/特殊字符（UTF-8）。
- ✅ **支持异构体**：每行一变体（100k人类行），is_canonical 标记代表性。
- ✅ **数据血缘追溯**：统一 source_version 列追踪多源，辅以 fetch_date/audit 日志。
- ✅ **多源整合**：优先 Swiss-Prot > TrEMBL；补 GO/STRING/PhosphoSite（功能/PTM/动态）；外键 on uniprot_id，outer join 补漏（空<1%）。
- ✅ **动态预备**：主表静态实体，辅表 edges 起步超图交互（STRING score>0.7，Reactome 免疫通路）。

### 1.2 范围界定
- **物种**：人类（Homo sapiens, Taxonomy ID: 9606），filter 非蛋白编码。
- **数据截止日期**：2025-11-01（月增，cron 查新版）。
- **优先级**：UniProt 主干；结构 PDB (resolution<3Å) > AlphaFold (pLDDT>70) > None；功能 UniProt > GO。
- **完整率目标**：>99%（QA: 空字段 API 补，status=OK>95%）。
- **存储估**：主表 100k 行，200MB（序列 GZIP 压缩 50%）。

***

## 2. 主表设计：`protein_master`

主表静态实体底座（20 核心字段），每行一蛋白变体，支持 KG 节点（属性 TEXT/JSON）。[3][4]

### 2.1 基本标识字段

| 字段名 | 数据类型 | 约束 | 说明 | 示例 |
|--------|---------|------|------|------|
| `uniprot_id` | VARCHAR(20) | **PRIMARY KEY** | UniProt Accession（含-数字异构体） | P04637-2 |
| `uniprot_entry_name` | VARCHAR(50) | NOT NULL | UniProt Entry Name | TP53_HUMAN |
| `entry_type` | ENUM('Swiss-Prot', 'TrEMBL') | NOT NULL | 条目类型（reviewed > automatic） | Swiss-Prot |
| `protein_name_recommended` | VARCHAR(500) | NOT NULL | 推荐名称 | Cellular tumor antigen p53 |
| `protein_name_alternative` | TEXT | NULLABLE | 其他名称（JSON 数组） | ["Tumor suppressor p53", "Antigen NY-CO-13"] |
| `protein_name_short` | VARCHAR(100) | NULLABLE | 简称 | p53 |
| `is_canonical` | BOOLEAN | NOT NULL | 代表性异构体 | TRUE |
| `isoform_name` | VARCHAR(100) | NULLABLE | 异构体名称 | Isoform 2 |

**字段说明**：uniprot_id 唯一（explode UniProt isoforms array）；entry_type 优先 Swiss-Prot（50k 高质）。[5][1]

### 2.2 基因关联字段

| 字段名 | 数据类型 | 约束 | 说明 | 示例 |
|--------|---------|------|------|------|
| `hgnc_id` | VARCHAR(20) | NULLABLE | HGNC ID（权威命名） | HGNC:11998 |
| `hgnc_symbol` | VARCHAR(50) | NULLABLE | 官方基因符号 | TP53 |
| `gene_name_official` | VARCHAR(500) | NULLABLE | 基因全称 | tumor protein p53 |
| `gene_synonyms` | TEXT | NULLABLE | 同义词（JSON 数组，去重 HGNC+UniProt） | ["P53", "TRP53", "BCC7"] |
| `ncbi_gene_id` | VARCHAR(20) | NULLABLE | NCBI Gene ID（外键） | 7157 |
| `ensembl_gene_id` | VARCHAR(20) | NULLABLE | Ensembl 基因 ID | ENSG00000141510 |
| `ensembl_transcript_id` | VARCHAR(20) | NULLABLE | Ensembl 转录本 ID（异构体源） | ENST00000269305 |

**字段说明**：一 uniprot_id 一基因，多基因 fallback Ensembl（不匹配 1% audit）；synonyms 优化检索。[4][3]

### 2.3 序列信息字段

| 字段名 | 数据类型 | 约束 | 说明 | 示例 |
|--------|---------|------|------|------|
| `sequence` | TEXT | NOT NULL | 完整氨基酸序列（单字母 raw，无 FASTA 头） | MEEPQSDPSVEPPLSQ... |
| `sequence_length` | INTEGER | NOT NULL | 长度（残基数） | 393 |
| `sequence_mass` | INTEGER | NOT NULL | 分子量 (Da) | 43653 |
| `sequence_crc64` | VARCHAR(16) | NULLABLE | 校验和（去重） | 9A3C5C8E4C3C8B1A |
| `sequence_md5` | VARCHAR(32) | NULLABLE | MD5 哈希 | 7C8F0E... |
| `sequence_version` | INTEGER | NOT NULL | 序列版本 | 4 |
| `sequence_modified_date` | DATE | NOT NULL | 最后修改日期 | 2025-10-15 |

**字段说明**：sequence 初级结构核心，校验防重复；长度 = len(sequence) 校验。[1][5]

### 2.4 功能注释字段（完整原文）

| 字段名 | 数据类型 | 约束 | 说明 | 数据源 |
|--------|---------|------|------|--------|
| `function_text` | TEXT | NULLABLE | 功能描述原文 | UniProt (cc_function) + GO 补 |
| `catalytic_activity` | TEXT | NULLABLE | 催化活性原文 | UniProt (cc_catalytic_activity) |
| `cofactor` | TEXT | NULLABLE | 辅因子原文 | UniProt (cc_cofactor) |
| `activity_regulation` | TEXT | NULLABLE | 调控机制原文 | UniProt (cc_activity_regulation) |
| `subunit_structure` | TEXT | NULLABLE | 亚基组成原文 | UniProt (cc_subunit) |
| `pathway` | TEXT | NULLABLE | 通路原文（动态 hint） | UniProt (cc_pathway) + Reactome |
| `interaction_partners` | TEXT | NULLABLE | 互作蛋白（JSON，预 edges） | UniProt (cc_interaction) + STRING |

**重要提示**：TEXT 无限制，保留文献引用；空 function API 拉 GO definition (evidence=IDA) 补 10%。[6][4]

### 2.5 亚细胞定位字段

| 字段名 | 数据类型 | 约束 | 说明 | 数据源 |
|--------|---------|------|------|--------|
| `subcellular_location` | TEXT | NULLABLE | 亚细胞定位原文 | UniProt (cc_subcellular_location) |
| `tissue_specificity` | TEXT | NULLABLE | 组织表达原文 | UniProt (cc_tissue_specificity) + HPA 补 |
| `developmental_stage` | TEXT | NULLABLE | 发育表达原文 | UniProt (cc_developmental_stage) |
| `induction` | TEXT | NULLABLE | 诱导条件原文 | UniProt (cc_induction) |

### 2.6 结构信息字段

| 字段名 | 数据类型 | 约束 | 说明 | 数据源 |
|--------|---------|------|------|--------|
| `structure_available` | ENUM('PDB', 'AlphaFold', 'None') | NOT NULL | 结构可用性（优先级） | PDB > AlphaFold |
| `pdb_ids` | TEXT | NULLABLE | PDB ID 列表（JSON 数组） | UniProt cross + PDB API |
| `pdb_best_resolution` | FLOAT | NULLABLE | 最佳分辨率 (Å, <3 高质) | PDB (min 值) |
| `alphafold_model_url` | VARCHAR(500) | NULLABLE | AlphaFold 下载链接 | AlphaFold DB |
| `alphafold_mean_plddt` | FLOAT | NULLABLE | 平均 pLDDT (0-100, >70 高可信) | AlphaFold (np.mean) |
| `alphafold_version` | VARCHAR(20) | NULLABLE | 版本 | v4 |

**字段说明**：pLDDT >90 Very high, 70-90 High；None 时 audit '补 AlphaFold'。[3][1]

### 2.7 疾病关联字段

| 字段名 | 数据类型 | 约束 | 说明 | 数据源 |
|--------|---------|------|------|--------|
| `disease_involvement` | TEXT | NULLABLE | 疾病关联原文 | UniProt (cc_disease) |
| `disease_ids` | TEXT | NULLABLE | 疾病 ID（JSON） | {"OMIM": ["191170"], "MeSH": ["D002277"]} |
| `polymorphism` | TEXT | NULLABLE | 多态性描述 | UniProt (cc_polymorphism) |
| `mutagenesis_summary` | TEXT | NULLABLE | 致病突变总结 | UniProt (cc_mutagenesis) |

### 2.8 翻译后修饰字段（PTM）

| 字段名 | 数据类型 | 约束 | 说明 | 数据源 |
|--------|---------|------|------|--------|
| `ptm_processing` | TEXT | NULLABLE | 加工原文 | UniProt (cc_ptm) + PhosphoSite |
| `ptm_sites` | TEXT | NULLABLE | 修饰位点（JSON，磷/糖/泛等） | UniProt features + PhosphoSite |

**JSON 示例**：[{"position": 15, "type": "Phosphoserine", "evidence": "IDA", "source": "UniProt"}]。主表概览，辅表详展开。[4][6]

### 2.9 数据来源与时间戳字段（数据血缘）

| 字段名 | 数据类型 | 约束 | 说明 | 示例 |
|--------|---------|------|------|------|
| `primary_source` | VARCHAR(50) | NOT NULL | 主要源 | UniProt |
| `source_version` | VARCHAR(20) | NOT NULL | 版本（多源逗号分隔） | "UniProt:2025_04, HGNC:2025-10" |
| `last_modified` | DATE | NOT NULL | 最后修改 | 2025-10-15 |
| `integration_date` | TIMESTAMP | NOT NULL | 整合时间 | 2025-10-26 14:30:00 |
| `audit_status` | ENUM('OK', 'MULTI', 'MISSING') | NOT NULL | 质量状态 | OK |
| `audit_comments` | TEXT | NULLABLE | 补漏日志 | "补 GO function" |

**字段说明**：source_version 统一追踪（e.g., UniProt 2025_04）；audit 捕获不匹配/空，支持审稿。[3][4]

***

## 3. 辅助表设计

### 3.1 蛋白质别名表：`protein_aliases`

| 字段名 | 数据类型 | 约束 | 说明 |
|--------|---------|------|------|
| `alias_id` | SERIAL | PRIMARY KEY | 自增 |
| `uniprot_id` | VARCHAR(20) | FOREIGN KEY (protein_master) | 关联主表 |
| `alias_type` | ENUM('protein_alt', 'gene_alt', 'gene_prev', 'short') | NOT NULL | 类型 |
| `alias_value` | VARCHAR(500) | NOT NULL | 值 |
| `source` | VARCHAR(50) | NOT NULL | 来源 (HGNC/UniProt) |
| `fetch_date` | DATE | NOT NULL | 获取日期 |

**示例**：uniprot_id=P04637, alias_type=gene_alt, alias_value=P53, source=HGNC。[4]

### 3.2 结构特征表：`protein_features`

存储二级/三级结构域/位点（UniProt features + PDB）。

| 字段名 | 数据类型 | 约束 | 说明 |
|--------|---------|------|------|
| `feature_id` | SERIAL | PRIMARY KEY | 自增 |
| `uniprot_id` | VARCHAR(20) | FOREIGN KEY | 关联 |
| `feature_type` | VARCHAR(50) | NOT NULL | 类型 (Domain/Region/Binding/Helix 等) |
| `description` | TEXT | NULLABLE | 描述 |
| `start_pos` | INTEGER | NOT NULL | 起始 |
| `end_pos` | INTEGER | NOT NULL | 终止 |
| `evidence` | VARCHAR(20) | NULLABLE | 证据 (IDA/预测) |
| `source` | VARCHAR(50) | NOT NULL | 来源 |

### 3.3 PDB 结构详情表：`pdb_structures`

| 字段名 | 数据类型 | 约束 | 说明 |
|--------|---------|------|------|
| `pdb_id` | VARCHAR(4) | PRIMARY KEY | PDB ID |
| `uniprot_id` | VARCHAR(20) | FOREIGN KEY | 关联 |
| `title` | TEXT | NULLABLE | 标题 |
| `method` | VARCHAR(50) | NOT NULL | 方法 (X-RAY/NMR) |
| `resolution` | FLOAT | NULLABLE | 分辨率 |
| `r_factor` | FLOAT | NULLABLE | R-factor |
| `ligands` | TEXT | NULLABLE | 配体 (JSON) |
| `release_date` | DATE | NOT NULL | 发布 |
| `citation` | TEXT | NULLABLE | 文献 |
| `fetch_date` | DATE | NOT NULL | 获取 |

### 3.4 疾病变异表：`disease_variants`

| 字段名 | 数据类型 | 约束 | 说明 |
|--------|---------|------|------|
| `variant_id` | SERIAL | PRIMARY KEY | 自增 |
| `uniprot_id` | VARCHAR(20) | FOREIGN KEY | 关联 |
| `position` | INTEGER | NOT NULL | 位置 |
| `original_aa` | CHAR(1) | NOT NULL | 原 AA |
| `variant_aa` | CHAR(1) | NOT NULL | 变异 AA |
| `disease` | VARCHAR(500) | NULLABLE | 疾病 |
| `dbsnp_id` | VARCHAR(20) | NULLABLE | dbSNP |
| `clinvar_id` | VARCHAR(20) | NULLABLE | ClinVar |
| `significance` | ENUM('Pathogenic', 'Benign', 'Uncertain') | NULLABLE | 临床 |
| `evidence` | VARCHAR(20) | NOT NULL | 证据 |

### 3.5 新增：交互边表：`protein_edges`（动态预备）

| 字段名 | 数据类型 | 约束 | 说明 |
|--------|---------|------|------|
| `edge_id` | SERIAL | PRIMARY KEY | 自增 |
| `source_uniprot_id` | VARCHAR(20) | FOREIGN KEY | 源蛋白 |
| `target_id` | VARCHAR(50) | NOT NULL | 目标 (uniprot_id 或 pathway: R-HSA-...) |
| `edge_type` | VARCHAR(50) | NOT NULL | 类型 (binding/pathway) |
| `score` | FLOAT | NULLABLE | 置信 (STRING 0-1) |
| `source` | VARCHAR(50) | NOT NULL | 来源 (STRING/Reactome) |
| `fetch_date` | DATE | NOT NULL | 获取 |

**说明**：起步 30k 免疫边（Reactome filter），支持超图（weight=score）。[6][3]

***

## 4. 索引设计

### 4.1 主表索引
```sql
CREATE INDEX idx_uniprot_id ON protein_master(uniprot_id);  -- PK 已隐含
CREATE INDEX idx_hgnc_symbol ON protein_master(hgnc_symbol);
CREATE INDEX idx_ncbi_gene_id ON protein_master(ncbi_gene_id);
CREATE INDEX idx_ensembl_gene_id ON protein_master(ensembl_gene_id);
CREATE INDEX idx_entry_type ON protein_master(entry_type);
CREATE INDEX idx_is_canonical ON protein_master(is_canonical);
CREATE INDEX idx_structure_available ON protein_master(structure_available);
CREATE INDEX idx_function_text ON protein_master USING GIN(to_tsvector('english', function_text));  -- 全文搜索
```

### 4.2 辅助表索引
```sql
-- aliases
CREATE INDEX idx_alias_uniprot ON protein_aliases(uniprot_id);
CREATE INDEX idx_alias_value ON protein_aliases USING GIN(to_tsvector('english', alias_value));

-- features/edges 等类似 uniprot_id + type
CREATE INDEX idx_feature_uniprot ON protein_features(uniprot_id);
CREATE INDEX idx_edge_source ON protein_edges(source_uniprot_id);
```

**说明**：GIN/FTS 优化模糊/全文查询（e.g., "p53 function"）。[4]

***

## 5. 数据校验规则

### 5.1 必填字段检查
- ✅ `uniprot_id`：格式 [A-Z0-9]{6,10}(-[1-9]\d*)?，唯一 >99%。
- ✅ `sequence`：仅 20 标准 AA (ACDEFGHIKLMNPQRSTVWY)，len=sequence_length。
- ✅ `entry_type`：仅 ENUM 值。
- ✅ `audit_status`：OK 若无 missing/muti。

### 5.2 数据一致性检查
- ✅ `pdb_ids` 非空 → structure_available='PDB'，resolution >0 & <4。
- ✅ `alphafold_mean_plddt` 在 0-100，>70 优先 AlphaFold。
- ✅ `gene_synonyms` 去重，无 HGNC 冲突（权威）。
- ✅ 空率 <1%（function/PTM API 补，日志 audit_comments）。
- ✅ 引用完整：辅表 uniprot_id EXISTS 主表。

### 5.3 质量审计
- 运行 QA：空字段 %、ID 唯一、版本一致（source_version ≤2025-11-01）。
- 示例脚本：Pandas df.isnull().sum() / len(df) <0.01。[3][4]

***

## 6. 示例记录

### TP53（P04637）完整记录示例

```json
{
  "uniprot_id": "P04637",
  "uniprot_entry_name": "TP53_HUMAN",
  "entry_type": "Swiss-Prot",
  "protein_name_recommended": "Cellular tumor antigen p53",
  "protein_name_alternative": ["Tumor suppressor p53", "Phosphoprotein p53", "Antigen NY-CO-13"],
  "protein_name_short": "p53",
  "is_canonical": true,
  "isoform_name": null,
  "hgnc_id": "HGNC:11998",
  "hgnc_symbol": "TP53",
  "gene_name_official": "tumor protein p53",
  "gene_synonyms": ["P53", "TRP53", "BCC7", "LFS1"],
  "ncbi_gene_id": "7157",
  "ensembl_gene_id": "ENSG00000141510",
  "ensembl_transcript_id": "ENST00000269305",
  "sequence": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD",
  "sequence_length": 393,
  "sequence_mass": 43653,
  "sequence_crc64": "9A3C5C8E4C3C8B1A",
  "sequence_md5": "7C8F0E...",
  "sequence_version": 4,
  "sequence_modified_date": "2025-10-15",
  "function_text": "Acts as a tumor suppressor... (完整 UniProt + GO apoptotic process 定义)",
  "subcellular_location": "Cytoplasm. Nucleus... (完整)",
  "structure_available": "PDB",
  "pdb_ids": ["1TUP", "1TSR", "2OCJ"],
  "pdb_best_resolution": 1.8,
  "alphafold_model_url": "https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v4.pdb",
  "alphafold_mean_plddt": 85.3,
  "alphafold_version": "v4",
  "disease_involvement": "Li-Fraumeni syndrome (LFS)... (完整)",
  "primary_source": "UniProt",
  "source_version": "UniProt:2025_04, HGNC:2025-10, GO:2025-10-01",
  "last_modified": "2025-10-15",
  "integration_date": "2025-10-26 14:30:00",
  "audit_status": "OK",
  "audit_comments": null
}
```

***

## 7. 与1017/v2.0版本的关键差异

| 设计维度 | 1017（错误） | v2.0 | v3.0（优化） |
|----------|--------------|------|--------------|
| **主键** | ncbi_gene_id | uniprot_id | uniprot_id（+异构体展开） |
| **异构体处理** | 丢失 | 完整（is_canonical） | +ensembl_transcript 源 |
| **文本字段** | VARCHAR(400) | TEXT | TEXT + JSONB (PostgreSQL 查询优) |
| **结构信息** | 仅 flag | 分辨率/pLDDT | +阈值 (>70/ <3Å) + ligands JSON |
| **数据血缘** | 无 | 每源 date | 统一 source_version + audit_status |
| **多源整合** | 单一 | 5 源 | 8+ 源 (加 GO/STRING/PhosphoSite) |
| **辅助表** | 仅别名 | 4 表 | 5 表 (+ edges 动态) |
| **校验** | 无 | 基本 | +空率<1%/唯一>99%/全文索引 |

[5][6]

***

## 8. 实施建议

### 8.1 数据库选择
- **PostgreSQL** (推荐)：JSONB/GIN 支持，TEXT TOAST 压缩长文；100k 行高效。
- **MySQL 8.0+**：JSON 类型，需 innodb_large_prefix。
- **SQLite**：测试用，<10k 行。

### 8.2 存储优化
- sequence：GZIP 压缩（PostgreSQL toast），减 50%。
- function_text 等：全文索引 GIN，提升检索（e.g., KG 模糊）。
- 总 200MB，备份 raw/processed。

### 8.3 增量更新策略
```sql
-- 查需更新
SELECT uniprot_id FROM protein_master WHERE last_modified < '2025-10-15' OR audit_status != 'OK';
-- Pandas ETL: df_new.merge(df_old, on='uniprot_id', how='outer') → append diff
```
Cron 月跑，API 分批 <25/秒。[4]

### 8.4 ETL 脚本提示
- Python + Pandas/SQLAlchemy：HGNC 种子 → UniProt API 展开 → outer join 多源 → to_sql('protein_master')。
- 示例：notebooks/schema_etl.py（一周内免疫 TP53 示例）。[3]

***

## 9. 下一步工作

1. ✅ Schema 设计（v3.0 优化）。
2. ⏳ 建表 SQL：docs/create_schema.sql (PostgreSQL)。
3. ⏳ ETL 流程：notebooks/etl_protein.py (HGNC→UniProt→补漏→QA)。
4. ⏳ 数据获取：API 批量 + raw 缓存 (data/raw/)。
5. ⏳ 质控：QA_Checklist.md (空<1%, 版本一致, audit>95%)。
6. ⏳ 动态扩展：edges 表起步免疫交互 (30k 行)。

**文档版本**：v3.0  
**创建时间**：2025-10-26  
**更新时间**：2025-10-26  
**状态**：对齐矫正/整合方案，待实施审核。  
**下次审核**：2025-10-27（跑小规模测试）。[1][6]

[1](https://cloud.tencent.com/developer/article/1168324)
[2](https://www.cnblogs.com/miyuanbiotech/p/12198840.html)
[3](https://blog.csdn.net/u011559552/article/details/142875765)
[4](https://www.fanruan.com/blog/article/315761/)
[5](http://ibi.zju.edu.cn/bioinplant/courses/chap6.pdf)
[6](https://blog.csdn.net/2401_83941020/article/details/140603565)
[7](https://patents.google.com/patent/CN106605228B/zh)
[8](https://www.pibb.ac.cn/pibbcn/article/html/20240402)
