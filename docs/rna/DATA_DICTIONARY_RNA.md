Identity Rule:

主键 (rna_id) 必须使用 RNAcentral URS 格式，且必须包含物种后缀（如 URS00000478B7_9606）。

rna_id 是全局唯一的，不允许重复。

字段定义表 (Tier 1 - v1 必填)
v1 阶段，如果不填满这些字段，就不能上线。
| 字段名 (Field)  | 类型            | 示例                 | 规则与来源                                                |
| ------------ | ------------- | ------------------ | ---------------------------------------------------- |
| rna_id       | String (PK)   | URS00000478B7_9606 | 主键。来自 RNAcentral。必须以 _9606 结尾​。 |
| rna_type     | String (Enum) | mirna, mrna        | 统一转小写。v1 仅限: mirna, mrna, transcript, lncrna。        |
| rna_name     | String        | hsa-miR-21-5p      | miRNA 用 miRBase 名；mRNA 用 Symbol 或 Transcript Name。   |
| sequence     | String        | UAGCUUAUC...       | 必须非空。只允许 A/C/G/U/N 字符。                               |
| sequence_len | Integer       | 22                 | 由 sequence 计算得出。                                     |
| taxon_id     | Integer       | 9606               | 固定值。方便未来扩展多物种，v1 全填 9606。                            |
| symbol       | String        | MIR21, TP53        | 对应的 Gene Symbol。用于和蛋白层对齐。                            |
| source       | String        | RNAcentral;miRBase | 数据来源，多个用分号 ; 分隔。                                     |
| fetch_date   | Date          | 2025-12-18         | ETL 运行日期。                                            |
source_version	String	RNAcentral:25;miRBase:22.1	多源用分号分隔。格式：源名:版本号。必须非空。
hgnc_id	String	HGNC:7097	格式 HGNC:xxxxx。mRNA 覆盖率目标 ≥70%，miRNA 可空。
ensembl_gene_id	String	ENSG00000141510	格式 ENSGxxxxxxxxxxx。mRNA 覆盖率目标 ≥80%，miRNA 可空。
ncbi_gene_id	Integer	7157	NCBI Gene ID。mRNA 覆盖率目标 ≥60%，miRNA 可空。


### Gene 映射字段优先级规则

**对于 mRNA/transcript**：
- **必须满足**：`hgnc_id`, `ensembl_gene_id`, `ncbi_gene_id` 三者**至少有一个非空**
- **优先级**：`ensembl_gene_id` > `hgnc_id` > `ncbi_gene_id`
- **覆盖率要求**：mRNA 整体的 gene ID 非空率 ≥80%

**对于 miRNA**：
- `symbol` 能填则填（如 MIR21），无则空
- `hgnc_id` 能填则填（部分 miRNA 有对应的 MIR 基因），无则空
- **不强制**要求 miRNA 有 gene ID（因为部分成熟 miRNA 无基因注释）



预留字段 (Tier 2 - v1 允许为空)
这些字段你现在就要把列建好，ETL 里如果方便拿就填，拿不到就填 null / NA。

字段名	备注
mirbase_id	专门存 MIMAT0000076 这种 ID，方便精确检索。
ensembl_transcript_id	专门存 ENST00000xxx。
hgnc_id	基因 ID，如 HGNC:7097。
rfam_id	Rfam 家族 ID (Tier 2) ​。
secondary_structure	二级结构 Dot-bracket 字符串 (Tier 2)。
pdb_ids	相关 PDB 结构 ID 列表，如 1A34;2KOC。

划定 v1 边界 (Define Scope)
## v1 排除范围 (Out of Scope)
以下内容属于 L2/L3 或后续版本，v1 阶段**严禁**引入：

1. **RPI 互作数据**：所有 RNA-Protein, RNA-RNA 互作（RNAInter, NPInter）属于 L2 关系表，不放在本实体表中 [file:1][file:2]。
2. **复杂功能注释**：GO 注释、ENCODE 表达量数据暂不抓取。
3. **非 Human 数据**：严禁混入小鼠或其他模式生物数据。
4. **其他小 RNA**：tRNA, rRNA, snRNA 除非顺手抓取，否则不专门清洗，v1 集中精力在 miRNA 和 mRNA。

