  我分三层说：

- A. 总体阶段和依赖关系（先看大图）
- B. 每个阶段的具体任务（清晰可执行）
- C. 后续 R2 扩展留口（先标方向，不现在做）
  
  下面只讲策略和步骤，不写具体代码。
  
  ———
  
A. 总体阶段 & 依赖关系（大图）

  我们先把“必须按顺序做的”和“可以并行做的”拆出来。

阶段划分（只针对 RNA）

- R0：设计定盘星（Schema & 文档，不拉数据）
  - 输出：RNA Data Dictionary + 版本锚点说明
  - 必须先做，后面所有 ETL 都依赖它
- R1：RNA v1 种子集（真正落数据）
  - 拆成两条管线：
    - R1-miR：人类 miRNA 节点（完整）
    - R1-long：人类长链 RNA 节点（先选 mRNA/transcript 管线）
  - 这两条可以并行做，但都依赖 R0 的字段定义
  - 最后合并成：rna_master_v1.tsv
- R1-QA：质量检查 + README
  - 在 R1-miR & R1-long 都出结果后做
  - 输出：rna_master_v1_qc.txt + docs/RNA_README.md（或补充到总 README）
    
    > 顺序关系可以理解为：
    > R0（Schema） → R1-miR / R1-long（可并行） → 合并 → R1-QA & README
    
    ———
    
B. 每个阶段的具体任务

R0：设计定盘星（不写 ETL，只写“表头 & 规则”）

  目标：把永远不会后悔的东西一次性写清楚，后续不返工。

  R0-1：确定主键和 ID 风格

- 决定：
  - rna_id 用 RNAcentral 的人类 URS，全写形式：URSxxxxxxxx_9606
  - taxon_id 可以有一列，全填 9606，也可以暂时不加（你说只做人类，这里不纠结多物种了）
- 规则写进文档：
  - “RNA master 表中的主键 rna_id 与 RNAcentral URS 完全一致（人类：末尾 _9606），所有 xref 和扩展都围绕它。”
    
    R0-2：写 RNA master 的字段表（Data Dictionary）
    
    在 docs/ 下写一个类似 DATA_DICTIONARY_RNA.md 的文档，按 Tier 分层（v1 必做 / 可以先空 / 未来扩展）：
    
- Tier 1（v1 必做字段）
  - rna_id：主键（RNAcentral URS）
  - rna_type / biotype：mirna / mrna / lncrna 等
  - rna_name：比如 hsa-miR-21-5p, TP53 transcript variant 1
  - sequence：RNA 序列（A/C/G/U/N）
  - sequence_len：长度（按 sequence 算）
  - hgnc_id, symbol：能填则填，不强求 100%
  - ncbi_gene_id, ensembl_gene_id：长链 RNA 尽量填
  - mirbase_id：仅 miRNA 有，v1 要求尽量填
  - source：例如 RNAcentral;miRBase / RNAcentral;Ensembl
  - source_version：例如 RNAcentral:25;miRBase:22.1
  - fetch_date：抓取日期（整表常量也可以）
- Tier 2（v1 先预留列，可以大量为空）
  - genbank_id, refseq_id
  - ensembl_transcript_id
  - rfam_id, secondary_structure_dotbracket
  - pdb_ids, structure_source（Rfam/PDB/RhoFold/None）
  - sequence_md5
- Tier 3（只写在文档里，不在 v1 填）
  - ENCODE 表达相关字段（组织/细胞类型表达）
  - 功能注释文本（lncRNA 机制等）
  - 审计字段（质量标签）
    
    R0-3：写清版本锚点
    
    在一个小节里写清：
    
- RNAcentral Release 25 (2025-06-18)
- Rfam 15.0 (2024/2025)
- miRBase 22.1 (或 2025 更新版)
- 说明：RNA master v1 基于这些版本构建，后续更新将记录版本变更。
  
  > R0 做完，你有一份“RNA master 的表头 + 版本说明”，以后改 ETL 都不用改设计。
  
  ———
  
R1-miR：miRNA 管线（v1 必做，和 long RNA 可并行）

  目的：
  做出一批“人类 miRNA 节点”，主键是 URS，能映射回 miRBase ID 和（有则）基因。

  R1-miR-1：确定输入源

  你可以先不关心实际文件在哪，只在文档里钉规则：

- 源 1：miRBase 人类 miRNA 列表（ID + 序列）
- 源 2：RNAcentral 的 miRNA cross-ref（miRBase ↔ URS）
  
  R1-miR-2：构造“miRNA URS 表”
  
  逻辑步骤（以后写代码时照着做）：
  
1. 从 miRBase 拿到 human miRNA 列表：
  - 拿到：mirbase_id（hsa-xxx）、mirna_name、sequence
2. 用 RNAcentral 的 mapping 表把每个 mirbase_id 映射成一个或多个 rna_id (URS…_9606)：
  - 如果一个 miRBase ID 对应多个 URS：
    - 先选“人类 + 主 canonical 条目”，规则写清楚（例如选最长、或 RNAcentral 标记为代表的那个）
  - 如果某个 miRNA 暂时没有 URS：
    - 可以先跳过，或记录到一个 “未映射列表”（方便之后补）
3. 得到一张中间表：
  - 每行：rna_id, mirbase_id, rna_name, sequence, sequence_len, rna_type='mirna'
    
    R1-miR-3：补 gene 映射（能补多少补多少）
    
- 有些 miRNA 会有 MIR 基因映射（HGNC 上的 MIRxxx），你可以：
  - 用 RNAcentral/miRBase 提供的 gene mapping（如果有）
  - 没有的话，v1 可以接受基因映射为空（重点是 ID + 序列先有）
    
    填字段：
    
- hgnc_id（有则填）
- symbol（MIR 基因符号，如 MIR21）
- ncbi_gene_id, ensembl_gene_id（有就填，没有空）
  
  R1-miR-4：输出 miRNA 子表
  
  得到一个独立表，比如：
  
- 文件：rna_master_mirna_v1.tsv
- 行：所有 mapping 成功的人类 miRNA URS
- 字段：按 Tier 1 拿能填的都填上；Tier 2 留空也没关系。
  
  ———
  
R1-long：长链 RNA 管线（建议先做 mRNA/transcript）

  目的：
  给“蛋白–基因–转录本”这条链补上 RNA 实体——mRNA/transcript，是最自然的第一选择。

  你已经有蛋白 v1 + HGNC + 基因 ID，所以可以这样设计：

  R1-long-1：定义 gene seed 集合

- Seed 来自：
  - protein_master_v6_clean.tsv 中所有有 hgnc_id/symbol 的 gene
  - 或者来自 hgnc_core.tsv 中筛出的 protein-coding gene（本地用，不必上传）
    
    输出：一个 gene 列表（hgnc_id/symbol/ensembl_gene_id），作为“我们关心的基因”。
    
    R1-long-2：从 RNAcentral / Ensembl 抽 transcript URS
    
    目标：只抽这部分：
    
- 物种 = human
- biotype = protein_coding transcript（mRNA 等）
- gene_id 在 seed 集合内
  
  流程设想：
  
1. 用 RNAcentral 或 Ensembl 的 transcript 列表：
  - 对每条 transcript 拿到：URS, ensembl_transcript_id, ensembl_gene_id, sequence 等
2. 用 seed gene 集过滤：
  - 只保留属于 seed gene 的那些 transcript
    
    R1-long-3：构造长链 RNA 子表
    
    对每一个满足条件的 URS：
    
- 填 Tier 1 字段：
  - rna_id = URS
  - rna_type = 'mrna' 或 'transcript'（你可以在 schema 里统一定义）
  - rna_name（用 Ensembl 描述，或简名）
  - sequence, sequence_len
  - hgnc_id, symbol
  - ensembl_gene_id, ensembl_transcript_id
  - source = 'RNAcentral;Ensembl', source_version, fetch_date
- Tier 2：
  - 能拿到的 xref（RefSeq 等）可以填一部分；
  - Rfam/PDB/结构先不管。
    
    输出：
    
- 文件：rna_master_mrna_v1.tsv
- 行：所有 seed gene 对应的 transcript URS 节点
  
  ———
  
R1-merge：合并为统一的 rna_master_v1.tsv

  R1-miR 和 R1-long 做完，可以合并：

1. 合并规则：
  
- 直接按行拼接：
  - miRNA 子表 + mRNA 子表 → rna_master_v1.tsv
- 保证：
  - rna_id 在全表中是唯一的（如果两边有重复 URS，要合并行）
  - 列集合 = R0 里定义的 full schema（Tier 1/2 都在）
    
2. 输出：
  
- rna_master_v1.tsv：
  - 行：miRNA URS + mRNA/transcript URS
  - 列：所有 Tier 1 + Tier 2 字段
    
    ———
    
R1-QA：质量验收 + README

  和蛋白 v6 的 QA 类似：

  R1-QA-1：写一个 QA 清单（将来用脚本实现）

  检查指标：

- 主键：
  - rna_id 非空率 = 100%
  - rna_id uniqueness = 100%
- 序列：
  - sequence 非空率 ≈ 100%
  - 只含 A/C/G/U/N
  - sequence_len > 0，记录 min/max
- 类型：
  - rna_type 非空率 = 100%
  - rna_type 值集合包含至少 mirna 和 mrna
  - 各类型的行数统计一下（方便写到 README）
- 对齐能力：
  - 在 mRNA 部分：hgnc_id 或 ensembl_gene_id 非空率 > 某个位点（比如 >80%）
  - 在 miRNA 部分：mirbase_id 非空率接近 100%
- 版本与日期：
  - source_version, fetch_date 非空率 = 100%
    
    R1-QA-2：RNA README（写清楚 v1 的边界）
    
    在 docs/ 里写一个 RNA_README.md，内容包括：
    
- v1 覆盖：
  - 全体人类 miRNA（miRBase+RNAcentral）
  - 与当前 protein-coding gene 集相关的 mRNA/transcript
- v1 未覆盖：
  - lncRNA / tRNA / rRNA / 其他小 RNA
  - Rfam 二级结构 / PDB / RhoFold / ENCODE / RPI 边
- 质量概览（用 QA 结果）：
  - 行数、类型分布、序列完整率、gene 对齐率等
    
    ———
    
C. 可并行 vs 必须串行的总结

  最后再帮你提炼一下“哪些可以并行、哪些必须按顺序”：

- 必须先做：
  - R0：Schema + Data Dictionary + 版本锚点
- 在 R0 完成之后，可以并行：
  - R1-miR：miRNA 子表构建
  - R1-long：mRNA/transcript 子表构建
- 必须在两条管线之后做：
  - R1-merge：合并成 rna_master_v1.tsv
  - R1-QA：对 rna_master_v1.tsv 进行 QA + 写 README
    
    这样，你的执行顺序就很清晰：
    
1. 先写清楚表长什么样（R0）
2. 然后你可以今天摸 miRNA、明天摸 transcript，两条互不堵
3. 最后合并 + QA + 文档 → RNA v1 收尾
