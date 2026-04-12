# RNA 组件核查矩阵（2026-04-11）

> 核查口径：仅做“现状证据化 + 缺口分级”，不新增大数据抓取。
> 
> 说明：当前工作区缺少 RNA 原始输入与 RNA 产物文件，以下“当前状态”区分为：
> - **代码层已实现**：脚本/contract 可证
> - **产物层未验证**：本地 `data/output/rna*` 缺失

## 0) 本地数据可用性快照（阻塞项）

- 缺失：`data/raw/rna/rnacentral/id_mapping.tsv(.gz)`
- 缺失：`data/raw/rna/mirbase/mature.fa.gz`
- 缺失：`data/raw/rna/ensembl/Homo_sapiens.GRCh38.cdna.all.fa`
- 缺失：`data/output/rna_master_v1.tsv`、`data/output/rna_xref_mrna_enst_urs_v1.tsv` 等 RNA 产物

证据：`pipelines/rna/reports/rna_component_check_2026-04-11.json` 的 `local_file_checks`

---

## 1) 统一模板核查表

| 组件 | 属性字段 | 数据类型 | 获取方式 | 存储形式 | 示例 | 用途/指标 | 指标 | 层次结构 | 语义关系 | 数据来源 | 维护机构 | 更新频率 | 当前状态 | 证据路径 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 序列数据 | RNA 序列 (`sequence`) | String | miRNA: `mature.fa` + RNAcentral 映射；mRNA: Ensembl cDNA FASTA 解析 | `rna_master_v1.tsv` 列 | `UAGCUUAUC...` | 序列完整性与下游比对 | 非空率=100%，字符集=ACGUN | L1 实体主表 | `rna_id`-has_sequence->`sequence` | RNAcentral/miRBase/Ensembl | 文档未显式登记 | 随版本重跑 | **[P0][阻塞级] 代码层已实现，产物层未验证（本地无 RNA 输出）** | `pipelines/rna/scripts/04_finalize_mirna_table_v3.py:96-99,123-127`; `pipelines/rna/scripts/08_finalize_mrna_table.py:103-106,129-134`; `pipelines/rna/contracts/rna_master_v1.json:51-61` |
| 序列数据 | URS 主键一致性 (`rna_id`) | String(PK) | miRNA 用 URS；mRNA 当前由 ENST 派生 | 主表 `rna_id` | `URS00000478B7_9606` / `ENST..._9606` | 全局主键统一与跨库 join | 目标：主表统一 canonical URS | L1 主键层 | legacy_id -> canonical_id | RNAcentral + Ensembl | 文档未显式登记 | 随版本重跑 | **[P0][阻塞级] 规范要求 URS，但 mRNA 实现仍为 ENST_9606 临时键** | `docs/rna/DATA_DICTIONARY_RNA.md:3-5,11`; `pipelines/rna/scripts/07_extract_mrna_from_ensembl.py:173-180`; `pipelines/rna/contracts/rna_master_v1.json:91-96`; `pipelines/rna/scripts/11_mrna_enst_to_urs_xref.py:96-109` |
| 序列数据 | RNA 类型 (`rna_type`) | Enum | 由脚本写入固定类型 | 主表 `rna_type` | `mirna`,`mrna` | 类型分层统计 | 枚举约束通过 | L1 分类层 | RNA subclass-of RNA type | RNAcentral/Ensembl | 文档未显式登记 | 随版本重跑 | **[P1][可延期] 已实现 mirna/mrna；lncrna 未落地** | `pipelines/rna/contracts/rna_master_v1.json:44-49`; `docs/rna/DATA_DICTIONARY_RNA.md:12` |
| 序列数据 | 物种来源 (`taxon_id`) | Integer | 脚本硬编码 9606 + contract 校验 | 主表 `taxon_id` | `9606` | 物种一致性 | 人类占比=100% | L1 元数据层 | RNA belongs_to taxon | RNAcentral/Ensembl | 文档未显式登记 | 随版本重跑 | **[P0][阻塞级] 代码层就绪；需产物验证 100% 覆盖** | `pipelines/rna/scripts/04_finalize_mirna_table_v3.py:100`; `pipelines/rna/scripts/08_finalize_mrna_table.py:106`; `pipelines/rna/contracts/rna_master_v1.json:37-42` |
| 结构数据 | 二级结构 (`secondary_structure`) | Dot-bracket/String | 当前未抓取，仅预留 | 主表预留空列 | `(((...)))` | 结构层分析 | 非空率目标（v2） | L1 扩展列 | RNA-has_secondary_structure | Rfam（计划） | 文档未显式登记 | 随 Rfam 扩展 | **[P1][可延期] 字段存在但脚本统一置空** | `pipelines/rna/scripts/04_finalize_mirna_table_v3.py:111-113`; `pipelines/rna/scripts/08_finalize_mrna_table.py:117-119`; `docs/rna/DATA_DICTIONARY_RNA.md:47-49` |
| 结构数据 | Rfam 家族 (`rfam_id`) | String | 当前未抓取，仅预留 | 主表预留空列 | `RF00001` | 家族归类/结构推断 | 格式校验 `RF\d+`（目标） | L1 扩展列 | RNA-member_of-Rfam | Rfam | 文档未显式登记 | 随 Rfam 扩展 | **[P1][可延期] 字段存在但置空** | `pipelines/rna/scripts/04_finalize_mirna_table_v3.py:111`; `pipelines/rna/scripts/08_finalize_mrna_table.py:117`; `docs/rna/RNA_SOURCES_AND_VERSIONS.md:10,23` |
| 结构数据 | 协方差模型 (`covariance_model`) | String/ModelRef | 未实现 | 无 | `RFxxxxx.cm` | 家族级结构比对 | 覆盖率/命中阈值（待定义） | L2 结构模型层 | RNA-scored_by-CM | Rfam CM（计划） | 文档未显式登记 | 未定义 | **[P2][仅文档追踪] 仓库无字段、无脚本、无 contract** | `pipelines/rna/contracts/rna_master_v1.json`（无该字段）；`pipelines/rna/scripts`（无 CM 步骤） |
| 结构数据 | 三维结构 (`pdb_ids`) | String(list) | 当前未抓取，仅预留 | 主表预留空列 | `1A34;2KOC` | 3D 结构引用 | 格式校验（目标） | L1 扩展列 | RNA-mapped_to-PDB | PDB/wwPDB（计划） | 文档未显式登记 | 未定义 | **[P1][可延期] 字段存在但置空** | `pipelines/rna/scripts/04_finalize_mirna_table_v3.py:113`; `pipelines/rna/scripts/08_finalize_mrna_table.py:119`; `docs/rna/DATA_DICTIONARY_RNA.md:49` |
| 结构数据 | 预测结构 (`predicted_structure_ref`) | String/URI | 未实现 | 无 | `rhofold://...` | 结构补全 | 覆盖率（待定义） | L2 预测层 | RNA-predicted_structure->model | 预测模型（计划） | 文档未显式登记 | 未定义 | **[P2][仅文档追踪] 无字段、无脚本、无规范** | `pipelines/rna/contracts/rna_master_v1.json`（无该字段）；`docs/rna/*`（无预测结构规范） |
| 功能注释 | miRNA 注释 (`mirbase_id`,`mimat_id/xref`) | String | `mature.fa` + `id_mapping`；附加 xref 脚本 | 主表列 + xref 表 | `hsa-miR-21-5p`,`MIMAT0000076` | miRNA 精确检索 | `mirbase_id` 非空率（目标） | L1 功能标签层 | miRNA-linked_to-miRBase ID | RNAcentral/miRBase | 文档未显式登记 | 随版本重跑 | **[P0][阻塞级] 主表字段已实现；xref 仅有脚本未纳入主 run/contract** | `pipelines/rna/scripts/04_finalize_mirna_table_v3.py:97,109`; `pipelines/rna/scripts/10_mirna_xref_and_conflicts.py:50-63`; `pipelines/rna/run.sh`（未调用 script10） |
| 功能注释 | lncRNA 注释 | String/Enum | 未实现 | 无 | `lncrna` | 功能层扩展 | 类型覆盖率（待定义） | L1 类型扩展 | RNA-instance_of-lncRNA | Ensembl/RNAcentral（计划） | 文档未显式登记 | 未定义 | **[P1][可延期] 文档提及但 contract 与 run 未落地** | `docs/rna/DATA_DICTIONARY_RNA.md:12`; `pipelines/rna/contracts/rna_master_v1.json:44-49`; `pipelines/rna/run.sh:8-18` |
| 功能注释 | 表达与证据 | Numeric/Text | 未抓取（文档明确暂不抓） | 无 | `TPM`,`evidence_level` | 可解释性与组织特异性 | 覆盖率/证据分级（待定义） | L2 证据层 | RNA-supported_by-expression evidence | ENCODE 等（计划） | 文档未显式登记 | 未定义 | **[P2][仅文档追踪] 当前明确 out-of-scope** | `docs/rna/DATA_DICTIONARY_RNA.md:56` |
| 类型特征 | tRNA 反密码子 | String | 未实现 | 无 | `GAA` | 小 RNA 类型特征 | 非空率（待定义） | L2 特征层 | tRNA-has_anticodon | RNAcentral/Rfam（计划） | 文档未显式登记 | 未定义 | **[P2][仅文档追踪] 当前 out-of-scope** | `docs/rna/DATA_DICTIONARY_RNA.md:58`; `pipelines/rna/run.sh:8-18` |
| 类型特征 | rRNA 位点 | String/Coordinate | 未实现 | 无 | `18S:...` | rRNA 定位特征 | 覆盖率（待定义） | L2 特征层 | rRNA-has_site | RNAcentral/Rfam（计划） | 文档未显式登记 | 未定义 | **[P2][仅文档追踪] 当前 out-of-scope** | `docs/rna/DATA_DICTIONARY_RNA.md:58`; `pipelines/rna/run.sh:8-18` |
| 互作锚点 | RBP 位点 | BED-like/Coordinate | 未实现 | 无 | `ELAVL1@chr...` | RNA-蛋白互作锚点 | 位点覆盖率（待定义） | L2 关系锚点层 | RNA-binds-RBP@site | ENCODE/eCLIP（计划） | 文档未显式登记 | 未定义 | **[P2][仅文档追踪] 当前 out-of-scope（RPI 不进主表）** | `docs/rna/DATA_DICTIONARY_RNA.md:55` |
| 标识映射 | miRNA 外部 ID 映射（miRBase） | TSV relation | 专门脚本生成 xref + 冲突审计 | `rna_xref_mirna_v1.tsv`（计划） | `URS... ↔ MIMAT...` | 跨库检索与冲突审计 | 一对多冲突可审计 | L1.5 映射层 | RNA-xref-miRBase | RNAcentral/miRBase | 文档未显式登记 | 随版本重跑 | **[P1][可延期] 脚本存在但未并入 run 与 contract；本地无产物** | `pipelines/rna/scripts/10_mirna_xref_and_conflicts.py:8-14,50-89`; `pipelines/rna/run.sh:8-21`; `pipelines/rna/contracts/`（无 xref contract） |
| 标识映射 | mRNA ENST↔URS 映射 | TSV relation | 可选脚本扫描 `id_mapping` 生成 xref/coverage 文本 | `rna_xref_mrna_enst_urs_v1.tsv` | `ENST..._9606 ↔ URS..._9606` | 主键统一与兼容层 | 覆盖率目标建议 ≥70% | L1.5 映射层 | legacy ENST -> canonical URS | RNAcentral | 文档未显式登记 | 随版本重跑 | **[P0][阻塞级] run 中为 optional 且无 contract/validation，主键统一风险高** | `pipelines/rna/run.sh:20-21`; `pipelines/rna/scripts/11_mrna_enst_to_urs_xref.py:8-18,96-127`; `pipelines/rna/contracts/`（无 xref contract） |

---

## 2) Top10 缺口清单（含分级与工作量）

| 排名 | 缺口 | 缺失项分级 | 优先级 | 预计工作量 | 说明 | 证据 |
|---|---|---|---|---|---|---|
| 1 | mRNA 主键仍为 ENST 派生，未统一 URS canonical | 阻塞级 | P0 | M | 与 URS 规范冲突，影响跨库实体一致性 | `docs/rna/DATA_DICTIONARY_RNA.md:3-5`; `pipelines/rna/scripts/07_extract_mrna_from_ensembl.py:177-180` |
| 2 | ENST↔URS 映射链路仅 optional，缺 contract+validation | 阻塞级 | P0 | M | 关键映射未纳入强制门禁 | `pipelines/rna/run.sh:20-21`; `pipelines/rna/scripts/11_mrna_enst_to_urs_xref.py:96-127` |
| 3 | 当前工作区缺 RNA 原始输入，无法重跑/验证 | 阻塞级 | P0 | S | 缺 `id_mapping`/`mature.fa`/`cdna.fa` | `pipelines/rna/reports/rna_component_check_2026-04-11.json` `local_file_checks` |
| 4 | 当前工作区缺 RNA 输出产物，覆盖率不可审计 | 阻塞级 | P0 | S | 无 `rna_master_v1.tsv` 等，现状仅能代码层审计 | 同上 |
| 5 | `rfam_id`、`secondary_structure` 仅预留空列 | 可延期 | P1 | M | 结构层空白，影响家族/结构分析 | `pipelines/rna/scripts/04_finalize_mirna_table_v3.py:111-113`; `pipelines/rna/scripts/08_finalize_mrna_table.py:117-119` |
| 6 | 协方差模型（CM）无字段无流程 | 可延期 | P1 | M | 无法支撑 Rfam/CM 证据链 | `pipelines/rna/contracts/rna_master_v1.json`（无 CM 字段） |
| 7 | 3D 结构映射（PDB/mmCIF）未实现 | 可延期 | P1 | M | `pdb_ids` 仅占位 | `docs/rna/DATA_DICTIONARY_RNA.md:49`; 脚本置空证据同上 |
| 8 | lncRNA 类型未落地（仅文档提及） | 可延期 | P1 | M | 当前 contract 只允许 mirna/mrna | `docs/rna/DATA_DICTIONARY_RNA.md:12`; `pipelines/rna/contracts/rna_master_v1.json:44-49` |
| 9 | 表达证据层未建设 | 仅文档追踪 | P2 | L | 文档明确暂不抓取 | `docs/rna/DATA_DICTIONARY_RNA.md:56` |
| 10 | RBP/tRNA/rRNA 特征与锚点未建设 | 仅文档追踪 | P2 | L | 当前版本边界外，需后续 pipeline | `docs/rna/DATA_DICTIONARY_RNA.md:55,58` |

---

## 3) 缺外部数据时的下载清单（按协作约束）

> 当前已触发“缺外部数据”条件，先给下载清单，不进行大数据抓取。

| 文件 | 建议版本 | 放置路径 | 校验方式（sha256） |
|---|---|---|---|
| RNAcentral `id_mapping.tsv.gz` | Release 25 | `data/raw/rna/rnacentral/id_mapping.tsv.gz` | `sha256sum data/raw/rna/rnacentral/id_mapping.tsv.gz` |
| miRBase `mature.fa.gz` | 22.1（通过 RNAcentral 口径） | `data/raw/rna/mirbase/mature.fa.gz` | `sha256sum data/raw/rna/mirbase/mature.fa.gz` |
| Ensembl `Homo_sapiens.GRCh38.cdna.all.fa.gz` | 与当前 Ensembl release 对齐 | `data/raw/rna/ensembl/Homo_sapiens.GRCh38.cdna.all.fa.gz` | `sha256sum data/raw/rna/ensembl/Homo_sapiens.GRCh38.cdna.all.fa.gz` |

来源锚点：`docs/rna/RNA_SOURCES_AND_VERSIONS.md:7-10,35-46`
