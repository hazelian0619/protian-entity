# RNA 组件核查矩阵（2026-04-12）

> 核查口径：
> - 方案 B 已生效：**mRNA canonical=ENST，miRNA canonical=URS**。
> - A 的 70%/90% 为 **KPI（默认非阻塞）**。
> - 本次只做现状重刷与缺口复盘，不新增大规模抓取。

## 0) 本地文件快照（实查）

| 文件 | 存在 | 大小(bytes) | 数据行数 | 证据 |
|---|---:|---:|---:|---|
| `data/output/rna_master_v1.tsv` | ✅ | 744,099,663 | 292,274 | `data/output/rna_master_v1.tsv` |
| `data/output/rna_xref_mrna_enst_urs_v2.tsv` | ✅ | 1,470 | 10 | `data/output/rna_xref_mrna_enst_urs_v2.tsv` |
| `data/output/rna_id_canonical_map_v1.tsv` | ✅ | 54,033,961 | 292,274 | `data/output/rna_id_canonical_map_v1.tsv` |
| `data/output/rna_structure_rfam_v1.tsv` | ✅ | 224,787 | 1,096 | `data/output/rna_structure_rfam_v1.tsv` |
| `data/output/rna_expression_evidence_v1.tsv` | ✅ | 51,857,745 | 289,614 | `data/output/rna_expression_evidence_v1.tsv` |
| `data/output/rna_rbp_sites_v1.tsv` | ✅ | 4,399,904 | 31,084 | `data/output/rna_rbp_sites_v1.tsv` |
| `data/output/rna_pdb_structures_v1.tsv` | ❌ | 0 | - | `data/output/rna_pdb_structures_v1.tsv` |
| `pipelines/rna_xref_enst_urs/reports/mrna_enst_urs_coverage_v2.json` | ✅ | 2,160 | - | `pipelines/rna_xref_enst_urs/reports/mrna_enst_urs_coverage_v2.json` |
| `pipelines/rna_xref_enst_urs/run.sh` | ✅ | 3,400 | - | `pipelines/rna_xref_enst_urs/run.sh` |

> 实查结论：`data/output/rna_pdb_structures_v1.tsv` 当前不存在；其余本次主线要求的 6 个 RNA 产物文件均存在。

## 1) 统一字段模板重刷（组件/字段/获取方式/状态/证据路径）

| 组件 | 字段 | 获取方式 | 状态 | 证据路径 |
|---|---|---|---|---|
| 序列主表 | rna_master_v1 (rna_id/rna_type/sequence/taxon_id 等) | RNA 主流程产物（merge 后主表） | 已落地；主表存在，且 row_count=292,274 | `data/output/rna_master_v1.tsv`<br>`pipelines/rna/contracts/rna_master_v1.json` |
| 方案B主键口径 | mRNA canonical=ENST；miRNA canonical=URS（主表 rna_id） | 按 rna_type 对 rna_id 正则核查 | 已满足：mRNA ENST 命中率=100%；miRNA URS 命中率=100% | `data/output/rna_master_v1.tsv`<br>`pipelines/rna/contracts/rna_master_v1.json` |
| mRNA ENST↔URS 映射 | rna_xref_mrna_enst_urs_v2 | rna_xref_enst_urs pipeline（id_mapping 流式 join） | 已落地；覆盖率0.0035%，未达70%/90% KPI；默认非阻塞 | `data/output/rna_xref_mrna_enst_urs_v2.tsv`<br>`pipelines/rna_xref_enst_urs/reports/mrna_enst_urs_coverage_v2.json`<br>`pipelines/rna_xref_enst_urs/run.sh` |
| RNA canonical 桥接 | rna_id_canonical_map_v1 | rna_master + xref 规则映射 | 已落地；非空/唯一校验通过；mRNA 中有10条 canonical_changed | `data/output/rna_id_canonical_map_v1.tsv`<br>`pipelines/rna_id_canonical/reports/rna_id_canonical_map_v1.build.json`<br>`pipelines/rna_id_canonical/reports/rna_id_canonical_map_v1.validation.json` |
| 结构补充（Rfam） | rna_structure_rfam_v1 (rfam_id, secondary_structure) | RNAcentral/Rfam 映射 + 启发式补全 | 已落地；row_count=1,096；contract通过 | `data/output/rna_structure_rfam_v1.tsv`<br>`pipelines/rna_rfam_structure/reports/rna_structure_rfam_v1.metrics.json`<br>`pipelines/rna_rfam_structure/reports/rna_structure_rfam_v1.validation.json` |
| 表达证据 | rna_expression_evidence_v1 | ENCODE expression 证据清洗并 join 到 rna_id | 已落地；row_count=289,614；join_rate=1.0 | `data/output/rna_expression_evidence_v1.tsv`<br>`pipelines/rna_expression_rbp/reports/rna_expression_rbp_v1.metrics.json`<br>`pipelines/rna_expression_rbp/reports/rna_expression_evidence_v1.validation.json` |
| RBP 位点证据 | rna_rbp_sites_v1 | ENCODE eCLIP 位点清洗并 join 到 rna_id | 已落地；row_count=31,084；join_rate=1.0 | `data/output/rna_rbp_sites_v1.tsv`<br>`pipelines/rna_expression_rbp/reports/rna_expression_rbp_v1.metrics.json`<br>`pipelines/rna_expression_rbp/reports/rna_rbp_sites_v1.validation.json` |
| RNA 3D 结构 | rna_pdb_structures_v1 | 应由 RNA PDB 子流程产出（当前仓库未见该流程） | 缺失：目标文件不存在 | `data/output/rna_pdb_structures_v1.tsv`<br>`docs/rna/DATA_DICTIONARY_RNA.md` |
| 主表 Tier2 结构列填充 | rna_master_v1.{rfam_id,secondary_structure,pdb_ids} | 主流程脚本当前保留空列 | 待策略化：主表中3列当前均为空（依赖补充表外连） | `data/output/rna_master_v1.tsv`<br>`pipelines/rna/scripts/04_finalize_mirna_table_v3.py`<br>`pipelines/rna/scripts/08_finalize_mrna_table.py` |
| 文档口径一致性 | DATA_DICTIONARY 主键口径 vs 方案B | 对比数据字典与现行 pipeline 口径 | 存在偏差：字典仍写“rna_id 必须 URS”，与 mRNA ENST 主线不一致 | `docs/rna/DATA_DICTIONARY_RNA.md`<br>`pipelines/rna/contracts/rna_master_v1.json`<br>`data/output/rna_master_v1.tsv` |

## 2) Top10 缺口（分级 + 估时）

| 排名 | 缺口 | 分级 | 优先级 | 工时 | 说明 | 证据 |
|---:|---|---|---|---|---|---|
| 1 | 缺少 RNA 3D 结构产物文件 data/output/rna_pdb_structures_v1.tsv | 可延期 | P1 | M | 当前无 RNA PDB 子流程与产物，3D 结构仍为空白能力。 | `data/output/rna_pdb_structures_v1.tsv`<br>`docs/rna/DATA_DICTIONARY_RNA.md` |
| 2 | mRNA ENST↔URS 覆盖率仅 0.0035%，未达 70%/90% KPI | 可延期 | P1 | M | 目前仅10条映射命中；默认非阻塞，仅在 STRICT_COVERAGE_GATE=1 时阻塞。 | `pipelines/rna_xref_enst_urs/reports/mrna_enst_urs_coverage_v2.json`<br>`pipelines/rna_xref_enst_urs/run.sh` |
| 3 | canonical map 中有 10 条 mRNA 被切换为 URS canonical | 可延期 | P1 | S | 与“mRNA canonical=ENST”协作口径存在边界歧义，需要明确策略。 | `data/output/rna_id_canonical_map_v1.tsv`<br>`pipelines/rna_id_canonical/reports/rna_id_canonical_map_v1.build.json` |
| 4 | DATA_DICTIONARY 仍要求 rna_id 全量 URS，与方案B不一致 | 可延期 | P1 | S | 文档与实际产物口径不统一，易造成后续实现误判。 | `docs/rna/DATA_DICTIONARY_RNA.md`<br>`data/output/rna_master_v1.tsv` |
| 5 | 主表 Tier2 结构列（rfam_id/secondary_structure/pdb_ids）仍全空 | 可延期 | P1 | M | 当前采用“主表空列+外部证据表”模式，需要明确消费约定。 | `data/output/rna_master_v1.tsv`<br>`data/output/rna_structure_rfam_v1.tsv` |
| 6 | miRNA xref 独立产物（rna_xref_mirna_v1.tsv）未物化 | 仅文档追踪 | P2 | M | 已有脚本历史实现，但当前主线未产出该文件。 | `data/output/rna_xref_mirna_v1.tsv`<br>`pipelines/rna/scripts/10_mirna_xref_and_conflicts.py` |
| 7 | lncRNA/transcript 类型尚未纳入 v1 主流程 | 仅文档追踪 | P2 | M | 当前主线以 mirna/mrna 为边界。 | `pipelines/rna/contracts/rna_master_v1.json`<br>`docs/rna/DATA_DICTIONARY_RNA.md` |
| 8 | covariance_model 字段与流程未落地 | 仅文档追踪 | P2 | L | 当前无 contract 字段与构建脚本。 | `pipelines/rna/contracts/rna_master_v1.json` |
| 9 | predicted_structure_ref 字段与流程未落地 | 仅文档追踪 | P2 | L | 预测结构仍处于后续规划。 | `pipelines/rna/contracts/rna_master_v1.json`<br>`docs/rna/DATA_DICTIONARY_RNA.md` |
| 10 | 旧文件名/旧报告口径仍在部分文档中残留（v1 xref/12-12 缺失结论） | 仅文档追踪 | P2 | S | 当前实际主线已切到 v2/v1 新产物，但旧描述未完全清理。 | `docs/rna/RNA_COMPONENT_CHECKLIST_2026-04-11.md`<br>`pipelines/rna/reports/rna_component_check_2026-04-11.json`<br>`data/output/rna_xref_mrna_enst_urs_v2.tsv` |

## 3) 结论（主线口径）

- **主线是否已解卡**：是（已解卡）。
- **当前唯一主线阻塞项（如果有）**：无。
- **A 低覆盖是否阻塞**：否（默认 `STRICT_COVERAGE_GATE=0` 时仅告警；仅 strict mode 才阻塞）。

结论证据：
- `data/output/rna_master_v1.tsv`
- `data/output/rna_xref_mrna_enst_urs_v2.tsv`
- `data/output/rna_id_canonical_map_v1.tsv`
- `data/output/rna_structure_rfam_v1.tsv`
- `data/output/rna_expression_evidence_v1.tsv`
- `data/output/rna_rbp_sites_v1.tsv`
- `pipelines/rna_xref_enst_urs/reports/mrna_enst_urs_coverage_v2.json`
- `pipelines/rna_xref_enst_urs/run.sh`
