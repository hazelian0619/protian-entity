# Protein L1（v6_clean）复盘 + 决策说明（工业化闭环 v1）

> 目标：Protein L1 的数据本体已经成熟，本阶段的主线是把它包装成“工业级交付”：
> **可复现（别人能跑）+ 可追溯（这份数据是哪份）+ 可回归（指标不莫名其妙退化）**。

## 1) 背景与主线对齐（为什么要做这件事）

在知识图谱的执行流里（见 `docs/KG_EXECUTION_FLOW.md`），L1 实体表是所有 join 的地基：

- `gene_master` 会从 Protein/RNA 汇总 gene 映射
- `edges_ppi` 会用 Protein 做 join coverage 门禁

Protein 主表本身即使“很好”，也会面临一个现实问题：

- 团队协作中，**每个人机器上拿到的文件可能不同**（版本、被替换、被改写）
- 下游出现差异时，如果没有“指纹 + QA”，很难定位到底是数据变了、代码变了还是环境变了

因此需要把 Protein L1 做成与 RNA 同范式的交付闭环：

1) contract（定义什么叫合格）
2) QA report（客观指标与 PASS/FAIL）
3) manifest（文件指纹 + 追溯）

## 2) 输入与产物（做什么）

### 输入（现成产物，不重建）

- `data/processed/protein_master_v6_clean.tsv`

### 输出（小文件，可提交 git）

- `pipelines/protein/reports/protein_master_v6.validation.json`
- `pipelines/protein/reports/protein_master_v6.manifest.json`

### 核心原则：不做什么

- 不联网、不下载任何数据
- 不修改 `protein_master_v6_clean.tsv`
- 不引入新的大文件进 git

> 这份 pipeline 的定位是“交付标准化与自证”，不是“再跑一次蛋白 ETL”。

## 3) v1 的口径（抓大放小的取舍）

本阶段只做“能挡住回归”的质量门禁，避免陷入过度工程：

### 3.1 强门禁（必须 100%）

- 主键：`uniprot_id` 非空率=100%，且唯一
- 基础有效性：`sequence` 非空率=100%，字符集合法
- 溯源锚点：`fetch_date` 非空且格式合法（`YYYY-MM-DD`）

### 3.2 覆盖率门禁（目标 >= 0.99）

- `alphafold_mean_plddt` 非空率
- `ensembl_gene_id` 非空率（并做格式 sanity）
- `ncbi_gene_id` 非空率（并做格式 sanity）

这些门禁的目标是“防退化”而不是“追求理论 100%”：

- 上游版本变化、清洗策略调整可能造成极小波动
- 如果门槛过严，会导致 pipeline 永远 FAIL，团队会倾向于绕开 QA（反而失去主线意义）

## 4) 执行闭环（如何复现）

最短命令：

```bash
cd /Users/pluviophile/graph/1218
bash pipelines/protein/run.sh
```

它内部做两件事：

1) `tools/kg_validate_table.py` 执行 contract，生成 validation.json
2) `tools/kg_make_manifest.py` 对输入表生成 manifest.json

可选：通过环境变量指定 manifest 的 data_version：

```bash
DATA_VERSION=protein-l1-v6-clean bash pipelines/protein/run.sh
```

## 5) 结果如何被使用（它如何服务主线）

这套交付的价值不在于“又多了两个 JSON”，而在于：

- **同事/下游**拿到一份 protein 表，先跑一键 QA：确认“能用”
- **对齐问题**出现时，用 manifest 的 sha256 快速确认“是不是同一份文件”
- **版本升级**时，用 QA 指标对比确认“质量没有悄悄退化”（coverage、重复、格式）

这会显著减少后续 gene_master/ppi 任务的沟通成本与排查成本。

## 6) 已知边界与后续 v2 方向

### 6.1 已知边界（v1 刻意不做）

- 跨列一致性（例如 `sequence_len == len(sequence)`）
- 数值范围约束（例如 `alphafold_mean_plddt` 0–100）
- 更强的溯源（例如强制 `source_version`）

原因：当前通用 validator 的 rule 类型有限，v1 的主线目标是“先把交付闭环跑通”。

### 6.2 v2 方向（如果要更强）

- 扩展 `tools/kg_validate_table.py` 支持更多 rule（range、cross-field check）
- 或在 `pipelines/protein/run.sh` 增加轻量附加检查（失败即 FAIL，并写入报告）
- 如果 Protein 未来版本频繁迭代且体积继续增长，考虑引入 Release 分发（像 RNA 一样），减少 repo 历史膨胀

