# Protein L1（v6_clean）：工业化入口（QA + Manifest）

这条 pipeline 面向的不是“重建 protein 主表”，而是把现成产物标准化成可复现、可追溯的数据交付。

**输入（已存在）**：`data/processed/protein_master_v6_clean.tsv`

**输出（小文件，可进 git）**：
- `pipelines/protein/reports/protein_master_v6.validation.json`：质量门禁（PASS/FAIL + 指标）
- `pipelines/protein/reports/protein_master_v6.manifest.json`：指纹与追溯（sha256/size/mtime/git commit）

> 目的：后续的 `gene_master`、`edges_ppi` 等下游任务能“用同一套方式”验证 Protein L1 是否合格、是否是同一份数据。

## 最短复现（推荐）

在仓库根目录运行：

```bash
cd /Users/pluviophile/graph/1218
bash pipelines/protein/run.sh
```

成功标志：终端出现 `[PASS] protein_master_v6 ...`，并生成上述 `reports/*.json`。

## Contract（数据承诺）是什么？

- Contract 文件：`pipelines/protein/contracts/protein_master_v6.json`
- 执行工具：`python3 tools/kg_validate_table.py`

Contract 的定位是“可执行的数据承诺”（schema + 质量门槛），用于防回归，而不是追求理论最严。

## 常见卡点 / 隐藏点（先看这里）

1) **我以为会生成新的 protein 主表**：不会。
   - 本 pipeline 不做 ETL，不会写 `data/output/`。
   - 它只对现有 `data/processed/protein_master_v6_clean.tsv` 做 QA + manifest。

2) **我跑了但找不到输出**：请确认你从 repo root 跑，且查看目录：
   - `pipelines/protein/reports/`

3) **FAIL 了怎么定位**：
   - 先打开 `pipelines/protein/reports/protein_master_v6.validation.json`
   - 找 `passed=false` 的 rule id；它就是失败原因。

4) **为什么没有 `taxon_id`/`source_version` 之类字段**：
   - 以当前主表的现实字段为准（contract 不强行要求不存在的列）。
   - 溯源字段以 `source` + `fetch_date` 为主。

5) **manifest 里 `git_commit` 为空**：
   - 少数环境下可能无法读取 git commit（例如不在 git 仓库中运行）。
   - 不影响 sha256/size/mtime 作为文件指纹的用途。

## 进一步阅读（维护者/传承）

- 复盘与决策记录：`pipelines/protein/POSTMORTEM_protein_master_v6.md`

