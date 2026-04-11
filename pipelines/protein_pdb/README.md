# protein_pdb：PDB 结构详情补洗（v1）

把 `protein_master_v6_clean.tsv` 里的 `pdb_ids` 展开为结构化详情表：

- `pdb_id`
- `uniprot_id`
- `experimental_method`
- `resolution`
- `release_date`
- `ligand_count`
- `source`
- `fetch_date`

## 输入

- `data/processed/protein_master_v6_clean.tsv`（必需列：`uniprot_id`, `pdb_ids`）

## 输出

- 数据表：`data/output/protein/pdb_structures_v1.tsv`
- Contract：`pipelines/protein_pdb/contracts/pdb_structures_v1.json`
- 报告：`pipelines/protein_pdb/reports/*.json`

## 一键运行（先 smoke，再 full）

在仓库根目录执行：

```bash
cd /Users/pluviophile/graph/1218
bash pipelines/protein_pdb/run.sh
```

默认流程：
1. 先跑 smoke（`SMOKE_MAX_UNIQUE_PDB=200`）
2. 再跑全量
3. contract 校验 + QA + manifest

可选环境变量：

- `SMOKE_MAX_UNIQUE_PDB`：smoke 阶段唯一 PDB 数量（默认 `200`）
- `PDB_API_FAIL_THRESHOLD`：API 失败率中断阈值（默认 `0.05`）
- `PDB_BATCH_SIZE`：GraphQL 批次大小（默认 `10000`）
- `PDB_TIMEOUT`：单次 API 请求超时秒数（默认 `40`）
- `PDB_RETRIES`：API 重试次数（默认 `3`）
- `PDB_SLEEP_BETWEEN_BATCHES`：批次间 sleep（默认 `0.0`）
- `DATA_VERSION`：manifest 记录用数据版本（默认 `kg-data-local`）

## 验收口径

- **可追溯**：`(uniprot_id, pdb_id)` 必须能回溯到主表 `pdb_ids` 展开结果
- **分辨率合法范围**：QA gate 默认检查 `0.1 <= resolution <= 25.0`
- **行数一致性**：QA 依据 `build_report.selected_pairs` 校验（smoke 校验样本，full 校验全量）
- **API 审计**：`pdb_structures_v1.api_audit.json` 明确记录
  - missing（RCSB 返回未命中）
  - api_error（网络/限流/批次失败）

## API 批量失败中断机制

默认先走 RCSB API（GraphQL 批量，失败时回退 core/entry 单条接口）。

若 `api_error_unique_pdb / selected_unique_pdb > PDB_API_FAIL_THRESHOLD`，会：

1. 立即中断构建（退出码 3）
2. 在报告里写入手动下载清单、放置路径、SHA256 校验命令：
   - `pipelines/protein_pdb/reports/pdb_structures_v1.api_audit.json`
   - `pipelines/protein_pdb/reports/pdb_structures_v1.build.json`

## 文件说明

- `scripts/01_build_pdb_structures_v1.py`：展开主表并抓取 PDB 详情
- `scripts/02_qa_pdb_structures_v1.py`：追溯性 + 分辨率范围 + API 失败率 QA
- `contracts/pdb_structures_v1.json`：表结构与基础规则
- `run.sh`：标准入口（smoke → full → validate → QA → manifest）
