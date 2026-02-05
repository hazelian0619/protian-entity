# 人类蛋白质与 RNA 知识图谱（工业级数据产品）

本仓库以工业级数据产品的实践组织：可审计的 ETL 代码、数据契约（contracts）、以及 QA 报告；大型 L1 数据以 GitHub Releases 发布并附带 `manifest.json` 与校验报告。

快速索引
- 核心蛋白质表（L1）：`data/processed/protein_master_v6_clean.tsv`（v6 快照）
- RNA（L1）数据以 Release 发布：`rna-l1-v1`（包含 `.tsv.gz`、`manifest.json`、QA 报告）
- 验证工具：`tools/kg_validate_table.py`

快速开始
```bash
# 验证 protein 主表并生成报告
python3 tools/kg_validate_table.py \
  --contract pipelines/protein/contracts/protein_master_v6.json \
  --table data/processed/protein_master_v6_clean.tsv \
  --out build/validate/protein_master_v6_report.json
```

CI 与非交互运行
- CI 作业应为非交互模式（`--yes` / `--ci` 或 `CI=true`）。
- 管理员可参考 `docs/GITHUB_ACTIONS_SETTINGS.md` 配置仓库 Actions 权限（放行 fork 工作流会带来安全风险，请谨慎）。

数据发布流程
1. 使用 `pipelines/` 生成 L1 数据与 `manifest.json`（含 row counts、checksums、commit SHA、build timestamp）。
2. 运行全部验证合同并保存报告。
3. 建立 GitHub Release，上传数据包与 `manifest.json`、QA 报告。

贡献与治理
- 请参阅 `CONTRIBUTING.md` 了解分支、PR、验证与发布要求。提交修改涉及数据表或合同时请附带验证报告。

许可
- MIT — 详见 `LICENSE`

联系方式
- 仓库所有者：`@hazelian0619`，或在仓库中打开 issue / PR
