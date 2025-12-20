# 仓库使用说明（这仓库到底怎么用）

这不是“纯代码仓库”，而是一个 **工业化数据产品仓库**：

- **代码 / 规范 / QA 报告 / manifest** 进入 git（可审计、可复现、可回归）
- **大体积数据产物**通过 GitHub Releases 分发（可下载、可校验、可回滚）

项目大纲（`check/项目大纲.pdf`）的核心目标是：

- 三类实体（节点）：Protein / RNA / Molecules（小分子）
- 三类关系（边）：PPI / PSI / RPI
- 分层推进：L1（实体标准化）→ L2（关系证据化）→ L3（本体语义）→ L4（路径推理）

本仓库主要覆盖的是 **L1/L2 的“数据交付”层**：把表做成可复现的标准件。

---

## 1) 你是谁：按角色选择最短路径

### A. 我只是要用数据（同事/下游）

- Protein（L1）：直接读 `data/processed/protein_master_v6_clean.tsv`
  - 可选自证：`bash pipelines/protein/run.sh`（生成/更新 QA + manifest）
- RNA（L1, v1）：优先下载 Release `rna-l1-v1`（根 `README.md` 有入口）
- Molecules（L1+PSI L2）：优先下载 Release（根 `README.md` 有入口）

这类使用者的关键需求：

- “数据在哪”
- “怎么确认我拿到的是合格版本、且和别人一致”

因此每条 pipeline 都会提供 QA 报告与 manifest。

### B. 我要本地重跑 / 复现（开发者）

你需要进入 `pipelines/<name>/README.md`，按该 pipeline 的说明准备 `data/raw/` 输入，然后运行：

- `bash pipelines/rna/run.sh`
- `bash pipelines/molecules/run.sh`

注意：本地跑出来的大文件通常落在 `data/output/`，默认不提交到 git。

### C. 我要发布一个版本（维护者/发版）

你需要确保一个版本具备完整闭环：

1) contract（合约）
2) QA report（门禁报告）
3) manifest（产物指纹）
4) postmortem（复盘：口径/取舍/已知边界）

大文件产物（TSV/SQLite）走 GitHub Releases；仓库里保留可审计的“交付说明 + 报告 + 指纹”。

---

## 2) 目录结构：只记住这 5 个目录

### `data/`（数据分层）

- `data/raw/`：外部下载的原始输入（大多不进 git）
- `data/output/`：本地跑 pipeline 的重建产物（大文件，不进 git）
- `data/processed/`：仓库内直接可用的交付表（Protein 当前主要在这里）

### `pipelines/`（每条数据产品的统一入口）

每条 pipeline 都按同一范式组织：

- `run.sh`：一键入口（能复现/能生成报告）
- `contracts/*.json`：数据合约（机器可读规则）
- `reports/*.json`：QA/manifest（小文件，可进 git，用于回归）
- `README.md`：同事手册（怎么用/怎么跑/怎么排错）
- `POSTMORTEM_*.md`：维护者复盘（为什么这么做、卡点、v2方向）

### `tools/`（通用 QA/manifest 工具）

- `tools/kg_validate_table.py`：按 contract 验证 TSV，输出 `validation.json`
- `tools/kg_make_manifest.py`：生成 `manifest.json`（sha256/size/mtime/git_commit）

### `docs/`（规范与字典）

- `docs/KG_EXECUTION_FLOW.md`：为什么要 contract→pipeline→QA→manifest
- `docs/QUALITY_GATES.md`：质量门禁原则
- `docs/*DATA_DICTIONARY*`：字段字典

### `scripts/`（历史 ETL/一次性脚本）

这里更多是早期构建/排查脚本，不作为同事入口。
稳定入口应当以 `pipelines/*` 为准。

---

## 3) 项目大纲如何映射到仓库实现（L1–L4）

| 层级 | 大纲目标 | 本仓库对应交付物 | 当前状态 |
|---|---|---|---|
| L1 实体标准化 | Protein/RNA/Molecules 节点表 | `data/processed/*` 或 Release + `pipelines/*` | Protein ✅ / RNA ✅ / Molecules ✅ |
| L2 关系证据化 | PPI/PSI/RPI（边 + 证据 + 评分） | `pipelines/<edges>/` 产出 edges+evidence + QA/manifest | PSI ✅（molecules/M3）；PPI/RPI ⏳ |
| L3 语义对齐 | GO/KEGG/ChEBI 等本体挂接与约束 | 术语表/映射表 + contract/QA | ⏳ |
| L4 路径推理 | 多跳路径、查询模板 | 图数据库/服务层 | ⏳ |

说明：本仓库当前的“主线价值”是把 L1/L2 的表交付做成工业化闭环；L3/L4 往往需要另一个 serving 层（Neo4j/RDF/API）。

---

## 4) 目前已实现的三条主线交付

### Protein（L1）

- 数据：`data/processed/protein_master_v6_clean.tsv`
- 工业化入口：`pipelines/protein/README.md`
- 复盘：`pipelines/protein/POSTMORTEM_protein_master_v6.md`

### RNA（L1, v1）

- 入口：`pipelines/rna/README.md`
- 规范：`docs/rna/README.md`
- 复盘：`pipelines/rna/POSTMORTEM_rna_master_v1.md`

### Molecules（L1 + PSI L2）

- 入口：`pipelines/molecules/README.md`
- 规范：`docs/molecules/README.md`
- 复盘：`pipelines/molecules/POSTMORTEM_molecules_m1_v1.md` / `pipelines/molecules/POSTMORTEM_molecules_m2_v1.md` / `pipelines/molecules/POSTMORTEM_molecules_m3_v1.md`

---

## 5) 下一步还缺什么（按大纲主线）

最关键的缺口不是“再多一张表”，而是把剩下的边也做成同样的工业化交付：

1) **PPI（L2）标准交付**：把 `protein_edges.tsv` 变成 edges 表 + evidence 表（含评分与溯源），并配 contract/QA/manifest/postmortem
2) **RPI（L2）标准交付**：RNA–Protein interaction（边 + 证据）
3) **Gene 底座表（对齐层）**：`gene_master`（跨 Protein/RNA 的 xrefs 底座表），用于统一 gene_key 与覆盖率回归
4) **L3 本体对齐**：GO/KEGG/ChEBI 的术语表与映射表（先做最小可用）

