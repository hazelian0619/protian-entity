# 1218 数据集主线进度复盘（给同事也能看懂）

更新时间：2025-12-20

这份文档的目标只有一个：

> **用“人话”解释清楚：这个 GitHub 仓库是干什么的、怎么用、我们做到哪了、还差什么。**

它对齐 `check/项目大纲.pdf` 的主线：

- 三类实体（节点）：Protein / RNA / Molecules（小分子）
- 三类关系（边）：PPI / PSI / RPI
- 分层推进：L1（实体标准化）→ L2（关系证据化）→ L3（本体语义）→ L4（路径推理）

本仓库当前主要交付 **L1/L2 的数据产品层**：把表做成可复现的标准件。

---

## 0) 先讲清楚：这仓库的“使用心智模型”

你可以把这个仓库理解为一句话：

> **仓库里放“规则与证明”，大数据放 Release。**

具体来说：

- 仓库（git）保存：pipeline（怎么跑）+ contract（什么叫合格）+ QA reports（是否合格）+ manifest（文件指纹）+ postmortem（复盘/口径/卡点）
- GitHub Releases 保存：真正的大体积产物（`.tsv.gz` / `.sqlite` 等），并附带 manifest/QA 便于校验与回滚

如果你只记住一个入口：先看 `docs/REPO_USAGE.md`。

---

## 1) 3 分钟怎么用（按角色给最短路径）

### 1.1 我只是想用数据（不跑代码）

- Protein：直接读仓库内 `data/processed/protein_master_v6_clean.tsv`
- RNA：直接下载 Release `rna-l1-v1`（根 `README.md` 有链接）
- Molecules：直接下载 Release（根 `README.md` 有链接）

你只要会做一件事：用 manifest 校验“我拿到的是同一份文件”。

### 1.2 我要验证数据是否合格（最推荐的自检方式）

- Protein：`bash pipelines/protein/run.sh`（只生成/更新 `reports/*.json` 小文件）

验证失败时不要猜：直接看 `validation.json` 里哪条 rule `passed=false`。

### 1.3 我是开发者，需要本地重跑（会用到 raw 输入）

- RNA：准备 `data/raw/rna/`（见 `docs/rna/RNA_SOURCES_AND_VERSIONS.md`）→ `bash pipelines/rna/run.sh`
- Molecules：准备 `data/raw/molecules/`（见 `pipelines/molecules/README.md`）→ `bash pipelines/molecules/run.sh`

注意：本地重跑生成的大文件一般写到 `data/output/`，默认不提交 git。

---

## 2) 5 分钟汇报版（抓主线）

**主线目标（对齐 check 文档）**

- 知识图谱 = 三类实体（Protein / RNA / Small molecule）+ 三类关系（PPI / PSI / RPI）+ 边的证据层（provenance + score + method）。
- 工程落地采用“工业级数据产品”方式：**pipeline + contract + QA reports + manifest + Release**。

**当前进度（按 L 层）**

- L1（实体）：Protein ✅、RNA ✅、Molecules ✅（交付形态齐全；重跑时依赖 raw 输入）
- L2（关系+证据）：PSI ✅（小分子/ChEMBL M3 已有 evidence + gates），PPI/RPI ⏳（需要补成标准交付）
- L3（本体语义）：⏳（GO/KEGG/ChEBI 的术语表/层级/约束尚未产品化）
- L4（路径推理）：⏳（通常是图数据库/服务层，不在当前仓库主线范围内）

**一句话总结“我是如何洗这些数据集的”**

- 我们不把“清洗”当成一次性的脚本，而是把每张表当成可发布的数据产品：
  - 固定 schema（contract）
  - 自动化质量门禁（QA reports：PASS/FAIL + 指标）
  - 可追溯指纹（manifest：sha256/size/mtime/git commit）
  - 大产物走 GitHub Releases（避免仓库膨胀、可回滚）

**统一清洗/交付流程（所有数据集都按这个模板走）**

```
raw inputs → pipeline scripts → output tables → contract validation (PASS/FAIL) → manifest → GitHub Release
             (可复现)            (大文件)        (小 JSON 可审计)           (可追溯)
```

---

## 3) 从哪里查看（GitHub 与本地的对应关系）

### 1.1 GitHub 仓库入口（链接）

（链接由本地 `git remote show origin` 推导；本环境不直接联网抓取 GitHub 页面内容）

- Repo 首页：`https://github.com/hazelian0619/protian-entity`
- Releases：`https://github.com/hazelian0619/protian-entity/releases`
- RNA L1 v1 Release（已在仓库 README 中引用）：`https://github.com/hazelian0619/protian-entity/releases/tag/rna-l1-v1`

### 1.1.1 直接可点的“汇报用”链接（建议复制到飞书/Notion）

- 总入口：`https://github.com/hazelian0619/protian-entity/blob/main/README.md`
- 执行流（L1–L6）：`https://github.com/hazelian0619/protian-entity/blob/main/docs/KG_EXECUTION_FLOW.md`
- 质量门禁口径：`https://github.com/hazelian0619/protian-entity/blob/main/docs/QUALITY_GATES.md`
- 发布原则：`https://github.com/hazelian0619/protian-entity/blob/main/docs/DATA_RELEASES.md`

- 仓库怎么用（建议同事第一篇就看这个）：`https://github.com/hazelian0619/protian-entity/blob/main/docs/REPO_USAGE.md`

- Protein L1 工业化入口（分支）：`https://github.com/hazelian0619/protian-entity/tree/feat/protein-l1-industrial/pipelines/protein`
- Protein L1 QA 报告（PASS/FAIL + 指标）：`https://github.com/hazelian0619/protian-entity/blob/feat/protein-l1-industrial/pipelines/protein/reports/protein_master_v6.validation.json`
- Protein L1 manifest（sha256 + git commit）：`https://github.com/hazelian0619/protian-entity/blob/feat/protein-l1-industrial/pipelines/protein/reports/protein_master_v6.manifest.json`
- RNA L1 pipeline（已在 main）：`https://github.com/hazelian0619/protian-entity/tree/main/pipelines/rna`
- RNA L1 Release：`https://github.com/hazelian0619/protian-entity/releases/tag/rna-l1-v1`

- Molecules pipeline（分支）：`https://github.com/hazelian0619/protian-entity/tree/feat/protein-l1-industrial/pipelines/molecules`
- Molecules docs（分支）：`https://github.com/hazelian0619/protian-entity/tree/feat/protein-l1-industrial/docs/molecules`
- Molecules Releases：
  - L1（M1+M2）：`https://github.com/hazelian0619/protian-entity/releases/tag/molecules-l1-v1`
  - PSI L2（M3）：`https://github.com/hazelian0619/protian-entity/releases/tag/molecules-psi-l2-v1`

### 1.2 GitHub 上建议的浏览路径（给同事）

1) 先看 `README.md`（总入口，告诉你数据在哪里、哪些走 Release）
2) 再看 `docs/`（执行流、质量门禁、发布原则）
3) 最后看 `pipelines/<name>/`（可复现的具体实现：run.sh、contracts、reports、postmortem）

### 1.3 分支与“为什么我在 GitHub 看不到某些东西”

根据本地仓库信息：

- 远端默认分支：`main`
- 远端已存在分支：`feat/protein-l1-industrial`、`feat/rna-l1-v1`

说明：

- `main` 是“同事默认看到的内容”。
- 新增的 Protein 工业化与 Molecules pipeline 当前在 `feat/protein-l1-industrial` 分支（PR 合并后才会进入 `main`）。

如果同事在 GitHub 只看 `main`：

- 只会看到已经合并进 `main` 的内容
- Feature 分支里的 pipeline/文档需要切换 branch 才能看到

---

## 4) 数据产品地图（我们到底产出了哪些“可用东西”）

下面按“表/产物”列出：它属于哪一层、来源是什么、产出在哪里、同事怎么验证。

### 4.1 L1：Protein 实体（v6_clean）

**定位**：Protein 实体主表（L1），为下游 Gene spine、PPI、PSI 等提供端点与属性。

- 主表（已在仓库内）：`data/processed/protein_master_v6_clean.tsv`
- 工业化入口（QA + manifest，一键验证）：
  - `pipelines/protein/run.sh`
  - `pipelines/protein/contracts/protein_master_v6.json`
  - `pipelines/protein/reports/protein_master_v6.validation.json`
  - `pipelines/protein/reports/protein_master_v6.manifest.json`

**“怎么理解清洗”**：

- 这条 pipeline 不重建主表（主表已存在且体积适中）；它把“这份表是否合格、是否发生回归、是否是同一份数据”用 contract+manifest 固化下来。
- 典型门禁：`uniprot_id` 唯一/非空、`sequence` 非空与字符集、溯源字段 `source/fetch_date`、关键覆盖率（如 `alphafold_mean_plddt`）。

**关键规模（用于汇报一句话）**：`protein_master_v6_clean.tsv` = 19,135 rows（见 `pipelines/protein/reports/protein_master_v6.validation.json`）。

**复现命令（同事）**：

```bash
bash pipelines/protein/run.sh
```

### 4.2 L1：RNA 实体（v1）

**定位**：RNA 实体主表（L1，miRNA + mRNA/transcript）。产物体积大，因此“数据走 Release，代码/门禁进仓库”。

- Pipeline：`pipelines/rna/run.sh`
- Contract：`pipelines/rna/contracts/rna_master_v1.json`
- QA/manifest（小文件，便于审计）：
  - `pipelines/rna/reports/rna_master_v1.validation.json`
  - `pipelines/rna/reports/rna_v1.manifest.json`
- 大产物（本地/Release）：`data/output/rna_master_v1.tsv`（或 Release 里的 `.tsv.gz`）

**“怎么理解清洗”**（按脚本流水线）：

- 从 raw（miRBase/RNAcentral/Ensembl）抽取 → 映射 ID → 补 gene mapping → 修复/合并子表 → 生成统一 RNA 主表。
- Contract 门禁覆盖：主键唯一、`taxon_id==9606`、`rna_type` 枚举、序列字符集（ACGUN）、mRNA 的 gene mapping 覆盖率等。

**同事怎么获取数据**：优先使用 Release `rna-l1-v1`。

**关键规模（用于汇报一句话）**：`rna_master_v1.tsv` = 292,274 rows（见 `pipelines/rna/reports/rna_master_v1.validation.json`）。

### 4.3 L1.5：Gene spine（gene_master + xrefs）

**定位**：跨实体对齐底座：把 Protein 与 RNA 统一对齐到 `gene_key`，用于后续跨模态查询与边构建。

- 产物（当前本地已存在）：
  - `data/output/gene_master_v1.tsv`
  - `data/output/gene_xref_protein_v1.tsv`
  - `data/output/gene_xref_rna_v1.tsv`
- Pipeline（本地已存在，但尚未提交到 GitHub）：`pipelines/gene/`

**“怎么理解清洗”**：

- 使用确定性的 `gene_key` 选择规则（优先 ENSG，其次 HGNC，再次 NCBI Gene，再次 symbol），避免“同一个基因多 key”导致的漂移。
- 输出同时给出 xref 表，明确 Protein/RNA 各自如何落到同一个 gene_key（便于追溯与修复）。

**关键规模（用于汇报一句话）**：`gene_master_v1.tsv` ≈ 20,874 rows（本地跑 `pipelines/gene/run.sh` 后可生成验证报告）。

### 4.4 L1 + L2：Molecules（小分子）+ PSI（ChEMBL）

**定位**：

- L1：小分子实体（结构锚点 + 理化属性）
- L2：PSI（Protein–Small molecule activity）证据抽取与边聚合

当前实现主要在：`pipelines/molecules/`（已纳入 `feat/protein-l1-industrial` 分支；合并后进入 `main`）。

- 入口：`pipelines/molecules/run.sh`
- 发布辅助：`pipelines/molecules/package_release.sh`
- QA/门禁与 manifest：`pipelines/molecules/reports/*.json`

**卡点（与“清洗”直接相关）**：该 pipeline 依赖外部 raw 输入（ChEMBL DB/chemreps），默认不进 git。

---

## 5) L2 关系层（PPI/PSI/RPI）进度与缺口

### 3.1 PPI（STRING）

现状：我们已经在仓库内有 STRING 原始边表：

- `data/processed/protein_edges.tsv`（约 88 万条，字段：source/target_uniprot、combined_score、source、fetch_date）

**关键规模（用于汇报一句话）**：`protein_edges.tsv` = 884,555 rows（当前为 STRING_v12.0）。

缺口：还没有把它做成 L2 标准交付（edges + evidence + contract + QA + manifest）。

这意味着：

- 图谱“有边数据”，但缺少“可审计的产品化版本”
- 后续引入更多证据类型（PSI-MI、文献、实验等）时，缺少统一的 evidence 框架

### 3.2 PSI（ChEMBL）

现状：`pipelines/molecules/` 已覆盖 PSI 的 evidence 抽取与门禁，但依赖 ChEMBL raw 输入就位。

### 3.3 RPI

现状：未看到实现/产物（需要决定数据源与证据口径）。

---

## 6) 下一步缺口（对齐项目大纲：我们还要补什么）

按 `check/项目大纲.pdf` 的主线，最关键的缺口不是“多写点脚本”，而是把剩下的关系也做成同样的工业化交付：

1) **PPI（L2）标准交付**
   - 现状：只有 `data/processed/protein_edges.tsv`（原始边表）
   - 目标：`edges_ppi_v1.tsv`（边主表）+ `ppi_evidence_v1.tsv`（证据表）+ contract/QA/manifest/postmortem

2) **RPI（L2）标准交付**（RNA–Protein interaction）
   - 需要先确定 v1 数据源与证据口径（实验/数据库/文本挖掘？）

3) **Gene 对齐底座（gene_master）纳入仓库**
   - 现状：本地已有 `pipelines/gene/`（可复现产物与 QA），但尚未合并到 git 主线
   - 价值：这是跨 Protein/RNA 的统一对齐 spine，下游边构建与多跳查询都需要它

4) **L3 本体语义（GO/KEGG/ChEBI）最小产品化**
   - 建议先做最小可用：术语表 + 映射表 + contract/QA

---

## 7) 同事复现与自检（建议作为汇报收尾页）

在仓库根目录运行：

```bash
# Protein L1: QA + manifest（不会生成大文件，只生成 reports）
bash pipelines/protein/run.sh

# RNA L1: 需要 data/raw/rna 输入就位
bash pipelines/rna/run.sh

# Gene spine: 需要 RNA 输出已生成（data/output/rna_master_v1.tsv）
bash pipelines/gene/run.sh

# Molecules: 需要 data/raw/molecules 输入就位（ChEMBL）
bash pipelines/molecules/run.sh
```

---

## 8) 当前最重要的卡点与下一步（用于项目推进讨论）

1) **把 PPI(STRING) 补成 L2 标准交付**（建议优先做）
   - 新增 `pipelines/edges_ppi/*`：从 `protein_edges.tsv` 生成 `edges_ppi_v1.tsv` + `ppi_evidence_v1.tsv`，并配齐 contract/QA/manifest。

2) **把 Gene spine pipeline 纳入 git 并对齐忽略策略**
   - 代码/contract/小报告需要进 git；大 TSV 产物不进 git。

3) **明确 molecules 的“输入就位/发布策略”**
   - 同事要复现，必须知道 ChEMBL DB 怎么获取、放哪、版本如何锁定。
