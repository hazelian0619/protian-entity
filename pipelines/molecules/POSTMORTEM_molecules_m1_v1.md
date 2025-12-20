# M1（小分子实体底座）全流程复盘 + 卡点反思

> 目标：把 ChEMBL chemreps 变成“可复现、可追溯、可回归”的小分子实体底座（L1 静态层），为后续 M2（理化属性/RDKit）与 M3（PSI 活性证据边）提供统一锚点。

## 1. 背景与目标对齐（为什么做 M1）

从 `check/组件细节.docx` 的定义看，小分子实体至少需要：
- **结构标识**：SMILES / InChI / InChIKey（用于唯一性、去重、跨库 join）
- **数据库标识**：PubChem CID / ChEMBL ID / DrugBank ID / ZINC ID（用于跨库对齐）
- 后续增强：理化属性、3D、活性数据、分类（ChEBI/ATC）、ADME/Tox

而 `check/项目大纲.pdf` 的路线是：
- L1：先把实体（节点）做稳（“它是什么”）
- L2：再做关系（边）并附证据（“为什么相关”）
- L3：做本体对齐（术语层级/语义检索）
- L4：做多跳推理

**因此 M1 的定位**：
- 不追求“全字段/全来源”；
- 先把“结构唯一锚点 + 版本化 + QC 回归”做成可交付产品。

## 2. 输入数据（抓什么）

- 输入文件：`data/raw/molecules/chembl_36_chemreps.txt`
- 行数：2,854,816（含表头；数据行 2,854,815）
- 字段（表头）：
  - `chembl_id`
  - `canonical_smiles`
  - `standard_inchi`
  - `standard_inchi_key`

> 选择该数据源的原因：它稳定、可离线、结构锚点齐全，适合做“实体底座”。

## 3. 输出产物（产出什么）

- SQLite 数据库：`data/output/molecules/molecules_m1.sqlite`
- QC 报告：
  - 人可读：`pipelines/molecules/reports/m1/qc_report.md`
  - 机器可读：`pipelines/molecules/reports/m1/qc_report.json`
- 可复跑脚本：`pipelines/molecules/scripts/m1_build_sqlite.py`

## 4. M1 流程（怎么抓：一步步做了什么）

### Step 0：定义“底座=结构锚点 + 质量标记 + 可追溯”

M1 **不做**：去盐、标准化、3D、活性数据抽取。
M1 **只做**：
- 结构字段落库（SMILES/InChI/InChIKey/ChEMBL ID）
- 质量标记（has_dot, smiles_len）
- 去重汇总为实体表（以 InChIKey 为锚）
- 输出 QC 与元数据（可回归）

### Step 1：流式解析 TSV → 写入 `molecule_raw_chembl36`

- 使用 Python 标准库 `csv.DictReader`（制表符分隔）逐行读取
- 对 `\N` 与空字符串规范化为 NULL
- 生成行号 `ingest_rownum` 便于定位与复盘

Reject 规则（进入 `molecule_rejects_chembl36`）：
- `smiles` 缺失
- `inchikey` 缺失
- `inchikey` 长度异常（不等于 27）

> 本次数据质量极好：rejects=0。

### Step 2：构建实体表（按 InChIKey 聚合）

- `molecule_entity_chembl36`：
  - 主键：`inchikey`
  - 衍生：`inchikey_connectivity = inchikey[:14]`
  - 代表值：`rep_chembl_id`、`smiles`、`inchi`
  - 统计：`n_rows`（该 InChIKey 在 raw 中出现次数）

- `molecule_entity_strict_chembl36`（严格版，方便后续“抓大放小”）：
  - 默认规则：`smiles_len <= 200` 且 `has_dot = 0`
  - 目的：先排除明显的多组分/大分子长尾，快速进入可计算属性/可做建模的主体集合

### Step 3：生成 QC 报告

QC 重点覆盖：
- 规模与损耗：raw/entity/strict/rejects
- 完整性：缺失率
- 多组分比例（SMILES 含 `.`）
- 过长结构比例（SMILES 长度长尾）
- 最长 SMILES 示例（便于确认长尾是什么类型）

## 5. 运行方式（如何复现）

```bash
python3 pipelines/molecules/scripts/m1_build_sqlite.py \
  --input data/raw/molecules/chembl_36_chemreps.txt \
  --db data/output/molecules/molecules_m1.sqlite \
  --outdir pipelines/molecules/reports/m1 \
  --batch-size 20000 \
  --progress-every 500000
```

## 6. 结果汇总（抓取效果如何）

来自 `pipelines/molecules/reports/m1/qc_report.md`：
- 数据行数：2,854,815
- `raw`：2,854,815
- `rejects`：0
- `entity`：2,854,815
- `entity_strict`：2,699,403
- 多组分（has_dot=1）：118,896（4.1648%）
- 过长（smiles_len>200）：38,009（1.3314%）
- 同时满足两者：1,493

严格集排除构成：
- 排除总计：155,412
- 仅多组分：117,403
- 仅过长：36,516
- 两者都有：1,493

SMILES 长度分布（采样/统计）：
- p50≈51
- p95≈106
- p99≈229
- max=2093（属于长尾；QC 报告里展示了前 20 个前缀样例）

性能：全量运行约 28 秒完成（本机）。

## 7. 自查方法（如何证明结果可信）

### 7.1 表规模核对

```bash
sqlite3 data/output/molecules/molecules_m1.sqlite \
  "SELECT 'raw',count(*) FROM molecule_raw_chembl36
   UNION ALL SELECT 'rejects',count(*) FROM molecule_rejects_chembl36
   UNION ALL SELECT 'entity',count(*) FROM molecule_entity_chembl36
   UNION ALL SELECT 'entity_strict',count(*) FROM molecule_entity_strict_chembl36;"
```

### 7.2 多组分与长尾核对

```bash
sqlite3 data/output/molecules/molecules_m1.sqlite \
  "SELECT ROUND(AVG(has_dot),6) AS dot_rate FROM molecule_raw_chembl36;"
```

## 8. 卡点反思（为什么 M1 之前会卡住）

### 8.1 工程复现性卡点（根因）

- 之前脚本依赖 `pandas/pubchempy/tqdm`，但环境未安装，导致无法跑通。
- 文件路径与命名不统一：需要明确输入文件在 `data/raw/molecules/`。
- 数据源与表名混用（ChEMBL 数据写到 `molecules_pubchem`）会导致溯源混乱。

**M1 的解决办法**：
- 纯标准库、离线、流式处理（无第三方依赖）
- 固定输入/输出路径
- 产物 + QC + 元数据一次性闭环

### 8.2 口径卡点（需要尽快定，但不应阻塞 M1）

- **多组分（SMILES 含 `.`）怎么处理**：去盐？取最大组分？全部保留？
- **小分子边界**：是否排除肽/寡核苷酸等大分子长尾？阈值怎么定？

**M1 的策略**：先标记（has_dot/smiles_len），并提供 strict 子集让你快速推进后续。

## 9. 下一步建议（怎么抓大放小，尽快推到 M2/M3）

### 9.1 M2（理化属性）建议路线

RDKit 计算的 2D 理化描述符（MW/LogP/TPSA/HBD/HBA/RotB 等）本质是 **CPU 任务**。
- 不需要 GPU（RDKit 不用 GPU；除非后续要做 3D 构象生成/深度模型训练）。

抓大放小策略：
1) 先对 `molecule_entity_strict_chembl36` 抽样（如 10 万）跑通 RDKit 属性计算与 QC。
2) 再决定全量：
   - 如果目标是“通用小分子属性库”，就对 strict 全量计算。
   - 如果目标是“贴合 PSI/药物发现”，优先对 **有活性证据的化合物子集** 计算（需要 M3 先产出 PSI 边或从 ChEMBL activities 得到 compound 子集）。

### 9.2 M3（PSI 证据边）建议路线

- 不建议走 PubChem API 临时抓 10k（不可复现 + 偏差大）。
- 建议直接使用 ChEMBL release（activities/assays/targets），并把 assay 上下文、单位、置信度纳入边属性。

## 10. 你需要做的最小决策（不影响 M1，但决定 M2/M3 的方向）

1) 小分子边界是否采用 strict（默认：`has_dot=0` 且 `smiles_len<=200`）作为主入口？
2) 多组分策略：
   - A：先排除（用于快速建模/属性计算）
   - B：保留，但在 M2/M3 进一步标准化（去盐/主组分提取）
