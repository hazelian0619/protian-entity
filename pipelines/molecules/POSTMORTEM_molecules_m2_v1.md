# M2（RDKit 理化属性）全流程复盘 + 决策说明

> 目标：在 M1 结构锚点（InChIKey + SMILES）基础上，计算小分子的 2D 理化描述符（MW/LogP/TPSA/HBD/HBA/RotB 等），作为小分子节点的多模态属性之一，支持筛选、排序、建模与后续 PSI 边的特征增强。

## 1. 背景与目标对齐（为什么做 M2）

在 `check/组件细节.docx` 的“小分子实体”中，理化属性是核心属性模块之一：
- 分子量（MW）
- LogP
- TPSA
- HBD/HBA
- 旋转键数等

这些属性的价值：
- 为“成药性/筛选/粗过滤”提供可解释的数值依据
- 为后续模型（图嵌入/预测）提供基础特征
- 在 PSI 边（蛋白-小分子）上，可用于排序、分层与质量控制（例如极端值/异常值识别）

## 2. 输入与输出（做什么）

### 输入
- 来自 M1 的 SQLite：`data/output/molecules/molecules_m1.sqlite`
- 使用的源表：`molecule_entity_strict_chembl36`
  - strict 口径默认：`has_dot=0 AND smiles_len<=200`
  - 目的：先排除明显多组分/超长长尾，快速获得“更像小分子”的主体集合

### 输出
- 仍写回同一 SQLite（便于联表）：`data/output/molecules/molecules_m1.sqlite`
- 输出表：`molecule_props_rdkit_v1`
- QC 输出目录：`pipelines/molecules/reports/m2/`
  - 最新：`pipelines/molecules/reports/m2/qc_m2_rdkit_latest.md`
  - 每次运行快照：`pipelines/molecules/reports/m2/qc_m2_rdkit_<run_id>.md`

## 3. M2 的范围界定（不做什么）

M2 **不做**：
- 化学标准化（去盐/标准化 tautomer 等）
- 3D 构象生成/量化计算
- 生物活性数据抽取（这是 M3）

M2 **只做**：
- 用 RDKit 从 SMILES 构建分子
- 计算一组稳定的 2D 描述符

## 4. 计算字段（输出 schema）

输出表 `molecule_props_rdkit_v1` 字段：
- `inchikey`：主键
- `rdkit_canonical_smiles`：RDKit canonical SMILES（用于标准化对照）
- `mol_wt`：分子量
- `mol_logp`：Crippen LogP
- `tpsa`：拓扑极性表面积
- `hbd`：氢键供体数
- `hba`：氢键受体数
- `rot_bonds`：可旋转键数
- `rings`：环数
- `heavy_atoms`：重原子数
- `frac_csp3`：sp3 碳比例
- `rdkit_sanitize_ok`：是否通过标准 sanitize（1=正常，0=使用 fallback）
- `rdkit_error`：fallback 的原因/警告（用于排查）
- `computed_at_utc`：计算时间

## 5. 执行流程（一步步怎么做）

### Step 0：环境准备（为什么需要 venv）

- 系统默认 `python` 是 2.7，不适合 RDKit。
- 使用 Python 3.x 建立虚拟环境（例如 `python3 -m venv .venv`）

安装依赖（需要联网）：
- `rdkit-pypi`
- `tqdm`
- **关键兼容点**：`rdkit-pypi (2022.9.x)` 在 macOS 上需要 `numpy<2`

### Step 1：先跑样本（验证正确性与稳定性）

- 先跑 10 万条样本（`--limit 100000`）
- 目标：验证解析成功率、速度、输出 schema 与 QC 机制

样本运行结果（首次验证）：
- ok=100000, fail=0
- 吞吐约 2.3k rows/s

### Step 2：全量补齐（skip-existing 增量模式）

- 全量运行时启用 skip-existing（默认）
  - 逻辑：只选择 `source_table` 里还没有出现在输出表的 inchikey
  - 好处：可断点续跑、可多次执行补齐，不怕中断

### Step 3：失败处理策略（保质保量的关键）

实际全量过程中出现极少数 RDKit sanitize 问题（例如 kekulize 失败）。

处理策略：
1) 先尝试 `MolFromSmiles(smiles, sanitize=True)`
2) 若失败：
   - 用 `sanitize=False` 读取
   - 再尝试 `SanitizeMol`
   - 若仍失败：尝试跳过 kekulize 的 sanitize（保留分子用于计算多数 2D 描述符）
3) 将是否 fallback 记录到 `rdkit_sanitize_ok` 与 `rdkit_error`

> 目标是“覆盖率优先 + 可追溯”，避免因为极少数结构导致全量结果缺口。

### Step 4：生成 QC

QC 输出为：
- `qc_m2_rdkit_latest.md/json`
- `qc_m2_rdkit_<run_id>.md/json`

包含：existing_before、selected、ok、fail、速度等。

## 6. 最终结果（覆盖率与质量）

最终覆盖（严格集）：
- strict：2,699,403
- props：2,699,403
- missing_props：0
- rejects：0
- `rdkit_sanitize_ok=0`：1（极少数 fallback 解析）

自查命令：

```bash
sqlite3 data/output/molecules/molecules_m1.sqlite \
  "SELECT 'strict', COUNT(*) FROM molecule_entity_strict_chembl36;
   SELECT 'props', COUNT(*) FROM molecule_props_rdkit_v1;
   SELECT 'missing_props', (
      SELECT COUNT(*)
      FROM molecule_entity_strict_chembl36 s
      LEFT JOIN molecule_props_rdkit_v1 p ON s.inchikey=p.inchikey
      WHERE p.inchikey IS NULL
   );
   SELECT 'sanitize_ok_0', COUNT(*) FROM molecule_props_rdkit_v1 WHERE rdkit_sanitize_ok=0;"
```

属性范围 sanity check（示例）：
- `mol_wt`：min≈4.003, mean≈418.85, max≈2631.14
- `mol_logp`：min≈-16.89, mean≈3.52, max≈38.90
- `tpsa`：min=0.0, mean≈86.30, max≈863.01

> 注：极端值存在是正常的（长尾仍可能包含一些“偏大/偏极性”的结构），后续可在 M3 子集或下游任务中做更严格过滤。

## 7. 卡点与反思（为什么会出现“卡很久/不动”的体验）

### 7.1 计算本质是 CPU + I/O 密集
- 2.7M 分子全量计算包含：SMILES 解析、descriptor 计算、SQLite 写入与索引维护。
- 即使吞吐 2k+/s，也需要 20–45 分钟量级。

### 7.2 进度输出与缓冲
- 某些终端环境会对 stderr/stdout 缓冲，导致“看起来没输出”。
- 脚本已增加 `flush=True` 保证进度实时输出。

### 7.3 二进制依赖兼容问题（NumPy 2）
- 初次安装时 `numpy==2.x` 会导致 rdkit-pypi 导入失败（binary ABI 不兼容）。
- 已固定为 `numpy<2` 解决。

## 8. 关键决策（贴合目的：保质保量怎么做）

- **为什么先 strict**：
  - M2 是“属性层”，不应被多组分/超长长尾拖累，先把主体集合做稳。
  - strict 仍覆盖 269 万级别，足以支撑大部分分析。

- **为什么要 fallback 而不是直接丢弃**：
  - 丢弃会造成缺口与后续 join 的不一致；fallback 让覆盖率更完整。
  - 通过 `rdkit_sanitize_ok/rdkit_error` 保留可追溯性。

- **是否需要 GPU**：不需要。

## 9. 下一步建议（如何尽快推到 M3）

- M3 的核心是 PSI 边（蛋白-小分子活性证据）。建议使用 ChEMBL release 的 activities/assays/targets。
- 下载 ChEMBL SQLite dump 后，运行 `pipelines/molecules/scripts/m3_extract_chembl_psi.py` 先小跑 `--limit 100000` 验证 schema/QC，再按需全量。

相关计划文档：`pipelines/molecules/plan.md`
