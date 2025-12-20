# M3 决策清单（v1，已根据 full QC 填充）

> 目的：把“需要专业生信背景才能定的点”先降维为可观测指标。
> 当 M3 smoke-run 的 QC 出来后，只需要在本页填数字/勾选，就能确定 v1 默认口径。

## 1) 输入版本

- ChEMBL release：36
- 输入 DB 路径：`data/raw/molecules/chembl_36/chembl_36.db`
- run_id（full QC）：`m3_20251219T170848Z`

## 2) 目标：先把数据准备全，不强行拍板生物学口径

策略：
- 证据表（activity/assay/doc）尽量保留上下文字段（便于专家后续筛选/重聚合）。
- 聚合边表可先按一个“默认可解释规则”生成，后续可在证据表基础上重建。

## 3) v1 默认过滤（可调参数）

- `standard_type` allowlist：IC50 / Ki / Kd / EC50
- `min_confidence`：6
- 是否要求 `compound_inchikey`：是（`--require-inchikey`，本次 missing=0）
- 是否要求 `target_uniprot_accession`：是（`--require-uniprot`，本次 missing=0）

## 4) 关键卡点（需要 QC 驱动决策）

### 4.1 target 映射是否“爆炸”

指标（由 smoke-run 填）：
- evidence 行数：3,113,487
- distinct `target_tid`：9,210
- distinct `target_chembl_id`：9,210
- distinct `target_uniprot_accession`（代表 accession）：8,795
- target→component 的平均映射数（`target_uniprot_map_v1`，按 target_tid）：约 1.29（max=78）

决策：
- v1 边端点：
  - 默认用于 KG 的端点：`target_uniprot_accession`（missing=0，且更易与蛋白实体对齐）
  - 同时保留：`target_chembl_id`（用于审计/回溯/必要时改用 ChEMBL target 口径）

### 4.2 活性可比性（units/type/pchembl）

指标：
- `standard_units` Top10（full QC）：
  - nM: 3,104,738
  - ug.mL-1: 6,389
  - (NULL): 539
  - 10^2 uM: 162
  - %: 154
  - /s: 106
  - 10'-3/s: 92
  - 10'-1/s: 77
  - 10'-2/s: 75
  - 10^3nM: 58
- `pchembl_value` 缺失率：约 26.65%
- `pchembl_value_eff`（补齐后）缺失率：约 0.32%

决策：
- v1 排序主指标：`pchembl_value_eff`（推荐；优先用 pchembl，缺失时用单位换算补齐）

### 4.3 与 L1/M2 的 join 覆盖率（工程主线指标）

指标：
- evidence 中 `compound_inchikey` 能 join 到 `data/output/molecules/molecules_m1.sqlite` 的 `molecule_props_rdkit_v1` 覆盖率：0.958236

决策：
- 若覆盖率低：优先排查 `compound_structures` 缺失、或是否需要更宽松/更严格的 `require-inchikey` 口径。

## 5) v1 聚合规则（边表）

- 聚合 key（默认）：(compound_inchikey, target_key, standard_type)
- 代表证据选择：按 `pchembl_value_eff` desc → `evidence_score_v1` desc → `confidence_score` desc → `activity_id` asc

其中 `target_key` 在库中体现为：

- `psi_edges_v1.target_uniprot_accession`（用于实体对齐）
- `psi_edges_v1.target_chembl_id`（用于审计/回溯）
