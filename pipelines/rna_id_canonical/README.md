# rna_id_canonical pipeline（RNA 统一标识桥接 v1）

该 pipeline 产出 **legacy_rna_id → canonical_rna_id** 的桥接表，不修改 `rna_master_v1.tsv` 主键。

核心原则：

- 兼容旧键：`legacy_rna_id` 全量保留
- 对齐新键：优先使用 RNAcentral URS（来自 `rna_xref_mrna_enst_urs_v2.tsv`）
- 不可映射时回退：保留 legacy ENST（保证下游可用性）

## Inputs

- `data/output/rna_master_v1.tsv`
- `data/output/rna_xref_mrna_enst_urs_v2.tsv`（助手 A 产物）
- `data/output/gene_xref_rna_v1.tsv`（可选，用于下游 join 审计）

## Outputs

- 主输出（本地保留）  
  `data/output/rna_id_canonical_map_v1.tsv`
- 合约  
  `pipelines/rna_id_canonical/contracts/rna_id_canonical_map_v1.json`
- 报告  
  `pipelines/rna_id_canonical/reports/rna_id_canonical_map_v1.build.json`  
  `pipelines/rna_id_canonical/reports/rna_id_canonical_map_v1.validation.json`  
  `pipelines/rna_id_canonical/reports/rna_id_canonical_map_v1.gene_xref_join_audit.json`

另含 sample 报告（`*.sample.*.json`），满足“先最小样本后全量”。

## Run

在仓库根目录运行：

```bash
bash pipelines/rna_id_canonical/run.sh
```

可选环境变量：

- `SAMPLE_ROWS`（默认 `5000`）
- `SAMPLE_GENE_ROWS`（默认 `5000`）
- `REQUIRE_GENE_XREF`（默认 `0`；设为 `1` 时缺少 `gene_xref_rna_v1.tsv` 会导致 gate 失败）

## Canonicalization 规则

1. **xref 命中（mRNA 少量）**  
   `legacy ENST -> canonical URS`，`canonical_system=RNAcentral`，`confidence=high`
2. **legacy 本身是 URS（miRNA）**  
   `canonical = legacy`，`canonical_system=RNAcentral`，`confidence=high`
3. **其余 ENST**  
   保持 `canonical = legacy`，`canonical_system=EnsemblTranscript`，`confidence=medium`

## QA Gates

- `legacy_rna_id` 覆盖率必须 100%
- `canonical_rna_id` 非空率必须 100%
- 当提供 `gene_xref_rna_v1.tsv` 时：legacy key join 成功率必须 100%
- 若未提供且 `REQUIRE_GENE_XREF=0`：该 gate 作为可选审计项

## 统一核查表（RNA）

| 组件 | 属性字段 | 数据类型 | 获取方式 | 存储形式 | 示例 | 用途/指标 | 指标 | 层次结构 | 语义关系 | 数据来源 | 维护机构 | 更新频率 | 当前状态 | 证据路径 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| RNA canonical map v1 | `legacy_rna_id`,`canonical_rna_id`,`canonical_system`,`confidence` 等 10 列 | TSV | `rna_master_v1` 与 `rna_xref_mrna_enst_urs_v2` 规则映射 | `data/output/rna_id_canonical_map_v1.tsv` | `ENST00000000233_9606 -> ENST00000000233_9606` / `ENST00000514884_9606 -> URS00002EBCD8_9606` | 统一标识桥接，兼容 legacy、对齐 canonical | 覆盖率、非空率、join 成功率、分层变化率 | RNA 主表上层桥接 | legacy key 到 canonical key | RNA master + RNA xref + gene_xref | 本仓库数据装配流程 | 随补洗重跑 | 已落地 | `pipelines/rna_id_canonical/reports/*.json` |
