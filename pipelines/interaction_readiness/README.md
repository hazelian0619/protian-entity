# interaction_readiness（互作层 L2 本地就绪度核查）

目标：不改业务语义，仅对 **PPI / PSI / RPI** 做同口径现状量化、差距分级、收口建议。

## 核查对象

- PPI：`edges_ppi_v1` / `ppi_evidence_v1`
- PSI：`drug_target_*` + `molecules M3 psi_*`
- RPI：`rna_protein_edges_v1` / `rna_protein_evidence_v1`

## 输出

- `docs/interaction/INTERACTION_L2_LOCAL_READINESS_2026-04-12.md`
- `pipelines/interaction_readiness/reports/interaction_l2_readiness_2026-04-12.json`
- `pipelines/interaction_readiness/reports/interaction_l2_readiness_2026-04-12.validation.json`
- `pipelines/interaction_readiness/reports/interaction_l2_readiness_2026-04-12.gates.json`
- `pipelines/interaction_readiness/reports/interaction_l2_readiness_2026-04-12.manifest.json`

## 维度（统一口径）

- 组件、属性字段、数据类型、获取方式、存储形式、示例
- 用途/指标、指标、层次结构、语义关系
- 数据来源、维护机构、更新频率、当前状态、证据路径

## 运行

```bash
bash pipelines/interaction_readiness/run.sh
```

流程：
1. 互作产物工业化落盘（sample）  
   - 从历史目录（默认 `../1218`）抽取并标准化到当前仓 canonical 路径  
2. 互作产物工业化落盘（full）  
3. PPI/PSI 标准化表 contract 校验  
4. readiness 最小样本核查（sample）  
5. readiness 全量核查（full）  
6. readiness report contract 校验  
7. gates 生成  
8. manifest 生成

## 可选参数

默认会尝试扫描 `../1218` 作为仓外候选证据（仅用于 gap 建议，不计入本仓就绪）。
如需指定其他候选目录：

```bash
INTERACTION_OFFREPO_ROOT=/path/to/legacy_workspace \
  bash pipelines/interaction_readiness/run.sh
```

如需指定“工业化落盘”的源目录（默认 `../1218`）：

```bash
INTERACTION_SOURCE_ROOT=/path/to/legacy_workspace \
  bash pipelines/interaction_readiness/run.sh
```

## 验收口径

- 三类互作都有同口径状态表：`已完成 / 部分 / 缺失`
- Top10 缺口按 `阻塞级/可延期` + `S/M/L`
- 明确“本地通过后再上传”的前置门槛列表
