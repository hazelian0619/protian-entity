# RNA Rfam + 二级结构补洗管线（rna_rfam_structure_v1）

将 `rfam_id` 与 `secondary_structure` 从预留空列推进到可用结构层，输出：

- `data/output/rna_structure_rfam_v1.tsv`

## 输入

- `data/output/rna_master_v1.tsv`
- `data/raw/rna/rnacentral/id_mapping.tsv(.gz)`
- `data/raw/rna/rnacentral/rfam_annotations.tsv.gz`
- `data/raw/rna/rfam/family.txt.gz`
- `data/raw/rna/rfam/Rfam.seed.gz`
- （当 `rna_master_v1.tsv` 含 ENST 主键时）
  - `data/output/rna_xref_mrna_enst_urs_v2.tsv`（或 `v1`）

## 输出字段

- `rna_id`
- `rfam_id`
- `rfam_desc`
- `secondary_structure`
- `structure_source`
- `score`
- `evalue`
- `source_version`
- `fetch_date`

## 目录结构

- `pipelines/rna_rfam_structure/run.sh`
- `pipelines/rna_rfam_structure/scripts/build_rna_rfam_structure.py`
- `pipelines/rna_rfam_structure/contracts/rna_structure_rfam_v1.json`
- `pipelines/rna_rfam_structure/reports/*.json`

## 运行

```bash
bash pipelines/rna_rfam_structure/run.sh
```

该脚本会：
1. 先做输入检查（缺数据时立刻中断并写阻塞报告）
2. 再跑最小样本（`--limit 200`）
3. 再跑全量
4. 最后按 contract 输出 validation 报告

映射策略（按优先级）：
- RNAcentral `rfam_annotations.tsv.gz`（URS -> Rfam）
- 对 `miRNA` 在无法命中 RNAcentral 映射时，回退到 `Rfam family name` 启发式映射（如 `hsa-miR-660-3p -> mir-660`）
- 二级结构来自 `Rfam.seed.gz` 中 `SS_cons`

## 报告文件

- `rna_structure_rfam_v1.blocked_or_ready.json`：输入检查结果（阻塞/就绪）
- `rna_structure_rfam_v1.sample.metrics.json`：最小样本统计
- `rna_structure_rfam_v1.metrics.json`：全量统计（含按 `rna_type` 分层覆盖率）
- `rna_structure_rfam_v1.validation.json`：contract 校验

## 手动下载中断条件与清单

当输入缺失或线上接口不稳定时，请改用 Rfam 离线包，并放置到：

- `data/raw/rna/rfam/family.txt.gz`
- `data/raw/rna/rfam/Rfam.seed.gz`

同时准备：

- `data/raw/rna/rnacentral/id_mapping.tsv.gz`
- `data/raw/rna/rnacentral/rfam_annotations.tsv.gz`
- `data/output/rna_master_v1.tsv`

校验方式（示例）：

```bash
sha256sum data/raw/rna/rfam/family.txt.gz
sha256sum data/raw/rna/rfam/Rfam.seed.gz
sha256sum data/raw/rna/rnacentral/id_mapping.tsv.gz
sha256sum data/raw/rna/rnacentral/rfam_annotations.tsv.gz
sha256sum data/output/rna_master_v1.tsv
```

## 统一核查模板（本管线）

| 组件 | 属性字段 | 数据类型 | 获取方式 | 存储形式 | 示例 | 用途/指标 | 指标 | 层次结构 | 语义关系 | 数据来源 | 维护机构 | 更新频率 | 当前状态 | 证据路径 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 结构数据 | `rfam_id`,`secondary_structure` | `String`,`Dot-bracket` | RNAcentral `rfam_annotations` + Rfam `family/seed`（miRNA 回退启发式） | `data/output/rna_structure_rfam_v1.tsv` | `RF00027`, `(((...)))` | 结构层回链与家族注释 | `rna_id`回链率、`rfam_id`格式、dot-bracket合法率、分层覆盖率 | L1.5 结构补充层 | RNA `has_family` Rfam；RNA `has_secondary_structure` 结构字符串 | RNAcentral, Rfam | RNAcentral Consortium, Rfam/EMBL-EBI | 随 RNA/Rfam 版本重跑 | `inputs_ready`（可运行） | `pipelines/rna_rfam_structure/reports/rna_structure_rfam_v1.metrics.json` |
