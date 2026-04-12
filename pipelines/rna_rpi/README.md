# RNA-Protein Interaction 最小可用管线（rna_rpi_v1）

目标：产出 L2 最小可用 `edges + evidence + QA + manifest`。

## 输入

- `data/output/rna_master_v1.tsv`
- `data/processed/protein_master_v6_clean.tsv`
- RPI 离线快照（RNAInter / NPInter / starBase 任一）
  - 推荐放置：`data/raw/rpi/`
  - 推荐命名：
    - `data/raw/rpi/rnainter_human.tsv.gz`
    - `data/raw/rpi/npinter_human.tsv.gz`
    - `data/raw/rpi/starbase_human.tsv.gz`

> 输入支持 `TSV/CSV/TXT`，支持 `.gz`。

## 输出

- `data/output/edges/rna_protein_edges_v1.tsv`
- `data/output/evidence/rna_protein_evidence_v1.tsv`
- `pipelines/rna_rpi/contracts/*.json`
- `pipelines/rna_rpi/reports/*.json`

## 字段

### 1) edges：`rna_protein_edges_v1.tsv`

`edge_id,src_type,src_id,dst_type,dst_id,predicate,directed,best_score,source,source_version,fetch_date`

### 2) evidence：`rna_protein_evidence_v1.tsv`

`evidence_id,edge_id,evidence_type,method,score,reference,source,source_version,fetch_date`

## 运行

```bash
bash pipelines/rna_rpi/run.sh
```

### 可选：自动下载 starBase/ENCORI 快照（本次 A 已采用）

```bash
mkdir -p data/raw/rpi
curl -L -A 'Mozilla/5.0' \
  'https://rnasysu.com/encori/api/RBPTarget/?assembly=hg38&geneType=mRNA&RBP=all&clipExpNum=1&pancancerNum=0&target=all&cellType=all' \
  -o data/raw/rpi/starbase_human.tsv
sha256sum data/raw/rpi/starbase_human.tsv
```

本次产物使用：`data/raw/rpi/starbase_human.tsv`（ENCORI API）。

流程：
1. preflight 检查（缺输入立即中断）
2. 最小样本 (`--limit 200`)
3. 全量构建
4. edges contract 校验
5. evidence contract 校验
6. 生成 manifest

## 报告

- `rna_rpi_v1.blocked_or_ready.json`
- `rna_rpi_v1.sample.metrics.json`
- `rna_rpi_v1.metrics.json`
- `rna_rpi_v1.gates.json`
- `rna_protein_edges_v1.validation.json`
- `rna_protein_evidence_v1.validation.json`
- `rna_rpi_v1.manifest.json`

## 验收门槛

- `src_id`（RNA）join coverage >= 0.99
- `dst_id`（Protein）join coverage >= 0.99
- `edge_id` 唯一，且 self-loop = 0
- evidence 中 `method` 或 `reference` 至少一项非空（脚本输出覆盖率）

## 手动下载中断条件（缺 RPI 原始输入）

若 preflight 报 `blocked_missing_inputs`，按下列方式补齐并重跑：

- 下载文件清单与版本（任一可用）：
  - RNAInter human RNA-protein interactions snapshot
  - NPInter human ncRNA-protein interactions snapshot
  - starBase CLIP-supported RNA-protein interactions snapshot
- 放置路径：`data/raw/rpi/`
- 完整性校验命令：

```bash
sha256sum data/raw/rpi/rnainter_human.tsv.gz
sha256sum data/raw/rpi/npinter_human.tsv.gz
sha256sum data/raw/rpi/starbase_human.tsv.gz
```
