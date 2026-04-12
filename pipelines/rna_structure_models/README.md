# rna_structure_models：结构模型补齐（v1）

产出两张结构模型相关表：

- `data/output/rna_covariance_models_index_v1.tsv`
- `data/output/rna_predicted_structures_v1.tsv`

默认流程：**preflight → sample → full → validation → 文件引用 QA → manifest**。

---

## 输入

必需：

- `data/raw/rna/rfam/family.txt.gz`
- CM 来源二选一：
  - `data/raw/rna/rfam/Rfam.cm.gz`（bundle）
  - `data/raw/rna/rfam/cm/`（按 RFxxxxx 拆分 CM 文件）
- 预测结构来源二选一：
  - `data/raw/rna/predicted_structures/manifest.tsv`
  - `data/raw/rna/predicted_structures/`（扫描 RNAfold/RhoFold 产物）

支持扫描的结构文件后缀：

- RNAfold：`.dbn`, `.ct`, `.bpseq`
- RhoFold：`.pdb`, `.cif`, `.mmcif`

---

## 输出字段

### 1) 协方差模型索引

- `rfam_id`
- `cm_file`
- `ga_threshold`
- `source_version`
- `fetch_date`

### 2) 预测结构登记

- `rna_id`
- `model_tool`
- `structure_file`
- `confidence`
- `energy`
- `source_version`
- `fetch_date`

---

## 合约与报告

- 合约：
  - `pipelines/rna_structure_models/contracts/rna_covariance_models_index_v1.json`
  - `pipelines/rna_structure_models/contracts/rna_predicted_structures_v1.json`
- 报告：
  - `pipelines/rna_structure_models/reports/rna_structure_models_v1.blocked_or_ready.json`
  - `pipelines/rna_structure_models/reports/rna_structure_models_v1.sample.metrics.json`
  - `pipelines/rna_structure_models/reports/rna_structure_models_v1.metrics.json`
  - `pipelines/rna_structure_models/reports/rna_covariance_models_index_v1.validation.json`
  - `pipelines/rna_structure_models/reports/rna_predicted_structures_v1.validation.json`
  - `pipelines/rna_structure_models/reports/rna_structure_models_v1.file_existence.json`
  - `pipelines/rna_structure_models/reports/rna_structure_models_v1.manifest.json`

---

## 运行

```bash
bash pipelines/rna_structure_models/run.sh
```

可选环境变量：

- `DATA_VERSION`（manifest 用）
- `SOURCE_VERSION_RFAM`（默认 `Rfam:current`）
- `SOURCE_VERSION_PRED`（默认 `RNAfold/RhoFold:local`）
- `SAMPLE_LIMIT_COV`（默认 `200`）
- `SAMPLE_LIMIT_PRED`（默认 `200`）
- `INPUT_RFAM_FAMILY`（默认 `data/raw/rna/rfam/family.txt.gz`）
- `INPUT_RFAM_CM_BUNDLE`（默认 `data/raw/rna/rfam/Rfam.cm.gz`）
- `INPUT_RFAM_CM_DIR`（默认 `data/raw/rna/rfam/cm`）
- `INPUT_PREDICTED_ROOT`（默认 `data/raw/rna/predicted_structures`）
- `INPUT_PREDICTED_MANIFEST`（默认 `data/raw/rna/predicted_structures/manifest.tsv`）

---

## 缺输入中断策略

若缺输入，preflight 会直接中断并在 `blocked_or_ready` 报告输出下载清单（含路径与 sha256 命令），不会 silently skip。

示例校验命令：

```bash
sha256sum data/raw/rna/rfam/family.txt.gz
sha256sum data/raw/rna/rfam/Rfam.cm.gz
find data/raw/rna/rfam/cm -type f -name '*.cm*' -print0 | xargs -0 sha256sum
find data/raw/rna/predicted_structures -type f -print0 | xargs -0 sha256sum
```
