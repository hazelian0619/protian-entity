# Molecules（small molecules）L1+L2（v1）

这条管线把 **小分子实体（L1）** 与 **蛋白-小分子活性证据（L2 / PSI）** 做成可复现数据产物。

项目主线口径：先把 **可复现的数据集与证据层** 做充分，把“专业生信/生物学判断”尽量后置到 evidence/query 层。

## 同事怎么用（最短路径）

优先直接下载 Release 数据（不需要跑代码）：

- Molecules L1（M1+M2）：Release `molecules-l1-v1`
  - https://github.com/hazelian0619/protian-entity/releases/tag/molecules-l1-v1
- Molecules PSI L2（M3）：Release `molecules-psi-l2-v1`
  - https://github.com/hazelian0619/protian-entity/releases/tag/molecules-psi-l2-v1

Release 里会包含：

- 大体积 SQLite 数据库（分卷压缩）
- QA 报告（JSON）
- manifest（JSON，sha256/大小/时间/commit）
- postmortem（复盘）

## 本地重跑（开发者）

### 1) 准备输入数据（不入 git）

放到以下路径：

- `data/raw/molecules/chembl_36_chemreps.txt`
- `data/raw/molecules/chembl_36/chembl_36.db`

说明：

- `chemreps` 只用于 M1（结构锚点）
- `chembl_36.db` 是 M3 必需的 evidence 来源（activities/assays/targets/docs）

### 2) 运行

从仓库根目录：

```bash
bash pipelines/molecules/run.sh
```

输出：

- `data/output/molecules/molecules_m1.sqlite`（M1+M2）
- `data/output/molecules/chembl_m3.sqlite`（M3）
- `pipelines/molecules/reports/**`（QC/gates 小文件）

## QA（门禁）

- `pipelines/molecules/scripts/m3_quality_gates.py` 会输出 `molecules_m3_v1.gates.json`，用于“可不可以入库/发版”。

## 复盘（维护者）

- `pipelines/molecules/POSTMORTEM_molecules_m1_v1.md`
- `pipelines/molecules/POSTMORTEM_molecules_m2_v1.md`
- `pipelines/molecules/POSTMORTEM_molecules_m3_v1.md`
解压示例（zstd）：

```bash
zstd -d molecules_m1.sqlite.zst -o molecules_m1.sqlite
zstd -d chembl_m3.sqlite.zst -o chembl_m3.sqlite
```

