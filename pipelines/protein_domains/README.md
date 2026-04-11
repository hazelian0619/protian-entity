# Protein 结构域管线（protein_domains_interpro_v1）

将 `data/processed/protein_master_v6_clean.tsv` 中的 `uniprot_id, domains` 升级为 InterPro/Pfam 结构化映射表。

## 输出

- 主表：`data/output/protein/protein_domains_interpro_v1.tsv`
- 字段：`uniprot_id, interpro_id, pfam_id, entry_name, start, end, source_version, fetch_date`

## 目录结构

- `pipelines/protein_domains/run.sh`
- `pipelines/protein_domains/scripts/build_protein_domains_interpro.py`
- `pipelines/protein_domains/contracts/protein_domains_interpro_v1.json`
- `pipelines/protein_domains/reports/*.json`

## 依赖

```bash
python3 -m pip install requests
```

## 运行

在仓库根目录执行：

```bash
bash pipelines/protein_domains/run.sh
```

脚本会：
1. 先跑最小样本（`--limit 30`）
2. 再跑全量（InterPro API）
3. 用 contract 做结构校验

## 报告说明

- `protein_domains_interpro_v1.sample.metrics.json`：样本验收报告
- `protein_domains_interpro_v1.metrics.json`：全量覆盖率/缺失/格式校验/API吞吐
- `protein_domains_interpro_v1.audit.json`：完整审计（失败ID、缺失ID、非可回连ID）
- `protein_domains_interpro_v1.validation.json`：contract 校验结果

## API 中断与手动下载策略

当 **全量** 运行出现以下任一情况时，脚本会中断并在审计中给出 bulk 下载方案：

- API 失败率 `> 5%`
- 吞吐 `uniprot/sec < --min-throughput`（默认 5）

建议本地目录：`data/raw/protein_domains/`。
审计里会给出需下载文件、预期大小与 `sha256sum` 校验命令。
