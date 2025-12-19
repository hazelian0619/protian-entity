# RNA Master 数据源与版本锚点

## 1. 版本锚点（v1 锁定）

| 数据源 | 版本号 | 发布日期 | 用途 | FTP/API |
|--------|--------|----------|------|---------|
| **RNAcentral** | Release 25 | 2025-06-18 | 主 ID 体系（URS）+ 跨库映射 | https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/ |
| **miRBase** | 22.1 | 2018 | miRNA 序列 + ID（通过 RNAcentral 汇聚） | 通过 RNAcentral xref 获取 |
| **Ensembl** | 当前 Release | 动态更新 | 转录本注释 + gene 映射 | https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/ |
| **Rfam** | 15.0 | 2024-2025 | RNA 家族（v1 暂不用） | https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/ |

**v1 锁定原则**：
- 所有数据基于 RNAcentral Release 25 (2025-06-18)
- source_version 字段必须填写：`RNAcentral:25;miRBase:22.1`
- fetch_date 填写 ETL 运行日期（ISO 8601 格式）

---

## 2. 更新策略

### v1 → v2 触发条件（满足任一即升级）
1. RNAcentral 发布 Release 26
2. Rfam 数据开始实际导入（从 v1 的"预留"变为"有数据"）
3. 种子基因集扩展（从 protein-coding 扩展到全基因组）

### 版本快照保留
- 每个版本独立存档：`rna_master_v1.tsv`, `rna_master_v2.tsv`
- 版本间变更记录在 `CHANGELOG_RNA.md`
- 旧版本保留至少 1 年

---

## 3. 数据文件清单（v1 使用）

**miRNA 管线输入**：
- `id_mapping.tsv.gz`（RNAcentral xref，包含 URS ↔ miRBase 映射）
- `mature.fa.gz`（miRBase 成熟 miRNA 序列）

**mRNA 管线输入**：
- `id_mapping.tsv.gz`（RNAcentral xref，包含 URS ↔ Ensembl 映射）
- `Homo_sapiens.GRCh38.cdna.all.fa.gz`（Ensembl cDNA 序列）
- `protein_master_v6_clean.tsv`（种子基因集来源）

**文件位置**：
- 数据文件统一存放在 `data/raw/rna/` 目录
- ETL 脚本统一存放在 `scripts/rna/` 目录
