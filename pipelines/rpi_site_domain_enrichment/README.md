# RPI 位点/结构域/功能增强管线（rpi_site_domain_enrichment_v2）

目标：在 `rna_protein_edges_v1 + rna_protein_evidence_v1` 基础上补齐 RPI V2 上下文三张增强表：

- `data/output/evidence/rpi_site_context_v2.tsv`
- `data/output/evidence/rpi_domain_context_v2.tsv`
- `data/output/evidence/rpi_function_context_v2.tsv`

## 输入

必需输入：

- `data/output/edges/rna_protein_edges_v1.tsv`
- `data/output/evidence/rna_protein_evidence_v1.tsv`
- `data/output/protein/protein_domains_interpro_v1.tsv`
- `data/processed/protein_master_v6_clean.tsv`

可选输入（用于位点与方法补强，自动扫描 `data/raw/rpi/`）：

- `starBase/ENCORI` 快照（推荐：`data/raw/rpi/starbase_human.tsv`）
- `RNAInter` 快照（推荐：`data/raw/rpi/rnainter_human.tsv.gz`）
- `NPInter` 快照（推荐：`data/raw/rpi/npinter_human.tsv.gz`）

## 输出字段

### 1) 位点上下文：`rpi_site_context_v2.tsv`

`context_id,edge_id,reference,method_type,method_subtype,site_type,chromosome,start,end,strand,transcript_id,gene_id,cell_context,support_count,source,source_version,fetch_date`

### 2) 结构域上下文：`rpi_domain_context_v2.tsv`

`context_id,edge_id,protein_id,domain_class,interpro_id,pfam_id,domain_name,domain_start,domain_end,source,source_version,fetch_date`

### 3) 功能上下文：`rpi_function_context_v2.tsv`

`context_id,edge_id,protein_id,function_relation,inference_basis,evidence_snippet,source,source_version,fetch_date`

## 方法说明

- **RPI 方法类型标准化**：归一到 `CLIP / RIP / CHIRP / RAP / PARIS / OTHER`。
- **结合位点**：优先从 `reference` 回连原始快照（如 starBase clusterID）提取基因组/转录本坐标。
- **蛋白结构域**：使用 InterPro/Pfam 映射并归类为 `RRM / KH / CCCH / PFAM`。
- **功能关系**：基于 UniProt 注释关键词推断 `transcription_regulation / translation_regulation / splicing_regulation / post_transcriptional_regulation`。

## 运行

```bash
bash pipelines/rpi_site_domain_enrichment/run.sh
```

流程：
1. preflight 检查（缺输入即中断并输出 `blocked_missing_inputs`）
2. 最小样本 (`--limit 200`)
3. 全量构建
4. 三张表合同校验
5. 生成 manifest
6. 汇总 gates（覆盖率 + contract）

## 报告

- `rpi_site_domain_enrichment_v2.blocked_or_ready.json`
- `rpi_site_domain_enrichment_v2.sample.metrics.json`
- `rpi_site_domain_enrichment_v2.metrics.json`
- `rpi_site_context_v2.validation.json`
- `rpi_domain_context_v2.validation.json`
- `rpi_function_context_v2.validation.json`
- `rpi_site_domain_enrichment_v2.manifest.json`
- `rpi_site_domain_enrichment_v2.gates.json`

## 验收门槛（任务 C）

- `edge_id` 连接率 >= 0.99
- 方法字段覆盖率 >= 0.90
- 位点/结构域任一覆盖率 >= 0.60
- 合同校验 PASS

## 外部数据与校验

若缺输入阻塞，查看 `blocked_or_ready` 报告中的下载清单与 `sha256sum` 命令。
