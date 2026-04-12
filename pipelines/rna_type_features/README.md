# RNA Type Feature Pack (Assistant B)

产物：
- `data/output/rna_lnc_entries_v1.tsv`
- `data/output/rna_trna_features_v1.tsv`
- `data/output/rna_rrna_loci_v1.tsv`

输入：
- `data/output/rna_master_v1.tsv`
- `data/raw/rna/ensembl/Homo_sapiens.GRCh38.115.chr.gtf.gz`
- `data/raw/rna/rnacentral/id_mapping.tsv.gz`

## Run (v1 full pack)

```bash
bash pipelines/rna_type_features/run.sh
```

## tRNA anticodon optimization v2 (balanced)

新增不覆盖 v1 的优化流程：
- 输出：`data/output/rna_trna_features_v2.tsv`
- 冲突审计：`pipelines/rna_type_features/reports/rna_trna_anticodon_conflicts_v2.tsv`
- 指标：`pipelines/rna_type_features/reports/rna_trna_features_v2.metrics.json`

运行：

```bash
bash pipelines/rna_type_features/run_trna_v2.sh
```

## Reports

- `pipelines/rna_type_features/reports/rna_type_features_v1.metrics.json`
- `pipelines/rna_type_features/reports/rna_lnc_entries_v1.validation.json`
- `pipelines/rna_type_features/reports/rna_trna_features_v1.validation.json`
- `pipelines/rna_type_features/reports/rna_rrna_loci_v1.validation.json`
- `pipelines/rna_type_features/reports/*.manifest.json`
