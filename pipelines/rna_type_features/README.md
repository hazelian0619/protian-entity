# RNA Type Feature Pack (Assistant B)

产物：
- `data/output/rna_lnc_entries_v1.tsv`
- `data/output/rna_trna_features_v1.tsv`
- `data/output/rna_rrna_loci_v1.tsv`

输入：
- `data/output/rna_master_v1.tsv`
- `data/raw/rna/ensembl/Homo_sapiens.GRCh38.115.chr.gtf.gz`
- `data/raw/rna/rnacentral/id_mapping.tsv.gz`

## Run

```bash
bash pipelines/rna_type_features/run.sh
```

## Reports

- `pipelines/rna_type_features/reports/rna_type_features_v1.metrics.json`
- `pipelines/rna_type_features/reports/rna_lnc_entries_v1.validation.json`
- `pipelines/rna_type_features/reports/rna_trna_features_v1.validation.json`
- `pipelines/rna_type_features/reports/rna_rrna_loci_v1.validation.json`
- `pipelines/rna_type_features/reports/*.manifest.json`
