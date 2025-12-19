# RNA pipeline (v1)

This pipeline builds the **RNA entity tables** for the KG (miRNA + mRNA transcripts), following the project-level RNA spec.

- Canonical RNA docs: `docs/rna/DATA_DICTIONARY_RNA.md`, `docs/rna/RNA_SOURCES_AND_VERSIONS.md`
- Main outputs (not committed):
  - `data/output/rna_master_v1.tsv`
  - `data/output/rna_master_mirna_v1.tsv`
  - `data/output/rna_master_mrna_v1.tsv`

## Inputs

Required local files (large files are intentionally not tracked in git):

- `data/raw/rna/mirbase/mature.fa`
- `data/raw/rna/rnacentral/id_mapping.tsv`
- `data/raw/rna/ensembl/Homo_sapiens.GRCh38.cdna.all.fa`

Seed genes source (any one):
- `data/raw/protein_master_v6_clean.tsv` (optional), or
- `data/processed/protein_master_v6_clean.tsv` (tracked in this repo)

## Run

From repo root:

```bash
bash pipelines/rna/run.sh
```

## QA

After running, generate QA + manifest:

```bash
python3 tools/kg_validate_table.py \
  --contract pipelines/rna/contracts/rna_master_v1.json \
  --table data/output/rna_master_v1.tsv \
  --out pipelines/rna/reports/rna_master_v1.validation.json

python3 tools/kg_make_manifest.py \
  --data-version kg-data-local \
  --out pipelines/rna/reports/rna_master_v1.manifest.json \
  data/output/rna_master_v1.tsv
```

## Data release

- GitHub Releases (artifacts): https://github.com/hazelian0619/protian-entity/releases/tag/rna-l1-v1
