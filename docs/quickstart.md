# Quickstart: Download Datasets

## 1) Build or refresh release index

```bash
python3 scripts/build_release_index.py
```

## 2) Download one product (latest)

```bash
# RNA
python3 scripts/download_dataset.py --product rna --version latest --merge-chunks

# Molecule
python3 scripts/download_dataset.py --product molecule --version latest

# Interaction (chunked assets recommended)
python3 scripts/download_dataset.py --product interaction --version latest --merge-chunks

# Protein (repository snapshot URLs)
python3 scripts/download_dataset.py --product protein --version latest
```

## 3) Optional decompression

```bash
python3 scripts/download_dataset.py --product interaction --version latest --merge-chunks --decompress
```

> `.zst` decompression requires `zstd` command installed locally.

## 4) Validate release consistency locally

```bash
python3 scripts/validate_release_index.py \
  --index release/index.json \
  --schema release/schema/index.schema.json \
  --repo-root .

python3 scripts/check_release_consistency.py --index release/index.json --out release/consistency_report.json
```

## 5) Run release metadata regression tests

```bash
pytest -q tests/release
```
