# Protian Entity — Bio-Entity KG Foundation

This repository is a **multi-product data foundation** for:

- Protein Entity
- RNA Entity
- Small Molecule Entity
- Cross-Entity Interaction (PPI / PSI / RPI)

It is organized as a data product platform with contracts, QA reports, manifests, and release metadata.

## Unified dataset entrypoint

- Product metadata: `products/*/current.json`
- Global index: `release/index.json`

Generate/rebuild the index:

```bash
python3 scripts/build_release_index.py
```

## Download datasets (public-consumable flow)

Use one command per product:

```bash
# RNA
python3 scripts/download_dataset.py --product rna --version latest --merge-chunks

# Molecule
python3 scripts/download_dataset.py --product molecule --version latest

# Interaction
python3 scripts/download_dataset.py --product interaction --version latest --merge-chunks

# Protein (repository snapshot mode)
python3 scripts/download_dataset.py --product protein --version latest
```

## Verify release consistency

```bash
python3 scripts/validate_release_index.py \
  --index release/index.json \
  --schema release/schema/index.schema.json \
  --repo-root .

python3 scripts/check_release_consistency.py \
  --index release/index.json \
  --out release/consistency_report.json
```

This checks manifest/table/validation row-count consistency where available.

Optional local regression tests:

```bash
pytest -q tests/release
```

## Key directories

- `products/` — product metadata and current release pointers
- `release/` — unified release index + consistency outputs
- `pipelines/` — ETL implementations, contracts, reports
- `scripts/` — release index build, dataset download, consistency checks
- `docs/` — architecture, release policy, quickstart
- `tools/` — validation and manifest helpers

## Documentation

- Architecture: `docs/architecture.md`
- Release policy: `docs/release-policy.md`
- Download quickstart: `docs/quickstart.md`

## License

MIT — see `LICENSE`.
