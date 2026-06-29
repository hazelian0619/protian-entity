# Bio-Entity KG Foundation

The development engine for 4 bio-entity product lines:

- **Protein Entity**
- **RNA Entity**
- **Small Molecule Entity**
- **Cross-Entity Interactions** (PPI / PSI / RPI)

> **Looking for downloadable datasets?** Published data lives in [graph-rag](https://github.com/hazelian0619/graph-rag), the sole public-facing repository. This repo contains pipelines, raw data, and QA evidence used to produce those datasets.

---

## 1) What to use first (external users)

If you only care about downloading datasets, start here:

1. **Unified product index**: `release/index.json`
2. **Product pointers**: `products/*/current.json`
3. **One-command downloader**: `scripts/download_dataset.py`

### Quick download

```bash
# Protein
python3 scripts/download_dataset.py --product protein --version latest

# RNA
python3 scripts/download_dataset.py --product rna --version latest --merge-chunks

# Molecule
python3 scripts/download_dataset.py --product molecule --version latest

# Interaction
python3 scripts/download_dataset.py --product interaction --version latest --merge-chunks
```

---

## 2) Dataset catalog (clearly separated by product)

| Product | Scope | Latest pointer | Distribution mode | Primary manifest |
|---|---|---|---|---|
| Protein | Human protein L1 entity + provenance | `products/protein/current.json` | `repository_snapshot` | `pipelines/protein/reports/protein_master_v6.manifest.json` |
| RNA | Human RNA L1/L2 | `products/rna/current.json` | `github_release` | `build/releases/rna-l1l2-v2/manifest.json` |
| Molecule | Small-molecule L1/L2 | `products/molecule/current.json` | `github_release` | `dist/molecule-l1l2-v2/manifest.json` |
| Interaction | Cross-entity L2 (PPI/PSI/RPI) | `products/interaction/current.json` | `github_release` | `pipelines/interaction_release_local/reports/interaction_l2_v1.release_assets_manifest.json` |

> Source of truth for latest versions: `release/index.json`

---

## 3) Expected local download layout

After running downloader:

```text
downloads/
  protein/<version>/
  rna/<version>/
  molecule/<version>/
  interaction/<version>/
```

For large assets:
- use `--merge-chunks` to merge `*.part.000` files
- use `--decompress` to decompress `.gz` / `.zst`

---

## 4) Quality & release integrity checks

Run before trusting a release snapshot:

```bash
# rebuild unified index
python3 scripts/build_release_index.py

# validate release metadata schema + path integrity
python3 scripts/validate_release_index.py \
  --index release/index.json \
  --schema release/schema/index.schema.json \
  --repo-root .

# consistency checks (manifest / validation / asset integrity)
python3 scripts/check_release_consistency.py \
  --index release/index.json \
  --out release/consistency_report.json
```

Optional regression tests:

```bash
pytest -q tests/release
```

---

## 5) Repository architecture (for maintainers)

```text
products/   # product metadata and latest release pointers
release/    # unified index + schema + consistency report
pipelines/  # per-domain ETL contracts/reports/scripts
docs/       # architecture, release policy, quickstart
scripts/    # index build / validation / consistency / download tooling
tests/      # release metadata regression tests
```

Key docs:
- `docs/architecture.md`
- `docs/release-policy.md`
- `docs/quickstart.md`

---

## 6) Notes on legacy vs public entrypoint

This repository contains historical pipeline modules. For external consumption, **do not start from random pipeline folders**.

Always start from:
1. `release/index.json`
2. `products/*/current.json`
3. `scripts/download_dataset.py`

This is the stable public contract.

---

## License

MIT — see `LICENSE`.
