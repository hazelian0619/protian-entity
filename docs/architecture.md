# Bio-Entity KG Foundation — Architecture

This repository is organized as a **data product platform** with four public product lines:

- Protein Entity
- RNA Entity
- Small Molecule Entity
- Cross-Entity Interaction (PPI/PSI/RPI)

## Control Plane vs Data Plane

### Control Plane (in git)
- `products/*/current.json`: product-level release pointers
- `release/index.json`: global release index (machine-readable)
- `release/schema/index.schema.json`: release index schema contract
- `contracts/` + `pipelines/**/contracts/*.json`: schema and quality contracts
- `pipelines/**/reports/*.json`: validation/QA/gates/manifest evidence
- `scripts/`: release index build, schema/path validation, consistency check, dataset download tooling

### Data Plane (release assets)
Large `.tsv.gz` / `.tsv.zst` artifacts are distributed via GitHub Releases, not committed to git by default.

## Product Model

Each product is modeled with:

1. `products/<product>/product.json` (ownership, SLA, contract roots)
2. `products/<product>/current.json` (latest distribution metadata)

`release/index.json` is generated from these files and acts as the unified public entry point.

## Release Integrity Model

Release integrity is enforced by:

1. **Manifest checks** (rows, sha256, paths)
2. **Validation checks** (`*.validation.json`)
3. **Cross-check script**: `scripts/check_release_consistency.py`

If manifest rows, actual rows, and validation rows diverge, the consistency check fails.
