# Release Policy

## Public release channels

- `stable`: latest verified product tag per line
- `snapshot`: repository-contained artifacts when no GitHub release is published

## Product release contract

Every public release should include:

1. Data artifacts (`*.tsv.gz` / `*.tsv.zst`)
2. `manifest.json` with checksums + row counts
3. Validation and QA reports (`*.validation.json`, `*.qa.json`, `*.gates.json`)
4. Optional `SHA256SUMS.txt`

## Required checks before announcing release

Run:

```bash
python3 scripts/build_release_index.py
python3 scripts/validate_release_index.py --index release/index.json --schema release/schema/index.schema.json --repo-root .
python3 scripts/check_release_consistency.py --index release/index.json --out release/consistency_report.json
```

A release is considered externally consumable only when consistency status is PASS (or documented exceptions are approved).

## Consistency modes

Products may declare a consistency mode in `products/<product>/current.json`:

- `strict_local` (default): verify local table row counts against manifest + validation reports.
- `release_assets`: verify release-asset metadata integrity (asset names/sha256/release availability), without requiring local full TSV presence.

`release_assets` is recommended for products where canonical full artifacts are distributed via GitHub Releases instead of checked-in local tables.

## Asset rules

- Large data artifacts should be distributed via Release assets.
- Git should keep contracts, reports, manifests, and reproducibility metadata.
- Chunked assets (`*.part.000`) are allowed for very large files; reconstruction instructions must be provided.
