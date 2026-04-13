# Interaction L2 Dataset Delivery (2026-04-13)

## Scope
Interaction L2 result tables for:
- PPI semantic enrichment
- PSI activity/structure enrichment
- RPI site/domain/function enrichment
- Cross-validation & aggregate score
- Ontology mapping

## Delivery policy
- Keep pipeline code/contracts/QA reports in git (reviewable + reproducible).
- Keep large TSV result artifacts out of git (`data/output/**`).
- Deliver full result artifacts via GitHub Release assets (`*.tsv.zst`) with sha256 manifest.

## Why not commit all TSV directly to git
GitHub blob size/history constraints make direct commit of multi-GB tables non-sustainable.
Release assets are used as the canonical distribution channel for full dataset artifacts.

## Canonical reports
- `pipelines/interaction_release_local/reports/interaction_l2_v1.release_assets_manifest.json`
- `pipelines/interaction_release_local/reports/interaction_l2_v1.release_assets_delivery.json`

## Chunked assets
Most tables are single-file assets. Very large compressed assets are published in chunks:
- `psi_activity_context_v2.tsv.zst.part.000 ... part.009`
- `psi_structure_evidence_v2.tsv.zst.part.000 ... part.002`
- `interaction_aggregate_score_v2.tsv.zst.part.000 ... part.003`

## Download & reconstruct example
```bash
# 1) download all assets for this release
gh release download interaction-l2-v1 -D dist/interaction-l2-v1

# 2) reconstruct one chunked asset
cat dist/interaction-l2-v1/psi_activity_context_v2.tsv.zst.part.0* > \
  dist/interaction-l2-v1/psi_activity_context_v2.tsv.zst

# 3) decompress to raw TSV
zstd -d dist/interaction-l2-v1/psi_activity_context_v2.tsv.zst \
  -o data/output/evidence/psi_activity_context_v2.tsv
```
