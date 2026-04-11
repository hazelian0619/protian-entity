# Quality gates (self-check metrics)

This document defines *objective* acceptance criteria for each table.
Pipelines must generate validation reports (JSON) and fail fast if required thresholds are not met.

## Common gates (all tables)

- **Schema**: required columns present; no unexpected missing columns.
- **Primary key**: non-empty rate = 100%; uniqueness = 100%.
- **Type sanity**: numeric columns parse; date columns ISO format.
- **Provenance**: `source`, `source_version`, `fetch_date` non-empty.

## RNA master v1 (entity)

Applied to `data/output/rna_master_v1.tsv`:

- **rna_type**: only `mirna` or `mrna`.
- **taxon_id**: 100% == 9606.
- **sequence**: 100% non-empty; characters subset of {A,C,G,U,N}.
- **ID formats**:
  - for `mirna`: `rna_id` starts with `URS` and ends with `_9606`.
  - for `mrna`: `rna_id` starts with `ENST` and ends with `_9606`.
- **Gene mapping coverage (mrna)**:
  - `ensembl_gene_id` non-empty rate >= 99%.
  - `symbol` non-empty rate >= 99%.

## Protein master v6 (entity)

Applied to `data/processed/protein_master_v6_clean.tsv`:

- `uniprot_id` unique and non-empty.
- `taxon_id` must be `9606` for v1 (human-only release).
- `sequence` non-empty; AA alphabet validity.
- `alphafold_id` coverage target (project-defined; current README claims ~99.7%).

## Edge + evidence (L2)

(Template; implement once edge tables are standardized)

- edges: `(src_id, dst_id, predicate)` uniqueness within a release.
- evidence: each row must reference an existing edge id.
- evidence method controlled vocabulary (PSI-MI for PPI, etc.).
