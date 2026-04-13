# Interaction L2 PR Scope (2026-04-13)

This PR is **interaction-only** and follows the same review pattern as prior Protein/RNA submissions:
- commit pipeline code/contracts/reports/manifests to git;
- keep large evidence TSV outputs out of git (`data/output/**`);
- use reports + manifests for reproducibility and review.

## Included pipelines
- `pipelines/ppi_semantic_enrichment`
- `pipelines/psi_activity_structure_enrichment`
- `pipelines/rpi_site_domain_enrichment`
- `pipelines/interaction_cross_validation`
- `pipelines/interaction_ontology_mapping`
- `pipelines/interaction_release_local`

## Included quality evidence (JSON)
- A: PPI semantic enrichment QA/gates
- B: PSI activity+structure QA/gates
- C: RPI site/domain/function gates
- D: cross-validation + aggregate-score QA/gates
- E: ontology mapping gates
- local release package summary/gates/validation/manifest

## Explicitly excluded
- `data/output/**` large TSV artifacts
- downloaded raw sources under pipeline `data/raw/**`
- smoke/sample TSV files in reports

## Acceptance snapshot
- A PASS; B PASS; C PASS; D PASS; E PASS
- release-local package PASS
