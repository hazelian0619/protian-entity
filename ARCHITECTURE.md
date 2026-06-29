# Architecture

## Repository Role

**protian-entity is the development engine** for the Bio-Entity KG project. All pipelines are built here, data is produced and QA-tested here, and release packages are staged here before promotion to graph-rag.

This repository is NOT the external download entry point. To download published datasets, use [graph-rag](https://github.com/hazelian0619/graph-rag).

## Sister Repository

[graph-rag](https://github.com/hazelian0619/graph-rag) is the sole public-facing repository. Only curated, QA-passed results from this repository are promoted there for public delivery.

## Data Promotion Path

```
protian-entity (this repo)               graph-rag (public)
  pipeline produces                        curated results committed to
  QA gates pass             ------>        <domain>/release/external/<version>/
  manifest generated                       product pointer updated
                                           release index rebuilt
```

After QA gates pass, curated result tables are copied to graph-rag's `release/external/` directory and committed there. This repository retains the original pipeline outputs and evidence.

## Product Lines

Four product lines, all with pipelines in this repository:

| Product | Scope | Pipelines |
|---|---|---|
| Protein | Human protein L1 entity + provenance | protein, protein_domains, protein_isoform, protein_kegg, protein_pdb, protein_physchem, protein_physchem_extended, protein_provenance |
| RNA | Human RNA L1/L2 | rna, rna_expression_rbp, rna_external_xref, rna_id_canonical, rna_pdb, rna_rfam_structure, rna_rpi, rna_structure_models, rna_type_features, rna_xref_enst_urs |
| Molecule | Small molecule L1/L2 | molecule_3d_experimental_linker, molecule_3d_registry, molecule_activity_fusion, molecule_physchem_descriptors, molecule_pk_tox_layer, molecule_pk_tox_v2, molecule_pubchem_direct_xref, molecule_semantic_layer, molecule_structure_identifiers, molecule_xref_enrichment_v2, molecule_zinc_direct_xref, molecule_zinc_direct_xref_v2 |
| Interaction | Cross-entity L2 (PPI/PSI/RPI) | edges_ppi, ppi_semantic_enrichment, psi_activity_structure_enrichment, psi_condition_enrichment, rpi_site_domain_enrichment, interaction_cross_validation, interaction_ontology_mapping, interaction_readiness, interaction_release_local |

## Repository Layout

```
data/           # raw downloads, intermediate products, QA outputs
pipelines/      # 39 pipelines (contracts, scripts, reports) - the core engine
scripts/        # release tooling and operational scripts
  legacy/       # archived exploratory scripts (step1-11 era)
products/       # product definitions and version pointers
release/        # release index, consistency reports, schema
build/          # staging area for release packages (gitignored)
dist/           # packaged distributions (gitignored)
docs/           # development documentation
tools/          # table validation utilities
tests/          # test suite
```

## Download URLs

All `download_url` fields in `products/*/current.json` point to graph-rag raw URLs. This repository does not serve as a download source for published data.
