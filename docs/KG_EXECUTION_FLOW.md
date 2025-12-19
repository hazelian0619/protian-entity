# KG execution flow (industrialized)

This repository is evolving from a protein-centric dataset into a multi-entity, evidence-backed Knowledge Graph.

## Why the execution flow matters

Your project design explicitly requires:

- **L1 entity standardization** (protein / RNA / small molecule) with multi-modal attributes.
- **L2 relationship discovery** with **evidence provenance + confidence scoring**.
- Frequent upstream updates (UniProt, PDB, RNAcentral, STRING, ChEMBL, etc.), which implies:
  - version locking
  - reproducible rebuilds
  - regression checks and rollback

Therefore, the industrial-grade approach is:

1) **Define data contracts** (schemas + constraints + thresholds) as code.
2) **Build ETL pipelines** that always emit the same set of tables.
3) **Validate outputs** with quality gates.
4) **Publish artifacts** with a manifest (versions, checksums, counts).

## Layers mapping (L1–L6)

- **L1 (Now)**: entity tables (protein_master, rna_master, molecule_master) + core attributes.
- **L2 (Next)**: edge tables (ppi/psi/rpi) + evidence table (many-to-one to edges).
- **L3 (Later)**: ontology alignment tables (GO/KEGG/ChEBI) and semantic constraints.
- **L4–L6 (Later)**: query templates, reasoning, embeddings, link prediction.

## Table-first, graph-second

Treat TSV/Parquet tables as the **source of truth** and load them into Neo4j/RDF/etc as a *serving layer*.
This keeps rebuilds reproducible and makes QA measurable.
