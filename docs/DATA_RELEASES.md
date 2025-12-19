# Data releases

A **data release** is a snapshot of output tables + a manifest that makes the build reproducible.

## What gets released

- Entity tables: protein, rna, molecule
- Edge tables: ppi, psi, rpi, pathway membership, etc.
- Evidence table: per-edge provenance and scoring
- `manifest.json`: versions, checksums, row counts, build timestamp, git commit

## What never goes into git

- large raw inputs (RNAcentral id_mapping, Ensembl FASTA, etc.)
- large intermediate files
- large output artifacts

Use GitHub Releases or object storage for artifacts.
