# Molecules sources & license notes

This repository ships **derived** small-molecule datasets and PSI evidence datasets.

## Primary source (ChEMBL)

The molecules + PSI datasets in this repo are derived from **ChEMBL** (release 36 in v1).

What we use:

- `chembl_36_chemreps.txt` (chemreps): used for molecule structure identifiers (M1)
- `chembl_36.db` (SQLite dump): used for PSI evidence extraction (M3)

What we publish:

- `molecules_m1.sqlite` (derived entity table + mappings)
- `molecule_props_rdkit_v1` (derived RDKit descriptors written into the M1 sqlite)
- `chembl_m3.sqlite` (derived PSI evidence + aggregated edges)

## License / attribution

ChEMBL is an EMBL-EBI resource. Redistribution and reuse may be subject to the ChEMBL license terms.

Before redistributing release assets, confirm your compliance with ChEMBL’s current licensing and attribution requirements (and include the recommended citations).

This repository’s **code** is independent of ChEMBL; the **data assets** published in GitHub Releases are derived artifacts.

## Reproducibility promise

- We publish QA JSON + manifest JSON alongside release assets.
- Users can verify file integrity via sha256 in the manifest.
- The pipeline can rebuild the same table schema from locally downloaded upstream inputs.

