#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

python3 pipelines/rna/scripts/01_extract_human_mirna.py
python3 pipelines/rna/scripts/02_map_mirna_to_urs.py
python3 pipelines/rna/scripts/03_add_gene_mapping.py
python3 pipelines/rna/scripts/04_finalize_mirna_table_v3.py

python3 pipelines/rna/scripts/05_extract_seed_genes.py
python3 pipelines/rna/scripts/07_extract_mrna_from_ensembl.py
python3 pipelines/rna/scripts/08_finalize_mrna_table.py

python3 pipelines/rna/scripts/09a_repair_rna_subtables.py
python3 pipelines/rna/scripts/09b_merge_rna_fixed.py

# optional xref (maps ENST_9606 -> URS_9606)
python3 pipelines/rna/scripts/11_mrna_enst_to_urs_xref.py || true

echo "[DONE] RNA v1 pipeline finished."
