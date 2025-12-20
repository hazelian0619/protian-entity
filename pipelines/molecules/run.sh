#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

# Inputs (download separately; see docs)
CHEMREPS_TSV="data/raw/molecules/chembl_36_chemreps.txt"
CHEMBL_DB="data/raw/molecules/chembl_36/chembl_36.db"

# Outputs (large; do not commit to git)
OUT_M1_DB="data/output/molecules/molecules_m1.sqlite"
OUT_M3_DB="data/output/molecules/chembl_m3.sqlite"

mkdir -p data/output/molecules
mkdir -p pipelines/molecules/reports/m1 pipelines/molecules/reports/m2 pipelines/molecules/reports/m3

echo "[INFO] M1: build molecule entity sqlite"
python3 pipelines/molecules/scripts/m1_build_sqlite.py \
  --input "$CHEMREPS_TSV" \
  --db "$OUT_M1_DB" \
  --outdir pipelines/molecules/reports/m1

echo "[INFO] M2: compute RDKit descriptors (writes into M1 sqlite)"
python3 pipelines/molecules/scripts/m2_compute_rdkit_props.py \
  --db "$OUT_M1_DB" \
  --source-table molecule_entity_strict_chembl36 \
  --outdir pipelines/molecules/reports/m2 \
  --table molecule_props_rdkit_v1 \
  --limit 0 \
  --skip-existing

echo "[INFO] M3: extract ChEMBL PSI evidence + aggregate edges"
python3 pipelines/molecules/scripts/m3_check_chembl_db.py \
  --chembl-db "$CHEMBL_DB"

python3 pipelines/molecules/scripts/m3_extract_chembl_psi.py \
  --chembl-db "$CHEMBL_DB" \
  --out-db "$OUT_M3_DB" \
  --outdir pipelines/molecules/reports/m3 \
  --qc-prefix qc_m3_full \
  --types IC50,Ki,Kd,EC50 \
  --min-confidence 6 \
  --require-inchikey \
  --require-uniprot \
  --limit 0

echo "[INFO] M3: quality gates (join coverage against M1/M2)"
python3 pipelines/molecules/scripts/m3_quality_gates.py \
  --m3-db "$OUT_M3_DB" \
  --m1-db "$OUT_M1_DB" \
  --json-out pipelines/molecules/reports/m3/molecules_m3_v1.gates.json \
  --refs-sample 1000

echo "[DONE] molecules M1/M2/M3 pipeline finished."

