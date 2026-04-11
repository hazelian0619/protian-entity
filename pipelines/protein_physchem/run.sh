#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

mkdir -p data/output/protein pipelines/protein_physchem/reports

echo "[1/4] Dry run on minimal sample (n=20)..."
python3 pipelines/protein_physchem/scripts/build_protein_physchem.py \
  --input data/processed/protein_master_v6_clean.tsv \
  --output data/output/protein/protein_physchem_v1.sample.tsv \
  --report pipelines/protein_physchem/reports/protein_physchem_v1.sample.metrics.json \
  --limit 20

echo "[2/4] Full run on all rows..."
python3 pipelines/protein_physchem/scripts/build_protein_physchem.py \
  --input data/processed/protein_master_v6_clean.tsv \
  --output data/output/protein/protein_physchem_v1.tsv \
  --report pipelines/protein_physchem/reports/protein_physchem_v1.metrics.json

echo "[3/4] Contract validation..."
python3 tools/kg_validate_table.py \
  --contract pipelines/protein_physchem/contracts/protein_physchem_v1.json \
  --table data/output/protein/protein_physchem_v1.tsv \
  --out pipelines/protein_physchem/reports/protein_physchem_v1.validation.json

echo "[4/4] Done."
