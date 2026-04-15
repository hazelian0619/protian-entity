#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

mkdir -p data/output/protein pipelines/protein_physchem_extended/reports

DATA_VERSION="${DATA_VERSION:-kg-data-local}"

OUTPUT_TABLE="data/output/protein/protein_physchem_extended_v1.tsv"
CONTRACT="pipelines/protein_physchem_extended/contracts/protein_physchem_extended_v1.json"
QA_REPORT="pipelines/protein_physchem_extended/reports/protein_physchem_extended_v1.qa.json"
VALIDATION_REPORT="pipelines/protein_physchem_extended/reports/protein_physchem_extended_v1.validation.json"
MANIFEST_REPORT="pipelines/protein_physchem_extended/reports/protein_physchem_extended_v1.manifest.json"

echo "[1/3] Build extended physchem table + QA..."
python3 pipelines/protein_physchem_extended/scripts/build_protein_physchem_extended.py \
  --input data/processed/protein_master_v6_clean.tsv \
  --source-versions data/output/protein/protein_source_versions_v1.tsv \
  --core-table data/output/protein/protein_physchem_v1.tsv \
  --output "$OUTPUT_TABLE" \
  --qa "$QA_REPORT"

echo "[2/3] Validate contract..."
python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUTPUT_TABLE" \
  --out "$VALIDATION_REPORT"

echo "[3/3] Generate manifest..."
python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$MANIFEST_REPORT" \
  "$OUTPUT_TABLE" \
  "$CONTRACT" \
  "$VALIDATION_REPORT" \
  "$QA_REPORT"

echo "[DONE] protein_physchem_extended_v1"
