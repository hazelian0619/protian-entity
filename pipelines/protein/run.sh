#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

INPUT_TABLE="data/processed/protein_master_v6_clean.tsv"
CONTRACT="pipelines/protein/contracts/protein_master_v6.json"
REPORT_VALIDATION="pipelines/protein/reports/protein_master_v6.validation.json"
REPORT_MANIFEST="pipelines/protein/reports/protein_master_v6.manifest.json"

DATA_VERSION="${DATA_VERSION:-kg-data-local}"

if [[ ! -f "$INPUT_TABLE" ]]; then
  echo "[ERROR] missing input table: $INPUT_TABLE" >&2
  exit 2
fi

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$INPUT_TABLE" \
  --out "$REPORT_VALIDATION"

python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$INPUT_TABLE"

echo "[DONE] Protein v6 QA + manifest finished."
