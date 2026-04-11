#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

INPUT_TABLE="data/processed/protein_master_v6_clean.tsv"
OUTPUT_TABLE="data/output/protein/pdb_structures_v1.tsv"

PIPELINE_DIR="pipelines/protein_pdb"
CONTRACT="$PIPELINE_DIR/contracts/pdb_structures_v1.json"

REPORT_BUILD="$PIPELINE_DIR/reports/pdb_structures_v1.build.json"
REPORT_AUDIT="$PIPELINE_DIR/reports/pdb_structures_v1.api_audit.json"
REPORT_VALIDATION="$PIPELINE_DIR/reports/pdb_structures_v1.validation.json"
REPORT_QA="$PIPELINE_DIR/reports/pdb_structures_v1.qa.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/pdb_structures_v1.manifest.json"

SMOKE_TABLE="$PIPELINE_DIR/.cache/pdb_structures_v1.smoke.tsv"
SMOKE_BUILD="$PIPELINE_DIR/reports/pdb_structures_v1.smoke.build.json"
SMOKE_AUDIT="$PIPELINE_DIR/reports/pdb_structures_v1.smoke.api_audit.json"
SMOKE_VALIDATION="$PIPELINE_DIR/reports/pdb_structures_v1.smoke.validation.json"
SMOKE_QA="$PIPELINE_DIR/reports/pdb_structures_v1.smoke.qa.json"

DATA_VERSION="${DATA_VERSION:-kg-data-local}"
SMOKE_MAX_UNIQUE_PDB="${SMOKE_MAX_UNIQUE_PDB:-200}"
PDB_API_FAIL_THRESHOLD="${PDB_API_FAIL_THRESHOLD:-0.05}"
PDB_BATCH_SIZE="${PDB_BATCH_SIZE:-10000}"
PDB_TIMEOUT="${PDB_TIMEOUT:-40}"
PDB_RETRIES="${PDB_RETRIES:-3}"
PDB_SLEEP_BETWEEN_BATCHES="${PDB_SLEEP_BETWEEN_BATCHES:-0.0}"

if [[ ! -f "$INPUT_TABLE" ]]; then
  echo "[ERROR] missing input table: $INPUT_TABLE" >&2
  exit 2
fi

mkdir -p "$(dirname "$OUTPUT_TABLE")" "$PIPELINE_DIR/reports" "$PIPELINE_DIR/.cache"

echo "[STEP 1/5] smoke run (max_unique_pdb=${SMOKE_MAX_UNIQUE_PDB})"
python3 "$PIPELINE_DIR/scripts/01_build_pdb_structures_v1.py" \
  --input "$INPUT_TABLE" \
  --output "$SMOKE_TABLE" \
  --report-build "$SMOKE_BUILD" \
  --report-audit "$SMOKE_AUDIT" \
  --max-unique-pdb "$SMOKE_MAX_UNIQUE_PDB" \
  --batch-size "$PDB_BATCH_SIZE" \
  --timeout "$PDB_TIMEOUT" \
  --retries "$PDB_RETRIES" \
  --sleep-between-batches "$PDB_SLEEP_BETWEEN_BATCHES" \
  --fail-threshold "$PDB_API_FAIL_THRESHOLD"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$SMOKE_TABLE" \
  --out "$SMOKE_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_pdb_structures_v1.py" \
  --master "$INPUT_TABLE" \
  --table "$SMOKE_TABLE" \
  --build-report "$SMOKE_BUILD" \
  --out "$SMOKE_QA"

echo "[STEP 2/5] full run"
python3 "$PIPELINE_DIR/scripts/01_build_pdb_structures_v1.py" \
  --input "$INPUT_TABLE" \
  --output "$OUTPUT_TABLE" \
  --report-build "$REPORT_BUILD" \
  --report-audit "$REPORT_AUDIT" \
  --batch-size "$PDB_BATCH_SIZE" \
  --timeout "$PDB_TIMEOUT" \
  --retries "$PDB_RETRIES" \
  --sleep-between-batches "$PDB_SLEEP_BETWEEN_BATCHES" \
  --fail-threshold "$PDB_API_FAIL_THRESHOLD"

echo "[STEP 3/5] contract validation"
python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUTPUT_TABLE" \
  --out "$REPORT_VALIDATION"

echo "[STEP 4/5] QA checks"
python3 "$PIPELINE_DIR/scripts/02_qa_pdb_structures_v1.py" \
  --master "$INPUT_TABLE" \
  --table "$OUTPUT_TABLE" \
  --build-report "$REPORT_BUILD" \
  --out "$REPORT_QA"

echo "[STEP 5/5] manifest"
python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$OUTPUT_TABLE" \
  "$INPUT_TABLE" \
  "$CONTRACT" \
  "$REPORT_BUILD" \
  "$REPORT_VALIDATION" \
  "$REPORT_QA" \
  "$REPORT_AUDIT"

echo "[DONE] protein_pdb v1: smoke + full + validation + QA + manifest finished."
