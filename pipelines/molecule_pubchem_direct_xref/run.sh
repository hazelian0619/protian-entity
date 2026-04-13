#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

PIPELINE_DIR="pipelines/molecule_pubchem_direct_xref"

IN_CORE="${CORE_TSV_PATH:-data/output/molecules/molecule_xref_core_v2.tsv}"
IN_M1_DB="${M1_DB_PATH:-../12182/out/m1/molecules_m1.sqlite}"
OUT_TABLE="data/output/molecules/molecule_xref_pubchem_enhanced_v1.tsv"

CONTRACT="$PIPELINE_DIR/contracts/molecule_xref_pubchem_enhanced_v1.json"
REPORT_BUILD="$PIPELINE_DIR/reports/molecule_xref_pubchem_enhanced_v1.build.json"
REPORT_CONFLICT="$PIPELINE_DIR/reports/molecule_xref_pubchem_enhanced_v1.conflict_audit.json"
REPORT_VALIDATION="$PIPELINE_DIR/reports/molecule_xref_pubchem_enhanced_v1.validation.json"
REPORT_QA="$PIPELINE_DIR/reports/molecule_xref_pubchem_enhanced_v1.qa.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/molecule_xref_pubchem_enhanced_v1.manifest.json"

SMOKE_ROWS="${SMOKE_MAX_ROWS:-1500}"
SMOKE_TABLE="$PIPELINE_DIR/.cache/molecule_xref_pubchem_enhanced_v1.smoke.tsv"
SMOKE_BUILD="$PIPELINE_DIR/reports/molecule_xref_pubchem_enhanced_v1.smoke.build.json"
SMOKE_CONFLICT="$PIPELINE_DIR/reports/molecule_xref_pubchem_enhanced_v1.smoke.conflict_audit.json"
SMOKE_VALIDATION="$PIPELINE_DIR/reports/molecule_xref_pubchem_enhanced_v1.smoke.validation.json"
SMOKE_QA="$PIPELINE_DIR/reports/molecule_xref_pubchem_enhanced_v1.smoke.qa.json"

WORKERS="${PUBCHEM_WORKERS:-8}"
QPS="${PUBCHEM_QPS:-5.0}"
TIMEOUT="${PUBCHEM_TIMEOUT:-20}"
RETRIES="${PUBCHEM_RETRIES:-3}"

DATA_VERSION="${DATA_VERSION:-kg-data-local}"

mkdir -p "$(dirname "$OUT_TABLE")" "$PIPELINE_DIR/reports" "$PIPELINE_DIR/.cache"

missing=()
[[ -f "$IN_CORE" ]] || missing+=("$IN_CORE")
[[ -f "$IN_M1_DB" ]] || missing+=("$IN_M1_DB")
if (( ${#missing[@]} > 0 )); then
  echo "[ERROR] missing required inputs:" >&2
  printf '  - %s\n' "${missing[@]}" >&2
  exit 2
fi

echo "[INFO] core_tsv=$IN_CORE"
echo "[INFO] m1_db=$IN_M1_DB"
echo "[INFO] workers=$WORKERS qps=$QPS timeout=$TIMEOUT retries=$RETRIES"

echo "[STEP] smoke build"
python3 "$PIPELINE_DIR/scripts/build_molecule_xref_pubchem_enhanced_v1.py" \
  --core-tsv "$IN_CORE" \
  --m1-db "$IN_M1_DB" \
  --out "$SMOKE_TABLE" \
  --build-report "$SMOKE_BUILD" \
  --conflict-audit "$SMOKE_CONFLICT" \
  --max-rows "$SMOKE_ROWS" \
  --workers "$WORKERS" \
  --qps "$QPS" \
  --timeout "$TIMEOUT" \
  --retries "$RETRIES"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$SMOKE_TABLE" \
  --out "$SMOKE_VALIDATION"

python3 "$PIPELINE_DIR/scripts/qa_molecule_xref_pubchem_enhanced_v1.py" \
  --core-tsv "$IN_CORE" \
  --enhanced-tsv "$SMOKE_TABLE" \
  --build-report "$SMOKE_BUILD" \
  --conflict-audit "$SMOKE_CONFLICT" \
  --out "$SMOKE_QA" \
  --allow-row-subset 1

echo "[STEP] full build"
python3 "$PIPELINE_DIR/scripts/build_molecule_xref_pubchem_enhanced_v1.py" \
  --core-tsv "$IN_CORE" \
  --m1-db "$IN_M1_DB" \
  --out "$OUT_TABLE" \
  --build-report "$REPORT_BUILD" \
  --conflict-audit "$REPORT_CONFLICT" \
  --workers "$WORKERS" \
  --qps "$QPS" \
  --timeout "$TIMEOUT" \
  --retries "$RETRIES"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUT_TABLE" \
  --out "$REPORT_VALIDATION"

python3 "$PIPELINE_DIR/scripts/qa_molecule_xref_pubchem_enhanced_v1.py" \
  --core-tsv "$IN_CORE" \
  --enhanced-tsv "$OUT_TABLE" \
  --build-report "$REPORT_BUILD" \
  --conflict-audit "$REPORT_CONFLICT" \
  --out "$REPORT_QA"

python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$IN_CORE" \
  "$IN_M1_DB" \
  "$OUT_TABLE" \
  "$CONTRACT" \
  "$REPORT_BUILD" \
  "$REPORT_CONFLICT" \
  "$REPORT_VALIDATION" \
  "$REPORT_QA"

echo "[DONE] molecule_pubchem_direct_xref finished (smoke + full + validation + qa + manifest)."
