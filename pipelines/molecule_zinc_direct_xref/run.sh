#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

PIPELINE_DIR="pipelines/molecule_zinc_direct_xref"

resolve_input() {
  local primary="$1"
  shift
  if [[ -n "$primary" && -f "$primary" ]]; then
    echo "$primary"
    return 0
  fi
  for candidate in "$@"; do
    if [[ -f "$candidate" ]]; then
      echo "$candidate"
      return 0
    fi
  done
  echo ""
}

REQ_XREF_V2="${XREF_V2_PATH:-data/output/molecules/molecule_xref_core_v2.tsv}"
REQ_XREF_V1="${XREF_V1_PATH:-data/output/molecules/molecule_xref_core_v1.tsv}"
BASELINE_COVERAGE="${BASELINE_COVERAGE_PATH:-pipelines/molecule_3d_registry/reports/molecule_3d_registry_v1.coverage.json}"

IN_XREF="$(resolve_input "$REQ_XREF_V2" "$REQ_XREF_V1" "../1218/data/output/molecules/molecule_xref_core_v2.tsv" "../1218/data/output/molecules/molecule_xref_core_v1.tsv")"

OUT_TABLE="data/output/molecules/molecule_zinc_xref_v1.tsv"
CONTRACT="$PIPELINE_DIR/contracts/molecule_zinc_xref_v1.json"

SMOKE_MAX_ROWS="${SMOKE_MAX_ROWS:-3000}"
SMOKE_MAX_VENDORS="${SMOKE_MAX_VENDORS:-8}"
FULL_MAX_VENDORS="${FULL_MAX_VENDORS:-}"

ZINC_BASE_URL="${ZINC_BASE_URL:-http://files.docking.org/catalogs/}"
ZINC_TIERS="${ZINC_TIERS:-10 20 30 40 50}"
HTTP_TIMEOUT="${ZINC_HTTP_TIMEOUT:-35}"
HTTP_RETRIES="${ZINC_HTTP_RETRIES:-2}"

DATA_VERSION="${DATA_VERSION:-kg-data-local}"

SMOKE_TABLE="$PIPELINE_DIR/.cache/molecule_zinc_xref_v1.smoke.tsv"
SMOKE_BUILD="$PIPELINE_DIR/reports/molecule_zinc_xref_v1.smoke.build.json"
SMOKE_COVERAGE="$PIPELINE_DIR/reports/molecule_zinc_xref_v1.smoke.coverage.json"
SMOKE_CONFLICT="$PIPELINE_DIR/reports/molecule_zinc_xref_v1.smoke.conflict_audit.json"
SMOKE_VALIDATION="$PIPELINE_DIR/reports/molecule_zinc_xref_v1.smoke.validation.json"
SMOKE_QA="$PIPELINE_DIR/reports/molecule_zinc_xref_v1.smoke.qa.json"
SMOKE_MANUAL="$PIPELINE_DIR/reports/molecule_zinc_xref_v1.smoke.manual_download.json"

REPORT_BUILD="$PIPELINE_DIR/reports/molecule_zinc_xref_v1.build.json"
REPORT_COVERAGE="$PIPELINE_DIR/reports/molecule_zinc_xref_v1.coverage.json"
REPORT_CONFLICT="$PIPELINE_DIR/reports/molecule_zinc_xref_v1.conflict_audit.json"
REPORT_VALIDATION="$PIPELINE_DIR/reports/molecule_zinc_xref_v1.validation.json"
REPORT_QA="$PIPELINE_DIR/reports/molecule_zinc_xref_v1.qa.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/molecule_zinc_xref_v1.manifest.json"
REPORT_MANUAL="$PIPELINE_DIR/reports/molecule_zinc_xref_v1.manual_download.json"

mkdir -p "$(dirname "$OUT_TABLE")" "$PIPELINE_DIR/reports" "$PIPELINE_DIR/.cache"

missing=()
[[ -n "$IN_XREF" ]] || missing+=("$REQ_XREF_V2 or $REQ_XREF_V1")
if (( ${#missing[@]} > 0 )); then
  echo "[ERROR] missing required inputs for molecule_zinc_xref_v1:" >&2
  for f in "${missing[@]}"; do
    echo "  - $f" >&2
  done
  exit 2
fi

echo "[INFO] resolved inputs:"
echo "  xref             = $IN_XREF"
echo "  zinc_base_url    = $ZINC_BASE_URL"
echo "  zinc_tiers       = $ZINC_TIERS"
echo "  baseline_coverage= $BASELINE_COVERAGE"

echo "[STEP] smoke build"
python3 "$PIPELINE_DIR/scripts/01_build_molecule_zinc_xref_v1.py" \
  --xref "$IN_XREF" \
  --out "$SMOKE_TABLE" \
  --report "$SMOKE_BUILD" \
  --coverage-report "$SMOKE_COVERAGE" \
  --conflict-report "$SMOKE_CONFLICT" \
  --manual-report "$SMOKE_MANUAL" \
  --max-rows "$SMOKE_MAX_ROWS" \
  --max-vendors "$SMOKE_MAX_VENDORS" \
  --zinc-base-url "$ZINC_BASE_URL" \
  --tiers $ZINC_TIERS \
  --http-timeout "$HTTP_TIMEOUT" \
  --http-retries "$HTTP_RETRIES" \
  --baseline-coverage-report "$BASELINE_COVERAGE"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$SMOKE_TABLE" \
  --out "$SMOKE_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_zinc_xref_v1.py" \
  --xref "$IN_XREF" \
  --table "$SMOKE_TABLE" \
  --build-report "$SMOKE_BUILD" \
  --coverage-report "$SMOKE_COVERAGE" \
  --conflict-report "$SMOKE_CONFLICT" \
  --out "$SMOKE_QA"

echo "[STEP] full build"
FULL_MAX_VENDORS_ARG=()
if [[ -n "$FULL_MAX_VENDORS" ]]; then
  FULL_MAX_VENDORS_ARG=(--max-vendors "$FULL_MAX_VENDORS")
fi

python3 "$PIPELINE_DIR/scripts/01_build_molecule_zinc_xref_v1.py" \
  --xref "$IN_XREF" \
  --out "$OUT_TABLE" \
  --report "$REPORT_BUILD" \
  --coverage-report "$REPORT_COVERAGE" \
  --conflict-report "$REPORT_CONFLICT" \
  --manual-report "$REPORT_MANUAL" \
  "${FULL_MAX_VENDORS_ARG[@]}" \
  --zinc-base-url "$ZINC_BASE_URL" \
  --tiers $ZINC_TIERS \
  --http-timeout "$HTTP_TIMEOUT" \
  --http-retries "$HTTP_RETRIES" \
  --baseline-coverage-report "$BASELINE_COVERAGE"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUT_TABLE" \
  --out "$REPORT_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_zinc_xref_v1.py" \
  --xref "$IN_XREF" \
  --table "$OUT_TABLE" \
  --build-report "$REPORT_BUILD" \
  --coverage-report "$REPORT_COVERAGE" \
  --conflict-report "$REPORT_CONFLICT" \
  --out "$REPORT_QA"

python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$OUT_TABLE" \
  "$IN_XREF" \
  "$CONTRACT" \
  "$REPORT_BUILD" \
  "$REPORT_COVERAGE" \
  "$REPORT_CONFLICT" \
  "$REPORT_VALIDATION" \
  "$REPORT_QA"

echo "[DONE] molecule_zinc_xref_v1 finished (smoke + full + validation + QA + manifest)."
