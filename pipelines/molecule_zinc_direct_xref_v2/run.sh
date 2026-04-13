#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

PIPELINE_DIR="pipelines/molecule_zinc_direct_xref_v2"

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

REQ_ZINC_V1="${ZINC_V1_PATH:-data/output/molecules/molecule_zinc_xref_v1.tsv}"
REQ_STRUCT_V1="${STRUCT_IDENTIFIERS_PATH:-data/output/molecules/molecule_structure_identifiers_v1.tsv}"

IN_ZINC_V1="$(resolve_input "$REQ_ZINC_V1" "../1218/data/output/molecules/molecule_zinc_xref_v1.tsv" "../protian-entity/data/output/molecules/molecule_zinc_xref_v1.tsv")"
IN_STRUCT_V1="$(resolve_input "$REQ_STRUCT_V1" "../1218/data/output/molecules/molecule_structure_identifiers_v1.tsv" "../protian-entity/data/output/molecules/molecule_structure_identifiers_v1.tsv")"

OUT_TABLE="data/output/molecules/molecule_zinc_xref_v2.tsv"
CONTRACT="$PIPELINE_DIR/contracts/molecule_zinc_xref_v2.json"

SMOKE_MAX_ROWS="${SMOKE_MAX_ROWS:-5000}"
MIN_SMILES_DELTA="${MIN_SMILES_DELTA:-0.15}"
DATA_VERSION="${DATA_VERSION:-kg-data-local}"

SMOKE_TABLE="$PIPELINE_DIR/.cache/molecule_zinc_xref_v2.smoke.tsv"
SMOKE_BUILD="$PIPELINE_DIR/reports/molecule_zinc_xref_v2.smoke.build.json"
SMOKE_COVERAGE="$PIPELINE_DIR/reports/molecule_zinc_xref_v2.smoke.coverage.json"
SMOKE_CONFLICT="$PIPELINE_DIR/reports/molecule_zinc_xref_v2.smoke.conflict_audit.json"
SMOKE_VALIDATION="$PIPELINE_DIR/reports/molecule_zinc_xref_v2.smoke.validation.json"
SMOKE_QA="$PIPELINE_DIR/reports/molecule_zinc_xref_v2.smoke.qa.json"

REPORT_BUILD="$PIPELINE_DIR/reports/molecule_zinc_xref_v2.build.json"
REPORT_COVERAGE="$PIPELINE_DIR/reports/molecule_zinc_xref_v2.coverage.json"
REPORT_CONFLICT="$PIPELINE_DIR/reports/molecule_zinc_xref_v2.conflict_audit.json"
REPORT_VALIDATION="$PIPELINE_DIR/reports/molecule_zinc_xref_v2.validation.json"
REPORT_QA="$PIPELINE_DIR/reports/molecule_zinc_xref_v2.qa.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/molecule_zinc_xref_v2.manifest.json"

mkdir -p "$(dirname "$OUT_TABLE")" "$PIPELINE_DIR/reports" "$PIPELINE_DIR/.cache"

missing=()
[[ -n "$IN_ZINC_V1" ]] || missing+=("$REQ_ZINC_V1")
[[ -n "$IN_STRUCT_V1" ]] || missing+=("$REQ_STRUCT_V1")
if (( ${#missing[@]} > 0 )); then
  echo "[ERROR] missing required inputs for molecule_zinc_xref_v2:" >&2
  printf '  - %s\n' "${missing[@]}" >&2
  exit 2
fi

echo "[INFO] resolved inputs:"
echo "  zinc_v1   = $IN_ZINC_V1"
echo "  structure = $IN_STRUCT_V1"
echo "  min_smiles_delta = $MIN_SMILES_DELTA"

echo "[STEP] smoke build (max_rows=$SMOKE_MAX_ROWS)"
python3 "$PIPELINE_DIR/scripts/01_build_molecule_zinc_xref_v2.py" \
  --zinc-v1 "$IN_ZINC_V1" \
  --structure-identifiers "$IN_STRUCT_V1" \
  --out "$SMOKE_TABLE" \
  --report "$SMOKE_BUILD" \
  --coverage-report "$SMOKE_COVERAGE" \
  --conflict-report "$SMOKE_CONFLICT" \
  --max-rows "$SMOKE_MAX_ROWS"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$SMOKE_TABLE" \
  --out "$SMOKE_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_zinc_xref_v2.py" \
  --zinc-v1 "$IN_ZINC_V1" \
  --table "$SMOKE_TABLE" \
  --build-report "$SMOKE_BUILD" \
  --coverage-report "$SMOKE_COVERAGE" \
  --conflict-report "$SMOKE_CONFLICT" \
  --out "$SMOKE_QA" \
  --min-smiles-delta "$MIN_SMILES_DELTA" \
  --allow-row-subset 1

echo "[STEP] full build"
python3 "$PIPELINE_DIR/scripts/01_build_molecule_zinc_xref_v2.py" \
  --zinc-v1 "$IN_ZINC_V1" \
  --structure-identifiers "$IN_STRUCT_V1" \
  --out "$OUT_TABLE" \
  --report "$REPORT_BUILD" \
  --coverage-report "$REPORT_COVERAGE" \
  --conflict-report "$REPORT_CONFLICT"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUT_TABLE" \
  --out "$REPORT_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_zinc_xref_v2.py" \
  --zinc-v1 "$IN_ZINC_V1" \
  --table "$OUT_TABLE" \
  --build-report "$REPORT_BUILD" \
  --coverage-report "$REPORT_COVERAGE" \
  --conflict-report "$REPORT_CONFLICT" \
  --out "$REPORT_QA" \
  --min-smiles-delta "$MIN_SMILES_DELTA"

python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$OUT_TABLE" \
  "$IN_ZINC_V1" \
  "$IN_STRUCT_V1" \
  "$CONTRACT" \
  "$REPORT_BUILD" \
  "$REPORT_COVERAGE" \
  "$REPORT_CONFLICT" \
  "$REPORT_VALIDATION" \
  "$REPORT_QA"

echo "[DONE] molecule_zinc_xref_v2 finished (smoke + full + validation + QA + manifest)."
