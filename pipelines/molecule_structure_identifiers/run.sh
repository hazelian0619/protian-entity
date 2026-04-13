#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

PIPELINE_DIR="pipelines/molecule_structure_identifiers"

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

REQ_CORE_TSV="${CORE_TSV_PATH:-data/output/molecules/molecule_xref_core_v2.tsv}"
REQ_CHEMBL_DB="${CHEMBL_DB_PATH:-data/raw/molecules/chembl_36/chembl_36.db}"

IN_CORE_TSV="$(resolve_input "$REQ_CORE_TSV" "../1218/data/output/molecules/molecule_xref_core_v2.tsv")"
IN_CHEMBL_DB="$(resolve_input "$REQ_CHEMBL_DB" "../1218/data/raw/molecules/chembl_36/chembl_36.db" "../12182/data/chembl_36/chembl_36.db")"

OUT_TABLE="data/output/molecules/molecule_structure_identifiers_v1.tsv"
CONTRACT="$PIPELINE_DIR/contracts/molecule_structure_identifiers_v1.json"

SMOKE_MAX_ROWS="${SMOKE_MAX_ROWS:-2500}"
DATA_VERSION="${DATA_VERSION:-kg-data-local}"

PUBCHEM_WORKERS="${PUBCHEM_WORKERS:-8}"
PUBCHEM_QPS="${PUBCHEM_QPS:-6.0}"
PUBCHEM_TIMEOUT="${PUBCHEM_TIMEOUT:-20}"
PUBCHEM_RETRIES="${PUBCHEM_RETRIES:-3}"

SMOKE_TABLE="$PIPELINE_DIR/.cache/molecule_structure_identifiers_v1.smoke.tsv"
SMOKE_BUILD="$PIPELINE_DIR/reports/molecule_structure_identifiers_v1.smoke.build.json"
SMOKE_VALIDATION="$PIPELINE_DIR/reports/molecule_structure_identifiers_v1.smoke.validation.json"
SMOKE_QA="$PIPELINE_DIR/reports/molecule_structure_identifiers_v1.smoke.qa.json"

REPORT_BUILD="$PIPELINE_DIR/reports/molecule_structure_identifiers_v1.build.json"
REPORT_VALIDATION="$PIPELINE_DIR/reports/molecule_structure_identifiers_v1.validation.json"
REPORT_QA="$PIPELINE_DIR/reports/molecule_structure_identifiers_v1.qa.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/molecule_structure_identifiers_v1.manifest.json"
REPORT_MANUAL="$PIPELINE_DIR/reports/molecule_structure_identifiers_v1.manual_download.json"

mkdir -p "$(dirname "$OUT_TABLE")" "$PIPELINE_DIR/reports" "$PIPELINE_DIR/.cache"

missing=()
[[ -n "$IN_CORE_TSV" ]] || missing+=("$REQ_CORE_TSV")
[[ -n "$IN_CHEMBL_DB" ]] || missing+=("$REQ_CHEMBL_DB")

if (( ${#missing[@]} > 0 )); then
  echo "[ERROR] missing required inputs for molecule_structure_identifiers:" >&2
  printf '  - %s\n' "${missing[@]}" >&2

  python3 - <<'PY'
import json
from pathlib import Path

report = {
  "name": "molecule_structure_identifiers_v1.manual_download",
  "reason": "required inputs missing",
  "download_checklist": [
    {
      "name": "molecule_xref_core_v2.tsv",
      "expected_path": "data/output/molecules/molecule_xref_core_v2.tsv",
      "sha256_cmd": "sha256sum data/output/molecules/molecule_xref_core_v2.tsv"
    },
    {
      "name": "chembl_36.db",
      "expected_path": "data/raw/molecules/chembl_36/chembl_36.db",
      "sha256_cmd": "sha256sum data/raw/molecules/chembl_36/chembl_36.db"
    }
  ]
}
path = Path("pipelines/molecule_structure_identifiers/reports/molecule_structure_identifiers_v1.manual_download.json")
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(f"[INFO] manual checklist -> {path}")
PY
  exit 2
fi

echo "[INFO] resolved inputs:"
echo "  core_tsv  = $IN_CORE_TSV"
echo "  chembl_db = $IN_CHEMBL_DB"
echo "  workers=$PUBCHEM_WORKERS qps=$PUBCHEM_QPS timeout=$PUBCHEM_TIMEOUT retries=$PUBCHEM_RETRIES"

echo "[STEP] smoke build (max_rows=$SMOKE_MAX_ROWS)"
python3 "$PIPELINE_DIR/scripts/01_build_molecule_structure_identifiers_v1.py" \
  --core-tsv "$IN_CORE_TSV" \
  --chembl-db "$IN_CHEMBL_DB" \
  --out "$SMOKE_TABLE" \
  --build-report "$SMOKE_BUILD" \
  --max-rows "$SMOKE_MAX_ROWS" \
  --workers "$PUBCHEM_WORKERS" \
  --qps "$PUBCHEM_QPS" \
  --timeout "$PUBCHEM_TIMEOUT" \
  --retries "$PUBCHEM_RETRIES"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$SMOKE_TABLE" \
  --out "$SMOKE_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_structure_identifiers_v1.py" \
  --core-tsv "$IN_CORE_TSV" \
  --table "$SMOKE_TABLE" \
  --build-report "$SMOKE_BUILD" \
  --out "$SMOKE_QA" \
  --allow-row-subset 1

echo "[STEP] full build"
python3 "$PIPELINE_DIR/scripts/01_build_molecule_structure_identifiers_v1.py" \
  --core-tsv "$IN_CORE_TSV" \
  --chembl-db "$IN_CHEMBL_DB" \
  --out "$OUT_TABLE" \
  --build-report "$REPORT_BUILD" \
  --workers "$PUBCHEM_WORKERS" \
  --qps "$PUBCHEM_QPS" \
  --timeout "$PUBCHEM_TIMEOUT" \
  --retries "$PUBCHEM_RETRIES"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUT_TABLE" \
  --out "$REPORT_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_structure_identifiers_v1.py" \
  --core-tsv "$IN_CORE_TSV" \
  --table "$OUT_TABLE" \
  --build-report "$REPORT_BUILD" \
  --out "$REPORT_QA"

python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$IN_CORE_TSV" \
  "$IN_CHEMBL_DB" \
  "$OUT_TABLE" \
  "$CONTRACT" \
  "$REPORT_BUILD" \
  "$REPORT_VALIDATION" \
  "$REPORT_QA"

echo "[DONE] molecule_structure_identifiers finished (smoke + full + validation + QA + manifest)."
