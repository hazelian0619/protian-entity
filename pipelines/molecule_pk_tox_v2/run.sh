#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

PIPELINE_DIR="pipelines/molecule_pk_tox_v2"

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

REQ_XREF="${XREF_V2_PATH:-data/output/molecules/molecule_xref_core_v2.tsv}"
REQ_DRUGBANK_XML="${DRUGBANK_XML_PATH:-../1218/data/raw/drugbank/drugbank_complete_database_2026-01-04.xml}"

IN_XREF="$(resolve_input "$REQ_XREF" "../1218/data/output/molecules/molecule_xref_core_v2.tsv")"
IN_DRUGBANK_XML="$(resolve_input "$REQ_DRUGBANK_XML" "../1218/data/raw/drugbank/drugbank_complete_database_2026-01-04.xml")"

CHEMBL_CACHE="${CHEMBL_CACHE_PATH:-$PIPELINE_DIR/.cache/chembl_pk_tox_v2.cache.tsv}"

OUT_TABLE="data/output/molecules/molecule_pk_tox_v2.tsv"
CONTRACT="$PIPELINE_DIR/contracts/molecule_pk_tox_v2.json"

SMOKE_MAX_ROWS="${SMOKE_MAX_ROWS:-3000}"
DATA_VERSION="${DATA_VERSION:-kg-data-local}"
MIN_NORMALIZED_LEGAL_RATE="${MIN_NORMALIZED_LEGAL_RATE:-0.98}"
BASELINE_ROWS_WITH_ANY="${BASELINE_ROWS_WITH_ANY:-2370}"

SMOKE_TABLE="$PIPELINE_DIR/.cache/molecule_pk_tox_v2.smoke.tsv"
SMOKE_BUILD="$PIPELINE_DIR/reports/molecule_pk_tox_v2.smoke.build.json"
SMOKE_VALIDATION="$PIPELINE_DIR/reports/molecule_pk_tox_v2.smoke.validation.json"
SMOKE_QA="$PIPELINE_DIR/reports/molecule_pk_tox_v2.smoke.qa.json"

REPORT_BUILD="$PIPELINE_DIR/reports/molecule_pk_tox_v2.build.json"
REPORT_VALIDATION="$PIPELINE_DIR/reports/molecule_pk_tox_v2.validation.json"
REPORT_QA="$PIPELINE_DIR/reports/molecule_pk_tox_v2.qa.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/molecule_pk_tox_v2.manifest.json"
REPORT_MANUAL="$PIPELINE_DIR/reports/molecule_pk_tox_v2.manual_download.json"

mkdir -p "$(dirname "$OUT_TABLE")" "$PIPELINE_DIR/reports" "$PIPELINE_DIR/.cache"

missing=()
[[ -n "$IN_XREF" ]] || missing+=("$REQ_XREF")
[[ -n "$IN_DRUGBANK_XML" ]] || missing+=("$REQ_DRUGBANK_XML")

if (( ${#missing[@]} > 0 )); then
  echo "[ERROR] missing required inputs for molecule_pk_tox_v2:" >&2
  printf '  - %s\n' "${missing[@]}" >&2

  python3 - <<'PY'
import json
from pathlib import Path

report = {
  "name": "molecule_pk_tox_v2.manual_download",
  "reason": "required inputs missing",
  "download_checklist": [
    {
      "name": "molecule_xref_core_v2.tsv",
      "expected_path": "data/output/molecules/molecule_xref_core_v2.tsv",
      "sha256_cmd": "sha256sum data/output/molecules/molecule_xref_core_v2.tsv"
    },
    {
      "name": "drugbank_complete_database_2026-01-04.xml",
      "expected_path": "data/raw/drugbank/drugbank_complete_database_2026-01-04.xml",
      "sha256_cmd": "sha256sum data/raw/drugbank/drugbank_complete_database_2026-01-04.xml"
    }
  ]
}
path = Path("pipelines/molecule_pk_tox_v2/reports/molecule_pk_tox_v2.manual_download.json")
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(f"[INFO] manual checklist -> {path}")
PY
  exit 2
fi

echo "[INFO] resolved inputs:"
echo "  xref         = $IN_XREF"
echo "  drugbank_xml = $IN_DRUGBANK_XML"
echo "  chembl_cache = $CHEMBL_CACHE"

echo "[STEP] smoke build (max_rows=$SMOKE_MAX_ROWS)"
python3 "$PIPELINE_DIR/scripts/01_build_molecule_pk_tox_v2.py" \
  --xref "$IN_XREF" \
  --drugbank-xml "$IN_DRUGBANK_XML" \
  --chembl-cache "$CHEMBL_CACHE" \
  --out "$SMOKE_TABLE" \
  --report "$SMOKE_BUILD" \
  --max-rows "$SMOKE_MAX_ROWS"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$SMOKE_TABLE" \
  --out "$SMOKE_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_pk_tox_v2.py" \
  --xref "$IN_XREF" \
  --table "$SMOKE_TABLE" \
  --build-report "$SMOKE_BUILD" \
  --validation-report "$SMOKE_VALIDATION" \
  --out "$SMOKE_QA" \
  --baseline-rows-with-any 0 \
  --min-normalized-legal-rate "$MIN_NORMALIZED_LEGAL_RATE"

echo "[STEP] full build"
python3 "$PIPELINE_DIR/scripts/01_build_molecule_pk_tox_v2.py" \
  --xref "$IN_XREF" \
  --drugbank-xml "$IN_DRUGBANK_XML" \
  --chembl-cache "$CHEMBL_CACHE" \
  --out "$OUT_TABLE" \
  --report "$REPORT_BUILD"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUT_TABLE" \
  --out "$REPORT_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_pk_tox_v2.py" \
  --xref "$IN_XREF" \
  --table "$OUT_TABLE" \
  --build-report "$REPORT_BUILD" \
  --validation-report "$REPORT_VALIDATION" \
  --out "$REPORT_QA" \
  --baseline-rows-with-any "$BASELINE_ROWS_WITH_ANY" \
  --min-normalized-legal-rate "$MIN_NORMALIZED_LEGAL_RATE"

python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$OUT_TABLE" \
  "$IN_XREF" \
  "$IN_DRUGBANK_XML" \
  "$CONTRACT" \
  "$REPORT_BUILD" \
  "$REPORT_VALIDATION" \
  "$REPORT_QA"

echo "[DONE] molecule_pk_tox_v2 finished (smoke + full + validation + QA + manifest)."
