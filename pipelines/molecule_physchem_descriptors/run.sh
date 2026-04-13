#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

PIPELINE_DIR="pipelines/molecule_physchem_descriptors"

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
REQ_CHEMBL_DB="${CHEMBL_DB_PATH:-data/raw/molecules/chembl_36/chembl_36.db}"

IN_XREF="$(resolve_input "$REQ_XREF" "../1218/data/output/molecules/molecule_xref_core_v2.tsv")"
IN_CHEMBL_DB="$(resolve_input "$REQ_CHEMBL_DB" "../12182/data/chembl_36/chembl_36.db")"

OUT_TABLE="data/output/molecules/molecule_physchem_descriptors_v1.tsv"
CONTRACT="$PIPELINE_DIR/contracts/molecule_physchem_descriptors_v1.json"

SMOKE_MAX_ROWS="${SMOKE_MAX_ROWS:-5000}"
DATA_VERSION="${DATA_VERSION:-kg-data-local}"

SMOKE_TABLE="$PIPELINE_DIR/.cache/molecule_physchem_descriptors_v1.smoke.tsv"
SMOKE_BUILD="$PIPELINE_DIR/reports/molecule_physchem_descriptors_v1.smoke.build.json"
SMOKE_VALIDATION="$PIPELINE_DIR/reports/molecule_physchem_descriptors_v1.smoke.validation.json"
SMOKE_QA="$PIPELINE_DIR/reports/molecule_physchem_descriptors_v1.smoke.qa.json"

REPORT_BUILD="$PIPELINE_DIR/reports/molecule_physchem_descriptors_v1.build.json"
REPORT_VALIDATION="$PIPELINE_DIR/reports/molecule_physchem_descriptors_v1.validation.json"
REPORT_QA="$PIPELINE_DIR/reports/molecule_physchem_descriptors_v1.qa.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/molecule_physchem_descriptors_v1.manifest.json"
REPORT_MANUAL="$PIPELINE_DIR/reports/molecule_physchem_descriptors_v1.manual_download.json"

mkdir -p "$(dirname "$OUT_TABLE")" "$PIPELINE_DIR/reports" "$PIPELINE_DIR/.cache"

missing=()
[[ -n "$IN_XREF" ]] || missing+=("$REQ_XREF")
[[ -n "$IN_CHEMBL_DB" ]] || missing+=("$REQ_CHEMBL_DB")

if (( ${#missing[@]} > 0 )); then
  echo "[ERROR] missing required inputs for molecule_physchem_descriptors_v1:" >&2
  for f in "${missing[@]}"; do
    echo "  - $f" >&2
  done

  python3 - <<'PY'
import json
from pathlib import Path

report = {
    "name": "molecule_physchem_descriptors_v1.manual_download",
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
path = Path("pipelines/molecule_physchem_descriptors/reports/molecule_physchem_descriptors_v1.manual_download.json")
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(f"[INFO] manual download checklist -> {path}")
PY
  exit 2
fi

echo "[INFO] resolved inputs:"
echo "  xref      = $IN_XREF"
echo "  chembl_db = $IN_CHEMBL_DB"

echo "[STEP] smoke build (max_rows=$SMOKE_MAX_ROWS)"
python3 "$PIPELINE_DIR/scripts/build_molecule_physchem_descriptors_v1.py" \
  --xref "$IN_XREF" \
  --chembl-db "$IN_CHEMBL_DB" \
  --out "$SMOKE_TABLE" \
  --report "$SMOKE_BUILD" \
  --disable-pubchem-inchikey-fallback \
  --max-rows "$SMOKE_MAX_ROWS"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$SMOKE_TABLE" \
  --out "$SMOKE_VALIDATION"

python3 "$PIPELINE_DIR/scripts/qa_molecule_physchem_descriptors_v1.py" \
  --xref "$IN_XREF" \
  --table "$SMOKE_TABLE" \
  --build-report "$SMOKE_BUILD" \
  --out "$SMOKE_QA" \
  --min-backlink-rate 1.0 \
  --min-format-valid-rate 0.99 \
  --min-mappable-non-empty-rate 0.50 \
  --min-core-coverage 0.50

echo "[STEP] full build"
python3 "$PIPELINE_DIR/scripts/build_molecule_physchem_descriptors_v1.py" \
  --xref "$IN_XREF" \
  --chembl-db "$IN_CHEMBL_DB" \
  --out "$OUT_TABLE" \
  --report "$REPORT_BUILD" \
  --disable-pubchem-inchikey-fallback

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUT_TABLE" \
  --out "$REPORT_VALIDATION"

python3 "$PIPELINE_DIR/scripts/qa_molecule_physchem_descriptors_v1.py" \
  --xref "$IN_XREF" \
  --table "$OUT_TABLE" \
  --build-report "$REPORT_BUILD" \
  --out "$REPORT_QA" \
  --min-backlink-rate 1.0 \
  --min-format-valid-rate 0.99 \
  --min-mappable-non-empty-rate 0.90 \
  --min-core-coverage 0.70

python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$OUT_TABLE" \
  "$IN_XREF" \
  "$IN_CHEMBL_DB" \
  "$CONTRACT" \
  "$REPORT_BUILD" \
  "$REPORT_VALIDATION" \
  "$REPORT_QA"

echo "[DONE] molecule_physchem_descriptors_v1 finished (smoke + full + validation + QA + manifest)."
