#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

PIPELINE_DIR="pipelines/molecule_xref_enrichment_v2"

resolve_input() {
  local primary="$1"
  shift
  if [[ -f "$primary" ]]; then
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

# Canonical inputs (can be overridden by env)
REQ_V1_CORE="${V1_CORE_PATH:-data/output/molecules/molecule_xref_core_v1.tsv}"
REQ_DRUG_MASTER="${DRUG_MASTER_PATH:-data/output/drugbank/drug_master_v1.tsv}"
REQ_DRUG_XREF="${DRUG_XREF_PATH:-data/output/drugbank/drug_xref_molecules_v1.tsv}"
REQ_M1_DB="${M1_DB_PATH:-data/output/molecules/molecules_m1.sqlite}"
REQ_CHEMBL_DB="${CHEMBL_DB_PATH:-data/raw/molecules/chembl_36/chembl_36.db}"

# local fallback candidates (for shared workspace layout)
IN_V1_CORE="$(resolve_input "$REQ_V1_CORE" "../1218/data/output/molecules/molecule_xref_core_v1.tsv")"
IN_DRUG_MASTER="$(resolve_input "$REQ_DRUG_MASTER" "../1218/data/output/drugbank/drug_master_v1.tsv")"
IN_DRUG_XREF="$(resolve_input "$REQ_DRUG_XREF" "../1218/data/output/drugbank/drug_xref_molecules_v1.tsv")"
IN_M1_DB="$(resolve_input "$REQ_M1_DB" "../12182/out/m1/molecules_m1.sqlite")"
IN_CHEMBL_DB="$(resolve_input "$REQ_CHEMBL_DB" "../12182/data/chembl_36/chembl_36.db")"

OUT_TABLE="data/output/molecules/molecule_xref_core_v2.tsv"
CONTRACT="$PIPELINE_DIR/contracts/molecule_xref_core_v2.json"

SMOKE_MAX_ROWS="${SMOKE_MAX_ROWS:-5000}"
DATA_VERSION="${DATA_VERSION:-kg-data-local}"

SMOKE_TABLE="$PIPELINE_DIR/.cache/molecule_xref_core_v2.smoke.tsv"
SMOKE_BUILD="$PIPELINE_DIR/reports/molecule_xref_core_v2.smoke.build.json"
SMOKE_BACKFILL="$PIPELINE_DIR/reports/molecule_xref_core_v2.smoke.backfill_audit.json"
SMOKE_MISSING="$PIPELINE_DIR/reports/molecule_xref_core_v2.smoke.missing_audit.json"
SMOKE_VALIDATION="$PIPELINE_DIR/reports/molecule_xref_core_v2.smoke.validation.json"
SMOKE_QA="$PIPELINE_DIR/reports/molecule_xref_core_v2.smoke.qa.json"

REPORT_BUILD="$PIPELINE_DIR/reports/molecule_xref_core_v2.build.json"
REPORT_BACKFILL="$PIPELINE_DIR/reports/molecule_xref_core_v2.backfill_audit.json"
REPORT_MISSING="$PIPELINE_DIR/reports/molecule_xref_core_v2.missing_audit.json"
REPORT_VALIDATION="$PIPELINE_DIR/reports/molecule_xref_core_v2.validation.json"
REPORT_QA="$PIPELINE_DIR/reports/molecule_xref_core_v2.qa.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/molecule_xref_core_v2.manifest.json"
REPORT_MANUAL="$PIPELINE_DIR/reports/molecule_xref_core_v2.manual_download.json"

mkdir -p "$(dirname "$OUT_TABLE")" "$PIPELINE_DIR/reports" "$PIPELINE_DIR/.cache"

missing=()
[[ -n "$IN_V1_CORE" ]] || missing+=("$REQ_V1_CORE")
[[ -n "$IN_DRUG_MASTER" ]] || missing+=("$REQ_DRUG_MASTER")
[[ -n "$IN_DRUG_XREF" ]] || missing+=("$REQ_DRUG_XREF")
[[ -n "$IN_M1_DB" ]] || missing+=("$REQ_M1_DB")
[[ -n "$IN_CHEMBL_DB" ]] || missing+=("$REQ_CHEMBL_DB")

if (( ${#missing[@]} > 0 )); then
  echo "[ERROR] missing required inputs for molecule_xref_enrichment_v2:" >&2
  for f in "${missing[@]}"; do
    echo "  - $f" >&2
  done

  python3 - <<'PY'
import json
from pathlib import Path

report = {
    "name": "molecule_xref_core_v2.manual_download",
    "reason": "required inputs missing; cannot complete full backfill + enrichment gates",
    "download_checklist": [
        {
            "name": "molecule_xref_core_v1.tsv",
            "expected_path": "data/output/molecules/molecule_xref_core_v1.tsv",
            "sha256_cmd": "sha256sum data/output/molecules/molecule_xref_core_v1.tsv"
        },
        {
            "name": "drug_master_v1.tsv",
            "expected_path": "data/output/drugbank/drug_master_v1.tsv",
            "sha256_cmd": "sha256sum data/output/drugbank/drug_master_v1.tsv"
        },
        {
            "name": "drug_xref_molecules_v1.tsv",
            "expected_path": "data/output/drugbank/drug_xref_molecules_v1.tsv",
            "sha256_cmd": "sha256sum data/output/drugbank/drug_xref_molecules_v1.tsv"
        },
        {
            "name": "molecules_m1.sqlite",
            "expected_path": "data/output/molecules/molecules_m1.sqlite",
            "sha256_cmd": "sha256sum data/output/molecules/molecules_m1.sqlite"
        },
        {
            "name": "chembl_36.db",
            "expected_path": "data/raw/molecules/chembl_36/chembl_36.db",
            "sha256_cmd": "sha256sum data/raw/molecules/chembl_36/chembl_36.db"
        }
    ]
}
path = Path("pipelines/molecule_xref_enrichment_v2/reports/molecule_xref_core_v2.manual_download.json")
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(f"[INFO] manual download checklist -> {path}")
PY
  exit 2
fi

echo "[INFO] resolved inputs:"
echo "  v1_core      = $IN_V1_CORE"
echo "  drug_master  = $IN_DRUG_MASTER"
echo "  drug_xref    = $IN_DRUG_XREF"
echo "  m1_db        = $IN_M1_DB"
echo "  chembl_db    = $IN_CHEMBL_DB"

echo "[STEP] smoke build (max_rows=$SMOKE_MAX_ROWS)"
python3 "$PIPELINE_DIR/scripts/01_build_molecule_xref_core_v2.py" \
  --v1-core "$IN_V1_CORE" \
  --drug-master "$IN_DRUG_MASTER" \
  --drug-xref "$IN_DRUG_XREF" \
  --m1-db "$IN_M1_DB" \
  --chembl-db "$IN_CHEMBL_DB" \
  --out "$SMOKE_TABLE" \
  --report "$SMOKE_BUILD" \
  --backfill-audit "$SMOKE_BACKFILL" \
  --missing-audit "$SMOKE_MISSING" \
  --max-rows "$SMOKE_MAX_ROWS"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$SMOKE_TABLE" \
  --out "$SMOKE_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_xref_core_v2.py" \
  --v1-table "$IN_V1_CORE" \
  --v2-table "$SMOKE_TABLE" \
  --backfill-audit "$SMOKE_BACKFILL" \
  --build-report "$SMOKE_BUILD" \
  --out "$SMOKE_QA" \
  --require-positive-delta 0

echo "[STEP] full build"
python3 "$PIPELINE_DIR/scripts/01_build_molecule_xref_core_v2.py" \
  --v1-core "$IN_V1_CORE" \
  --drug-master "$IN_DRUG_MASTER" \
  --drug-xref "$IN_DRUG_XREF" \
  --m1-db "$IN_M1_DB" \
  --chembl-db "$IN_CHEMBL_DB" \
  --out "$OUT_TABLE" \
  --report "$REPORT_BUILD" \
  --backfill-audit "$REPORT_BACKFILL" \
  --missing-audit "$REPORT_MISSING"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUT_TABLE" \
  --out "$REPORT_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_xref_core_v2.py" \
  --v1-table "$IN_V1_CORE" \
  --v2-table "$OUT_TABLE" \
  --backfill-audit "$REPORT_BACKFILL" \
  --build-report "$REPORT_BUILD" \
  --out "$REPORT_QA" \
  --require-positive-delta 1

python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$OUT_TABLE" \
  "$IN_V1_CORE" \
  "$IN_DRUG_MASTER" \
  "$IN_DRUG_XREF" \
  "$IN_M1_DB" \
  "$IN_CHEMBL_DB" \
  "$CONTRACT" \
  "$REPORT_BUILD" \
  "$REPORT_BACKFILL" \
  "$REPORT_MISSING" \
  "$REPORT_VALIDATION" \
  "$REPORT_QA"

echo "[DONE] molecule_xref_enrichment_v2 finished (smoke + full + validation + QA + manifest)."
