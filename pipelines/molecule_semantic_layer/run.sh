#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

PIPELINE_DIR="pipelines/molecule_semantic_layer"

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
REQ_CHEBI_OBO="${CHEBI_OBO_PATH:-data/raw/molecules/chebi/chebi.obo}"
REQ_DRUG_MASTER="${DRUG_MASTER_PATH:-data/output/drugbank/drug_master_v1.tsv}"
REQ_CHEMBL_DB="${CHEMBL_DB_PATH:-data/raw/molecules/chembl_36/chembl_36.db}"
OPT_REGISTRY_3D="${REGISTRY_3D_PATH:-data/output/molecules/molecule_3d_registry_v1.tsv}"

IN_XREF="$(resolve_input "$REQ_XREF_V2" "$REQ_XREF_V1" "../1218/data/output/molecules/molecule_xref_core_v2.tsv" "../1218/data/output/molecules/molecule_xref_core_v1.tsv")"
IN_CHEBI_OBO="$(resolve_input "$REQ_CHEBI_OBO" "../1218/data/raw/molecules/chebi/chebi.obo")"
IN_DRUG_MASTER="$(resolve_input "$REQ_DRUG_MASTER" "../1218/data/output/drugbank/drug_master_v1.tsv")"
IN_CHEMBL_DB="$(resolve_input "$REQ_CHEMBL_DB" "../12182/data/chembl_36/chembl_36.db")"
IN_REGISTRY_3D="$(resolve_input "$OPT_REGISTRY_3D" "../protian-entity/data/output/molecules/molecule_3d_registry_v1.tsv")"

OUT_TABLE="data/output/molecules/molecule_semantic_tags_v1.tsv"
CONTRACT="$PIPELINE_DIR/contracts/molecule_semantic_tags_v1.json"

SMOKE_MAX_ROWS="${SMOKE_MAX_ROWS:-3000}"
DATA_VERSION="${DATA_VERSION:-kg-data-local}"

SMOKE_TABLE="$PIPELINE_DIR/.cache/molecule_semantic_tags_v1.smoke.tsv"
SMOKE_BUILD="$PIPELINE_DIR/reports/molecule_semantic_tags_v1.smoke.build.json"
SMOKE_HIERARCHY="$PIPELINE_DIR/reports/molecule_semantic_tags_v1.smoke.hierarchy.json"
SMOKE_COVERAGE="$PIPELINE_DIR/reports/molecule_semantic_tags_v1.smoke.coverage.json"
SMOKE_VALIDATION="$PIPELINE_DIR/reports/molecule_semantic_tags_v1.smoke.validation.json"
SMOKE_QA="$PIPELINE_DIR/reports/molecule_semantic_tags_v1.smoke.qa.json"

REPORT_BUILD="$PIPELINE_DIR/reports/molecule_semantic_tags_v1.build.json"
REPORT_HIERARCHY="$PIPELINE_DIR/reports/molecule_semantic_tags_v1.hierarchy.json"
REPORT_COVERAGE="$PIPELINE_DIR/reports/molecule_semantic_tags_v1.coverage.json"
REPORT_VALIDATION="$PIPELINE_DIR/reports/molecule_semantic_tags_v1.validation.json"
REPORT_QA="$PIPELINE_DIR/reports/molecule_semantic_tags_v1.qa.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/molecule_semantic_tags_v1.manifest.json"
REPORT_MANUAL="$PIPELINE_DIR/reports/molecule_semantic_tags_v1.manual_download.json"

mkdir -p "$(dirname "$OUT_TABLE")" "$PIPELINE_DIR/reports" "$PIPELINE_DIR/.cache"

missing=()
[[ -n "$IN_XREF" ]] || missing+=("$REQ_XREF_V2 or $REQ_XREF_V1")
[[ -n "$IN_CHEBI_OBO" ]] || missing+=("$REQ_CHEBI_OBO")
[[ -n "$IN_DRUG_MASTER" ]] || missing+=("$REQ_DRUG_MASTER")
[[ -n "$IN_CHEMBL_DB" ]] || missing+=("$REQ_CHEMBL_DB")

if (( ${#missing[@]} > 0 )); then
  echo "[ERROR] missing required inputs for molecule_semantic_layer:" >&2
  for f in "${missing[@]}"; do
    echo "  - $f" >&2
  done

  python3 - <<'PY'
import json
from pathlib import Path

report = {
    "name": "molecule_semantic_tags_v1.manual_download",
    "reason": "required inputs missing; cannot build semantic layer (ChEBI + ATC + DrugBank)",
    "download_checklist": [
      {
        "name": "molecule_xref_core_v2.tsv (preferred) or v1",
        "expected_path": "data/output/molecules/molecule_xref_core_v2.tsv",
        "sha256_cmd": "sha256sum data/output/molecules/molecule_xref_core_v2.tsv"
      },
      {
        "name": "ChEBI ontology (OBO)",
        "expected_path": "data/raw/molecules/chebi/chebi.obo",
        "sha256_cmd": "sha256sum data/raw/molecules/chebi/chebi.obo"
      },
      {
        "name": "drug_master_v1.tsv",
        "expected_path": "data/output/drugbank/drug_master_v1.tsv",
        "sha256_cmd": "sha256sum data/output/drugbank/drug_master_v1.tsv"
      },
      {
        "name": "chembl_36.db",
        "expected_path": "data/raw/molecules/chembl_36/chembl_36.db",
        "sha256_cmd": "sha256sum data/raw/molecules/chembl_36/chembl_36.db"
      }
    ]
}
path = Path("pipelines/molecule_semantic_layer/reports/molecule_semantic_tags_v1.manual_download.json")
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(f"[INFO] manual checklist -> {path}")
PY
  exit 2
fi

echo "[INFO] resolved inputs:"
echo "  xref        = $IN_XREF"
echo "  chebi_obo   = $IN_CHEBI_OBO"
echo "  drug_master = $IN_DRUG_MASTER"
echo "  chembl_db   = $IN_CHEMBL_DB"
if [[ -n "$IN_REGISTRY_3D" ]]; then
  echo "  registry_3d = $IN_REGISTRY_3D"
else
  echo "  registry_3d = <missing, optional>"
fi

REGISTRY_ARG=()
if [[ -n "$IN_REGISTRY_3D" ]]; then
  REGISTRY_ARG=(--registry-3d "$IN_REGISTRY_3D")
fi

echo "[STEP] smoke build (max_rows=$SMOKE_MAX_ROWS)"
python3 "$PIPELINE_DIR/scripts/01_build_molecule_semantic_tags_v1.py" \
  --xref "$IN_XREF" \
  --chebi-obo "$IN_CHEBI_OBO" \
  --drug-master "$IN_DRUG_MASTER" \
  --chembl-db "$IN_CHEMBL_DB" \
  "${REGISTRY_ARG[@]}" \
  --out "$SMOKE_TABLE" \
  --report "$SMOKE_BUILD" \
  --hierarchy-report "$SMOKE_HIERARCHY" \
  --coverage-report "$SMOKE_COVERAGE" \
  --max-rows "$SMOKE_MAX_ROWS"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$SMOKE_TABLE" \
  --out "$SMOKE_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_semantic_tags_v1.py" \
  --xref "$IN_XREF" \
  --table "$SMOKE_TABLE" \
  --build-report "$SMOKE_BUILD" \
  --hierarchy-report "$SMOKE_HIERARCHY" \
  --coverage-report "$SMOKE_COVERAGE" \
  --out "$SMOKE_QA"

echo "[STEP] full build"
python3 "$PIPELINE_DIR/scripts/01_build_molecule_semantic_tags_v1.py" \
  --xref "$IN_XREF" \
  --chebi-obo "$IN_CHEBI_OBO" \
  --drug-master "$IN_DRUG_MASTER" \
  --chembl-db "$IN_CHEMBL_DB" \
  "${REGISTRY_ARG[@]}" \
  --out "$OUT_TABLE" \
  --report "$REPORT_BUILD" \
  --hierarchy-report "$REPORT_HIERARCHY" \
  --coverage-report "$REPORT_COVERAGE"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUT_TABLE" \
  --out "$REPORT_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_semantic_tags_v1.py" \
  --xref "$IN_XREF" \
  --table "$OUT_TABLE" \
  --build-report "$REPORT_BUILD" \
  --hierarchy-report "$REPORT_HIERARCHY" \
  --coverage-report "$REPORT_COVERAGE" \
  --out "$REPORT_QA"

python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$OUT_TABLE" \
  "$IN_XREF" \
  "$IN_CHEBI_OBO" \
  "$IN_DRUG_MASTER" \
  "$IN_CHEMBL_DB" \
  "$CONTRACT" \
  "$REPORT_BUILD" \
  "$REPORT_HIERARCHY" \
  "$REPORT_COVERAGE" \
  "$REPORT_VALIDATION" \
  "$REPORT_QA"

echo "[DONE] molecule_semantic_layer finished (smoke + full + validation + QA + manifest)."
