#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

PIPELINE_DIR="pipelines/molecule_3d_experimental_linker"

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
REQ_PROTEIN_PDB="${PROTEIN_PDB_PATH:-data/output/protein/pdb_structures_v1.tsv}"

IN_XREF="$(resolve_input "$REQ_XREF_V2" "$REQ_XREF_V1" "../1218/data/output/molecules/molecule_xref_core_v2.tsv" "../1218/data/output/molecules/molecule_xref_core_v1.tsv")"
IN_PROTEIN_PDB="$(resolve_input "$REQ_PROTEIN_PDB" "../1218/data/output/protein/pdb_structures_v1.tsv")"

OUT_TABLE="data/output/molecules/molecule_3d_experimental_v1.tsv"
CONTRACT="$PIPELINE_DIR/contracts/molecule_3d_experimental_v1.json"

SMOKE_MAX_ROWS="${SMOKE_MAX_ROWS:-3000}"
SMOKE_MAX_PDB="${SMOKE_MAX_PDB:-1200}"
FULL_MAX_PDB="${FULL_MAX_PDB:-}"
DATA_VERSION="${DATA_VERSION:-kg-data-local}"

BATCH_SIZE="${RCSB_BATCH_SIZE:-120}"
HTTP_TIMEOUT="${RCSB_HTTP_TIMEOUT:-30}"
HTTP_RETRIES="${RCSB_HTTP_RETRIES:-2}"
HTTP_SLEEP="${RCSB_HTTP_SLEEP:-0.08}"

SMOKE_TABLE="$PIPELINE_DIR/.cache/molecule_3d_experimental_v1.smoke.tsv"
SMOKE_BUILD="$PIPELINE_DIR/reports/molecule_3d_experimental_v1.smoke.build.json"
SMOKE_COVERAGE="$PIPELINE_DIR/reports/molecule_3d_experimental_v1.smoke.coverage.json"
SMOKE_VALIDATION="$PIPELINE_DIR/reports/molecule_3d_experimental_v1.smoke.validation.json"
SMOKE_QA="$PIPELINE_DIR/reports/molecule_3d_experimental_v1.smoke.qa.json"

REPORT_BUILD="$PIPELINE_DIR/reports/molecule_3d_experimental_v1.build.json"
REPORT_COVERAGE="$PIPELINE_DIR/reports/molecule_3d_experimental_v1.coverage.json"
REPORT_VALIDATION="$PIPELINE_DIR/reports/molecule_3d_experimental_v1.validation.json"
REPORT_QA="$PIPELINE_DIR/reports/molecule_3d_experimental_v1.qa.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/molecule_3d_experimental_v1.manifest.json"
REPORT_MANUAL="$PIPELINE_DIR/reports/molecule_3d_experimental_v1.manual_download.json"

mkdir -p "$(dirname "$OUT_TABLE")" "$PIPELINE_DIR/reports" "$PIPELINE_DIR/.cache"

missing=()
[[ -n "$IN_XREF" ]] || missing+=("$REQ_XREF_V2 or $REQ_XREF_V1")
[[ -n "$IN_PROTEIN_PDB" ]] || missing+=("$REQ_PROTEIN_PDB")

if (( ${#missing[@]} > 0 )); then
  echo "[ERROR] missing required inputs for molecule_3d_experimental_v1:" >&2
  for f in "${missing[@]}"; do
    echo "  - $f" >&2
  done

  python3 - <<'PY'
import json
from pathlib import Path

report = {
    "name": "molecule_3d_experimental_v1.manual_download",
    "reason": "required inputs missing; cannot build experimental 3D mapping",
    "download_checklist": [
      {
        "name": "molecule_xref_core_v2.tsv (preferred) or v1",
        "expected_path": "data/output/molecules/molecule_xref_core_v2.tsv",
        "sha256_cmd": "sha256sum data/output/molecules/molecule_xref_core_v2.tsv"
      },
      {
        "name": "pdb_structures_v1.tsv",
        "expected_path": "data/output/protein/pdb_structures_v1.tsv",
        "sha256_cmd": "sha256sum data/output/protein/pdb_structures_v1.tsv"
      }
    ]
}
path = Path("pipelines/molecule_3d_experimental_linker/reports/molecule_3d_experimental_v1.manual_download.json")
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(f"[INFO] manual checklist -> {path}")
PY
  exit 2
fi

echo "[INFO] resolved inputs:"
echo "  xref        = $IN_XREF"
echo "  protein_pdb = $IN_PROTEIN_PDB"

echo "[STEP] smoke build (max_rows=$SMOKE_MAX_ROWS, max_pdb=$SMOKE_MAX_PDB)"
python3 "$PIPELINE_DIR/scripts/01_build_molecule_3d_experimental_v1.py" \
  --xref "$IN_XREF" \
  --protein-pdb "$IN_PROTEIN_PDB" \
  --out "$SMOKE_TABLE" \
  --report "$SMOKE_BUILD" \
  --coverage-report "$SMOKE_COVERAGE" \
  --max-rows "$SMOKE_MAX_ROWS" \
  --max-pdb "$SMOKE_MAX_PDB" \
  --batch-size "$BATCH_SIZE" \
  --http-timeout "$HTTP_TIMEOUT" \
  --http-retries "$HTTP_RETRIES" \
  --http-sleep "$HTTP_SLEEP"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$SMOKE_TABLE" \
  --out "$SMOKE_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_3d_experimental_v1.py" \
  --xref "$IN_XREF" \
  --table "$SMOKE_TABLE" \
  --build-report "$SMOKE_BUILD" \
  --coverage-report "$SMOKE_COVERAGE" \
  --out "$SMOKE_QA"

echo "[STEP] full build"
FULL_MAX_PDB_ARG=()
if [[ -n "$FULL_MAX_PDB" ]]; then
  FULL_MAX_PDB_ARG=(--max-pdb "$FULL_MAX_PDB")
fi

python3 "$PIPELINE_DIR/scripts/01_build_molecule_3d_experimental_v1.py" \
  --xref "$IN_XREF" \
  --protein-pdb "$IN_PROTEIN_PDB" \
  --out "$OUT_TABLE" \
  --report "$REPORT_BUILD" \
  --coverage-report "$REPORT_COVERAGE" \
  "${FULL_MAX_PDB_ARG[@]}" \
  --batch-size "$BATCH_SIZE" \
  --http-timeout "$HTTP_TIMEOUT" \
  --http-retries "$HTTP_RETRIES" \
  --http-sleep "$HTTP_SLEEP"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUT_TABLE" \
  --out "$REPORT_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_3d_experimental_v1.py" \
  --xref "$IN_XREF" \
  --table "$OUT_TABLE" \
  --build-report "$REPORT_BUILD" \
  --coverage-report "$REPORT_COVERAGE" \
  --out "$REPORT_QA"

python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$OUT_TABLE" \
  "$IN_XREF" \
  "$IN_PROTEIN_PDB" \
  "$CONTRACT" \
  "$REPORT_BUILD" \
  "$REPORT_COVERAGE" \
  "$REPORT_VALIDATION" \
  "$REPORT_QA"

echo "[DONE] molecule_3d_experimental_v1 finished (smoke + full + validation + QA + manifest)."
