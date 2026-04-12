#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

PIPELINE_DIR="pipelines/molecule_activity_fusion"

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

# required
REQ_CHEMBL_M3_DB="${CHEMBL_M3_DB_PATH:-data/output/molecules/chembl_m3.sqlite}"
IN_CHEMBL_M3_DB="$(resolve_input "$REQ_CHEMBL_M3_DB" "../12182/out/m3/chembl_m3.sqlite" "../1218/data/output/molecules/chembl_m3.sqlite")"

# optional externals (auto-download in build script if missing)
BINDINGDB_ZIP="${BINDINGDB_ZIP_PATH:-$PIPELINE_DIR/.cache/BindingDB_PubChem_202604_tsv.zip}"
PDBBIND_INDEX_TAR="${PDBBIND_INDEX_TAR_PATH:-$PIPELINE_DIR/.cache/PDBbind_v2020_plain_text_index.tar.gz}"

OUT_EVIDENCE="data/output/evidence/molecule_activity_evidence_v2.tsv"
OUT_EDGES="data/output/edges/molecule_target_edges_v2.tsv"

CONTRACT_EVIDENCE="$PIPELINE_DIR/contracts/molecule_activity_evidence_v2.json"
CONTRACT_EDGES="$PIPELINE_DIR/contracts/molecule_target_edges_v2.json"

SMOKE_ROWS="${SMOKE_ROWS:-50000}"
DATA_VERSION="${DATA_VERSION:-kg-data-local}"
NO_AUTO_DOWNLOAD="${NO_AUTO_DOWNLOAD:-0}"

SMOKE_EVIDENCE="$PIPELINE_DIR/.cache/molecule_activity_evidence_v2.smoke.tsv"
SMOKE_EDGES="$PIPELINE_DIR/.cache/molecule_target_edges_v2.smoke.tsv"
SMOKE_BUILD="$PIPELINE_DIR/reports/molecule_activity_fusion_v2.smoke.build.json"
SMOKE_CONFLICT_TSV="$PIPELINE_DIR/reports/molecule_activity_conflict_audit_v2.smoke.tsv"
SMOKE_CONFLICT_JSON="$PIPELINE_DIR/reports/molecule_activity_conflict_audit_v2.smoke.json"
SMOKE_MANUAL="$PIPELINE_DIR/reports/molecule_activity_fusion_v2.smoke.manual_download.json"
SMOKE_EVIDENCE_VALIDATION="$PIPELINE_DIR/reports/molecule_activity_evidence_v2.smoke.validation.json"
SMOKE_EDGES_VALIDATION="$PIPELINE_DIR/reports/molecule_target_edges_v2.smoke.validation.json"
SMOKE_QA="$PIPELINE_DIR/reports/molecule_activity_fusion_v2.smoke.qa.json"

REPORT_BUILD="$PIPELINE_DIR/reports/molecule_activity_fusion_v2.build.json"
REPORT_CONFLICT_TSV="$PIPELINE_DIR/reports/molecule_activity_conflict_audit_v2.tsv"
REPORT_CONFLICT_JSON="$PIPELINE_DIR/reports/molecule_activity_conflict_audit_v2.json"
REPORT_MANUAL="$PIPELINE_DIR/reports/molecule_activity_fusion_v2.manual_download.json"
REPORT_EVIDENCE_VALIDATION="$PIPELINE_DIR/reports/molecule_activity_evidence_v2.validation.json"
REPORT_EDGES_VALIDATION="$PIPELINE_DIR/reports/molecule_target_edges_v2.validation.json"
REPORT_QA="$PIPELINE_DIR/reports/molecule_activity_fusion_v2.qa.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/molecule_activity_fusion_v2.manifest.json"

mkdir -p "$(dirname "$OUT_EVIDENCE")" "$(dirname "$OUT_EDGES")" "$PIPELINE_DIR/reports" "$PIPELINE_DIR/.cache"

if [[ -z "$IN_CHEMBL_M3_DB" ]]; then
  echo "[ERROR] missing required input: $REQ_CHEMBL_M3_DB" >&2
  python3 - <<'PY'
import json
from pathlib import Path
report = {
  "name": "molecule_activity_fusion_v2.manual_download",
  "reason": "required chembl_m3.sqlite missing",
  "download_checklist": [
    {
      "name": "chembl_m3.sqlite",
      "expected_path": "data/output/molecules/chembl_m3.sqlite",
      "fallback_path": "../12182/out/m3/chembl_m3.sqlite",
      "sha256_cmd": "sha256sum data/output/molecules/chembl_m3.sqlite"
    }
  ]
}
p = Path("pipelines/molecule_activity_fusion/reports/molecule_activity_fusion_v2.manual_download.json")
p.parent.mkdir(parents=True, exist_ok=True)
p.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(f"[INFO] wrote checklist: {p}")
PY
  exit 2
fi

echo "[INFO] resolved inputs:"
echo "  chembl_m3_db     = $IN_CHEMBL_M3_DB"
echo "  bindingdb_zip    = $BINDINGDB_ZIP (optional, auto-download if missing)"
echo "  pdbbind_index    = $PDBBIND_INDEX_TAR (optional, auto-download if missing)"

AUTO_FLAG=""
if [[ "$NO_AUTO_DOWNLOAD" == "1" ]]; then
  AUTO_FLAG="--no-auto-download"
fi

echo "[STEP 1/6] smoke build (max_chembl_rows=$SMOKE_ROWS)"
python3 "$PIPELINE_DIR/scripts/01_build_molecule_activity_fusion_v2.py" \
  --chembl-m3-db "$IN_CHEMBL_M3_DB" \
  --out-evidence "$SMOKE_EVIDENCE" \
  --out-edges "$SMOKE_EDGES" \
  --report "$SMOKE_BUILD" \
  --conflict-audit-tsv "$SMOKE_CONFLICT_TSV" \
  --conflict-audit-json "$SMOKE_CONFLICT_JSON" \
  --manual-download-report "$SMOKE_MANUAL" \
  --bindingdb-zip "$BINDINGDB_ZIP" \
  --pdbbind-index-tar "$PDBBIND_INDEX_TAR" \
  --cache-dir "$PIPELINE_DIR/.cache" \
  --max-chembl-rows "$SMOKE_ROWS" \
  $AUTO_FLAG

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_EVIDENCE" \
  --table "$SMOKE_EVIDENCE" \
  --out "$SMOKE_EVIDENCE_VALIDATION"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_EDGES" \
  --table "$SMOKE_EDGES" \
  --out "$SMOKE_EDGES_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_activity_fusion_v2.py" \
  --chembl-m3-db "$IN_CHEMBL_M3_DB" \
  --evidence "$SMOKE_EVIDENCE" \
  --edges "$SMOKE_EDGES" \
  --conflict-audit-tsv "$SMOKE_CONFLICT_TSV" \
  --build-report "$SMOKE_BUILD" \
  --evidence-validation "$SMOKE_EVIDENCE_VALIDATION" \
  --edges-validation "$SMOKE_EDGES_VALIDATION" \
  --out "$SMOKE_QA" \
  --require-strict 0

echo "[STEP 2/6] full build"
python3 "$PIPELINE_DIR/scripts/01_build_molecule_activity_fusion_v2.py" \
  --chembl-m3-db "$IN_CHEMBL_M3_DB" \
  --out-evidence "$OUT_EVIDENCE" \
  --out-edges "$OUT_EDGES" \
  --report "$REPORT_BUILD" \
  --conflict-audit-tsv "$REPORT_CONFLICT_TSV" \
  --conflict-audit-json "$REPORT_CONFLICT_JSON" \
  --manual-download-report "$REPORT_MANUAL" \
  --bindingdb-zip "$BINDINGDB_ZIP" \
  --pdbbind-index-tar "$PDBBIND_INDEX_TAR" \
  --cache-dir "$PIPELINE_DIR/.cache" \
  $AUTO_FLAG

echo "[STEP 3/6] contract validation"
python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_EVIDENCE" \
  --table "$OUT_EVIDENCE" \
  --out "$REPORT_EVIDENCE_VALIDATION"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_EDGES" \
  --table "$OUT_EDGES" \
  --out "$REPORT_EDGES_VALIDATION"

echo "[STEP 4/6] QA gates"
python3 "$PIPELINE_DIR/scripts/02_qa_molecule_activity_fusion_v2.py" \
  --chembl-m3-db "$IN_CHEMBL_M3_DB" \
  --evidence "$OUT_EVIDENCE" \
  --edges "$OUT_EDGES" \
  --conflict-audit-tsv "$REPORT_CONFLICT_TSV" \
  --build-report "$REPORT_BUILD" \
  --evidence-validation "$REPORT_EVIDENCE_VALIDATION" \
  --edges-validation "$REPORT_EDGES_VALIDATION" \
  --out "$REPORT_QA" \
  --require-strict 1

echo "[STEP 5/6] manifest"
python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$OUT_EVIDENCE" \
  "$OUT_EDGES" \
  "$IN_CHEMBL_M3_DB" \
  "$CONTRACT_EVIDENCE" \
  "$CONTRACT_EDGES" \
  "$REPORT_BUILD" \
  "$REPORT_CONFLICT_TSV" \
  "$REPORT_CONFLICT_JSON" \
  "$REPORT_EVIDENCE_VALIDATION" \
  "$REPORT_EDGES_VALIDATION" \
  "$REPORT_QA" \
  "$REPORT_MANUAL"

echo "[STEP 6/6] final summary"
python3 - <<'PY'
import json
from pathlib import Path
qa = json.loads(Path("pipelines/molecule_activity_fusion/reports/molecule_activity_fusion_v2.qa.json").read_text(encoding="utf-8"))
status = "PASS" if qa.get("passed") else "FAIL"
print(f"[{status}] molecule_activity_fusion_v2")
for g in qa.get("gates", []):
    mark = "PASS" if g.get("passed") else "FAIL"
    print(f"  - {mark} {g.get('id')}")
PY

echo "[DONE] molecule_activity_fusion_v2 completed (smoke + full + validation + qa + manifest)."
