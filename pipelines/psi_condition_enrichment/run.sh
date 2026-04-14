#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

PIPELINE_DIR="pipelines/psi_condition_enrichment"

INPUT_V2="data/output/evidence/psi_activity_context_v2.tsv"
OUT_V3="data/output/evidence/psi_activity_context_v3.tsv"
OUT_AUDIT="data/output/evidence/psi_condition_parse_audit_v3.tsv"

SMOKE_V3="$PIPELINE_DIR/reports/psi_activity_context_v3.smoke.tsv"
SMOKE_AUDIT="$PIPELINE_DIR/reports/psi_condition_parse_audit_v3.smoke.tsv"

CONTRACT_V3="$PIPELINE_DIR/contracts/psi_activity_context_v3.json"
CONTRACT_AUDIT="$PIPELINE_DIR/contracts/psi_condition_parse_audit_v3.json"

REPORT_BLOCKED="$PIPELINE_DIR/reports/psi_condition_enrichment_v3.blocked_missing_inputs.json"
REPORT_SMOKE_BUILD="$PIPELINE_DIR/reports/psi_condition_enrichment_v3.smoke.build.json"
REPORT_SMOKE_V3_VALID="$PIPELINE_DIR/reports/psi_activity_context_v3.smoke.validation.json"
REPORT_SMOKE_AUDIT_VALID="$PIPELINE_DIR/reports/psi_condition_parse_audit_v3.smoke.validation.json"

REPORT_BUILD="$PIPELINE_DIR/reports/psi_condition_enrichment_v3.build.json"
REPORT_V3_VALID="$PIPELINE_DIR/reports/psi_activity_context_v3.validation.json"
REPORT_AUDIT_VALID="$PIPELINE_DIR/reports/psi_condition_parse_audit_v3.validation.json"
REPORT_QA="$PIPELINE_DIR/reports/psi_condition_enrichment_v3.qa.json"
REPORT_RESULT="$PIPELINE_DIR/reports/psi_condition_enrichment_v3.result_evidence.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/psi_condition_enrichment_v3.manifest.json"

REPORT_SAMPLE20="$PIPELINE_DIR/reports/psi_activity_context_v3.sample20.tsv"

SMOKE_ROWS="${SMOKE_ROWS:-50000}"
DATA_VERSION="${DATA_VERSION:-kg-psi-condition-v3-local}"
FETCH_DATE="${FETCH_DATE:-$(date -u +%F)}"

mkdir -p "$(dirname "$OUT_V3")" "$PIPELINE_DIR/reports"

CHEMBL_M3_DB="${CHEMBL_M3_DB:-}"
if [[ -z "$CHEMBL_M3_DB" ]]; then
  if [[ -s "data/output/molecules/chembl_m3.sqlite" ]]; then
    CHEMBL_M3_DB="data/output/molecules/chembl_m3.sqlite"
  elif [[ -s "/Users/pluviophile/graph/12182/out/m3/chembl_m3.sqlite" ]]; then
    CHEMBL_M3_DB="/Users/pluviophile/graph/12182/out/m3/chembl_m3.sqlite"
  fi
fi

OPTIONAL_BINDINGDB="data/input/bindingdb/bindingdb_conditions.tsv"

echo "[STEP 0/8] preflight + blocked_missing_inputs report"
python3 - <<PY
import json
from pathlib import Path

required = {
  "psi_activity_context_v2": "${INPUT_V2}",
}
optional = {
  "chembl_m3_db": "${CHEMBL_M3_DB}",
  "bindingdb_condition_source": "${OPTIONAL_BINDINGDB}",
}

missing_required = []
req_info = {}
for k, p in required.items():
    exists = bool(p) and Path(p).is_file()
    req_info[k] = {"path": p, "exists": exists}
    if not exists:
        missing_required.append(k)

opt_info = {}
for k, p in optional.items():
    exists = bool(p) and Path(p).is_file()
    size = Path(p).stat().st_size if exists else 0
    opt_info[k] = {"path": p, "exists": exists, "size_bytes": size, "used_by_pipeline": False}

obj = {
  "name": "psi_condition_enrichment_v3.blocked_missing_inputs",
  "required_inputs": req_info,
  "optional_inputs": opt_info,
  "missing_required": missing_required,
  "blocked": len(missing_required) > 0,
}
Path("${REPORT_BLOCKED}").write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"blocked": obj["blocked"], "missing_required": missing_required}))
PY

if [[ ! -f "$INPUT_V2" ]]; then
  echo "[BLOCKED] missing required input: $INPUT_V2" >&2
  exit 2
fi

if [[ ! -f "$CONTRACT_V3" || ! -f "$CONTRACT_AUDIT" ]]; then
  echo "[ERROR] missing contract(s) under $PIPELINE_DIR/contracts" >&2
  exit 2
fi

echo "[STEP 1/8] extractor unit tests"
python3 -m unittest pipelines/psi_condition_enrichment/tests/test_condition_extractors.py -v

echo "[STEP 2/8] smoke run (rows=${SMOKE_ROWS})"
python3 "$PIPELINE_DIR/scripts/01_build_psi_condition_enrichment_v3.py" \
  --in-v2 "$INPUT_V2" \
  --out-v3 "$SMOKE_V3" \
  --out-audit "$SMOKE_AUDIT" \
  --report "$REPORT_SMOKE_BUILD" \
  --fetch-date "$FETCH_DATE" \
  --max-rows "$SMOKE_ROWS"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_V3" \
  --table "$SMOKE_V3" \
  --out "$REPORT_SMOKE_V3_VALID"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_AUDIT" \
  --table "$SMOKE_AUDIT" \
  --out "$REPORT_SMOKE_AUDIT_VALID"


echo "[STEP 3/8] full run"
python3 "$PIPELINE_DIR/scripts/01_build_psi_condition_enrichment_v3.py" \
  --in-v2 "$INPUT_V2" \
  --out-v3 "$OUT_V3" \
  --out-audit "$OUT_AUDIT" \
  --report "$REPORT_BUILD" \
  --fetch-date "$FETCH_DATE"


echo "[STEP 4/8] contract validation"
python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_V3" \
  --table "$OUT_V3" \
  --out "$REPORT_V3_VALID"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_AUDIT" \
  --table "$OUT_AUDIT" \
  --out "$REPORT_AUDIT_VALID"


echo "[STEP 5/8] QA gates"
python3 "$PIPELINE_DIR/scripts/02_qa_psi_condition_enrichment_v3.py" \
  --v2 "$INPUT_V2" \
  --v3 "$OUT_V3" \
  --v3-validation "$REPORT_V3_VALID" \
  --audit-validation "$REPORT_AUDIT_VALID" \
  --build-report "$REPORT_BUILD" \
  --out "$REPORT_QA"


echo "[STEP 6/8] result evidence + sample20"
python3 "$PIPELINE_DIR/scripts/03_report_psi_condition_results_v3.py" \
  --v3 "$OUT_V3" \
  --audit "$OUT_AUDIT" \
  --build "$REPORT_BUILD" \
  --qa "$REPORT_QA" \
  --sample "$REPORT_SAMPLE20" \
  --out "$REPORT_RESULT" \
  --sample-size 20


echo "[STEP 7/8] manifest"
python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$INPUT_V2" \
  "$OUT_V3" \
  "$OUT_AUDIT" \
  "$CONTRACT_V3" \
  "$CONTRACT_AUDIT" \
  "$REPORT_BLOCKED" \
  "$REPORT_BUILD" \
  "$REPORT_V3_VALID" \
  "$REPORT_AUDIT_VALID" \
  "$REPORT_QA" \
  "$REPORT_RESULT" \
  "$REPORT_SAMPLE20"


echo "[STEP 8/8] summary"
python3 - <<'PY'
import json
from pathlib import Path
qa = json.loads(Path('pipelines/psi_condition_enrichment/reports/psi_condition_enrichment_v3.qa.json').read_text(encoding='utf-8'))
status = 'PASS' if qa.get('passed') else 'FAIL'
m2 = qa.get('metrics',{}).get('v2',{})
m3 = qa.get('metrics',{}).get('v3',{})
print(f'[{status}] psi_condition_enrichment_v3')
print(f"  condition_context: {m2.get('condition_context',0):.4f} -> {m3.get('condition_context',0):.4f}")
print(f"  condition_pH: {m2.get('condition_pH',0):.4f} -> {m3.get('condition_pH',0):.4f}")
print(f"  condition_temperature_c: {m2.get('condition_temperature_c',0):.4f} -> {m3.get('condition_temperature_c',0):.4f}")
print('  sample20: pipelines/psi_condition_enrichment/reports/psi_activity_context_v3.sample20.tsv')
PY

echo "[DONE] psi_condition_enrichment_v3 finished."
SH && chmod +x pipelines/psi_condition_enrichment/run.sh