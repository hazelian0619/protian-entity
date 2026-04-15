#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

DATE_TAG="2026-04-12"
REPORT_DIR="pipelines/interaction_readiness/reports"
DOC_OUT="docs/interaction/INTERACTION_L2_LOCAL_READINESS_${DATE_TAG}.md"
FINAL_REPORT="${REPORT_DIR}/interaction_l2_readiness_${DATE_TAG}.json"
SAMPLE_REPORT="${REPORT_DIR}/interaction_l2_readiness_${DATE_TAG}.sample.json"
VALIDATION_REPORT="${REPORT_DIR}/interaction_l2_readiness_${DATE_TAG}.validation.json"
GATES_REPORT="${REPORT_DIR}/interaction_l2_readiness_${DATE_TAG}.gates.json"
MANIFEST_REPORT="${REPORT_DIR}/interaction_l2_readiness_${DATE_TAG}.manifest.json"
MATERIALIZE_SAMPLE_REPORT="${REPORT_DIR}/interaction_artifact_materialization_${DATE_TAG}.sample.json"
MATERIALIZE_FULL_REPORT="${REPORT_DIR}/interaction_artifact_materialization_${DATE_TAG}.json"
PPI_EDGES_VALIDATION="pipelines/edges_ppi/reports/edges_ppi_v1.validation.json"
PPI_EVIDENCE_VALIDATION="pipelines/edges_ppi/reports/ppi_evidence_v1.validation.json"
PSI_EDGES_VALIDATION="pipelines/drugbank/reports/drug_target_edges_v1.validation.json"
PSI_EVIDENCE_VALIDATION="pipelines/drugbank/reports/drug_target_evidence_v1.validation.json"

mkdir -p "$REPORT_DIR" docs/interaction

OFFREPO_ARGS=()
SOURCE_ROOT="${INTERACTION_SOURCE_ROOT:-../1218}"
if [[ -n "${INTERACTION_OFFREPO_ROOT:-}" ]]; then
  OFFREPO_ARGS+=(--offrepo-root "${INTERACTION_OFFREPO_ROOT}")
elif [[ -d "../1218" ]]; then
  OFFREPO_ARGS+=(--offrepo-root "../1218")
fi

echo "[0/9] Materialize canonical artifacts (sample)..."
python3 pipelines/interaction_readiness/scripts/materialize_interaction_artifacts.py \
  --repo-root . \
  --source-root "$SOURCE_ROOT" \
  --mode sample \
  --sample-limit 200 \
  --report "$MATERIALIZE_SAMPLE_REPORT"

echo "[1/9] Materialize canonical artifacts (full)..."
python3 pipelines/interaction_readiness/scripts/materialize_interaction_artifacts.py \
  --repo-root . \
  --source-root "$SOURCE_ROOT" \
  --mode full \
  --sample-limit 200 \
  --report "$MATERIALIZE_FULL_REPORT"

echo "[2/9] Validate PPI/PSI standardized tables..."
python3 tools/kg_validate_table.py \
  --contract pipelines/edges_ppi/contracts/edges_ppi_v1.json \
  --table data/output/edges/edges_ppi_v1.tsv \
  --out "$PPI_EDGES_VALIDATION"
python3 tools/kg_validate_table.py \
  --contract pipelines/edges_ppi/contracts/ppi_evidence_v1.json \
  --table data/output/evidence/ppi_evidence_v1.tsv \
  --out "$PPI_EVIDENCE_VALIDATION"
python3 tools/kg_validate_table.py \
  --contract pipelines/drugbank/contracts/drug_target_edges_v1.json \
  --table data/output/edges/drug_target_edges_v1.tsv \
  --out "$PSI_EDGES_VALIDATION"
python3 tools/kg_validate_table.py \
  --contract pipelines/drugbank/contracts/drug_target_evidence_v1.json \
  --table data/output/evidence/drug_target_evidence_v1.tsv \
  --out "$PSI_EVIDENCE_VALIDATION"

echo "[3/9] Readiness sample run (minimal rows)..."
python3 pipelines/interaction_readiness/scripts/build_interaction_readiness.py \
  --repo-root . \
  --mode sample \
  --sample-rows 200 \
  --report "$SAMPLE_REPORT" \
  --doc-out "docs/interaction/INTERACTION_L2_LOCAL_READINESS_${DATE_TAG}.sample.md" \
  "${OFFREPO_ARGS[@]}"

echo "[4/9] Readiness full run..."
python3 pipelines/interaction_readiness/scripts/build_interaction_readiness.py \
  --repo-root . \
  --mode full \
  --sample-rows 200 \
  --report "$FINAL_REPORT" \
  --doc-out "$DOC_OUT" \
  "${OFFREPO_ARGS[@]}"

echo "[5/9] Validate readiness report contract..."
python3 pipelines/interaction_readiness/scripts/validate_interaction_readiness.py \
  --contract pipelines/interaction_readiness/contracts/interaction_l2_readiness_report_v1.json \
  --report "$FINAL_REPORT" \
  --out "$VALIDATION_REPORT"

echo "[6/9] Build readiness gates..."
python3 pipelines/interaction_readiness/scripts/build_readiness_gates.py \
  --report "$FINAL_REPORT" \
  --out "$GATES_REPORT"

echo "[7/9] Build manifest..."
python3 tools/kg_make_manifest.py \
  --data-version "kg-interaction-l2-readiness-${DATE_TAG}" \
  --out "$MANIFEST_REPORT" \
  "$MATERIALIZE_FULL_REPORT" \
  "$PPI_EDGES_VALIDATION" \
  "$PPI_EVIDENCE_VALIDATION" \
  "$PSI_EDGES_VALIDATION" \
  "$PSI_EVIDENCE_VALIDATION" \
  "$FINAL_REPORT" \
  "$VALIDATION_REPORT" \
  "$GATES_REPORT" \
  "$DOC_OUT"

echo "[8/9] Summary checks..."
python3 - <<'PY'
import json
path='pipelines/interaction_readiness/reports/interaction_l2_readiness_2026-04-12.gates.json'
j=json.load(open(path))
print('[SUMMARY] readiness_status=',j.get('status'),'gate_pass=',j.get('gate_pass_count'),'/',j.get('gate_total'),'blocking_gap_count=',j.get('blocking_gap_count'))
PY

echo "[9/9] DONE interaction_readiness pipeline complete."
