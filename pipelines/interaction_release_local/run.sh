#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

DATE_TAG="${DATE_TAG:-2026-04-13}"
PIPELINE_DIR="pipelines/interaction_release_local"
REPORT_DIR="$PIPELINE_DIR/reports"

SUMMARY_JSON="$REPORT_DIR/interaction_release_local_${DATE_TAG}.summary.json"
GATES_JSON="$REPORT_DIR/interaction_release_local_${DATE_TAG}.gates.json"
VALIDATION_JSON="$REPORT_DIR/interaction_release_local_${DATE_TAG}.validation.json"
MANIFEST_JSON="$REPORT_DIR/interaction_release_local_${DATE_TAG}.manifest.json"
CHECKLIST_MD="../_coord/互作补洗_本地收口Checklist_${DATE_TAG}.md"

mkdir -p "$REPORT_DIR"

echo "[1/4] Build release summary + checklist..."
python3 "$PIPELINE_DIR/scripts/build_interaction_release_local.py" \
  --release-tag "$DATE_TAG" \
  --summary-out "$SUMMARY_JSON" \
  --gates-out "$GATES_JSON" \
  --checklist-out "$CHECKLIST_MD"

echo "[2/4] Validate release summary contract..."
python3 "$PIPELINE_DIR/scripts/validate_release_summary.py" \
  --contract "$PIPELINE_DIR/contracts/interaction_release_local_summary_v1.json" \
  --summary "$SUMMARY_JSON" \
  --out "$VALIDATION_JSON"

echo "[3/4] Build manifest..."
python3 tools/kg_make_manifest.py \
  --data-version "kg-interaction-release-local-${DATE_TAG}" \
  --out "$MANIFEST_JSON" \
  "$SUMMARY_JSON" \
  "$GATES_JSON" \
  "$VALIDATION_JSON" \
  "$CHECKLIST_MD" \
  "pipelines/ppi_semantic_enrichment/reports/ppi_semantic_enrichment_v2.qa.json" \
  "pipelines/psi_activity_structure_enrichment/reports/psi_activity_structure_enrichment_v2.qa.json" \
  "pipelines/rpi_site_domain_enrichment/reports/rpi_site_domain_enrichment_v2.gates.json" \
  "pipelines/interaction_cross_validation/reports/interaction_cross_validation_v2.qa.json" \
  "pipelines/interaction_ontology_mapping/reports/interaction_ontology_mapping_v2.gates.json"

echo "[4/4] Done"
python3 - <<PY
import json
p="$GATES_JSON"
j=json.load(open(p))
print(f"[SUMMARY] release_status={j.get('status')} checks={j.get('checks')}")
print(f"[SUMMARY] checklist={"$CHECKLIST_MD"}")
PY
