#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

PIPELINE_DIR="pipelines/interaction_cross_validation"

resolve_first_existing() {
  for p in "$@"; do
    if [[ -n "$p" && -f "$p" ]]; then
      echo "$p"
      return 0
    fi
  done
  return 1
}

PPI_EDGES="$(resolve_first_existing \
  "data/output/edges/edges_ppi_v1.tsv" \
  "pipelines/edges_ppi/data/output/edges/edges_ppi_v1.tsv" || true)"
PPI_EVIDENCE="$(resolve_first_existing \
  "data/output/evidence/ppi_evidence_v1.tsv" \
  "pipelines/edges_ppi/data/output/evidence/ppi_evidence_v1.tsv" || true)"

PPI_METHOD_CONTEXT="$(resolve_first_existing "data/output/evidence/ppi_method_context_v2.tsv" || true)"
PPI_FUNCTION_CONTEXT="$(resolve_first_existing "data/output/evidence/ppi_function_context_v2.tsv" || true)"

PSI_EDGES="$(resolve_first_existing "data/output/edges/drug_target_edges_v1.tsv" || true)"
PSI_EVIDENCE="$(resolve_first_existing "data/output/evidence/drug_target_evidence_v1.tsv" || true)"
PSI_ACTIVITY_CONTEXT="$(resolve_first_existing "data/output/evidence/psi_activity_context_v2.tsv" || true)"
PSI_STRUCTURE_EVIDENCE="$(resolve_first_existing "data/output/evidence/psi_structure_evidence_v2.tsv" || true)"

RPI_EDGES="$(resolve_first_existing \
  "data/output/edges/rna_protein_edges_v1.tsv" \
  "../protian-entity/data/output/edges/rna_protein_edges_v1.tsv" \
  "data/output/edges/gene_to_rna_edges_v1.tsv" || true)"
RPI_EVIDENCE="$(resolve_first_existing \
  "data/output/evidence/rna_protein_evidence_v1.tsv" \
  "../protian-entity/data/output/evidence/rna_protein_evidence_v1.tsv" \
  "data/output/evidence/gene_to_rna_evidence_v1.tsv" || true)"

RPI_SITE_CONTEXT="$(resolve_first_existing \
  "data/output/evidence/rpi_site_context_v2.tsv" \
  "../protian-entity/data/output/evidence/rpi_site_context_v2.tsv" || true)"
RPI_DOMAIN_CONTEXT="$(resolve_first_existing \
  "data/output/evidence/rpi_domain_context_v2.tsv" \
  "../protian-entity/data/output/evidence/rpi_domain_context_v2.tsv" || true)"
RPI_FUNCTION_CONTEXT="$(resolve_first_existing \
  "data/output/evidence/rpi_function_context_v2.tsv" \
  "../protian-entity/data/output/evidence/rpi_function_context_v2.tsv" || true)"

MOLECULE_XREF="$(resolve_first_existing \
  "data/output/molecules/molecule_xref_core_v2.tsv" \
  "data/output/molecules/molecule_xref_core_v1.tsv" || true)"

OUT_CROSS="data/output/evidence/interaction_cross_validation_v2.tsv"
OUT_AGG="data/output/evidence/interaction_aggregate_score_v2.tsv"

SMOKE_CROSS="$PIPELINE_DIR/.cache/interaction_cross_validation_v2.smoke.tsv"
SMOKE_AGG="$PIPELINE_DIR/.cache/interaction_aggregate_score_v2.smoke.tsv"

REPORT_BLOCKED_OR_READY="$PIPELINE_DIR/reports/interaction_cross_validation_v2.blocked_or_ready.json"
REPORT_SMOKE_BUILD="$PIPELINE_DIR/reports/interaction_cross_validation_v2.smoke.build.json"
REPORT_SMOKE_CROSS_VALID="$PIPELINE_DIR/reports/interaction_cross_validation_v2.smoke.validation.json"
REPORT_SMOKE_AGG_VALID="$PIPELINE_DIR/reports/interaction_aggregate_score_v2.smoke.validation.json"
REPORT_SMOKE_QA="$PIPELINE_DIR/reports/interaction_cross_validation_v2.smoke.qa.json"

REPORT_BUILD="$PIPELINE_DIR/reports/interaction_cross_validation_v2.build.json"
REPORT_CROSS_VALID="$PIPELINE_DIR/reports/interaction_cross_validation_v2.validation.json"
REPORT_AGG_VALID="$PIPELINE_DIR/reports/interaction_aggregate_score_v2.validation.json"
REPORT_QA="$PIPELINE_DIR/reports/interaction_cross_validation_v2.qa.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/interaction_cross_validation_v2.manifest.json"

CONTRACT_CROSS="$PIPELINE_DIR/contracts/interaction_cross_validation_v2.json"
CONTRACT_AGG="$PIPELINE_DIR/contracts/interaction_aggregate_score_v2.json"

SMOKE_LIMIT_PER_TYPE="${SMOKE_LIMIT_PER_TYPE:-5000}"
DATA_VERSION="${DATA_VERSION:-kg-interaction-cross-validation-v2-local}"
FETCH_DATE="${FETCH_DATE:-$(date -u +%F)}"

mkdir -p "$(dirname "$OUT_CROSS")" "$PIPELINE_DIR/reports" "$PIPELINE_DIR/.cache"

common_args=(
  --ppi-edges "$PPI_EDGES"
  --ppi-evidence "$PPI_EVIDENCE"
  --psi-edges "$PSI_EDGES"
  --psi-evidence "$PSI_EVIDENCE"
  --psi-activity-context "$PSI_ACTIVITY_CONTEXT"
  --psi-structure-evidence "$PSI_STRUCTURE_EVIDENCE"
  --rpi-edges "$RPI_EDGES"
  --rpi-evidence "$RPI_EVIDENCE"
  --molecule-xref "$MOLECULE_XREF"
  --fetch-date "$FETCH_DATE"
)

if [[ -n "$PPI_METHOD_CONTEXT" ]]; then
  common_args+=(--ppi-method-context "$PPI_METHOD_CONTEXT")
fi
if [[ -n "$PPI_FUNCTION_CONTEXT" ]]; then
  common_args+=(--ppi-function-context "$PPI_FUNCTION_CONTEXT")
fi
if [[ -n "$RPI_SITE_CONTEXT" ]]; then
  common_args+=(--rpi-site-context "$RPI_SITE_CONTEXT")
fi
if [[ -n "$RPI_DOMAIN_CONTEXT" ]]; then
  common_args+=(--rpi-domain-context "$RPI_DOMAIN_CONTEXT")
fi
if [[ -n "$RPI_FUNCTION_CONTEXT" ]]; then
  common_args+=(--rpi-function-context "$RPI_FUNCTION_CONTEXT")
fi

echo "[STEP 0/7] preflight check"
set +e
python3 "$PIPELINE_DIR/scripts/01_build_interaction_cross_validation_v2.py" \
  "${common_args[@]}" \
  --cross-output "$SMOKE_CROSS" \
  --aggregate-output "$SMOKE_AGG" \
  --report "$REPORT_BLOCKED_OR_READY" \
  --check-inputs
preflight_rc=$?
set -e
if [[ $preflight_rc -ne 0 ]]; then
  echo "[BLOCKED] preflight failed, see $REPORT_BLOCKED_OR_READY" >&2
  exit $preflight_rc
fi

echo "[STEP 1/7] smoke run (limit_per_type=$SMOKE_LIMIT_PER_TYPE)"
python3 "$PIPELINE_DIR/scripts/01_build_interaction_cross_validation_v2.py" \
  "${common_args[@]}" \
  --cross-output "$SMOKE_CROSS" \
  --aggregate-output "$SMOKE_AGG" \
  --report "$REPORT_SMOKE_BUILD" \
  --limit-per-type "$SMOKE_LIMIT_PER_TYPE"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_CROSS" \
  --table "$SMOKE_CROSS" \
  --out "$REPORT_SMOKE_CROSS_VALID"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_AGG" \
  --table "$SMOKE_AGG" \
  --out "$REPORT_SMOKE_AGG_VALID"

python3 "$PIPELINE_DIR/scripts/02_qa_interaction_cross_validation_v2.py" \
  --cross "$SMOKE_CROSS" \
  --aggregate "$SMOKE_AGG" \
  --build-report "$REPORT_SMOKE_BUILD" \
  --cross-validation "$REPORT_SMOKE_CROSS_VALID" \
  --aggregate-validation "$REPORT_SMOKE_AGG_VALID" \
  --out "$REPORT_SMOKE_QA"

echo "[STEP 2/7] full run"
python3 "$PIPELINE_DIR/scripts/01_build_interaction_cross_validation_v2.py" \
  "${common_args[@]}" \
  --cross-output "$OUT_CROSS" \
  --aggregate-output "$OUT_AGG" \
  --report "$REPORT_BUILD"

echo "[STEP 3/7] contract validation"
python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_CROSS" \
  --table "$OUT_CROSS" \
  --out "$REPORT_CROSS_VALID"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_AGG" \
  --table "$OUT_AGG" \
  --out "$REPORT_AGG_VALID"

echo "[STEP 4/7] QA gates"
python3 "$PIPELINE_DIR/scripts/02_qa_interaction_cross_validation_v2.py" \
  --cross "$OUT_CROSS" \
  --aggregate "$OUT_AGG" \
  --build-report "$REPORT_BUILD" \
  --cross-validation "$REPORT_CROSS_VALID" \
  --aggregate-validation "$REPORT_AGG_VALID" \
  --out "$REPORT_QA"

echo "[STEP 5/7] manifest"
manifest_inputs=(
  "$OUT_CROSS"
  "$OUT_AGG"
  "$CONTRACT_CROSS"
  "$CONTRACT_AGG"
  "$REPORT_BUILD"
  "$REPORT_CROSS_VALID"
  "$REPORT_AGG_VALID"
  "$REPORT_QA"
  "$PPI_EDGES"
  "$PPI_EVIDENCE"
  "$PSI_EDGES"
  "$PSI_EVIDENCE"
  "$PSI_ACTIVITY_CONTEXT"
  "$PSI_STRUCTURE_EVIDENCE"
  "$RPI_EDGES"
  "$RPI_EVIDENCE"
  "$MOLECULE_XREF"
)

for f in "$PPI_METHOD_CONTEXT" "$PPI_FUNCTION_CONTEXT" "$RPI_SITE_CONTEXT" "$RPI_DOMAIN_CONTEXT" "$RPI_FUNCTION_CONTEXT"; do
  if [[ -n "$f" && -f "$f" ]]; then
    manifest_inputs+=("$f")
  fi
done

python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "${manifest_inputs[@]}"

echo "[STEP 6/7] summary"
python3 - <<'PY'
import json
from pathlib import Path
qa = json.loads(Path('pipelines/interaction_cross_validation/reports/interaction_cross_validation_v2.qa.json').read_text())
build = json.loads(Path('pipelines/interaction_cross_validation/reports/interaction_cross_validation_v2.build.json').read_text())
status = 'PASS' if qa.get('passed') else 'FAIL'
print(f'[{status}] interaction_cross_validation_v2')
print('  types_in_cross:', qa.get('metrics',{}).get('types_in_cross'))
print('  aggregate_score:', qa.get('metrics',{}).get('aggregate_score'))
print('  psi_b_update_integration:', build.get('psi_b_update_integration',{}))
PY

echo "[STEP 7/7] done"
