#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

PIPELINE_DIR="pipelines/psi_activity_structure_enrichment"

INPUT_PSI_DB="data/output/molecules/chembl_m3.sqlite"
INPUT_PROTEIN_PDB="data/output/protein/pdb_structures_v1.tsv"

OUT_ACTIVITY="data/output/evidence/psi_activity_context_v2.tsv"
OUT_STRUCTURE="data/output/evidence/psi_structure_evidence_v2.tsv"

CONTRACT_ACTIVITY="$PIPELINE_DIR/contracts/psi_activity_context_v2.json"
CONTRACT_STRUCTURE="$PIPELINE_DIR/contracts/psi_structure_evidence_v2.json"

REPORT_SMOKE_BUILD="$PIPELINE_DIR/reports/psi_activity_structure_enrichment_v2.smoke.build.json"
REPORT_SMOKE_ACTIVITY_VALIDATION="$PIPELINE_DIR/reports/psi_activity_context_v2.smoke.validation.json"
REPORT_SMOKE_STRUCTURE_VALIDATION="$PIPELINE_DIR/reports/psi_structure_evidence_v2.smoke.validation.json"

REPORT_BUILD="$PIPELINE_DIR/reports/psi_activity_structure_enrichment_v2.build.json"
REPORT_ACTIVITY_VALIDATION="$PIPELINE_DIR/reports/psi_activity_context_v2.validation.json"
REPORT_STRUCTURE_VALIDATION="$PIPELINE_DIR/reports/psi_structure_evidence_v2.validation.json"
REPORT_QA="$PIPELINE_DIR/reports/psi_activity_structure_enrichment_v2.qa.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/psi_activity_structure_enrichment_v2.manifest.json"

SMOKE_ACTIVITY_OUT="$PIPELINE_DIR/reports/psi_activity_context_v2.smoke.tsv"
SMOKE_STRUCTURE_OUT="$PIPELINE_DIR/reports/psi_structure_evidence_v2.smoke.tsv"

DATA_VERSION="${DATA_VERSION:-kg-data-local}"
FETCH_DATE="${FETCH_DATE:-$(date -u +%F)}"
SMOKE_ROWS="${SMOKE_ROWS:-50000}"
MIN_LIGAND_COUNT="${MIN_LIGAND_COUNT:-1}"

if [[ ! -f "$INPUT_PSI_DB" ]]; then
  echo "[ERROR] missing PSI DB: $INPUT_PSI_DB" >&2
  exit 2
fi
if [[ ! -f "$INPUT_PROTEIN_PDB" ]]; then
  echo "[ERROR] missing protein PDB table: $INPUT_PROTEIN_PDB" >&2
  exit 2
fi

mkdir -p "$(dirname "$OUT_ACTIVITY")" "$(dirname "$OUT_STRUCTURE")" "$PIPELINE_DIR/reports"

echo "[STEP 1/6] smoke run (minimal sample)"
python3 "$PIPELINE_DIR/scripts/01_build_psi_activity_structure_enrichment_v2.py" \
  --chembl-m3-db "$INPUT_PSI_DB" \
  --protein-pdb "$INPUT_PROTEIN_PDB" \
  --out-activity "$SMOKE_ACTIVITY_OUT" \
  --out-structure "$SMOKE_STRUCTURE_OUT" \
  --report "$REPORT_SMOKE_BUILD" \
  --fetch-date "$FETCH_DATE" \
  --max-rows "$SMOKE_ROWS" \
  --min-ligand-count "$MIN_LIGAND_COUNT"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_ACTIVITY" \
  --table "$SMOKE_ACTIVITY_OUT" \
  --out "$REPORT_SMOKE_ACTIVITY_VALIDATION"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_STRUCTURE" \
  --table "$SMOKE_STRUCTURE_OUT" \
  --out "$REPORT_SMOKE_STRUCTURE_VALIDATION"

echo "[STEP 2/6] full run"
python3 "$PIPELINE_DIR/scripts/01_build_psi_activity_structure_enrichment_v2.py" \
  --chembl-m3-db "$INPUT_PSI_DB" \
  --protein-pdb "$INPUT_PROTEIN_PDB" \
  --out-activity "$OUT_ACTIVITY" \
  --out-structure "$OUT_STRUCTURE" \
  --report "$REPORT_BUILD" \
  --fetch-date "$FETCH_DATE" \
  --min-ligand-count "$MIN_LIGAND_COUNT"

echo "[STEP 3/6] contract validation"
python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_ACTIVITY" \
  --table "$OUT_ACTIVITY" \
  --out "$REPORT_ACTIVITY_VALIDATION"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_STRUCTURE" \
  --table "$OUT_STRUCTURE" \
  --out "$REPORT_STRUCTURE_VALIDATION"

echo "[STEP 4/6] QA gates / acceptance"
python3 "$PIPELINE_DIR/scripts/02_qa_psi_activity_structure_enrichment_v2.py" \
  --psi-db "$INPUT_PSI_DB" \
  --activity-context "$OUT_ACTIVITY" \
  --structure-evidence "$OUT_STRUCTURE" \
  --activity-validation "$REPORT_ACTIVITY_VALIDATION" \
  --structure-validation "$REPORT_STRUCTURE_VALIDATION" \
  --out "$REPORT_QA"

echo "[STEP 5/6] manifest"
python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$OUT_ACTIVITY" \
  "$OUT_STRUCTURE" \
  "$INPUT_PSI_DB" \
  "$INPUT_PROTEIN_PDB" \
  "$CONTRACT_ACTIVITY" \
  "$CONTRACT_STRUCTURE" \
  "$REPORT_BUILD" \
  "$REPORT_ACTIVITY_VALIDATION" \
  "$REPORT_STRUCTURE_VALIDATION" \
  "$REPORT_QA"

echo "[STEP 6/6] acceptance summary"
python3 - <<'PY'
import json
from pathlib import Path
p = Path("pipelines/psi_activity_structure_enrichment/reports/psi_activity_structure_enrichment_v2.qa.json")
obj = json.loads(p.read_text(encoding="utf-8"))
status = "PASS" if obj.get("passed") else "FAIL"
print(f"[{status}] psi_activity_structure_enrichment_v2")
for g in obj.get("gates", []):
    mark = "PASS" if g.get("passed") else "FAIL"
    print(f"  - {mark} {g.get('id')}: {g.get('detail')}")
PY

echo "[DONE] psi_activity_structure_enrichment_v2: smoke + full + validation + QA + manifest finished."
