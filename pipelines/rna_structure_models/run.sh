#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

PIPELINE_DIR="pipelines/rna_structure_models"
REPORT_DIR="$PIPELINE_DIR/reports"
mkdir -p data/output "$REPORT_DIR"

COV_TABLE="data/output/rna_covariance_models_index_v1.tsv"
PRED_TABLE="data/output/rna_predicted_structures_v1.tsv"

CONTRACT_COV="$PIPELINE_DIR/contracts/rna_covariance_models_index_v1.json"
CONTRACT_PRED="$PIPELINE_DIR/contracts/rna_predicted_structures_v1.json"

REPORT_READY="$REPORT_DIR/rna_structure_models_v1.blocked_or_ready.json"
REPORT_SAMPLE="$REPORT_DIR/rna_structure_models_v1.sample.metrics.json"
REPORT_FULL="$REPORT_DIR/rna_structure_models_v1.metrics.json"
REPORT_VALID_COV="$REPORT_DIR/rna_covariance_models_index_v1.validation.json"
REPORT_VALID_PRED="$REPORT_DIR/rna_predicted_structures_v1.validation.json"
REPORT_FILE_QA="$REPORT_DIR/rna_structure_models_v1.file_existence.json"
REPORT_MANIFEST="$REPORT_DIR/rna_structure_models_v1.manifest.json"

DATA_VERSION="${DATA_VERSION:-kg-rna-structure-models-v1}"
SOURCE_VERSION_RFAM="${SOURCE_VERSION_RFAM:-Rfam:current}"
SOURCE_VERSION_PRED="${SOURCE_VERSION_PRED:-RNAfold/RhoFold:local}"
SAMPLE_LIMIT_COV="${SAMPLE_LIMIT_COV:-200}"
SAMPLE_LIMIT_PRED="${SAMPLE_LIMIT_PRED:-200}"
INPUT_RFAM_FAMILY="${INPUT_RFAM_FAMILY:-data/raw/rna/rfam/family.txt.gz}"
INPUT_RFAM_CM_BUNDLE="${INPUT_RFAM_CM_BUNDLE:-data/raw/rna/rfam/Rfam.cm.gz}"
INPUT_RFAM_CM_DIR="${INPUT_RFAM_CM_DIR:-data/raw/rna/rfam/cm}"
INPUT_PREDICTED_ROOT="${INPUT_PREDICTED_ROOT:-data/raw/rna/predicted_structures}"
INPUT_PREDICTED_MANIFEST="${INPUT_PREDICTED_MANIFEST:-data/raw/rna/predicted_structures/manifest.tsv}"

SCRIPT_BUILD="$PIPELINE_DIR/scripts/build_rna_structure_models.py"
SCRIPT_QA="$PIPELINE_DIR/scripts/qa_file_references.py"

echo "[0/7] Preflight check (inputs/download blockers)..."
python3 "$SCRIPT_BUILD" \
  --check-inputs \
  --rfam-family "$INPUT_RFAM_FAMILY" \
  --rfam-cm-bundle "$INPUT_RFAM_CM_BUNDLE" \
  --rfam-cm-dir "$INPUT_RFAM_CM_DIR" \
  --predicted-root "$INPUT_PREDICTED_ROOT" \
  --predicted-manifest "$INPUT_PREDICTED_MANIFEST" \
  --source-version-rfam "$SOURCE_VERSION_RFAM" \
  --source-version-pred "$SOURCE_VERSION_PRED" \
  --report "$REPORT_READY"

echo "[1/7] Dry run on minimal sample..."
python3 "$SCRIPT_BUILD" \
  --rfam-family "$INPUT_RFAM_FAMILY" \
  --rfam-cm-bundle "$INPUT_RFAM_CM_BUNDLE" \
  --rfam-cm-dir "$INPUT_RFAM_CM_DIR" \
  --predicted-root "$INPUT_PREDICTED_ROOT" \
  --predicted-manifest "$INPUT_PREDICTED_MANIFEST" \
  --limit-cov "$SAMPLE_LIMIT_COV" \
  --limit-pred "$SAMPLE_LIMIT_PRED" \
  --source-version-rfam "$SOURCE_VERSION_RFAM" \
  --source-version-pred "$SOURCE_VERSION_PRED" \
  --cov-output data/output/rna_covariance_models_index_v1.sample.tsv \
  --pred-output data/output/rna_predicted_structures_v1.sample.tsv \
  --report "$REPORT_SAMPLE"

echo "[2/7] Full run..."
python3 "$SCRIPT_BUILD" \
  --rfam-family "$INPUT_RFAM_FAMILY" \
  --rfam-cm-bundle "$INPUT_RFAM_CM_BUNDLE" \
  --rfam-cm-dir "$INPUT_RFAM_CM_DIR" \
  --predicted-root "$INPUT_PREDICTED_ROOT" \
  --predicted-manifest "$INPUT_PREDICTED_MANIFEST" \
  --source-version-rfam "$SOURCE_VERSION_RFAM" \
  --source-version-pred "$SOURCE_VERSION_PRED" \
  --cov-output "$COV_TABLE" \
  --pred-output "$PRED_TABLE" \
  --report "$REPORT_FULL"

echo "[3/7] Validate covariance table..."
python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_COV" \
  --table "$COV_TABLE" \
  --out "$REPORT_VALID_COV"

echo "[4/7] Validate predicted structures table..."
python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_PRED" \
  --table "$PRED_TABLE" \
  --out "$REPORT_VALID_PRED"

echo "[5/7] File reference QA..."
python3 "$SCRIPT_QA" \
  --cov-table "$COV_TABLE" \
  --pred-table "$PRED_TABLE" \
  --out "$REPORT_FILE_QA"

echo "[6/7] Build manifest..."
python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$COV_TABLE" \
  "$PRED_TABLE" \
  "$CONTRACT_COV" \
  "$CONTRACT_PRED" \
  "$REPORT_FULL" \
  "$REPORT_VALID_COV" \
  "$REPORT_VALID_PRED" \
  "$REPORT_FILE_QA"

echo "[7/7] Done."
