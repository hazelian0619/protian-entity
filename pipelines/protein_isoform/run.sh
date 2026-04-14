#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

mkdir -p data/output/protein pipelines/protein_isoform/reports

DATA_VERSION="${DATA_VERSION:-kg-data-local}"

ISOFORM_TABLE="data/output/protein/protein_isoform_v1.tsv"
MAP_TABLE="data/output/protein/protein_isoform_map_v1.tsv"

ISOFORM_CONTRACT="pipelines/protein_isoform/contracts/protein_isoform_v1.json"
MAP_CONTRACT="pipelines/protein_isoform/contracts/protein_isoform_map_v1.json"

ISOFORM_QA="pipelines/protein_isoform/reports/protein_isoform_v1.qa.json"
MAP_QA="pipelines/protein_isoform/reports/protein_isoform_map_v1.qa.json"

ISOFORM_VALIDATION="pipelines/protein_isoform/reports/protein_isoform_v1.validation.json"
MAP_VALIDATION="pipelines/protein_isoform/reports/protein_isoform_map_v1.validation.json"

ISOFORM_MANIFEST="pipelines/protein_isoform/reports/protein_isoform_v1.manifest.json"
MAP_MANIFEST="pipelines/protein_isoform/reports/protein_isoform_map_v1.manifest.json"

echo "[1/4] Build isoform layer + QA reports..."
python3 pipelines/protein_isoform/scripts/build_protein_isoform_layer.py \
  --input data/processed/protein_master_v6_clean.tsv \
  --source-versions data/output/protein/protein_source_versions_v1.tsv \
  --output-isoform "$ISOFORM_TABLE" \
  --output-map "$MAP_TABLE" \
  --qa-isoform "$ISOFORM_QA" \
  --qa-map "$MAP_QA"

echo "[2/4] Validate contracts..."
python3 tools/kg_validate_table.py \
  --contract "$ISOFORM_CONTRACT" \
  --table "$ISOFORM_TABLE" \
  --out "$ISOFORM_VALIDATION"

python3 tools/kg_validate_table.py \
  --contract "$MAP_CONTRACT" \
  --table "$MAP_TABLE" \
  --out "$MAP_VALIDATION"

echo "[3/4] Generate manifests..."
python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$ISOFORM_MANIFEST" \
  "$ISOFORM_TABLE" \
  "$ISOFORM_CONTRACT" \
  "$ISOFORM_VALIDATION" \
  "$ISOFORM_QA"

python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$MAP_MANIFEST" \
  "$MAP_TABLE" \
  "$MAP_CONTRACT" \
  "$MAP_VALIDATION" \
  "$MAP_QA"

echo "[4/4] Done."
