#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

mkdir -p data/output/evidence pipelines/interaction_ontology_mapping/reports

echo "[0/6] Preflight check (inputs/download blockers)..."
python3 pipelines/interaction_ontology_mapping/scripts/build_interaction_ontology_mapping.py \
  --check-inputs \
  --report pipelines/interaction_ontology_mapping/reports/interaction_ontology_mapping_v2.blocked_or_ready.json

echo "[1/6] Dry run on minimal sample (limit_per_type=2000)..."
python3 pipelines/interaction_ontology_mapping/scripts/build_interaction_ontology_mapping.py \
  --limit-per-type 2000 \
  --output data/output/evidence/interaction_ontology_mapping_v2.sample.tsv \
  --report pipelines/interaction_ontology_mapping/reports/interaction_ontology_mapping_v2.sample.metrics.json

echo "[2/6] Full run..."
python3 pipelines/interaction_ontology_mapping/scripts/build_interaction_ontology_mapping.py \
  --output data/output/evidence/interaction_ontology_mapping_v2.tsv \
  --report pipelines/interaction_ontology_mapping/reports/interaction_ontology_mapping_v2.metrics.json

echo "[3/6] Validate ontology-mapping table..."
python3 tools/kg_validate_table.py \
  --contract pipelines/interaction_ontology_mapping/contracts/interaction_ontology_mapping_v2.json \
  --table data/output/evidence/interaction_ontology_mapping_v2.tsv \
  --out pipelines/interaction_ontology_mapping/reports/interaction_ontology_mapping_v2.validation.json

echo "[4/6] Build manifest..."
python3 tools/kg_make_manifest.py \
  --data-version kg-interaction-ontology-mapping-v2 \
  --out pipelines/interaction_ontology_mapping/reports/interaction_ontology_mapping_v2.manifest.json \
  data/output/evidence/interaction_ontology_mapping_v2.tsv

echo "[5/6] Aggregate final gates (metrics + contract)..."
python3 pipelines/interaction_ontology_mapping/scripts/make_gates_report.py \
  --metrics pipelines/interaction_ontology_mapping/reports/interaction_ontology_mapping_v2.metrics.json \
  --validation pipelines/interaction_ontology_mapping/reports/interaction_ontology_mapping_v2.validation.json \
  --out pipelines/interaction_ontology_mapping/reports/interaction_ontology_mapping_v2.gates.json

echo "[6/6] Done."
