#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

mkdir -p data/output pipelines/rna_rfam_structure/reports

echo "[0/4] Preflight check (inputs/download blockers)..."
python3 pipelines/rna_rfam_structure/scripts/build_rna_rfam_structure.py \
  --check-inputs \
  --report pipelines/rna_rfam_structure/reports/rna_structure_rfam_v1.blocked_or_ready.json

echo "[1/4] Dry run on minimal sample (n=200)..."
python3 pipelines/rna_rfam_structure/scripts/build_rna_rfam_structure.py \
  --limit 200 \
  --output data/output/rna_structure_rfam_v1.sample.tsv \
  --report pipelines/rna_rfam_structure/reports/rna_structure_rfam_v1.sample.metrics.json || true

echo "[2/4] Full run..."
python3 pipelines/rna_rfam_structure/scripts/build_rna_rfam_structure.py \
  --output data/output/rna_structure_rfam_v1.tsv \
  --report pipelines/rna_rfam_structure/reports/rna_structure_rfam_v1.metrics.json

echo "[3/4] Contract validation..."
python3 tools/kg_validate_table.py \
  --contract pipelines/rna_rfam_structure/contracts/rna_structure_rfam_v1.json \
  --table data/output/rna_structure_rfam_v1.tsv \
  --out pipelines/rna_rfam_structure/reports/rna_structure_rfam_v1.validation.json

echo "[4/4] Done."
