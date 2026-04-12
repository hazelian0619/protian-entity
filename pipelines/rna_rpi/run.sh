#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

mkdir -p data/output/edges data/output/evidence pipelines/rna_rpi/reports
RPI_SOURCE_VERSION="${RPI_SOURCE_VERSION:-ENCORI_API_RBPTarget_hg38_mRNA_clipExpNum_ge_1_target_all_cellType_all}"

echo "[0/6] Preflight check (inputs/download blockers)..."
python3 pipelines/rna_rpi/scripts/build_rna_rpi.py \
  --check-inputs \
  --report pipelines/rna_rpi/reports/rna_rpi_v1.blocked_or_ready.json

echo "[1/6] Dry run on minimal sample (n=200)..."
python3 pipelines/rna_rpi/scripts/build_rna_rpi.py \
  --limit 200 \
  --source-version "$RPI_SOURCE_VERSION" \
  --edges-output data/output/edges/rna_protein_edges_v1.sample.tsv \
  --evidence-output data/output/evidence/rna_protein_evidence_v1.sample.tsv \
  --report pipelines/rna_rpi/reports/rna_rpi_v1.sample.metrics.json

echo "[2/6] Full run..."
python3 pipelines/rna_rpi/scripts/build_rna_rpi.py \
  --source-version "$RPI_SOURCE_VERSION" \
  --edges-output data/output/edges/rna_protein_edges_v1.tsv \
  --evidence-output data/output/evidence/rna_protein_evidence_v1.tsv \
  --report pipelines/rna_rpi/reports/rna_rpi_v1.metrics.json \
  --gates-report pipelines/rna_rpi/reports/rna_rpi_v1.gates.json

echo "[3/6] Validate edges table..."
python3 tools/kg_validate_table.py \
  --contract pipelines/rna_rpi/contracts/rna_protein_edges_v1.json \
  --table data/output/edges/rna_protein_edges_v1.tsv \
  --out pipelines/rna_rpi/reports/rna_protein_edges_v1.validation.json

echo "[4/6] Validate evidence table..."
python3 tools/kg_validate_table.py \
  --contract pipelines/rna_rpi/contracts/rna_protein_evidence_v1.json \
  --table data/output/evidence/rna_protein_evidence_v1.tsv \
  --out pipelines/rna_rpi/reports/rna_protein_evidence_v1.validation.json

echo "[5/6] Build manifest..."
python3 tools/kg_make_manifest.py \
  --data-version kg-rna-rpi-v1 \
  --out pipelines/rna_rpi/reports/rna_rpi_v1.manifest.json \
  data/output/edges/rna_protein_edges_v1.tsv \
  data/output/evidence/rna_protein_evidence_v1.tsv

echo "[6/6] Done."
