#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

mkdir -p data/output pipelines/rna_expression_rbp/reports

echo "[0/7] Prepare ENCODE minimal inputs (expression + RBP)..."
python3 pipelines/rna_expression_rbp/scripts/prepare_encode_minimal_inputs.py

echo "[1/7] Preflight check (inputs/download blockers)..."
python3 pipelines/rna_expression_rbp/scripts/build_rna_expression_rbp.py \
  --check-inputs \
  --report pipelines/rna_expression_rbp/reports/rna_expression_rbp_v1.blocked_or_ready.json

echo "[2/7] Dry run on minimal sample (n=200)..."
python3 pipelines/rna_expression_rbp/scripts/build_rna_expression_rbp.py \
  --limit 200 \
  --expression-output data/output/rna_expression_evidence_v1.sample.tsv \
  --rbp-output data/output/rna_rbp_sites_v1.sample.tsv \
  --report pipelines/rna_expression_rbp/reports/rna_expression_rbp_v1.sample.metrics.json

echo "[3/7] Full run..."
python3 pipelines/rna_expression_rbp/scripts/build_rna_expression_rbp.py \
  --expression-output data/output/rna_expression_evidence_v1.tsv \
  --rbp-output data/output/rna_rbp_sites_v1.tsv \
  --report pipelines/rna_expression_rbp/reports/rna_expression_rbp_v1.metrics.json

echo "[4/7] Validate expression table..."
python3 tools/kg_validate_table.py \
  --contract pipelines/rna_expression_rbp/contracts/rna_expression_evidence_v1.json \
  --table data/output/rna_expression_evidence_v1.tsv \
  --out pipelines/rna_expression_rbp/reports/rna_expression_evidence_v1.validation.json

echo "[5/7] Validate RBP table..."
python3 tools/kg_validate_table.py \
  --contract pipelines/rna_expression_rbp/contracts/rna_rbp_sites_v1.json \
  --table data/output/rna_rbp_sites_v1.tsv \
  --out pipelines/rna_expression_rbp/reports/rna_rbp_sites_v1.validation.json

echo "[6/7] Build manifest..."
python3 tools/kg_make_manifest.py \
  --data-version kg-rna-expression-rbp-v1 \
  --out pipelines/rna_expression_rbp/reports/rna_expression_rbp_v1.manifest.json \
  data/output/rna_expression_evidence_v1.tsv \
  data/output/rna_rbp_sites_v1.tsv

echo "[7/7] Done."
