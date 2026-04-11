#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

mkdir -p data/output/protein pipelines/protein_domains/reports

echo "[1/4] Dry run on minimal sample (n=30)..."
python3 pipelines/protein_domains/scripts/build_protein_domains_interpro.py \
  --input data/processed/protein_master_v6_clean.tsv \
  --output data/output/protein/protein_domains_interpro_v1.sample.tsv \
  --report pipelines/protein_domains/reports/protein_domains_interpro_v1.sample.metrics.json \
  --audit pipelines/protein_domains/reports/protein_domains_interpro_v1.sample.audit.json \
  --limit 30 \
  --workers 6

echo "[2/4] Full run on all proteins..."
python3 pipelines/protein_domains/scripts/build_protein_domains_interpro.py \
  --input data/processed/protein_master_v6_clean.tsv \
  --output data/output/protein/protein_domains_interpro_v1.tsv \
  --report pipelines/protein_domains/reports/protein_domains_interpro_v1.metrics.json \
  --audit pipelines/protein_domains/reports/protein_domains_interpro_v1.audit.json \
  --workers 8

echo "[3/4] Contract validation..."
python3 tools/kg_validate_table.py \
  --contract pipelines/protein_domains/contracts/protein_domains_interpro_v1.json \
  --table data/output/protein/protein_domains_interpro_v1.tsv \
  --out pipelines/protein_domains/reports/protein_domains_interpro_v1.validation.json

echo "[4/4] Done."
