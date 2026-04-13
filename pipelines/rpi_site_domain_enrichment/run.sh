#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

mkdir -p data/output/evidence pipelines/rpi_site_domain_enrichment/reports

echo "[0/8] Preflight check (inputs/download blockers)..."
python3 pipelines/rpi_site_domain_enrichment/scripts/build_rpi_site_domain_enrichment.py \
  --check-inputs \
  --report pipelines/rpi_site_domain_enrichment/reports/rpi_site_domain_enrichment_v2.blocked_or_ready.json

echo "[1/8] Dry run on minimal sample (n=200)..."
python3 pipelines/rpi_site_domain_enrichment/scripts/build_rpi_site_domain_enrichment.py \
  --limit 200 \
  --site-output data/output/evidence/rpi_site_context_v2.sample.tsv \
  --domain-output data/output/evidence/rpi_domain_context_v2.sample.tsv \
  --function-output data/output/evidence/rpi_function_context_v2.sample.tsv \
  --report pipelines/rpi_site_domain_enrichment/reports/rpi_site_domain_enrichment_v2.sample.metrics.json

echo "[2/8] Full run..."
python3 pipelines/rpi_site_domain_enrichment/scripts/build_rpi_site_domain_enrichment.py \
  --site-output data/output/evidence/rpi_site_context_v2.tsv \
  --domain-output data/output/evidence/rpi_domain_context_v2.tsv \
  --function-output data/output/evidence/rpi_function_context_v2.tsv \
  --report pipelines/rpi_site_domain_enrichment/reports/rpi_site_domain_enrichment_v2.metrics.json

echo "[3/8] Validate site context table..."
python3 tools/kg_validate_table.py \
  --contract pipelines/rpi_site_domain_enrichment/contracts/rpi_site_context_v2.json \
  --table data/output/evidence/rpi_site_context_v2.tsv \
  --out pipelines/rpi_site_domain_enrichment/reports/rpi_site_context_v2.validation.json

echo "[4/8] Validate domain context table..."
python3 tools/kg_validate_table.py \
  --contract pipelines/rpi_site_domain_enrichment/contracts/rpi_domain_context_v2.json \
  --table data/output/evidence/rpi_domain_context_v2.tsv \
  --out pipelines/rpi_site_domain_enrichment/reports/rpi_domain_context_v2.validation.json

echo "[5/8] Validate function context table..."
python3 tools/kg_validate_table.py \
  --contract pipelines/rpi_site_domain_enrichment/contracts/rpi_function_context_v2.json \
  --table data/output/evidence/rpi_function_context_v2.tsv \
  --out pipelines/rpi_site_domain_enrichment/reports/rpi_function_context_v2.validation.json

echo "[6/8] Build manifest..."
python3 tools/kg_make_manifest.py \
  --data-version kg-rpi-site-domain-enrichment-v2 \
  --out pipelines/rpi_site_domain_enrichment/reports/rpi_site_domain_enrichment_v2.manifest.json \
  data/output/evidence/rpi_site_context_v2.tsv \
  data/output/evidence/rpi_domain_context_v2.tsv \
  data/output/evidence/rpi_function_context_v2.tsv

echo "[7/8] Aggregate final gates (metrics + contracts)..."
python3 pipelines/rpi_site_domain_enrichment/scripts/make_gates_report.py \
  --metrics pipelines/rpi_site_domain_enrichment/reports/rpi_site_domain_enrichment_v2.metrics.json \
  --validations \
    pipelines/rpi_site_domain_enrichment/reports/rpi_site_context_v2.validation.json \
    pipelines/rpi_site_domain_enrichment/reports/rpi_domain_context_v2.validation.json \
    pipelines/rpi_site_domain_enrichment/reports/rpi_function_context_v2.validation.json \
  --out pipelines/rpi_site_domain_enrichment/reports/rpi_site_domain_enrichment_v2.gates.json

echo "[8/8] Done."
