#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

mkdir -p data/output pipelines/rna_type_features/reports

echo "[0/8] Preflight check..."
python3 pipelines/rna_type_features/scripts/build_rna_type_features.py \
  --check-inputs \
  --report pipelines/rna_type_features/reports/rna_type_features_v1.blocked_or_ready.json

echo "[1/8] Build type feature pack..."
python3 pipelines/rna_type_features/scripts/build_rna_type_features.py \
  --lnc-output data/output/rna_lnc_entries_v1.tsv \
  --trna-output data/output/rna_trna_features_v1.tsv \
  --rrna-output data/output/rna_rrna_loci_v1.tsv \
  --report pipelines/rna_type_features/reports/rna_type_features_v1.metrics.json

echo "[2/8] Validate lnc table..."
python3 tools/kg_validate_table.py \
  --contract pipelines/rna_type_features/contracts/rna_lnc_entries_v1.json \
  --table data/output/rna_lnc_entries_v1.tsv \
  --out pipelines/rna_type_features/reports/rna_lnc_entries_v1.validation.json

echo "[3/8] Validate tRNA table..."
python3 tools/kg_validate_table.py \
  --contract pipelines/rna_type_features/contracts/rna_trna_features_v1.json \
  --table data/output/rna_trna_features_v1.tsv \
  --out pipelines/rna_type_features/reports/rna_trna_features_v1.validation.json

echo "[4/8] Validate rRNA loci table..."
python3 tools/kg_validate_table.py \
  --contract pipelines/rna_type_features/contracts/rna_rrna_loci_v1.json \
  --table data/output/rna_rrna_loci_v1.tsv \
  --out pipelines/rna_type_features/reports/rna_rrna_loci_v1.validation.json

echo "[5/8] Manifest for lnc..."
python3 tools/kg_make_manifest.py \
  --data-version kg-rna-lnc-entries-v1 \
  --out pipelines/rna_type_features/reports/rna_lnc_entries_v1.manifest.json \
  data/output/rna_lnc_entries_v1.tsv

echo "[6/8] Manifest for tRNA..."
python3 tools/kg_make_manifest.py \
  --data-version kg-rna-trna-features-v1 \
  --out pipelines/rna_type_features/reports/rna_trna_features_v1.manifest.json \
  data/output/rna_trna_features_v1.tsv

echo "[7/8] Manifest for rRNA loci..."
python3 tools/kg_make_manifest.py \
  --data-version kg-rna-rrna-loci-v1 \
  --out pipelines/rna_type_features/reports/rna_rrna_loci_v1.manifest.json \
  data/output/rna_rrna_loci_v1.tsv

echo "[8/8] Done."
