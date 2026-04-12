#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

mkdir -p data/output pipelines/rna_type_features/reports

echo "[0/5] Preflight check..."
python3 pipelines/rna_type_features/scripts/build_rna_trna_features_v2.py \
  --check-inputs \
  --report pipelines/rna_type_features/reports/rna_trna_features_v2.blocked_or_ready.json

echo "[1/5] Build tRNA v2 (balanced)..."
python3 pipelines/rna_type_features/scripts/build_rna_trna_features_v2.py \
  --output data/output/rna_trna_features_v2.tsv \
  --conflict-audit pipelines/rna_type_features/reports/rna_trna_anticodon_conflicts_v2.tsv \
  --report pipelines/rna_type_features/reports/rna_trna_features_v2.metrics.json

echo "[2/5] Validate tRNA v2 table..."
python3 tools/kg_validate_table.py \
  --contract pipelines/rna_type_features/contracts/rna_trna_features_v2.json \
  --table data/output/rna_trna_features_v2.tsv \
  --out pipelines/rna_type_features/reports/rna_trna_features_v2.validation.json

echo "[3/5] Manifest..."
python3 tools/kg_make_manifest.py \
  --data-version kg-rna-trna-features-v2 \
  --out pipelines/rna_type_features/reports/rna_trna_features_v2.manifest.json \
  data/output/rna_trna_features_v2.tsv

echo "[4/5] Compare v1/v2 quick metrics..."
python3 - <<'PY'
import csv

def rate(path):
    rows = 0
    non = 0
    with open(path) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            rows += 1
            if (row.get('anticodon') or '').strip() != '':
                non += 1
    return rows, non, (non / rows if rows else 0.0)

r1 = rate('data/output/rna_trna_features_v1.tsv')
r2 = rate('data/output/rna_trna_features_v2.tsv')
print(f"v1 rows={r1[0]} anticodon_non_empty={r1[1]} rate={r1[2]:.6f}")
print(f"v2 rows={r2[0]} anticodon_non_empty={r2[1]} rate={r2[2]:.6f}")
PY

echo "[5/5] Done."
