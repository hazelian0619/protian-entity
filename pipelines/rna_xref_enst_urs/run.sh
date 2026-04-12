#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

OUT_XREF="data/output/rna_xref_mrna_enst_urs_v2.tsv"
CONTRACT="pipelines/rna_xref_enst_urs/contracts/rna_xref_mrna_enst_urs_v2.json"
MASTER_PATH="${MASTER_PATH:-data/output/rna_master_v1.tsv}"
ID_MAPPING_PATH="${ID_MAPPING_PATH:-data/raw/rna/rnacentral/id_mapping.tsv.gz}"

REPORT_COVERAGE="pipelines/rna_xref_enst_urs/reports/mrna_enst_urs_coverage_v2.json"
REPORT_CONFLICTS="pipelines/rna_xref_enst_urs/reports/mrna_enst_urs_conflicts_v2.json"
REPORT_MANUAL="pipelines/rna_xref_enst_urs/reports/mrna_enst_urs_manual_download_v2.json"
REPORT_VALIDATION="pipelines/rna_xref_enst_urs/reports/rna_xref_mrna_enst_urs_v2.validation.json"

REPORT_SAMPLE_COVERAGE="pipelines/rna_xref_enst_urs/reports/mrna_enst_urs_coverage_v2.sample.json"
REPORT_SAMPLE_CONFLICTS="pipelines/rna_xref_enst_urs/reports/mrna_enst_urs_conflicts_v2.sample.json"
REPORT_SAMPLE_MANUAL="pipelines/rna_xref_enst_urs/reports/mrna_enst_urs_manual_download_v2.sample.json"

SAMPLE_MASTER_ROWS="${SAMPLE_MASTER_ROWS:-5000}"
SAMPLE_IDMAP_LINES="${SAMPLE_IDMAP_LINES:-2000000}"
MIN_COVERAGE="${MIN_COVERAGE:-0.70}"
TARGET_COVERAGE="${TARGET_COVERAGE:-0.90}"
STRICT_COVERAGE_GATE="${STRICT_COVERAGE_GATE:-0}"

mkdir -p "$(dirname "$OUT_XREF")" "pipelines/rna_xref_enst_urs/reports"

# 1) minimal sample-first run
SAMPLE_XREF="$(mktemp -t rna_xref_mrna_enst_urs_v2.sample.XXXXXX.tsv)"
trap 'rm -f "$SAMPLE_XREF"' EXIT

python3 pipelines/rna_xref_enst_urs/scripts/01_build_rna_xref_mrna_enst_urs_v2.py \
  --master "$MASTER_PATH" \
  --id-mapping "$ID_MAPPING_PATH" \
  --out-xref "$SAMPLE_XREF" \
  --coverage-report "$REPORT_SAMPLE_COVERAGE" \
  --conflicts-report "$REPORT_SAMPLE_CONFLICTS" \
  --manual-download-report "$REPORT_SAMPLE_MANUAL" \
  --min-coverage "$MIN_COVERAGE" \
  --target-coverage "$TARGET_COVERAGE" \
  --max-master-rows "$SAMPLE_MASTER_ROWS" \
  --max-idmap-lines "$SAMPLE_IDMAP_LINES"

# 2) full run
python3 pipelines/rna_xref_enst_urs/scripts/01_build_rna_xref_mrna_enst_urs_v2.py \
  --master "$MASTER_PATH" \
  --id-mapping "$ID_MAPPING_PATH" \
  --out-xref "$OUT_XREF" \
  --coverage-report "$REPORT_COVERAGE" \
  --conflicts-report "$REPORT_CONFLICTS" \
  --manual-download-report "$REPORT_MANUAL" \
  --min-coverage "$MIN_COVERAGE" \
  --target-coverage "$TARGET_COVERAGE"

# 3) contract validation
python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUT_XREF" \
  --out "$REPORT_VALIDATION"

# 4) coverage gate (scheme-B: non-blocking KPI by default)
python3 - <<'PY'
import json
import os
from pathlib import Path

coverage = json.loads(Path("pipelines/rna_xref_enst_urs/reports/mrna_enst_urs_coverage_v2.json").read_text(encoding="utf-8"))
met_min = bool(coverage.get("thresholds", {}).get("met_minimum"))
strict_gate = os.environ.get("STRICT_COVERAGE_GATE", "0") == "1"
if met_min:
    print("[PASS] coverage gate met.")
    raise SystemExit(0)

print("[WARN] coverage KPI not met (non-blocking by default).")
manual = coverage.get("manual_download_report")
if manual:
    print(f"[WARN] manual checklist: {manual}")

if strict_gate:
    print("[BLOCKED] STRICT_COVERAGE_GATE=1, stop pipeline.")
    raise SystemExit(3)

print("[CONTINUE] STRICT_COVERAGE_GATE!=1, continue.")
raise SystemExit(0)
PY

echo "[DONE] rna_xref_mrna_enst_urs_v2 pipeline completed."
