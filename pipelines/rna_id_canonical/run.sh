#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

INPUT_MASTER="data/output/rna_master_v1.tsv"
INPUT_XREF="data/output/rna_xref_mrna_enst_urs_v2.tsv"
INPUT_GENE_XREF="data/output/gene_xref_rna_v1.tsv"
REQUIRE_GENE_XREF="${REQUIRE_GENE_XREF:-0}"

OUT_CANONICAL="data/output/rna_id_canonical_map_v1.tsv"
CONTRACT="pipelines/rna_id_canonical/contracts/rna_id_canonical_map_v1.json"

REPORT_BUILD="pipelines/rna_id_canonical/reports/rna_id_canonical_map_v1.build.json"
REPORT_VALIDATION="pipelines/rna_id_canonical/reports/rna_id_canonical_map_v1.validation.json"
REPORT_JOIN_AUDIT="pipelines/rna_id_canonical/reports/rna_id_canonical_map_v1.gene_xref_join_audit.json"

REPORT_SAMPLE_BUILD="pipelines/rna_id_canonical/reports/rna_id_canonical_map_v1.sample.build.json"
REPORT_SAMPLE_JOIN_AUDIT="pipelines/rna_id_canonical/reports/rna_id_canonical_map_v1.sample.gene_xref_join_audit.json"

SAMPLE_ROWS="${SAMPLE_ROWS:-5000}"
SAMPLE_GENE_ROWS="${SAMPLE_GENE_ROWS:-5000}"

if [[ ! -f "$INPUT_MASTER" ]]; then
  echo "[ERROR] missing input: $INPUT_MASTER" >&2
  exit 2
fi
if [[ ! -f "$INPUT_XREF" ]]; then
  echo "[ERROR] missing input: $INPUT_XREF" >&2
  echo "Please run assistant A pipeline first: bash pipelines/rna_xref_enst_urs/run.sh" >&2
  exit 2
fi
mkdir -p "$(dirname "$OUT_CANONICAL")" "pipelines/rna_id_canonical/reports"

# 1) sample-first run
SAMPLE_OUT="$(mktemp -t rna_id_canonical_map_v1.sample.XXXXXX.tsv)"
trap 'rm -f "$SAMPLE_OUT"' EXIT

python3 pipelines/rna_id_canonical/scripts/01_build_rna_id_canonical_map_v1.py \
  --master "$INPUT_MASTER" \
  --xref "$INPUT_XREF" \
  --out "$SAMPLE_OUT" \
  --report "$REPORT_SAMPLE_BUILD" \
  --max-rows "$SAMPLE_ROWS"

if [[ -f "$INPUT_GENE_XREF" ]]; then
  python3 pipelines/rna_id_canonical/scripts/02_audit_gene_xref_join_v1.py \
    --canonical-map "$SAMPLE_OUT" \
    --gene-xref-rna "$INPUT_GENE_XREF" \
    --out "$REPORT_SAMPLE_JOIN_AUDIT" \
    --max-rows "$SAMPLE_GENE_ROWS"
else
  echo "[WARN] missing optional input: $INPUT_GENE_XREF (skip sample join audit)"
fi

# 2) full run
python3 pipelines/rna_id_canonical/scripts/01_build_rna_id_canonical_map_v1.py \
  --master "$INPUT_MASTER" \
  --xref "$INPUT_XREF" \
  --out "$OUT_CANONICAL" \
  --report "$REPORT_BUILD"

# 3) contract validation
python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUT_CANONICAL" \
  --out "$REPORT_VALIDATION"

# 4) downstream join audit
if [[ -f "$INPUT_GENE_XREF" ]]; then
  python3 pipelines/rna_id_canonical/scripts/02_audit_gene_xref_join_v1.py \
    --canonical-map "$OUT_CANONICAL" \
    --gene-xref-rna "$INPUT_GENE_XREF" \
    --out "$REPORT_JOIN_AUDIT"
else
  echo "[WARN] missing optional input: $INPUT_GENE_XREF (skip full join audit)"
fi

# 5) gates
python3 - <<'PY'
import json
from pathlib import Path

build = json.loads(Path("pipelines/rna_id_canonical/reports/rna_id_canonical_map_v1.build.json").read_text(encoding="utf-8"))
join_path = Path("pipelines/rna_id_canonical/reports/rna_id_canonical_map_v1.gene_xref_join_audit.json")
import os
require_gene = os.environ.get("REQUIRE_GENE_XREF", "0") == "1"

ok_legacy = bool(build.get("gates", {}).get("legacy_coverage_100"))
ok_canonical = bool(build.get("gates", {}).get("canonical_non_empty_100"))
if join_path.exists():
    join = json.loads(join_path.read_text(encoding="utf-8"))
    ok_join = bool(join.get("gates", {}).get("legacy_join_success_100"))
else:
    ok_join = not require_gene

if ok_legacy and ok_canonical and ok_join:
    print("[PASS] canonical map gates passed.")
    raise SystemExit(0)

print("[FAIL] canonical map gates failed.")
print(f"  legacy_coverage_100={ok_legacy}")
print(f"  canonical_non_empty_100={ok_canonical}")
print(f"  legacy_join_success_100={ok_join} (REQUIRE_GENE_XREF={'1' if require_gene else '0'})")
raise SystemExit(3)
PY

echo "[DONE] rna_id_canonical_map_v1 pipeline completed."
