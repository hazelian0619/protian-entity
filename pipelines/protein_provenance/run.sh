#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

INPUT_MASTER="data/processed/protein_master_v6_clean.tsv"
INPUT_EDGES="data/processed/protein_edges.tsv"
INPUT_PTM="data/processed/ptm_sites.tsv"
INPUT_PATHWAY="data/processed/pathway_members.tsv"
INPUT_AF="data/processed/alphafold_quality.tsv"
INPUT_HGNC="data/processed/hgnc_core.tsv"

OUT_TABLE="data/output/protein/protein_source_versions_v1.tsv"
CONTRACT="pipelines/protein_provenance/contracts/protein_source_versions_v1.json"

REPORT_BUILD="pipelines/protein_provenance/reports/protein_source_versions_v1.build.json"
REPORT_VALIDATION="pipelines/protein_provenance/reports/protein_source_versions_v1.validation.json"
REPORT_QA="pipelines/protein_provenance/reports/protein_source_versions_v1.qa.json"
REPORT_MANIFEST="pipelines/protein_provenance/reports/protein_source_versions_v1.manifest.json"

SMOKE_ROWS="${SMOKE_ROWS:-2000}"
SMOKE_OUT="pipelines/protein_provenance/reports/protein_source_versions_v1.smoke.tsv"
SMOKE_BUILD="pipelines/protein_provenance/reports/protein_source_versions_v1.smoke.build.json"
SMOKE_VALIDATION="pipelines/protein_provenance/reports/protein_source_versions_v1.smoke.validation.json"
SMOKE_QA="pipelines/protein_provenance/reports/protein_source_versions_v1.smoke.qa.json"

DATA_VERSION="${DATA_VERSION:-kg-data-local}"

inputs=(
  "$INPUT_MASTER"
  "$INPUT_EDGES"
  "$INPUT_PTM"
  "$INPUT_PATHWAY"
  "$INPUT_AF"
  "$INPUT_HGNC"
)

missing=()
for f in "${inputs[@]}"; do
  [[ -f "$f" ]] || missing+=("$f")
done

if (( ${#missing[@]} > 0 )); then
  echo "[ERROR] missing required local inputs for protein provenance pipeline:" >&2
  for f in "${missing[@]}"; do
    echo "  - $f" >&2
  done
  echo >&2
  echo "下载清单（无需外部 API，使用本仓库标准落盘路径）：" >&2
  echo "  1) protein_master_v6_clean.tsv  -> data/processed/" >&2
  echo "  2) protein_edges.tsv            -> data/processed/" >&2
  echo "  3) ptm_sites.tsv                -> data/processed/" >&2
  echo "  4) pathway_members.tsv          -> data/processed/" >&2
  echo "  5) alphafold_quality.tsv        -> data/processed/" >&2
  echo "  6) hgnc_core.tsv                -> data/processed/" >&2
  echo "校验方式：" >&2
  echo "  - 列头检查: head -n 1 data/processed/<file>.tsv" >&2
  echo "  - 文件完整性: wc -l data/processed/<file>.tsv" >&2
  exit 2
fi

mkdir -p "$(dirname "$OUT_TABLE")" "pipelines/protein_provenance/reports"

echo "[STEP] smoke run (max rows per input = $SMOKE_ROWS)"
python3 pipelines/protein_provenance/scripts/01_build_protein_source_versions_v1.py \
  --repo-root "$REPO_ROOT" \
  --out "$SMOKE_OUT" \
  --out-build-report "$SMOKE_BUILD" \
  --max-rows-per-input "$SMOKE_ROWS"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$SMOKE_OUT" \
  --out "$SMOKE_VALIDATION"

python3 pipelines/protein_provenance/scripts/02_qa_protein_source_versions_v1.py \
  --table "$SMOKE_OUT" \
  --out "$SMOKE_QA"

echo "[STEP] full run"
python3 pipelines/protein_provenance/scripts/01_build_protein_source_versions_v1.py \
  --repo-root "$REPO_ROOT" \
  --out "$OUT_TABLE" \
  --out-build-report "$REPORT_BUILD"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUT_TABLE" \
  --out "$REPORT_VALIDATION"

python3 pipelines/protein_provenance/scripts/02_qa_protein_source_versions_v1.py \
  --table "$OUT_TABLE" \
  --out "$REPORT_QA"

python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$OUT_TABLE" \
  "$CONTRACT" \
  "$REPORT_BUILD" \
  "$REPORT_VALIDATION" \
  "$REPORT_QA" \
  "$INPUT_MASTER" \
  "$INPUT_EDGES" \
  "$INPUT_PTM" \
  "$INPUT_PATHWAY" \
  "$INPUT_AF" \
  "$INPUT_HGNC"

echo "[DONE] protein_provenance v1 finished."
