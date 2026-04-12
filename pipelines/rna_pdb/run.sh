#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

INPUT_MASTER="${INPUT_MASTER:-data/output/rna_master_v1.tsv}"
INPUT_ID_MAPPING="${INPUT_ID_MAPPING:-data/raw/rna/rnacentral/id_mapping.tsv.gz}"
INPUT_XREF="${INPUT_XREF:-data/output/rna_xref_mrna_enst_urs_v2.tsv}"

OUTPUT_TABLE="data/output/rna_pdb_structures_v1.tsv"

PIPELINE_DIR="pipelines/rna_pdb"
CONTRACT="$PIPELINE_DIR/contracts/rna_pdb_structures_v1.json"

REPORT_BUILD="$PIPELINE_DIR/reports/rna_pdb_structures_v1.build.json"
REPORT_AUDIT="$PIPELINE_DIR/reports/rna_pdb_structures_v1.api_audit.json"
REPORT_VALIDATION="$PIPELINE_DIR/reports/rna_pdb_structures_v1.validation.json"
REPORT_QA="$PIPELINE_DIR/reports/rna_pdb_structures_v1.qa.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/rna_pdb_structures_v1.manifest.json"

SAMPLE_TABLE="$PIPELINE_DIR/.cache/rna_pdb_structures_v1.sample.tsv"
SAMPLE_BUILD="$PIPELINE_DIR/reports/rna_pdb_structures_v1.sample.build.json"
SAMPLE_AUDIT="$PIPELINE_DIR/reports/rna_pdb_structures_v1.sample.api_audit.json"
SAMPLE_VALIDATION="$PIPELINE_DIR/reports/rna_pdb_structures_v1.sample.validation.json"
SAMPLE_QA="$PIPELINE_DIR/reports/rna_pdb_structures_v1.sample.qa.json"

DATA_VERSION="${DATA_VERSION:-kg-data-local}"
SOURCE_VERSION="${SOURCE_VERSION:-RNAcentral:25;RCSB:live}"
SAMPLE_MASTER_ROWS="${SAMPLE_MASTER_ROWS:-0}"
SAMPLE_IDMAP_LINES="${SAMPLE_IDMAP_LINES:-0}"
SMOKE_MAX_UNIQUE_PDB="${SMOKE_MAX_UNIQUE_PDB:-200}"
RCSB_API_FAIL_THRESHOLD="${RCSB_API_FAIL_THRESHOLD:-0.05}"
RCSB_BATCH_SIZE="${RCSB_BATCH_SIZE:-500}"
RCSB_TIMEOUT="${RCSB_TIMEOUT:-40}"
RCSB_RETRIES="${RCSB_RETRIES:-3}"
RCSB_SLEEP_BETWEEN_BATCHES="${RCSB_SLEEP_BETWEEN_BATCHES:-0.05}"

if [[ ! -f "$INPUT_MASTER" ]]; then
  echo "[ERROR] missing input master: $INPUT_MASTER" >&2
  exit 2
fi

if [[ ! -f "$INPUT_ID_MAPPING" ]]; then
  if [[ "$INPUT_ID_MAPPING" == *.gz && -f "${INPUT_ID_MAPPING%.gz}" ]]; then
    INPUT_ID_MAPPING="${INPUT_ID_MAPPING%.gz}"
  elif [[ "$INPUT_ID_MAPPING" != *.gz && -f "${INPUT_ID_MAPPING}.gz" ]]; then
    INPUT_ID_MAPPING="${INPUT_ID_MAPPING}.gz"
  else
    echo "[ERROR] missing id_mapping (.tsv/.gz): $INPUT_ID_MAPPING" >&2
    exit 2
  fi
fi

mkdir -p "$(dirname "$OUTPUT_TABLE")" "$PIPELINE_DIR/reports" "$PIPELINE_DIR/.cache"

XREF_ARGS=()
if [[ -f "$INPUT_XREF" ]]; then
  XREF_ARGS+=(--xref "$INPUT_XREF")
else
  echo "[WARN] xref file not found, continue without xref: $INPUT_XREF"
fi

echo "[STEP 1/5] sample run"
SAMPLE_ARGS=()
if [[ "$SAMPLE_MASTER_ROWS" =~ ^[0-9]+$ ]] && [[ "$SAMPLE_MASTER_ROWS" -gt 0 ]]; then
  SAMPLE_ARGS+=(--max-master-rows "$SAMPLE_MASTER_ROWS")
fi
if [[ "$SAMPLE_IDMAP_LINES" =~ ^[0-9]+$ ]] && [[ "$SAMPLE_IDMAP_LINES" -gt 0 ]]; then
  SAMPLE_ARGS+=(--max-idmap-lines "$SAMPLE_IDMAP_LINES")
fi

sample_cmd=(
  python3 "$PIPELINE_DIR/scripts/01_build_rna_pdb_structures_v1.py"
  --rna-master "$INPUT_MASTER"
  --id-mapping "$INPUT_ID_MAPPING"
  --output "$SAMPLE_TABLE"
  --report-build "$SAMPLE_BUILD"
  --report-audit "$SAMPLE_AUDIT"
  --source-version "$SOURCE_VERSION"
  --max-unique-pdb "$SMOKE_MAX_UNIQUE_PDB"
  --batch-size "$RCSB_BATCH_SIZE"
  --timeout "$RCSB_TIMEOUT"
  --retries "$RCSB_RETRIES"
  --sleep-between-batches "$RCSB_SLEEP_BETWEEN_BATCHES"
  --fail-threshold "$RCSB_API_FAIL_THRESHOLD"
)
if (( ${#SAMPLE_ARGS[@]} > 0 )); then
  sample_cmd+=("${SAMPLE_ARGS[@]}")
fi
if (( ${#XREF_ARGS[@]} > 0 )); then
  sample_cmd+=("${XREF_ARGS[@]}")
fi
"${sample_cmd[@]}"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$SAMPLE_TABLE" \
  --out "$SAMPLE_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_rna_pdb_structures_v1.py" \
  --master "$INPUT_MASTER" \
  --table "$SAMPLE_TABLE" \
  --build-report "$SAMPLE_BUILD" \
  --out "$SAMPLE_QA"

echo "[STEP 2/5] full run"
full_cmd=(
  python3 "$PIPELINE_DIR/scripts/01_build_rna_pdb_structures_v1.py"
  --rna-master "$INPUT_MASTER"
  --id-mapping "$INPUT_ID_MAPPING"
  --output "$OUTPUT_TABLE"
  --report-build "$REPORT_BUILD"
  --report-audit "$REPORT_AUDIT"
  --source-version "$SOURCE_VERSION"
  --batch-size "$RCSB_BATCH_SIZE"
  --timeout "$RCSB_TIMEOUT"
  --retries "$RCSB_RETRIES"
  --sleep-between-batches "$RCSB_SLEEP_BETWEEN_BATCHES"
  --fail-threshold "$RCSB_API_FAIL_THRESHOLD"
)
if (( ${#XREF_ARGS[@]} > 0 )); then
  full_cmd+=("${XREF_ARGS[@]}")
fi
"${full_cmd[@]}"

echo "[STEP 3/5] contract validation"
python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUTPUT_TABLE" \
  --out "$REPORT_VALIDATION"

echo "[STEP 4/5] QA checks"
python3 "$PIPELINE_DIR/scripts/02_qa_rna_pdb_structures_v1.py" \
  --master "$INPUT_MASTER" \
  --table "$OUTPUT_TABLE" \
  --build-report "$REPORT_BUILD" \
  --out "$REPORT_QA"

echo "[STEP 5/5] manifest"
python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$OUTPUT_TABLE" \
  "$INPUT_MASTER" \
  "$INPUT_ID_MAPPING" \
  "$CONTRACT" \
  "$REPORT_BUILD" \
  "$REPORT_VALIDATION" \
  "$REPORT_QA" \
  "$REPORT_AUDIT"

echo "[DONE] rna_pdb v1: sample + full + validation + QA + manifest finished."
