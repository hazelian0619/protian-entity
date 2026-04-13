#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

INPUT_EDGES="pipelines/edges_ppi/data/output/edges/edges_ppi_v1.tsv"
INPUT_MASTER="data/processed/protein_master_v6_clean.tsv"
INPUT_REACTOME="data/processed/pathway_members.tsv"
INPUT_KEGG="data/output/protein/protein_kegg_pathway_v1.tsv"

OUT_METHOD="data/output/evidence/ppi_method_context_v2.tsv"
OUT_FUNCTION="data/output/evidence/ppi_function_context_v2.tsv"

CONTRACT_METHOD="pipelines/ppi_semantic_enrichment/contracts/ppi_method_context_v2.json"
CONTRACT_FUNCTION="pipelines/ppi_semantic_enrichment/contracts/ppi_function_context_v2.json"

REPORT_SAMPLE="pipelines/ppi_semantic_enrichment/reports/ppi_semantic_enrichment_v2.sample.json"
REPORT_BUILD="pipelines/ppi_semantic_enrichment/reports/ppi_semantic_enrichment_v2.build.json"
REPORT_METHOD_VALIDATION="pipelines/ppi_semantic_enrichment/reports/ppi_method_context_v2.validation.json"
REPORT_FUNCTION_VALIDATION="pipelines/ppi_semantic_enrichment/reports/ppi_function_context_v2.validation.json"
REPORT_QA="pipelines/ppi_semantic_enrichment/reports/ppi_semantic_enrichment_v2.qa.json"
REPORT_MANIFEST="pipelines/ppi_semantic_enrichment/reports/ppi_semantic_enrichment_v2.manifest.json"
REPORT_RESULT_EVIDENCE="pipelines/ppi_semantic_enrichment/reports/ppi_semantic_enrichment_v2.result_evidence.json"
REPORT_BLOCKED_INPUTS="pipelines/ppi_semantic_enrichment/reports/ppi_semantic_enrichment_v2.blocked_missing_inputs.json"
SAMPLE_METHOD_20="pipelines/ppi_semantic_enrichment/reports/ppi_method_context_v2.sample20.tsv"
SAMPLE_FUNCTION_20="pipelines/ppi_semantic_enrichment/reports/ppi_function_context_v2.sample20.tsv"

CACHE_DIR="pipelines/ppi_semantic_enrichment/data/raw"
DATA_VERSION="${DATA_VERSION:-kg-data-local}"
SAMPLE_ROWS="${SAMPLE_ROWS:-5000}"

if [[ ! -f "$INPUT_EDGES" ]]; then
  echo "[ERROR] missing edges table: $INPUT_EDGES" >&2
  exit 2
fi
if [[ ! -f "$INPUT_MASTER" ]]; then
  echo "[ERROR] missing protein master table: $INPUT_MASTER" >&2
  exit 2
fi
if [[ ! -f "$INPUT_REACTOME" ]]; then
  echo "[ERROR] missing Reactome mapping: $INPUT_REACTOME" >&2
  exit 2
fi
if [[ ! -f "$INPUT_KEGG" ]]; then
  echo "[ERROR] missing KEGG mapping: $INPUT_KEGG" >&2
  exit 2
fi

mkdir -p "$(dirname "$OUT_METHOD")" "$(dirname "$REPORT_SAMPLE")" "$CACHE_DIR"

# 1) Required sample-first run (fast mode)
SAMPLE_METHOD_TMP="$(mktemp -t ppi_method_context_v2.sample.XXXXXX.tsv)"
SAMPLE_FUNCTION_TMP="$(mktemp -t ppi_function_context_v2.sample.XXXXXX.tsv)"
trap 'rm -f "$SAMPLE_METHOD_TMP" "$SAMPLE_FUNCTION_TMP"' EXIT

python3 pipelines/ppi_semantic_enrichment/scripts/01_build_ppi_semantic_enrichment_v2.py \
  --input-edges "$INPUT_EDGES" \
  --input-master "$INPUT_MASTER" \
  --input-reactome "$INPUT_REACTOME" \
  --input-kegg "$INPUT_KEGG" \
  --out-method "$SAMPLE_METHOD_TMP" \
  --out-function "$SAMPLE_FUNCTION_TMP" \
  --report "$REPORT_SAMPLE" \
  --cache-dir "$CACHE_DIR" \
  --max-rows "$SAMPLE_ROWS" \
  --sample-fast

# 2) Full run
python3 pipelines/ppi_semantic_enrichment/scripts/01_build_ppi_semantic_enrichment_v2.py \
  --input-edges "$INPUT_EDGES" \
  --input-master "$INPUT_MASTER" \
  --input-reactome "$INPUT_REACTOME" \
  --input-kegg "$INPUT_KEGG" \
  --out-method "$OUT_METHOD" \
  --out-function "$OUT_FUNCTION" \
  --report "$REPORT_BUILD" \
  --cache-dir "$CACHE_DIR"

# 3) Contract validation
python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_METHOD" \
  --table "$OUT_METHOD" \
  --out "$REPORT_METHOD_VALIDATION"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT_FUNCTION" \
  --table "$OUT_FUNCTION" \
  --out "$REPORT_FUNCTION_VALIDATION"

# 4) QA gates / acceptance metrics
python3 pipelines/ppi_semantic_enrichment/scripts/02_qa_ppi_semantic_enrichment_v2.py \
  --edges "$INPUT_EDGES" \
  --method-context "$OUT_METHOD" \
  --function-context "$OUT_FUNCTION" \
  --method-validation "$REPORT_METHOD_VALIDATION" \
  --function-validation "$REPORT_FUNCTION_VALIDATION" \
  --out "$REPORT_QA"

# 5) Manifest
python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$OUT_METHOD" \
  "$OUT_FUNCTION" \
  "$INPUT_EDGES" \
  "$INPUT_MASTER" \
  "$INPUT_REACTOME" \
  "$INPUT_KEGG"

# 6) Result-evidence report (coverage/source distribution/random sample20 + blocked inputs)
python3 pipelines/ppi_semantic_enrichment/scripts/03_emit_ppi_semantic_result_evidence_v2.py \
  --method-table "$OUT_METHOD" \
  --function-table "$OUT_FUNCTION" \
  --qa-report "$REPORT_QA" \
  --build-report "$REPORT_BUILD" \
  --sample-size 20 \
  --out-evidence-report "$REPORT_RESULT_EVIDENCE" \
  --out-blocked-report "$REPORT_BLOCKED_INPUTS" \
  --out-method-sample "$SAMPLE_METHOD_20" \
  --out-function-sample "$SAMPLE_FUNCTION_20"

echo "[DONE] ppi_semantic_enrichment_v2: sample + full + validation + QA + manifest + result evidence finished."
