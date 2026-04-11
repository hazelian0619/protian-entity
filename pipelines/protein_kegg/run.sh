#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

INPUT_MASTER="data/processed/protein_master_v6_clean.tsv"
INPUT_REACTOME="data/processed/pathway_members.tsv"

OUT_KEGG="data/output/protein/protein_kegg_pathway_v1.tsv"
CACHE_DIR="pipelines/protein_kegg/.cache"

CONTRACT="pipelines/protein_kegg/contracts/protein_kegg_pathway_v1.json"

REPORT_SAMPLE="pipelines/protein_kegg/reports/protein_kegg_pathway_v1.sample.json"
REPORT_BUILD="pipelines/protein_kegg/reports/protein_kegg_pathway_v1.build.json"
REPORT_VALIDATION="pipelines/protein_kegg/reports/protein_kegg_pathway_v1.validation.json"
REPORT_QA="pipelines/protein_kegg/reports/protein_kegg_pathway_v1.qa.json"
REPORT_MANIFEST="pipelines/protein_kegg/reports/protein_kegg_pathway_v1.manifest.json"

DATA_VERSION="${DATA_VERSION:-kg-data-local}"
SAMPLE_ROWS="${SAMPLE_ROWS:-200}"

if [[ ! -f "$INPUT_MASTER" ]]; then
  echo "[ERROR] missing input table: $INPUT_MASTER" >&2
  exit 2
fi

if [[ ! -f "$INPUT_REACTOME" ]]; then
  echo "[ERROR] missing reactome table: $INPUT_REACTOME" >&2
  exit 2
fi

mkdir -p "$(dirname "$OUT_KEGG")" "$(dirname "$REPORT_BUILD")" "$CACHE_DIR"

# 1) Minimal sample-first run (required by team convention)
SAMPLE_TMP="$(mktemp -t protein_kegg_sample.XXXXXX.tsv)"
trap 'rm -f "$SAMPLE_TMP"' EXIT

python3 pipelines/protein_kegg/scripts/01_build_protein_kegg_pathway_v1.py \
  --input-master "$INPUT_MASTER" \
  --out "$SAMPLE_TMP" \
  --report "$REPORT_SAMPLE" \
  --cache-dir "$CACHE_DIR" \
  --max-rows "$SAMPLE_ROWS"

# 2) Full run
python3 pipelines/protein_kegg/scripts/01_build_protein_kegg_pathway_v1.py \
  --input-master "$INPUT_MASTER" \
  --out "$OUT_KEGG" \
  --report "$REPORT_BUILD" \
  --cache-dir "$CACHE_DIR"

# 3) Contract validation
python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUT_KEGG" \
  --out "$REPORT_VALIDATION"

# 4) QA gates (join coverage + KEGG format + Reactome complementarity)
python3 pipelines/protein_kegg/scripts/02_qa_protein_kegg_pathway_v1.py \
  --kegg "$OUT_KEGG" \
  --master "$INPUT_MASTER" \
  --reactome "$INPUT_REACTOME" \
  --out "$REPORT_QA"

# 5) Manifest
python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$OUT_KEGG" \
  "$INPUT_MASTER" \
  "$INPUT_REACTOME"

echo "[DONE] protein_kegg_pathway_v1: sample + full + validation + QA + manifest finished."
