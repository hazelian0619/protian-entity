#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

PIPELINE_DIR="pipelines/rna_external_xref"

INPUT_MASTER="${INPUT_MASTER:-data/output/rna_master_v1.tsv}"
INPUT_ID_MAPPING="${INPUT_ID_MAPPING:-data/raw/rna/rnacentral/id_mapping.tsv.gz}"
INPUT_AUX_DIR="${INPUT_AUX_DIR:-data/raw/rna/aux_xref}"

OUT_TABLE="data/output/rna_external_xref_v1.tsv"
CONTRACT="$PIPELINE_DIR/contracts/rna_external_xref_v1.json"

REPORT_BUILD="$PIPELINE_DIR/reports/rna_external_xref_v1.build.json"
REPORT_COVERAGE="$PIPELINE_DIR/reports/rna_external_xref_v1.coverage.json"
REPORT_CONFLICTS="$PIPELINE_DIR/reports/rna_external_xref_v1.conflicts.json"
REPORT_VALIDATION="$PIPELINE_DIR/reports/rna_external_xref_v1.validation.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/rna_external_xref_v1.manifest.json"
REPORT_MANUAL="$PIPELINE_DIR/reports/rna_external_xref_v1.manual_download.json"

SAMPLE_TABLE="$PIPELINE_DIR/.cache/rna_external_xref_v1.sample.tsv"
SAMPLE_BUILD="$PIPELINE_DIR/reports/rna_external_xref_v1.sample.build.json"
SAMPLE_COVERAGE="$PIPELINE_DIR/reports/rna_external_xref_v1.sample.coverage.json"
SAMPLE_CONFLICTS="$PIPELINE_DIR/reports/rna_external_xref_v1.sample.conflicts.json"
SAMPLE_VALIDATION="$PIPELINE_DIR/reports/rna_external_xref_v1.sample.validation.json"

DATA_VERSION="${DATA_VERSION:-kg-rna-external-xref-v1}"
SOURCE_VERSION="${SOURCE_VERSION:-RNAcentral:25}"
FETCH_DATE="${FETCH_DATE:-$(date +%F)}"
MIN_BACKLINK_RATE="${MIN_BACKLINK_RATE:-0.99}"
SAMPLE_MASTER_ROWS="${SAMPLE_MASTER_ROWS:-20000}"
SAMPLE_IDMAP_LINES="${SAMPLE_IDMAP_LINES:-3000000}"

mkdir -p "$(dirname "$OUT_TABLE")" "$PIPELINE_DIR/reports" "$PIPELINE_DIR/.cache"
rm -f "$REPORT_MANUAL"

missing=()
[[ -f "$INPUT_MASTER" ]] || missing+=("$INPUT_MASTER")
[[ -f "$INPUT_ID_MAPPING" ]] || missing+=("$INPUT_ID_MAPPING")
[[ -d "$INPUT_AUX_DIR" ]] || missing+=("$INPUT_AUX_DIR (directory)")

aux_file_count=0
if [[ -d "$INPUT_AUX_DIR" ]]; then
  aux_file_count="$(find "$INPUT_AUX_DIR" -maxdepth 1 -type f \( -name '*.tsv' -o -name '*.tsv.gz' -o -name '*.txt' -o -name '*.txt.gz' \) | wc -l | tr -d ' ')"
  if [[ "$aux_file_count" == "0" ]]; then
    missing+=("$INPUT_AUX_DIR/*.{tsv,tsv.gz,txt,txt.gz}")
  fi
fi

if (( ${#missing[@]} > 0 )); then
  echo "[ERROR] missing required inputs for rna_external_xref_v1:" >&2
  for f in "${missing[@]}"; do
    echo "  - $f" >&2
  done

  python3 - <<'PY'
import json
from pathlib import Path

report = {
    "name": "rna_external_xref_v1.manual_download",
    "reason": "required inputs missing; cannot build external xref pack",
    "download_checklist": [
        {
            "name": "RNA master table",
            "expected_path": "data/output/rna_master_v1.tsv",
            "sha256_cmd": "shasum -a 256 data/output/rna_master_v1.tsv",
        },
        {
            "name": "RNAcentral id_mapping",
            "expected_path": "data/raw/rna/rnacentral/id_mapping.tsv.gz",
            "sha256_cmd": "shasum -a 256 data/raw/rna/rnacentral/id_mapping.tsv.gz",
        },
        {
            "name": "Aux ENST/RefSeq mapping pack",
            "expected_path": "data/raw/rna/aux_xref/",
            "required_columns": [
                "ensembl_transcript_id",
                "refseq_mrna OR rna_nucleotide_accession.version",
            ],
            "sha256_cmd": "cd data/raw/rna/aux_xref && shasum -a 256 * > SHA256SUMS.txt && shasum -a 256 -c SHA256SUMS.txt",
        },
    ],
}

out = Path("pipelines/rna_external_xref/reports/rna_external_xref_v1.manual_download.json")
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(f"[INFO] manual checklist -> {out}")
PY

  exit 2
fi

echo "[INFO] Inputs:"
echo "  master      = $INPUT_MASTER"
echo "  id_mapping  = $INPUT_ID_MAPPING"
echo "  aux_dir     = $INPUT_AUX_DIR (files=$aux_file_count)"

echo "[STEP 1/5] sample run"
python3 "$PIPELINE_DIR/scripts/01_build_rna_external_xref_v1.py" \
  --master "$INPUT_MASTER" \
  --id-mapping "$INPUT_ID_MAPPING" \
  --aux-dir "$INPUT_AUX_DIR" \
  --out-table "$SAMPLE_TABLE" \
  --build-report "$SAMPLE_BUILD" \
  --coverage-report "$SAMPLE_COVERAGE" \
  --conflicts-report "$SAMPLE_CONFLICTS" \
  --source-version "$SOURCE_VERSION" \
  --fetch-date "$FETCH_DATE" \
  --min-backlink-rate "$MIN_BACKLINK_RATE" \
  --max-master-rows "$SAMPLE_MASTER_ROWS" \
  --max-idmap-lines "$SAMPLE_IDMAP_LINES"

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$SAMPLE_TABLE" \
  --out "$SAMPLE_VALIDATION"

echo "[STEP 2/5] full run"
python3 "$PIPELINE_DIR/scripts/01_build_rna_external_xref_v1.py" \
  --master "$INPUT_MASTER" \
  --id-mapping "$INPUT_ID_MAPPING" \
  --aux-dir "$INPUT_AUX_DIR" \
  --out-table "$OUT_TABLE" \
  --build-report "$REPORT_BUILD" \
  --coverage-report "$REPORT_COVERAGE" \
  --conflicts-report "$REPORT_CONFLICTS" \
  --source-version "$SOURCE_VERSION" \
  --fetch-date "$FETCH_DATE" \
  --min-backlink-rate "$MIN_BACKLINK_RATE"

echo "[STEP 3/5] contract validation"
python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUT_TABLE" \
  --out "$REPORT_VALIDATION"

echo "[STEP 4/5] gates"
python3 - <<'PY'
import json
from pathlib import Path

coverage = json.loads(Path("pipelines/rna_external_xref/reports/rna_external_xref_v1.coverage.json").read_text(encoding="utf-8"))
validation = json.loads(Path("pipelines/rna_external_xref/reports/rna_external_xref_v1.validation.json").read_text(encoding="utf-8"))

backlink_rate = float(coverage.get("summary", {}).get("backlink_rate", 0.0))
required = float(coverage.get("summary", {}).get("min_backlink_rate_required", 0.99))
by_db = coverage.get("per_xref_db", {})
validation_passed = bool(validation.get("passed"))

ok = backlink_rate >= required and validation_passed and bool(by_db)
if ok:
    print("[PASS] gates passed.")
    raise SystemExit(0)

print("[FAIL] gates failed.")
print(f"  backlink_rate={backlink_rate:.6f}, required={required:.6f}")
print(f"  validation_passed={validation_passed}")
print(f"  per_xref_db_present={bool(by_db)}")
raise SystemExit(3)
PY

echo "[STEP 5/5] manifest"
python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$OUT_TABLE" \
  "$INPUT_MASTER" \
  "$INPUT_ID_MAPPING" \
  "$CONTRACT" \
  "$REPORT_BUILD" \
  "$REPORT_COVERAGE" \
  "$REPORT_CONFLICTS" \
  "$REPORT_VALIDATION"

echo "[DONE] rna_external_xref_v1 finished (sample + full + validation + coverage + conflicts + manifest)."
