#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

PIPELINE_DIR="pipelines/molecule_pk_tox_layer"

resolve_input() {
  local primary="$1"
  shift
  if [[ -n "$primary" && -f "$primary" ]]; then
    echo "$primary"
    return 0
  fi
  for candidate in "$@"; do
    if [[ -f "$candidate" ]]; then
      echo "$candidate"
      return 0
    fi
  done
  echo ""
}

REQ_XREF_V2="${XREF_V2_PATH:-data/output/molecules/molecule_xref_core_v2.tsv}"
REQ_DRUGBANK_XML="${DRUGBANK_XML_PATH:-data/raw/drugbank/drugbank_complete_database_2026-01-04.xml}"
OPT_PUBCHEM_TOX="${PUBCHEM_TOX_PATH:-data/raw/molecules/pubchem/pubchem_toxicity_summary_v1.tsv}"

IN_XREF="$(resolve_input "$REQ_XREF_V2" "../1218/data/output/molecules/molecule_xref_core_v2.tsv")"
IN_DRUGBANK_XML="$(resolve_input "$REQ_DRUGBANK_XML" "../1218/data/raw/drugbank/drugbank_complete_database_2026-01-04.xml" "/Users/pluviophile/graph/0108/full database.xml")"
IN_PUBCHEM_TOX="$(resolve_input "$OPT_PUBCHEM_TOX" "../1218/data/raw/molecules/pubchem/pubchem_toxicity_summary_v1.tsv")"

OUT_TABLE="data/output/molecules/molecule_pk_tox_v1.tsv"
CONTRACT="$PIPELINE_DIR/contracts/molecule_pk_tox_v1.json"

SMOKE_MAX_ROWS="${SMOKE_MAX_ROWS:-5000}"
DATA_VERSION="${DATA_VERSION:-kg-data-local}"

SMOKE_TABLE="$PIPELINE_DIR/.cache/molecule_pk_tox_v1.smoke.tsv"
SMOKE_BUILD="$PIPELINE_DIR/reports/molecule_pk_tox_v1.smoke.build.json"
SMOKE_VALIDATION="$PIPELINE_DIR/reports/molecule_pk_tox_v1.smoke.validation.json"
SMOKE_QA="$PIPELINE_DIR/reports/molecule_pk_tox_v1.smoke.qa.json"

REPORT_BUILD="$PIPELINE_DIR/reports/molecule_pk_tox_v1.build.json"
REPORT_VALIDATION="$PIPELINE_DIR/reports/molecule_pk_tox_v1.validation.json"
REPORT_QA="$PIPELINE_DIR/reports/molecule_pk_tox_v1.qa.json"
REPORT_MANIFEST="$PIPELINE_DIR/reports/molecule_pk_tox_v1.manifest.json"
REPORT_BLOCKED="$PIPELINE_DIR/reports/molecule_pk_tox_v1.blocked_or_manual_download.json"

mkdir -p "$(dirname "$OUT_TABLE")" "$PIPELINE_DIR/reports" "$PIPELINE_DIR/.cache"

if [[ -z "$IN_XREF" ]]; then
  echo "[ERROR] missing required input: $REQ_XREF_V2" >&2
  python3 - <<'PY'
import json
from pathlib import Path
report = {
  "name": "molecule_pk_tox_v1.blocked_or_manual_download",
  "status": "blocked_missing_required_input",
  "reason": "molecule_xref_core_v2.tsv is required",
  "download_checklist": [
    {
      "name": "molecule_xref_core_v2.tsv",
      "expected_path": "data/output/molecules/molecule_xref_core_v2.tsv",
      "sha256_cmd": "sha256sum data/output/molecules/molecule_xref_core_v2.tsv"
    }
  ]
}
p = Path("pipelines/molecule_pk_tox_layer/reports/molecule_pk_tox_v1.blocked_or_manual_download.json")
p.parent.mkdir(parents=True, exist_ok=True)
p.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(f"[INFO] blocked report -> {p}")
PY
  exit 2
fi

if [[ -z "$IN_DRUGBANK_XML" && -z "$IN_PUBCHEM_TOX" ]]; then
  echo "[ERROR] no usable PK/Tox sources (DrugBank extension or public tox source) were found." >&2
  python3 - <<'PY'
import json
from pathlib import Path
report = {
  "name": "molecule_pk_tox_v1.blocked_or_manual_download",
  "status": "blocked_missing_sources",
  "reason": "Neither DrugBank extension XML nor public toxicity source file is available.",
  "manual_download": {
    "required_any_of": [
      {
        "name": "DrugBank full XML (licensed)",
        "expected_path": "data/raw/drugbank/drugbank_complete_database_2026-01-04.xml",
        "authorization_required": True,
        "notes": "Requires institutional/commercial DrugBank license; place XML locally without committing."
      },
      {
        "name": "PubChem toxicity summary (local exported TSV/CSV)",
        "expected_path": "data/raw/molecules/pubchem/pubchem_toxicity_summary_v1.tsv",
        "authorization_required": False,
        "notes": "Any local curated export with CID + toxicity columns is accepted."
      }
    ]
  }
}
p = Path("pipelines/molecule_pk_tox_layer/reports/molecule_pk_tox_v1.blocked_or_manual_download.json")
p.parent.mkdir(parents=True, exist_ok=True)
p.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(f"[INFO] blocked report -> {p}")
PY
  exit 2
fi

python3 - <<PY
import json
from pathlib import Path
report = {
  "name": "molecule_pk_tox_v1.blocked_or_manual_download",
  "status": "ready",
  "reason": "At least one PK/Tox source is available.",
  "resolved_inputs": {
    "xref": "$IN_XREF",
    "drugbank_xml": "$IN_DRUGBANK_XML" if "$IN_DRUGBANK_XML" else None,
    "pubchem_tox": "$IN_PUBCHEM_TOX" if "$IN_PUBCHEM_TOX" else None
  },
  "missing_optional_inputs": [
    v for v in [
      None if "$IN_DRUGBANK_XML" else "$REQ_DRUGBANK_XML",
      None if "$IN_PUBCHEM_TOX" else "$OPT_PUBCHEM_TOX"
    ] if v
  ]
}
p = Path("$REPORT_BLOCKED")
p.parent.mkdir(parents=True, exist_ok=True)
p.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(f"[INFO] readiness report -> {p}")
PY

echo "[INFO] resolved inputs:"
echo "  xref         = $IN_XREF"
if [[ -n "$IN_DRUGBANK_XML" ]]; then
  echo "  drugbank_xml = $IN_DRUGBANK_XML"
else
  echo "  drugbank_xml = <missing, optional>"
fi
if [[ -n "$IN_PUBCHEM_TOX" ]]; then
  echo "  pubchem_tox  = $IN_PUBCHEM_TOX"
else
  echo "  pubchem_tox  = <missing, optional>"
fi

DRUGBANK_ARG=()
PUBCHEM_ARG=()
if [[ -n "$IN_DRUGBANK_XML" ]]; then
  DRUGBANK_ARG=(--drugbank-xml "$IN_DRUGBANK_XML")
fi
if [[ -n "$IN_PUBCHEM_TOX" ]]; then
  PUBCHEM_ARG=(--pubchem-tox "$IN_PUBCHEM_TOX")
fi

echo "[STEP] smoke build (max_rows=$SMOKE_MAX_ROWS)"
python3 "$PIPELINE_DIR/scripts/01_build_molecule_pk_tox_v1.py" \
  --xref "$IN_XREF" \
  --out "$SMOKE_TABLE" \
  --report "$SMOKE_BUILD" \
  --max-rows "$SMOKE_MAX_ROWS" \
  ${DRUGBANK_ARG[@]+"${DRUGBANK_ARG[@]}"} \
  ${PUBCHEM_ARG[@]+"${PUBCHEM_ARG[@]}"}

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$SMOKE_TABLE" \
  --out "$SMOKE_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_pk_tox_v1.py" \
  --xref "$IN_XREF" \
  --table "$SMOKE_TABLE" \
  --build-report "$SMOKE_BUILD" \
  --out "$SMOKE_QA" \
  --min-rows-with-data 0

echo "[STEP] full build"
python3 "$PIPELINE_DIR/scripts/01_build_molecule_pk_tox_v1.py" \
  --xref "$IN_XREF" \
  --out "$OUT_TABLE" \
  --report "$REPORT_BUILD" \
  ${DRUGBANK_ARG[@]+"${DRUGBANK_ARG[@]}"} \
  ${PUBCHEM_ARG[@]+"${PUBCHEM_ARG[@]}"}

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$OUT_TABLE" \
  --out "$REPORT_VALIDATION"

python3 "$PIPELINE_DIR/scripts/02_qa_molecule_pk_tox_v1.py" \
  --xref "$IN_XREF" \
  --table "$OUT_TABLE" \
  --build-report "$REPORT_BUILD" \
  --out "$REPORT_QA" \
  --min-rows-with-data 0

ROWS_WITH_DATA="$(python3 - <<'PY'
import json
from pathlib import Path
build = json.loads(Path("pipelines/molecule_pk_tox_layer/reports/molecule_pk_tox_v1.build.json").read_text(encoding="utf-8"))
print(int(build.get("metrics", {}).get("rows_with_any_pk_tox", 0)))
PY
)"

if [[ "$ROWS_WITH_DATA" -le 0 ]]; then
  echo "[ERROR] no PK/Tox values extracted from available sources." >&2
  python3 - <<'PY'
import json
from pathlib import Path
p = Path("pipelines/molecule_pk_tox_layer/reports/molecule_pk_tox_v1.blocked_or_manual_download.json")
report = {
  "name": "molecule_pk_tox_v1.blocked_or_manual_download",
  "status": "blocked_no_extracted_data",
  "reason": "Source files are present but no target PK/Tox fields were extracted.",
  "manual_download": {
    "suggestions": [
      "Verify DrugBank XML is full licensed export (not demo/sample).",
      "Provide a public toxicity table with CID and at least one of ld50/mutagenicity.",
      "Re-run with DRUGBANK_XML_PATH / PUBCHEM_TOX_PATH explicitly set."
    ]
  }
}
p.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(f"[INFO] blocked report -> {p}")
PY
  exit 2
fi

python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$OUT_TABLE" \
  "$IN_XREF" \
  "$CONTRACT" \
  "$REPORT_BUILD" \
  "$REPORT_VALIDATION" \
  "$REPORT_QA" \
  "$REPORT_BLOCKED"

echo "[DONE] molecule_pk_tox_layer finished (smoke + full + validation + qa + manifest)."
