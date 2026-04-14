#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Dict

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipelines.psi_condition_enrichment.scripts.condition_extractors import extract_condition_bundle, non_empty


NEW_COLUMNS = [
    "condition_context_json",
    "condition_extract_confidence",
    "condition_extract_source_field",
]

AUDIT_HEADER = [
    "audit_id",
    "activity_id",
    "edge_id",
    "conflict_flag",
    "conflict_fields",
    "selected_condition_pH",
    "selected_condition_temperature_c",
    "selected_condition_context",
    "condition_extract_source_field",
    "condition_extract_confidence",
    "resolution_strategy",
    "raw_condition_candidates_json",
    "raw_source_snapshot_json",
    "fetch_date",
]

ALLOWED_REL = {"=", "<", ">", "<=", ">=", "~"}


def utc_now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def safe_rate(num: int, den: int) -> float:
    if den <= 0:
        return 0.0
    return num / den


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-v2", type=Path, required=True)
    ap.add_argument("--out-v3", type=Path, required=True)
    ap.add_argument("--out-audit", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--fetch-date", default=dt.datetime.now(dt.timezone.utc).date().isoformat())
    ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("--source-version-suffix", default="condition_enrich_v3")
    args = ap.parse_args()

    if not args.in_v2.exists():
        raise SystemExit(f"[ERROR] missing input: {args.in_v2}")

    args.out_v3.parent.mkdir(parents=True, exist_ok=True)
    args.out_audit.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)

    tmp_v3 = args.out_v3.with_suffix(args.out_v3.suffix + ".tmp")
    tmp_audit = args.out_audit.with_suffix(args.out_audit.suffix + ".tmp")

    counts = Counter()
    rule_hits = Counter()
    source_field_hits = Counter()

    with (
        args.in_v2.open("r", encoding="utf-8", newline="") as fi,
        tmp_v3.open("w", encoding="utf-8", newline="") as fo,
        tmp_audit.open("w", encoding="utf-8", newline="") as fa,
    ):
        reader = csv.DictReader(fi, delimiter="\t")
        in_fields = list(reader.fieldnames or [])
        if not in_fields:
            raise SystemExit("[ERROR] empty input header")

        required_base = {
            "edge_id",
            "activity_id",
            "standard_type",
            "standard_relation",
            "assay_type",
            "condition_context",
            "condition_pH",
            "condition_temperature_c",
            "source_version",
        }
        miss = required_base - set(in_fields)
        if miss:
            raise SystemExit(f"[ERROR] input v2 missing columns: {sorted(miss)}")

        out_fields = list(in_fields)
        for c in NEW_COLUMNS:
            if c not in out_fields:
                out_fields.append(c)

        writer = csv.DictWriter(fo, fieldnames=out_fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()

        wa = csv.DictWriter(fa, fieldnames=AUDIT_HEADER, delimiter="\t", lineterminator="\n")
        wa.writeheader()

        for idx, row in enumerate(reader, start=1):
            if args.max_rows > 0 and idx > args.max_rows:
                break

            counts["rows_total"] += 1

            if non_empty(row.get("edge_id")):
                counts["v2_edge_id_non_empty"] += 1
            if non_empty(row.get("standard_type")):
                counts["v2_standard_type_non_empty"] += 1
            if str(row.get("standard_relation") or "").strip() in ALLOWED_REL:
                counts["v2_standard_relation_allowed"] += 1
            if non_empty(row.get("assay_type")):
                counts["v2_assay_type_non_empty"] += 1

            if non_empty(row.get("condition_context")):
                counts["v2_condition_context_non_empty"] += 1
            if non_empty(row.get("condition_pH")):
                counts["v2_condition_pH_non_empty"] += 1
            if non_empty(row.get("condition_temperature_c")):
                counts["v2_condition_temperature_non_empty"] += 1

            bundle = extract_condition_bundle(row)
            out_row: Dict[str, str] = dict(row)
            out_row["condition_pH"] = bundle["condition_pH"]
            out_row["condition_temperature_c"] = bundle["condition_temperature_c"]
            out_row["condition_system"] = bundle["condition_system"]
            out_row["condition_context"] = bundle["condition_context"]
            out_row["condition_context_json"] = bundle["condition_context_json"]
            out_row["condition_extract_confidence"] = bundle["condition_extract_confidence"]
            out_row["condition_extract_source_field"] = bundle["condition_extract_source_field"]

            source_version_old = str(row.get("source_version") or "").strip()
            if source_version_old:
                out_row["source_version"] = f"{source_version_old}|{args.source_version_suffix}"
            else:
                out_row["source_version"] = args.source_version_suffix

            writer.writerow(out_row)

            if non_empty(out_row.get("edge_id")):
                counts["v3_edge_id_non_empty"] += 1
            if non_empty(out_row.get("standard_type")):
                counts["v3_standard_type_non_empty"] += 1
            if str(out_row.get("standard_relation") or "").strip() in ALLOWED_REL:
                counts["v3_standard_relation_allowed"] += 1
            if non_empty(out_row.get("assay_type")):
                counts["v3_assay_type_non_empty"] += 1

            if non_empty(out_row.get("condition_context")):
                counts["v3_condition_context_non_empty"] += 1
            if non_empty(out_row.get("condition_pH")):
                counts["v3_condition_pH_non_empty"] += 1
            if non_empty(out_row.get("condition_temperature_c")):
                counts["v3_condition_temperature_non_empty"] += 1
            if non_empty(out_row.get("condition_context_json")):
                counts["v3_condition_context_json_non_empty"] += 1
            if non_empty(out_row.get("condition_extract_confidence")):
                counts["v3_condition_extract_confidence_non_empty"] += 1
            if non_empty(out_row.get("condition_extract_source_field")):
                counts["v3_condition_extract_source_field_non_empty"] += 1

            rh = str(bundle.get("_rule_hits") or "").strip()
            if rh:
                for hit in rh.split(";"):
                    if hit:
                        rule_hits[hit] += 1

            src_fields = str(out_row.get("condition_extract_source_field") or "").strip()
            if src_fields:
                for sf in src_fields.split(";"):
                    source_field_hits[sf] += 1

            if bundle.get("conflict_flag") == "true":
                counts["conflict_rows"] += 1
                wa.writerow(
                    {
                        "audit_id": f"psi-cond-audit-{idx}",
                        "activity_id": str(row.get("activity_id") or ""),
                        "edge_id": str(row.get("edge_id") or ""),
                        "conflict_flag": "true",
                        "conflict_fields": bundle.get("conflict_fields", ""),
                        "selected_condition_pH": out_row.get("condition_pH", ""),
                        "selected_condition_temperature_c": out_row.get("condition_temperature_c", ""),
                        "selected_condition_context": out_row.get("condition_context", ""),
                        "condition_extract_source_field": out_row.get("condition_extract_source_field", ""),
                        "condition_extract_confidence": out_row.get("condition_extract_confidence", ""),
                        "resolution_strategy": "priority: v2_existing > explicit_regex > inferred_buffer_default > inferred_cell_culture; keep raw candidates",
                        "raw_condition_candidates_json": bundle.get("_audit_candidates_json", ""),
                        "raw_source_snapshot_json": bundle.get("_source_snapshot_json", ""),
                        "fetch_date": args.fetch_date,
                    }
                )

    tmp_v3.replace(args.out_v3)
    tmp_audit.replace(args.out_audit)

    rows = counts["rows_total"]

    report = {
        "name": "psi_condition_enrichment_v3.build",
        "created_at": utc_now_iso(),
        "inputs": {
            "psi_activity_context_v2": str(args.in_v2),
            "max_rows": args.max_rows,
        },
        "outputs": {
            "psi_activity_context_v3": str(args.out_v3),
            "psi_condition_parse_audit_v3": str(args.out_audit),
        },
        "counts": dict(counts),
        "coverage": {
            "v2": {
                "edge_id": safe_rate(counts["v2_edge_id_non_empty"], rows),
                "standard_type": safe_rate(counts["v2_standard_type_non_empty"], rows),
                "standard_relation_allowed": safe_rate(counts["v2_standard_relation_allowed"], rows),
                "assay_type": safe_rate(counts["v2_assay_type_non_empty"], rows),
                "condition_context": safe_rate(counts["v2_condition_context_non_empty"], rows),
                "condition_pH": safe_rate(counts["v2_condition_pH_non_empty"], rows),
                "condition_temperature_c": safe_rate(counts["v2_condition_temperature_non_empty"], rows),
            },
            "v3": {
                "edge_id": safe_rate(counts["v3_edge_id_non_empty"], rows),
                "standard_type": safe_rate(counts["v3_standard_type_non_empty"], rows),
                "standard_relation_allowed": safe_rate(counts["v3_standard_relation_allowed"], rows),
                "assay_type": safe_rate(counts["v3_assay_type_non_empty"], rows),
                "condition_context": safe_rate(counts["v3_condition_context_non_empty"], rows),
                "condition_pH": safe_rate(counts["v3_condition_pH_non_empty"], rows),
                "condition_temperature_c": safe_rate(counts["v3_condition_temperature_non_empty"], rows),
                "condition_context_json": safe_rate(counts["v3_condition_context_json_non_empty"], rows),
                "condition_extract_confidence": safe_rate(counts["v3_condition_extract_confidence_non_empty"], rows),
                "condition_extract_source_field": safe_rate(counts["v3_condition_extract_source_field_non_empty"], rows),
            },
        },
        "rule_hits": dict(sorted(rule_hits.items())),
        "source_field_hits": dict(sorted(source_field_hits.items())),
        "conflict_summary": {
            "conflict_rows": counts["conflict_rows"],
            "conflict_rate": safe_rate(counts["conflict_rows"], rows),
        },
    }

    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(
        "[OK] psi_condition_enrichment_v3 build "
        f"rows={rows:,} ctx={report['coverage']['v3']['condition_context']:.4f} "
        f"pH={report['coverage']['v3']['condition_pH']:.4f} "
        f"temp={report['coverage']['v3']['condition_temperature_c']:.4f} "
        f"conflicts={counts['conflict_rows']:,}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
