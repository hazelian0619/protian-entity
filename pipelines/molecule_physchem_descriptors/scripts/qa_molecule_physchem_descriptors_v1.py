#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Set

INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
FLOAT_RE = re.compile(r"^-?\d+(?:\.\d+)?$")
INT_RE = re.compile(r"^-?\d+$")

FIELDS_FLOAT = ["molecular_weight", "logp", "tpsa"]
FIELDS_INT = ["hbd", "hba", "rotatable_bonds"]
FIELDS_ALL = FIELDS_FLOAT + FIELDS_INT
UNIT_STRATEGY_TEXT = "MW:Da|LogP:unitless|TPSA:Å²|HBD:count|HBA:count|RotatableBonds:count"


@dataclass
class Gate:
    gate_id: str
    passed: bool
    detail: Dict[str, Any]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize(x: object) -> str:
    return str(x or "").strip()


def split_multi(x: object) -> List[str]:
    text = normalize(x)
    if not text:
        return []
    return [t.strip() for t in text.split(";") if t.strip()]


def load_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        return list(r)


def is_valid_numeric(field: str, value: str) -> bool:
    text = normalize(value)
    if not text:
        return False
    if field in FIELDS_FLOAT:
        if not FLOAT_RE.match(text):
            return False
        try:
            float(text)
            return True
        except Exception:
            return False
    if field in FIELDS_INT:
        if not INT_RE.match(text):
            return False
        try:
            int(text)
            return True
        except Exception:
            return False
    return False


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xref", type=Path, required=True)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--build-report", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-backlink-rate", type=float, default=1.0)
    ap.add_argument("--min-format-valid-rate", type=float, default=0.99)
    ap.add_argument("--min-mappable-non-empty-rate", type=float, default=0.90)
    ap.add_argument("--min-core-coverage", type=float, default=0.70)
    args = ap.parse_args()

    for p in [args.xref, args.table, args.build_report]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    xref_rows = load_tsv(args.xref)
    table_rows = load_tsv(args.table)
    build = json.loads(args.build_report.read_text(encoding="utf-8"))

    xref_map: Dict[str, Dict[str, Any]] = {}
    for row in xref_rows:
        ik = normalize(row.get("inchikey", "")).upper()
        if not ik:
            continue
        has_mapping = bool(split_multi(row.get("chembl_id", "")) or split_multi(row.get("pubchem_cid", "")))
        xref_map[ik] = {"has_mapping": has_mapping}

    seen: Set[str] = set()
    duplicate_ik = 0
    bad_ik = 0
    missing_backlink = 0

    numeric_non_empty = 0
    numeric_valid = 0
    numeric_invalid_samples: List[Dict[str, Any]] = []

    field_non_empty = {f: 0 for f in FIELDS_ALL}
    field_valid = {f: 0 for f in FIELDS_ALL}
    mappable_non_empty = {f: 0 for f in FIELDS_ALL}

    mappable_rows = 0
    rows_with_core = 0
    unit_strategy_mismatch = 0
    source_empty = 0
    source_version_empty = 0
    fetch_date_empty = 0

    for row in table_rows:
        ik = normalize(row.get("inchikey", "")).upper()
        if not INCHIKEY_RE.match(ik):
            bad_ik += 1
        if ik in seen:
            duplicate_ik += 1
        seen.add(ik)

        xref_rec = xref_map.get(ik)
        if xref_rec is None:
            missing_backlink += 1
            has_mapping = False
        else:
            has_mapping = bool(xref_rec.get("has_mapping", False))

        if has_mapping:
            mappable_rows += 1

        if all(normalize(row.get(f, "")) for f in ["molecular_weight", "logp", "tpsa"]):
            rows_with_core += 1

        for f in FIELDS_ALL:
            value = normalize(row.get(f, ""))
            if not value:
                continue
            field_non_empty[f] += 1
            numeric_non_empty += 1
            ok = is_valid_numeric(f, value)
            if ok:
                field_valid[f] += 1
                numeric_valid += 1
            else:
                if len(numeric_invalid_samples) < 20:
                    numeric_invalid_samples.append({"inchikey": ik, "field": f, "value": value})
            if has_mapping:
                mappable_non_empty[f] += 1

        if normalize(row.get("unit_strategy", "")) != UNIT_STRATEGY_TEXT:
            unit_strategy_mismatch += 1

        if not normalize(row.get("source", "")):
            source_empty += 1
        if not normalize(row.get("source_version", "")):
            source_version_empty += 1
        if not normalize(row.get("fetch_date", "")):
            fetch_date_empty += 1

    rows = len(table_rows)
    backlink_rate = ((rows - missing_backlink) / rows) if rows else 0.0
    numeric_format_valid_rate = (numeric_valid / numeric_non_empty) if numeric_non_empty else 0.0
    core_coverage_rate = (rows_with_core / rows) if rows else 0.0

    mappable_field_rates = {
        f: ((mappable_non_empty[f] / mappable_rows) if mappable_rows else 0.0) for f in FIELDS_ALL
    }
    mappable_non_empty_total = int(sum(mappable_non_empty.values()))
    mappable_total_cells = int(mappable_rows * len(FIELDS_ALL))
    mappable_overall_non_empty_rate = (
        (mappable_non_empty_total / mappable_total_cells) if mappable_total_cells else 0.0
    )
    field_valid_rates = {
        f: ((field_valid[f] / field_non_empty[f]) if field_non_empty[f] else 1.0) for f in FIELDS_ALL
    }

    build_rows = int(build.get("metrics", {}).get("rows_written", rows))
    build_core_rate = float(build.get("metrics", {}).get("core_triple_coverage_rate", core_coverage_rate))

    gates = [
        Gate(
            gate_id="inchikey_unique_and_valid",
            passed=(bad_ik == 0 and duplicate_ik == 0),
            detail={"bad_inchikey": bad_ik, "duplicate_inchikey": duplicate_ik},
        ),
        Gate(
            gate_id="backlink_rate_1p0",
            passed=(backlink_rate >= args.min_backlink_rate),
            detail={
                "backlink_rate": backlink_rate,
                "min_required": args.min_backlink_rate,
                "missing_backlink_rows": missing_backlink,
            },
        ),
        Gate(
            gate_id="numeric_format_valid_rate",
            passed=(numeric_format_valid_rate >= args.min_format_valid_rate),
            detail={
                "numeric_non_empty_cells": numeric_non_empty,
                "numeric_valid_cells": numeric_valid,
                "valid_rate": numeric_format_valid_rate,
                "min_required": args.min_format_valid_rate,
                "invalid_samples_top20": numeric_invalid_samples,
                "field_valid_rates": field_valid_rates,
            },
        ),
        Gate(
            gate_id="mappable_subset_non_empty_rate",
            passed=(mappable_overall_non_empty_rate >= args.min_mappable_non_empty_rate),
            detail={
                "mappable_rows": mappable_rows,
                "overall_non_empty_rate": mappable_overall_non_empty_rate,
                "non_empty_cells": mappable_non_empty_total,
                "total_cells": mappable_total_cells,
                "field_rates": mappable_field_rates,
                "min_required": args.min_mappable_non_empty_rate,
            },
        ),
        Gate(
            gate_id="core_triple_coverage",
            passed=(core_coverage_rate >= args.min_core_coverage),
            detail={
                "core_fields": ["molecular_weight", "logp", "tpsa"],
                "coverage_rate": core_coverage_rate,
                "min_required": args.min_core_coverage,
                "rows_with_core": rows_with_core,
                "rows": rows,
            },
        ),
        Gate(
            gate_id="metadata_fields_non_empty_and_unit_strategy",
            passed=(
                source_empty == 0
                and source_version_empty == 0
                and fetch_date_empty == 0
                and unit_strategy_mismatch == 0
            ),
            detail={
                "source_empty_rows": source_empty,
                "source_version_empty_rows": source_version_empty,
                "fetch_date_empty_rows": fetch_date_empty,
                "unit_strategy_mismatch_rows": unit_strategy_mismatch,
            },
        ),
        Gate(
            gate_id="build_report_consistency",
            passed=(build_rows == rows and abs(build_core_rate - core_coverage_rate) < 1e-9),
            detail={
                "build_rows": build_rows,
                "actual_rows": rows,
                "build_core_triple_coverage_rate": build_core_rate,
                "actual_core_triple_coverage_rate": core_coverage_rate,
            },
        ),
    ]

    report = {
        "name": "molecule_physchem_descriptors_v1.qa",
        "created_at": utc_now(),
        "inputs": {
            "xref": str(args.xref),
            "table": str(args.table),
            "build_report": str(args.build_report),
        },
        "thresholds": {
            "min_backlink_rate": args.min_backlink_rate,
            "min_format_valid_rate": args.min_format_valid_rate,
            "min_mappable_non_empty_rate": args.min_mappable_non_empty_rate,
            "min_core_coverage": args.min_core_coverage,
        },
        "metrics": {
            "rows": rows,
            "xref_rows": len(xref_rows),
            "backlink_rate": backlink_rate,
            "mappable_rows": mappable_rows,
            "numeric_non_empty_cells": numeric_non_empty,
            "numeric_valid_cells": numeric_valid,
            "numeric_format_valid_rate": numeric_format_valid_rate,
            "field_non_empty_counts": field_non_empty,
            "field_valid_rates": field_valid_rates,
            "mappable_field_non_empty_rates": mappable_field_rates,
            "mappable_overall_non_empty_rate": mappable_overall_non_empty_rate,
            "core_triple_coverage_rate": core_coverage_rate,
        },
        "gates": [{"id": g.gate_id, "passed": g.passed, "detail": g.detail} for g in gates],
        "passed": all(g.passed for g in gates),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if report["passed"] else "FAIL"
    print(f"[{status}] molecule_physchem_descriptors_v1 QA -> {args.out}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
