#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List


DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")
EXPECTED_DATASETS = [
    "protein_master_v6_clean.tsv",
    "protein_edges.tsv",
    "ptm_sites.tsv",
    "pathway_members.tsv",
    "alphafold_quality.tsv",
    "hgnc_core.tsv",
]


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _non_empty(v: str) -> bool:
    return (v or "").strip() not in ("", "NA", "N/A", "None", "null")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    if not args.table.exists():
        raise SystemExit(f"[ERROR] missing table: {args.table}")

    dataset_counter: Counter[str] = Counter()
    rows: List[Dict[str, str]] = []

    with args.table.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {"dataset", "primary_source", "source_version", "fetch_date", "evidence_field", "notes"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] table missing required columns: {sorted(missing)}")
        for row in r:
            rows.append(row)
            dataset_counter[(row.get("dataset") or "").strip()] += 1

    total_rows = len(rows)
    source_version_non_empty = sum(1 for r in rows if _non_empty(r.get("source_version", "")))
    fetch_date_valid = sum(1 for r in rows if DATE_RE.match((r.get("fetch_date") or "").strip() or ""))

    source_version_non_empty_rate = 1.0 if total_rows == 0 else source_version_non_empty / total_rows
    fetch_date_valid_rate = 1.0 if total_rows == 0 else fetch_date_valid / total_rows

    missing_expected = [d for d in EXPECTED_DATASETS if dataset_counter[d] == 0]
    duplicated = {k: v for k, v in dataset_counter.items() if v > 1}
    unexpected = sorted([k for k in dataset_counter if k and k not in EXPECTED_DATASETS])

    master_row = next((r for r in rows if (r.get("dataset") or "").strip() == "protein_master_v6_clean.tsv"), None)
    master_has_multi_source_rule = False
    if master_row:
        sv = (master_row.get("source_version") or "").strip()
        notes = (master_row.get("notes") or "").strip().lower()
        master_has_multi_source_rule = ("|" in sv and "uniprot" in sv.lower()) or ("composite" in notes)

    gates = [
        {
            "id": "source_version_non_empty_100pct",
            "passed": abs(source_version_non_empty_rate - 1.0) < 1e-12,
            "detail": {"rate": source_version_non_empty_rate, "non_empty": source_version_non_empty, "total_rows": total_rows},
        },
        {
            "id": "fetch_date_iso_yyyy_mm_dd_100pct",
            "passed": abs(fetch_date_valid_rate - 1.0) < 1e-12,
            "detail": {"rate": fetch_date_valid_rate, "valid_rows": fetch_date_valid, "total_rows": total_rows},
        },
        {
            "id": "each_dataset_exactly_one_row",
            "passed": (not missing_expected) and (not duplicated) and (not unexpected),
            "detail": {
                "missing_expected": missing_expected,
                "duplicated": duplicated,
                "unexpected": unexpected,
                "dataset_counter": dict(dataset_counter),
            },
        },
        {
            "id": "protein_master_multi_source_rule_documented",
            "passed": master_has_multi_source_rule,
            "detail": {
                "dataset": "protein_master_v6_clean.tsv",
                "rule_detected": master_has_multi_source_rule,
            },
        },
    ]

    passed = all(g["passed"] for g in gates)

    report: Dict[str, Any] = {
        "name": "protein_source_versions_v1.qa",
        "created_at": _utc_now(),
        "input_table": str(args.table),
        "metrics": {
            "rows": total_rows,
            "source_version_non_empty_rate": source_version_non_empty_rate,
            "fetch_date_valid_rate": fetch_date_valid_rate,
            "expected_dataset_count": len(EXPECTED_DATASETS),
            "actual_dataset_count": len(dataset_counter),
        },
        "gates": gates,
        "passed": passed,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if passed else "FAIL"
    print(f"[{status}] protein_source_versions_v1 QA -> {args.out}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
