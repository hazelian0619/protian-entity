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


@dataclass
class Gate:
    gate_id: str
    passed: bool
    detail: Dict[str, Any]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize(x: str) -> str:
    return (x or "").strip()


def split_multi(x: str) -> List[str]:
    t = normalize(x)
    if not t:
        return []
    return [s.strip() for s in t.split(";") if s.strip()]


def load_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        return list(r)


def non_empty_rate(rows: List[Dict[str, str]], col: str) -> float:
    if not rows:
        return 0.0
    n = sum(1 for r in rows if normalize(r.get(col, "")) != "")
    return n / len(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xref", type=Path, required=True)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--build-report", type=Path, required=True)
    ap.add_argument("--coverage-report", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-non-empty-rate", type=float, default=0.99)
    args = ap.parse_args()

    for p in [args.xref, args.table, args.build_report, args.coverage_report]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    xref_rows = load_tsv(args.xref)
    reg_rows = load_tsv(args.table)
    build = json.loads(args.build_report.read_text(encoding="utf-8"))
    coverage = json.loads(args.coverage_report.read_text(encoding="utf-8"))

    xref_iks = {normalize(r.get("inchikey", "")).upper() for r in xref_rows if normalize(r.get("inchikey", ""))}

    bad_ik = 0
    dup_ik = 0
    seen: Set[str] = set()
    missing_in_xref: List[str] = []

    allowed_source_types = {"computational", "experimental", "mixed", "unknown"}
    bad_source_types = 0

    computational_rows = 0
    experimental_rows = 0

    for row in reg_rows:
        ik = normalize(row.get("inchikey", "")).upper()
        if not INCHIKEY_RE.match(ik):
            bad_ik += 1
        if ik in seen:
            dup_ik += 1
        else:
            seen.add(ik)

        if ik not in xref_iks:
            missing_in_xref.append(ik)

        sts = set(split_multi(row.get("source_types", "")))
        if not sts:
            bad_source_types += 1
        elif not sts.issubset(allowed_source_types):
            bad_source_types += 1

        if "computational" in sts:
            computational_rows += 1
        if "experimental" in sts:
            experimental_rows += 1

    source_rate = non_empty_rate(reg_rows, "source_providers")
    fmt_rate = non_empty_rate(reg_rows, "structure_formats")
    avail_rate = non_empty_rate(reg_rows, "availability_status")

    gates = [
        Gate(
            gate_id="inchikey_unique_and_valid",
            passed=(bad_ik == 0 and dup_ik == 0),
            detail={"bad_inchikey": bad_ik, "duplicate_inchikey": dup_ik},
        ),
        Gate(
            gate_id="inchikey_backlink_to_xref",
            passed=(len(missing_in_xref) == 0),
            detail={
                "missing_count": len(missing_in_xref),
                "missing_sample": sorted(missing_in_xref)[:50],
                "backlink_rate": ((len(reg_rows) - len(missing_in_xref)) / len(reg_rows)) if reg_rows else 0.0,
            },
        ),
        Gate(
            gate_id="source_field_non_empty_rate",
            passed=(source_rate >= args.min_non_empty_rate),
            detail={"rate": source_rate, "min_required": args.min_non_empty_rate},
        ),
        Gate(
            gate_id="format_field_non_empty_rate",
            passed=(fmt_rate >= args.min_non_empty_rate),
            detail={"rate": fmt_rate, "min_required": args.min_non_empty_rate},
        ),
        Gate(
            gate_id="availability_field_non_empty_rate",
            passed=(avail_rate >= args.min_non_empty_rate),
            detail={"rate": avail_rate, "min_required": args.min_non_empty_rate},
        ),
        Gate(
            gate_id="source_types_controlled_vocab",
            passed=(bad_source_types == 0),
            detail={"bad_source_type_rows": bad_source_types, "allowed": sorted(allowed_source_types)},
        ),
        Gate(
            gate_id="computational_type_present",
            passed=(computational_rows > 0),
            detail={"computational_rows": computational_rows, "experimental_rows": experimental_rows},
        ),
        Gate(
            gate_id="report_consistency",
            passed=(
                int(build.get("metrics", {}).get("rows_written", len(reg_rows))) == len(reg_rows)
                and int(coverage.get("rows", len(reg_rows))) == len(reg_rows)
            ),
            detail={
                "build_rows": int(build.get("metrics", {}).get("rows_written", len(reg_rows))),
                "coverage_rows": int(coverage.get("rows", len(reg_rows))),
                "actual_rows": len(reg_rows),
            },
        ),
    ]

    report = {
        "name": "molecule_3d_registry_v1.qa",
        "created_at": utc_now(),
        "inputs": {
            "xref": str(args.xref),
            "table": str(args.table),
            "build_report": str(args.build_report),
            "coverage_report": str(args.coverage_report),
        },
        "metrics": {
            "rows": len(reg_rows),
            "xref_rows": len(xref_rows),
            "backlink_missing": len(missing_in_xref),
            "source_non_empty_rate": source_rate,
            "format_non_empty_rate": fmt_rate,
            "availability_non_empty_rate": avail_rate,
            "computational_rows": computational_rows,
            "experimental_rows": experimental_rows,
        },
        "gates": [{"id": g.gate_id, "passed": g.passed, "detail": g.detail} for g in gates],
        "passed": all(g.passed for g in gates),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if report["passed"] else "FAIL"
    print(f"[{status}] molecule_3d_registry_v1 QA -> {args.out}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
