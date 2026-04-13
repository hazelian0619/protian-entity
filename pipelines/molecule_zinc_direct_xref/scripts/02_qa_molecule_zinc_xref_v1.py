#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
ZINC_ID_RE = re.compile(r"^ZINC\d+$")


@dataclass
class Gate:
    gate_id: str
    passed: bool
    detail: Dict[str, Any]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize(x: str) -> str:
    return (x or "").strip()


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
    ap.add_argument("--conflict-report", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-backlink-rate", type=float, default=0.99)
    ap.add_argument("--min-source-anchor-rate", type=float, default=1.0)
    ap.add_argument("--min-coverage-rate", type=float, default=0.05)
    ap.add_argument("--min-coverage-delta-vs-baseline", type=float, default=0.02)
    args = ap.parse_args()

    for p in [args.xref, args.table, args.build_report, args.coverage_report, args.conflict_report]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    xref_rows = load_tsv(args.xref)
    out_rows = load_tsv(args.table)

    build = json.loads(args.build_report.read_text(encoding="utf-8"))
    coverage = json.loads(args.coverage_report.read_text(encoding="utf-8"))
    conflict = json.loads(args.conflict_report.read_text(encoding="utf-8"))

    xref_ik = {normalize(r.get("inchikey", "")).upper() for r in xref_rows if normalize(r.get("inchikey", ""))}

    bad_ik = 0
    bad_zinc = 0
    missing_backlink: List[str] = []
    duplicate_rows = 0
    seen: Set[Tuple[str, str, str, str]] = set()

    for r in out_rows:
        ik = normalize(r.get("inchikey", "")).upper()
        zid = normalize(r.get("zinc_id", "")).upper()

        if not INCHIKEY_RE.match(ik):
            bad_ik += 1
        if not ZINC_ID_RE.match(zid):
            bad_zinc += 1

        key = (
            ik,
            zid,
            normalize(r.get("supplier_name", "")),
            normalize(r.get("supplier_code", "")),
            normalize(r.get("catalog_tier", "")),
        )
        if key in seen:
            duplicate_rows += 1
        else:
            seen.add(key)

        if ik not in xref_ik:
            missing_backlink.append(ik)

    backlink_rate = ((len(out_rows) - len(missing_backlink)) / len(out_rows)) if out_rows else 0.0

    source_rate = non_empty_rate(out_rows, "source")
    source_version_rate = non_empty_rate(out_rows, "source_version")
    fetch_date_rate = non_empty_rate(out_rows, "fetch_date")

    cov_rate = float(coverage.get("coverage_rate", 0.0))
    baseline_rate = float(coverage.get("baseline_rate", 0.0))
    cov_delta = float(coverage.get("coverage_delta_vs_baseline", cov_rate - baseline_rate))

    conflict_metrics = conflict.get("metrics", {})
    conflict_has_structure = isinstance(conflict.get("samples", {}), dict) and "inchikey_with_multi_zinc" in conflict.get("samples", {})

    gates = [
        Gate(
            gate_id="rows_non_empty_and_format_ok",
            passed=(len(out_rows) > 0 and bad_ik == 0 and bad_zinc == 0 and duplicate_rows == 0),
            detail={
                "rows": len(out_rows),
                "bad_inchikey": bad_ik,
                "bad_zinc_id": bad_zinc,
                "duplicate_rows": duplicate_rows,
            },
        ),
        Gate(
            gate_id="inchikey_backlink_rate",
            passed=(backlink_rate >= args.min_backlink_rate),
            detail={
                "rate": backlink_rate,
                "min_required": args.min_backlink_rate,
                "missing_count": len(missing_backlink),
                "missing_sample": sorted(set(missing_backlink))[:50],
            },
        ),
        Gate(
            gate_id="source_anchor_fields_non_empty",
            passed=(
                source_rate >= args.min_source_anchor_rate
                and source_version_rate >= args.min_source_anchor_rate
                and fetch_date_rate >= args.min_source_anchor_rate
            ),
            detail={
                "source_rate": source_rate,
                "source_version_rate": source_version_rate,
                "fetch_date_rate": fetch_date_rate,
                "min_required": args.min_source_anchor_rate,
            },
        ),
        Gate(
            gate_id="coverage_significantly_higher_than_baseline",
            passed=(cov_rate >= args.min_coverage_rate and cov_delta >= args.min_coverage_delta_vs_baseline),
            detail={
                "coverage_rate": cov_rate,
                "baseline_rate": baseline_rate,
                "coverage_delta": cov_delta,
                "min_coverage_rate": args.min_coverage_rate,
                "min_delta_vs_baseline": args.min_coverage_delta_vs_baseline,
            },
        ),
        Gate(
            gate_id="conflict_audit_present",
            passed=(conflict_has_structure and isinstance(conflict_metrics, dict)),
            detail={
                "conflict_metrics": conflict_metrics,
                "has_structured_samples": conflict_has_structure,
            },
        ),
        Gate(
            gate_id="report_consistency",
            passed=(
                int(build.get("metrics", {}).get("rows_written", len(out_rows))) == len(out_rows)
                and int(coverage.get("rows", len(out_rows))) == len(out_rows)
            ),
            detail={
                "build_rows": int(build.get("metrics", {}).get("rows_written", len(out_rows))),
                "coverage_rows": int(coverage.get("rows", len(out_rows))),
                "actual_rows": len(out_rows),
            },
        ),
    ]

    report = {
        "name": "molecule_zinc_xref_v1.qa",
        "created_at": utc_now(),
        "inputs": {
            "xref": str(args.xref),
            "table": str(args.table),
            "build_report": str(args.build_report),
            "coverage_report": str(args.coverage_report),
            "conflict_report": str(args.conflict_report),
        },
        "metrics": {
            "rows": len(out_rows),
            "xref_rows": len(xref_rows),
            "backlink_rate": backlink_rate,
            "source_rate": source_rate,
            "source_version_rate": source_version_rate,
            "fetch_date_rate": fetch_date_rate,
            "coverage_rate": cov_rate,
            "baseline_rate": baseline_rate,
            "coverage_delta_vs_baseline": cov_delta,
        },
        "gates": [{"id": g.gate_id, "passed": g.passed, "detail": g.detail} for g in gates],
        "passed": all(g.passed for g in gates),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if report["passed"] else "FAIL"
    print(f"[{status}] molecule_zinc_xref_v1 QA -> {args.out}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
