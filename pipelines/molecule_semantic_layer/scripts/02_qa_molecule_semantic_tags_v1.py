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
    ap.add_argument("--hierarchy-report", type=Path, required=True)
    ap.add_argument("--coverage-report", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-backlink-rate", type=float, default=0.99)
    ap.add_argument("--min-source-non-empty-rate", type=float, default=0.99)
    args = ap.parse_args()

    for p in [args.xref, args.table, args.build_report, args.hierarchy_report, args.coverage_report]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    xref_rows = load_tsv(args.xref)
    out_rows = load_tsv(args.table)

    build = json.loads(args.build_report.read_text(encoding="utf-8"))
    hierarchy = json.loads(args.hierarchy_report.read_text(encoding="utf-8"))
    coverage = json.loads(args.coverage_report.read_text(encoding="utf-8"))

    xref_ik = {normalize(r.get("inchikey", "")).upper() for r in xref_rows if normalize(r.get("inchikey", ""))}

    missing_backlink: List[str] = []
    bad_ik = 0
    dup_ik = 0
    seen: Set[str] = set()

    chebi_rows = 0
    chebi_expandable_rows = 0

    for r in out_rows:
        ik = normalize(r.get("inchikey", "")).upper()
        if not INCHIKEY_RE.match(ik):
            bad_ik += 1
        if ik in seen:
            dup_ik += 1
        else:
            seen.add(ik)

        if ik not in xref_ik:
            missing_backlink.append(ik)

        chebi_id = normalize(r.get("chebi_id", ""))
        chebi_ancestor = normalize(r.get("chebi_ancestor_ids", ""))
        if chebi_id and chebi_id != "NA":
            chebi_rows += 1
            if chebi_ancestor and chebi_ancestor != "NA":
                chebi_expandable_rows += 1

    backlink_rate = ((len(out_rows) - len(missing_backlink)) / len(out_rows)) if out_rows else 0.0
    chebi_expandable_rate = (chebi_expandable_rows / chebi_rows) if chebi_rows else 0.0

    source_cols = ["chebi_source", "atc_source", "drugbank_source", "zinc_source", "source_version"]
    source_rates = {c: non_empty_rate(out_rows, c) for c in source_cols}

    gates = [
        Gate(
            gate_id="inchikey_unique_and_valid",
            passed=(bad_ik == 0 and dup_ik == 0),
            detail={"bad_inchikey": bad_ik, "duplicate_inchikey": dup_ik},
        ),
        Gate(
            gate_id="inchikey_backlink_rate",
            passed=(backlink_rate >= args.min_backlink_rate),
            detail={
                "backlink_rate": backlink_rate,
                "min_required": args.min_backlink_rate,
                "missing_count": len(missing_backlink),
                "missing_sample": sorted(missing_backlink)[:50],
            },
        ),
        Gate(
            gate_id="chebi_hierarchy_expandable",
            passed=(chebi_rows == 0 or chebi_expandable_rate >= 0.99),
            detail={
                "chebi_rows": chebi_rows,
                "chebi_expandable_rows": chebi_expandable_rows,
                "expandable_rate": chebi_expandable_rate,
                "hierarchy_report_expandable_rate": hierarchy.get("expandable_rate_among_chebi_rows", 0.0),
            },
        ),
        Gate(
            gate_id="classification_source_anchor_non_empty",
            passed=all(rate >= args.min_source_non_empty_rate for rate in source_rates.values()),
            detail={
                "rates": source_rates,
                "min_required": args.min_source_non_empty_rate,
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
        "name": "molecule_semantic_tags_v1.qa",
        "created_at": utc_now(),
        "inputs": {
            "xref": str(args.xref),
            "table": str(args.table),
            "build_report": str(args.build_report),
            "hierarchy_report": str(args.hierarchy_report),
            "coverage_report": str(args.coverage_report),
        },
        "metrics": {
            "rows": len(out_rows),
            "xref_rows": len(xref_rows),
            "backlink_rate": backlink_rate,
            "chebi_rows": chebi_rows,
            "chebi_expandable_rows": chebi_expandable_rows,
            "chebi_expandable_rate": chebi_expandable_rate,
            "source_non_empty_rates": source_rates,
        },
        "gates": [{"id": g.gate_id, "passed": g.passed, "detail": g.detail} for g in gates],
        "passed": all(g.passed for g in gates),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if report["passed"] else "FAIL"
    print(f"[{status}] molecule_semantic_tags_v1 QA -> {args.out}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
