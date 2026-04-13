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
FIELDS = [
    "bioavailability",
    "clearance",
    "half_life",
    "volume_of_distribution",
    "ld50",
    "mutagenicity",
]
EVIDENCE_ALLOWED = {"structured/high", "text_mined/medium", "inferred/low", "none"}


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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xref", type=Path, required=True)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--build-report", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-backlink-rate", type=float, default=0.999)
    ap.add_argument("--min-rows-with-data", type=int, default=0)
    args = ap.parse_args()

    for p in [args.xref, args.table, args.build_report]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    xref_rows = load_tsv(args.xref)
    table_rows = load_tsv(args.table)
    build = json.loads(args.build_report.read_text(encoding="utf-8"))

    xref_ik = {normalize(r.get("inchikey", "")).upper() for r in xref_rows if normalize(r.get("inchikey", ""))}

    bad_ik = 0
    dup_ik = 0
    missing_backlink = 0
    seen: Set[str] = set()
    rows_with_any = 0
    anchor_violations = 0
    evidence_violations = 0
    field_non_empty = {f: 0 for f in FIELDS}

    for r in table_rows:
        ik = normalize(r.get("inchikey", "")).upper()
        if not INCHIKEY_RE.match(ik):
            bad_ik += 1
        if ik in seen:
            dup_ik += 1
        seen.add(ik)
        if ik not in xref_ik:
            missing_backlink += 1

        any_field = False
        for f in FIELDS:
            val = normalize(r.get(f, ""))
            ev = normalize(r.get(f"{f}_evidence", "none"))
            src = normalize(r.get(f"{f}_source", "none"))

            if ev not in EVIDENCE_ALLOWED:
                evidence_violations += 1

            if val:
                any_field = True
                field_non_empty[f] += 1
                if ev == "none" or src == "none" or not src:
                    anchor_violations += 1
        if any_field:
            rows_with_any += 1

    backlink_rate = ((len(table_rows) - missing_backlink) / len(table_rows)) if table_rows else 0.0
    build_rows = int(build.get("metrics", {}).get("rows_written", len(table_rows)))
    build_rows_with_any = int(build.get("metrics", {}).get("rows_with_any_pk_tox", rows_with_any))

    gates = [
        Gate(
            gate_id="inchikey_unique_and_valid",
            passed=(bad_ik == 0 and dup_ik == 0),
            detail={"bad_inchikey": bad_ik, "duplicate_inchikey": dup_ik},
        ),
        Gate(
            gate_id="xref_backlink_rate",
            passed=(backlink_rate >= args.min_backlink_rate),
            detail={
                "backlink_rate": backlink_rate,
                "min_required": args.min_backlink_rate,
                "missing_backlink_rows": missing_backlink,
            },
        ),
        Gate(
            gate_id="evidence_tier_domain",
            passed=(evidence_violations == 0),
            detail={"violations": evidence_violations, "allowed": sorted(EVIDENCE_ALLOWED)},
        ),
        Gate(
            gate_id="value_requires_anchor",
            passed=(anchor_violations == 0),
            detail={"anchor_violations": anchor_violations},
        ),
        Gate(
            gate_id="report_consistency",
            passed=(build_rows == len(table_rows) and build_rows_with_any == rows_with_any),
            detail={
                "build_rows": build_rows,
                "actual_rows": len(table_rows),
                "build_rows_with_any_pk_tox": build_rows_with_any,
                "actual_rows_with_any_pk_tox": rows_with_any,
            },
        ),
        Gate(
            gate_id="rows_with_data_gate",
            passed=(rows_with_any >= args.min_rows_with_data),
            detail={
                "rows_with_any_pk_tox": rows_with_any,
                "min_required": args.min_rows_with_data,
            },
        ),
    ]

    report = {
        "name": "molecule_pk_tox_v1.qa",
        "created_at": utc_now(),
        "inputs": {
            "xref": str(args.xref),
            "table": str(args.table),
            "build_report": str(args.build_report),
        },
        "metrics": {
            "rows": len(table_rows),
            "xref_rows": len(xref_rows),
            "rows_with_any_pk_tox": rows_with_any,
            "field_non_empty_counts": field_non_empty,
            "backlink_rate": backlink_rate,
            "anchor_violations": anchor_violations,
            "evidence_violations": evidence_violations,
        },
        "gates": [{"id": g.gate_id, "passed": g.passed, "detail": g.detail} for g in gates],
        "passed": all(g.passed for g in gates),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if report["passed"] else "FAIL"
    print(f"[{status}] molecule_pk_tox_v1 QA -> {args.out}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
