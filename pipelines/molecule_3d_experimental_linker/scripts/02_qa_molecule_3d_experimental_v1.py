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
PDB_RE = re.compile(r"^[0-9A-Z]{4}$")


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
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-backlink-rate", type=float, default=0.99)
    ap.add_argument("--min-pdb-non-empty-rate", type=float, default=0.99)
    ap.add_argument("--min-method-non-empty-rate", type=float, default=0.99)
    ap.add_argument("--min-resolution-non-empty-rate", type=float, default=0.99)
    args = ap.parse_args()

    for p in [args.xref, args.table, args.build_report, args.coverage_report]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    xref_rows = load_tsv(args.xref)
    out_rows = load_tsv(args.table)

    build = json.loads(args.build_report.read_text(encoding="utf-8"))
    coverage = json.loads(args.coverage_report.read_text(encoding="utf-8"))

    xref_ik = {normalize(r.get("inchikey", "")).upper() for r in xref_rows if normalize(r.get("inchikey", ""))}

    bad_ik = 0
    bad_pdb = 0
    dup_key = 0
    seen: Set[str] = set()
    missing_backlink: List[str] = []

    experimental_rows = 0
    computational_rows = 0
    bad_source_type = 0

    for r in out_rows:
        ik = normalize(r.get("inchikey", "")).upper()
        if not INCHIKEY_RE.match(ik):
            bad_ik += 1

        pdb_id = normalize(r.get("pdb_id", "")).upper()
        if not PDB_RE.match(pdb_id):
            bad_pdb += 1

        key = "|".join(
            [
                ik,
                normalize(r.get("uniprot_id", "")),
                pdb_id,
                normalize(r.get("ligand_ccd_id", "")),
                normalize(r.get("pdb_entity_id", "")),
            ]
        )
        if key in seen:
            dup_key += 1
        else:
            seen.add(key)

        if ik not in xref_ik:
            missing_backlink.append(ik)

        st = normalize(r.get("source_type", "")).lower()
        if st == "experimental":
            experimental_rows += 1
        elif st == "computational":
            computational_rows += 1
        else:
            bad_source_type += 1

    backlink_rate = ((len(out_rows) - len(missing_backlink)) / len(out_rows)) if out_rows else 0.0

    pdb_rate = non_empty_rate(out_rows, "pdb_id")
    method_rate = non_empty_rate(out_rows, "experimental_method")
    resolution_rate = non_empty_rate(out_rows, "resolution")

    gates = [
        Gate(
            gate_id="inchikey_and_pdb_format_with_unique_rows",
            passed=(bad_ik == 0 and bad_pdb == 0 and dup_key == 0),
            detail={
                "bad_inchikey": bad_ik,
                "bad_pdb_id": bad_pdb,
                "duplicate_rows": dup_key,
            },
        ),
        Gate(
            gate_id="inchikey_backlink_rate_mappable_subset",
            passed=(backlink_rate >= args.min_backlink_rate),
            detail={
                "rate": backlink_rate,
                "min_required": args.min_backlink_rate,
                "missing_count": len(missing_backlink),
                "missing_sample": sorted(set(missing_backlink))[:50],
            },
        ),
        Gate(
            gate_id="experimental_rows_present",
            passed=(experimental_rows > 0),
            detail={
                "rows": len(out_rows),
                "experimental_rows": experimental_rows,
            },
        ),
        Gate(
            gate_id="source_type_split_experimental_vs_computational",
            passed=(experimental_rows == len(out_rows) and computational_rows == 0 and bad_source_type == 0),
            detail={
                "experimental_rows": experimental_rows,
                "computational_rows": computational_rows,
                "bad_source_type_rows": bad_source_type,
            },
        ),
        Gate(
            gate_id="pdb_method_resolution_non_empty_rate",
            passed=(
                pdb_rate >= args.min_pdb_non_empty_rate
                and method_rate >= args.min_method_non_empty_rate
                and resolution_rate >= args.min_resolution_non_empty_rate
            ),
            detail={
                "pdb_non_empty_rate": pdb_rate,
                "method_non_empty_rate": method_rate,
                "resolution_non_empty_rate": resolution_rate,
                "min_required": {
                    "pdb": args.min_pdb_non_empty_rate,
                    "method": args.min_method_non_empty_rate,
                    "resolution": args.min_resolution_non_empty_rate,
                },
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
        "name": "molecule_3d_experimental_v1.qa",
        "created_at": utc_now(),
        "inputs": {
            "xref": str(args.xref),
            "table": str(args.table),
            "build_report": str(args.build_report),
            "coverage_report": str(args.coverage_report),
        },
        "metrics": {
            "rows": len(out_rows),
            "xref_rows": len(xref_rows),
            "experimental_rows": experimental_rows,
            "computational_rows": computational_rows,
            "backlink_rate": backlink_rate,
            "pdb_non_empty_rate": pdb_rate,
            "method_non_empty_rate": method_rate,
            "resolution_non_empty_rate": resolution_rate,
        },
        "gates": [{"id": g.gate_id, "passed": g.passed, "detail": g.detail} for g in gates],
        "passed": all(g.passed for g in gates),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if report["passed"] else "FAIL"
    print(f"[{status}] molecule_3d_experimental_v1 QA -> {args.out}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
