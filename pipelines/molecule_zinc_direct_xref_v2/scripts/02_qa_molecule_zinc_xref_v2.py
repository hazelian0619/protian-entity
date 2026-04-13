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
BACKFILL_SOURCE_TAG = "molecule_structure_identifiers_v1.canonical_smiles"


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
    return [p.strip() for p in t.split(";") if p.strip()]


def load_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        return list(r)


def non_empty_rate(rows: List[Dict[str, str]], col: str) -> float:
    if not rows:
        return 0.0
    return sum(1 for r in rows if normalize(r.get(col, ""))) / len(rows)


def map_key(row: Dict[str, str]) -> Tuple[str, str, str, str, str]:
    return (
        normalize(row.get("inchikey", "")).upper(),
        normalize(row.get("zinc_id", "")).upper(),
        normalize(row.get("supplier_name", "")),
        normalize(row.get("supplier_code", "")),
        normalize(row.get("catalog_tier", "")),
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--zinc-v1", type=Path, required=True)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--build-report", type=Path, required=True)
    ap.add_argument("--coverage-report", type=Path, required=True)
    ap.add_argument("--conflict-report", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-smiles-delta", type=float, default=0.15)
    ap.add_argument("--allow-row-subset", type=int, default=0)
    args = ap.parse_args()

    for p in [args.zinc_v1, args.table, args.build_report, args.coverage_report, args.conflict_report]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    allow_row_subset = bool(args.allow_row_subset)
    v1_rows = load_tsv(args.zinc_v1)
    v2_rows = load_tsv(args.table)
    build = json.loads(args.build_report.read_text(encoding="utf-8"))
    coverage = json.loads(args.coverage_report.read_text(encoding="utf-8"))
    conflict = json.loads(args.conflict_report.read_text(encoding="utf-8"))

    v1_by_key: Dict[Tuple[str, str, str, str, str], Dict[str, str]] = {map_key(r): r for r in v1_rows}
    v2_by_key: Dict[Tuple[str, str, str, str, str], Dict[str, str]] = {map_key(r): r for r in v2_rows}

    v1_key_set = set(v1_by_key.keys())
    v2_key_set = set(v2_by_key.keys())

    if allow_row_subset:
        mapping_equivalent = v2_key_set.issubset(v1_key_set)
        baseline_rows = [v1_by_key[k] for k in v2_key_set if k in v1_by_key]
    else:
        mapping_equivalent = (len(v2_rows) == len(v1_rows)) and (v2_key_set == v1_key_set)
        baseline_rows = v1_rows

    bad_inchikey = sum(1 for r in v2_rows if not INCHIKEY_RE.match(normalize(r.get("inchikey", "")).upper()))
    bad_zinc = sum(1 for r in v2_rows if not ZINC_ID_RE.match(normalize(r.get("zinc_id", "")).upper()))

    base_smiles_rate = non_empty_rate(baseline_rows, "smiles")
    v2_smiles_rate = non_empty_rate(v2_rows, "smiles")
    smiles_delta = v2_smiles_rate - base_smiles_rate

    source_rate = non_empty_rate(v2_rows, "source")
    source_version_rate = non_empty_rate(v2_rows, "source_version")
    fetch_date_rate = non_empty_rate(v2_rows, "fetch_date")

    backfilled_rows = 0
    bad_backfill_rows = 0
    changed_nonempty_smiles = 0
    for key, r2 in v2_by_key.items():
        r1 = v1_by_key.get(key)
        if r1 is None:
            continue

        s1 = normalize(r1.get("smiles", ""))
        s2 = normalize(r2.get("smiles", ""))
        flag = normalize(r2.get("smiles_backfill_flag", ""))
        src = normalize(r2.get("smiles_backfill_source", ""))
        source_chain = split_multi(r2.get("source", ""))

        if s1 and s2 and s1 != s2:
            changed_nonempty_smiles += 1

        if (not s1) and s2:
            backfilled_rows += 1
            if flag != "1" or src != BACKFILL_SOURCE_TAG or BACKFILL_SOURCE_TAG not in source_chain:
                bad_backfill_rows += 1

    coverage_delta_report = float(coverage.get("smiles_non_empty_delta", smiles_delta))

    conflict_metrics = conflict.get("metrics", {})
    conflict_has_structure = (
        isinstance(conflict.get("samples", {}), dict)
        and "inchikey_with_multi_zinc" in conflict.get("samples", {})
        and "zinc_with_multi_inchikey" in conflict.get("samples", {})
    )

    gates = [
        Gate(
            gate_id="zinc_mapping_equivalent_to_v1",
            passed=mapping_equivalent,
            detail={
                "allow_row_subset": allow_row_subset,
                "v1_rows": len(v1_rows),
                "v2_rows": len(v2_rows),
                "v1_unique_keys": len(v1_key_set),
                "v2_unique_keys": len(v2_key_set),
                "missing_in_v1": len(v2_key_set - v1_key_set),
                "missing_in_v2": len(v1_key_set - v2_key_set),
            },
        ),
        Gate(
            gate_id="format_integrity",
            passed=(bad_inchikey == 0 and bad_zinc == 0),
            detail={
                "bad_inchikey_rows": bad_inchikey,
                "bad_zinc_id_rows": bad_zinc,
            },
        ),
        Gate(
            gate_id="smiles_non_empty_rate_improved_ge_15pp",
            passed=(smiles_delta >= args.min_smiles_delta),
            detail={
                "baseline_smiles_rate": base_smiles_rate,
                "v2_smiles_rate": v2_smiles_rate,
                "delta": smiles_delta,
                "delta_pp": smiles_delta * 100.0,
                "min_required_delta": args.min_smiles_delta,
                "min_required_delta_pp": args.min_smiles_delta * 100.0,
            },
        ),
        Gate(
            gate_id="provenance_non_empty_rate_100pct",
            passed=(source_rate == 1.0 and source_version_rate == 1.0 and fetch_date_rate == 1.0),
            detail={
                "source_rate": source_rate,
                "source_version_rate": source_version_rate,
                "fetch_date_rate": fetch_date_rate,
            },
        ),
        Gate(
            gate_id="backfill_audit_integrity",
            passed=(bad_backfill_rows == 0 and changed_nonempty_smiles == 0),
            detail={
                "rows_backfilled": backfilled_rows,
                "bad_backfilled_rows": bad_backfill_rows,
                "changed_nonempty_smiles_rows": changed_nonempty_smiles,
            },
        ),
        Gate(
            gate_id="conflict_audit_preserved",
            passed=(isinstance(conflict_metrics, dict) and conflict_has_structure),
            detail={
                "conflict_metrics": conflict_metrics,
                "has_structured_samples": conflict_has_structure,
            },
        ),
        Gate(
            gate_id="report_consistency",
            passed=(
                int(build.get("metrics", {}).get("rows_written", len(v2_rows))) == len(v2_rows)
                and abs(float(build.get("metrics", {}).get("smiles_rate_delta", smiles_delta)) - smiles_delta) < 1e-9
                and abs(coverage_delta_report - smiles_delta) < 1e-9
            ),
            detail={
                "build_rows": int(build.get("metrics", {}).get("rows_written", len(v2_rows))),
                "actual_rows": len(v2_rows),
                "build_smiles_delta": float(build.get("metrics", {}).get("smiles_rate_delta", smiles_delta)),
                "qa_smiles_delta": smiles_delta,
                "coverage_smiles_delta": coverage_delta_report,
            },
        ),
    ]

    report = {
        "name": "molecule_zinc_xref_v2.qa",
        "created_at": utc_now(),
        "inputs": {
            "zinc_v1": str(args.zinc_v1),
            "table": str(args.table),
            "build_report": str(args.build_report),
            "coverage_report": str(args.coverage_report),
            "conflict_report": str(args.conflict_report),
        },
        "metrics": {
            "v1_rows": len(v1_rows),
            "v2_rows": len(v2_rows),
            "baseline_smiles_rate": base_smiles_rate,
            "v2_smiles_rate": v2_smiles_rate,
            "smiles_delta": smiles_delta,
            "smiles_delta_pp": smiles_delta * 100.0,
            "source_rate": source_rate,
            "source_version_rate": source_version_rate,
            "fetch_date_rate": fetch_date_rate,
            "rows_backfilled": backfilled_rows,
        },
        "gates": [{"id": g.gate_id, "passed": g.passed, "detail": g.detail} for g in gates],
        "passed": all(g.passed for g in gates),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if report["passed"] else "FAIL"
    print(f"[{status}] molecule_zinc_xref_v2 QA -> {args.out}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
