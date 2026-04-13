#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple

INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize(x: str) -> str:
    return (x or "").strip()


def load_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        return list(reader)


def non_empty_rate(rows: List[Dict[str, str]], column: str) -> float:
    if not rows:
        return 0.0
    return sum(1 for r in rows if normalize(r.get(column, ""))) / len(rows)


def unique_ik_stats(rows: List[Dict[str, str]]) -> Tuple[int, int]:
    seen = set()
    dup = 0
    for row in rows:
        ik = normalize(row.get("inchikey", "")).upper()
        if ik in seen:
            dup += 1
        else:
            seen.add(ik)
    return len(seen), dup


@dataclass
class Gate:
    gate_id: str
    passed: bool
    detail: Dict[str, Any]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--core-tsv", type=Path, required=True)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--build-report", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-inchi-rate", type=float, default=0.70)
    ap.add_argument("--min-inchi-mappable-rate", type=float, default=0.90)
    ap.add_argument("--min-smiles-rate", type=float, default=0.75)
    ap.add_argument("--allow-row-subset", type=int, default=0)
    args = ap.parse_args()

    for p in [args.core_tsv, args.table, args.build_report]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    allow_row_subset = bool(args.allow_row_subset)
    core_rows = load_tsv(args.core_tsv)
    table_rows = load_tsv(args.table)
    build = json.loads(args.build_report.read_text(encoding="utf-8"))

    core_iks = {normalize(r.get("inchikey", "")).upper() for r in core_rows if normalize(r.get("inchikey", ""))}
    table_iks = [normalize(r.get("inchikey", "")).upper() for r in table_rows if normalize(r.get("inchikey", ""))]
    table_ik_set = set(table_iks)

    valid_ik_rate = (
        sum(1 for ik in table_iks if bool(INCHIKEY_RE.match(ik))) / len(table_iks)
        if table_iks
        else 0.0
    )

    backlink_matched = len(table_ik_set.intersection(core_iks))
    backlink_rate = (backlink_matched / len(table_ik_set)) if table_ik_set else 0.0

    inchi_rate = non_empty_rate(table_rows, "inchi")
    smiles_rate = non_empty_rate(table_rows, "canonical_smiles")
    source_rate = non_empty_rate(table_rows, "source")
    source_version_rate = non_empty_rate(table_rows, "source_version")
    fetch_date_rate = non_empty_rate(table_rows, "fetch_date")

    mappable_rows = [r for r in table_rows if normalize(r.get("chembl_id", "")) or normalize(r.get("pubchem_cid", ""))]
    inchi_mappable_rate = non_empty_rate(mappable_rows, "inchi") if mappable_rows else 0.0

    uniq_count, dup_count = unique_ik_stats(table_rows)

    if allow_row_subset:
        backlink_gate = len(table_rows) <= len(core_rows) and table_ik_set.issubset(core_iks)
    else:
        backlink_gate = len(table_rows) == len(core_rows) and table_ik_set == core_iks

    inchi_gate = (inchi_rate >= args.min_inchi_rate) or (
        len(mappable_rows) > 0 and inchi_mappable_rate >= args.min_inchi_mappable_rate
    )

    build_inchi = float(build.get("metrics", {}).get("inchi_non_empty_rate", inchi_rate))
    build_smiles = float(build.get("metrics", {}).get("canonical_smiles_non_empty_rate", smiles_rate))

    gates = [
        Gate(
            gate_id="row_count_and_inchikey_backlink",
            passed=backlink_gate,
            detail={
                "core_rows": len(core_rows),
                "table_rows": len(table_rows),
                "allow_row_subset": allow_row_subset,
                "backlink_rate": backlink_rate,
                "table_keys_in_core": len(table_ik_set.intersection(core_iks)),
                "table_unique_keys": len(table_ik_set),
            },
        ),
        Gate(
            gate_id="inchikey_unique_and_valid_format",
            passed=(dup_count == 0 and valid_ik_rate == 1.0),
            detail={
                "unique_inchikey_count": uniq_count,
                "duplicate_inchikey_count": dup_count,
                "inchikey_format_rate": valid_ik_rate,
            },
        ),
        Gate(
            gate_id="inchi_coverage_acceptance",
            passed=inchi_gate,
            detail={
                "inchi_rate_all": inchi_rate,
                "min_inchi_rate_all": args.min_inchi_rate,
                "mappable_subset_rows": len(mappable_rows),
                "inchi_rate_mappable_subset": inchi_mappable_rate,
                "min_inchi_rate_mappable_subset": args.min_inchi_mappable_rate,
                "acceptance_logic": "all>=0.70 OR mappable_subset>=0.90",
            },
        ),
        Gate(
            gate_id="canonical_smiles_coverage_acceptance",
            passed=(smiles_rate >= args.min_smiles_rate),
            detail={
                "canonical_smiles_rate": smiles_rate,
                "min_canonical_smiles_rate": args.min_smiles_rate,
            },
        ),
        Gate(
            gate_id="provenance_columns_non_empty",
            passed=(source_rate == 1.0 and source_version_rate == 1.0 and fetch_date_rate == 1.0),
            detail={
                "source_non_empty_rate": source_rate,
                "source_version_non_empty_rate": source_version_rate,
                "fetch_date_non_empty_rate": fetch_date_rate,
            },
        ),
        Gate(
            gate_id="build_report_consistency",
            passed=(abs(build_inchi - inchi_rate) < 1e-9 and abs(build_smiles - smiles_rate) < 1e-9),
            detail={
                "qa_inchi_rate": inchi_rate,
                "build_inchi_rate": build_inchi,
                "qa_smiles_rate": smiles_rate,
                "build_smiles_rate": build_smiles,
            },
        ),
    ]

    report = {
        "name": "molecule_structure_identifiers_v1.qa",
        "created_at": utc_now(),
        "inputs": {
            "core_tsv": str(args.core_tsv),
            "table": str(args.table),
            "build_report": str(args.build_report),
        },
        "metrics": {
            "core_rows": len(core_rows),
            "table_rows": len(table_rows),
            "inchikey_backlink_rate": backlink_rate,
            "inchi_non_empty_rate": inchi_rate,
            "inchi_non_empty_rate_mappable_subset": inchi_mappable_rate,
            "canonical_smiles_non_empty_rate": smiles_rate,
            "source_non_empty_rate": source_rate,
            "source_version_non_empty_rate": source_version_rate,
            "fetch_date_non_empty_rate": fetch_date_rate,
            "mappable_subset_rows": len(mappable_rows),
        },
        "gates": [{"id": g.gate_id, "passed": g.passed, "detail": g.detail} for g in gates],
        "passed": all(g.passed for g in gates),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if report["passed"] else "FAIL"
    print(f"[{status}] molecule_structure_identifiers_v1 QA -> {args.out}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
