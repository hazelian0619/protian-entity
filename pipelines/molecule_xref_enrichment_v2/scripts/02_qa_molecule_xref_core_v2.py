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
CHEMBL_TOKEN_RE = re.compile(r"^CHEMBL\d+$")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize(x: str) -> str:
    return (x or "").strip()


def split_multi(x: str) -> List[str]:
    t = normalize(x)
    if not t:
        return []
    return [p.strip() for p in t.split(";") if p.strip()]


@dataclass
class Gate:
    gate_id: str
    passed: bool
    detail: Dict[str, Any]


def _load_table(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        return list(r)


def summarize(rows: List[Dict[str, str]]) -> Dict[str, Any]:
    n = len(rows)
    ik_seen: Set[str] = set()
    ik_dup = 0
    bad_ik = 0
    bad_chembl_token = 0

    chembl_rows = 0
    pubchem_rows = 0

    dbids: Set[str] = set()
    dbids_chembl: Set[str] = set()
    dbids_pubchem: Set[str] = set()

    backfill_rows = 0
    backfill_rows_missing_conf = 0

    for row in rows:
        ik = normalize(row.get("inchikey", "")).upper()
        if ik in ik_seen:
            ik_dup += 1
        ik_seen.add(ik)
        if not INCHIKEY_RE.match(ik):
            bad_ik += 1

        chembl_ids = [x.upper() for x in split_multi(row.get("chembl_id", ""))]
        pubchem_ids = split_multi(row.get("pubchem_cid", ""))
        dbid_tokens = split_multi(row.get("drugbank_id", ""))

        if chembl_ids:
            chembl_rows += 1
        if pubchem_ids:
            pubchem_rows += 1

        dbids.update(dbid_tokens)
        if chembl_ids:
            dbids_chembl.update(dbid_tokens)
        if pubchem_ids:
            dbids_pubchem.update(dbid_tokens)

        for c in chembl_ids:
            if not CHEMBL_TOKEN_RE.match(c):
                bad_chembl_token += 1

        strategies = set(split_multi(row.get("match_strategy", "")))
        conf = normalize(row.get("confidence", ""))
        if any(s.startswith("backfill_") for s in strategies):
            backfill_rows += 1
            if conf not in {"high", "medium", "low"}:
                backfill_rows_missing_conf += 1

    return {
        "rows": n,
        "inchikey_duplicates": ik_dup,
        "inchikey_bad_format": bad_ik,
        "chembl_bad_token": bad_chembl_token,
        "rows_with_chembl": chembl_rows,
        "rows_with_pubchem": pubchem_rows,
        "chembl_row_rate": (chembl_rows / n) if n else 0.0,
        "pubchem_row_rate": (pubchem_rows / n) if n else 0.0,
        "drugbank_ids": len(dbids),
        "drugbank_ids_with_chembl": len(dbids_chembl),
        "drugbank_ids_with_pubchem": len(dbids_pubchem),
        "drugbank_chembl_rate": (len(dbids_chembl) / len(dbids)) if dbids else 0.0,
        "drugbank_pubchem_rate": (len(dbids_pubchem) / len(dbids)) if dbids else 0.0,
        "backfill_rows": backfill_rows,
        "backfill_rows_missing_conf": backfill_rows_missing_conf,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--v1-table", type=Path, required=True)
    ap.add_argument("--v2-table", type=Path, required=True)
    ap.add_argument("--backfill-audit", type=Path, required=True)
    ap.add_argument("--build-report", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--require-positive-delta", type=int, default=1)
    args = ap.parse_args()

    for p in [args.v1_table, args.v2_table, args.backfill_audit, args.build_report]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    v1_rows = _load_table(args.v1_table)
    v2_rows = _load_table(args.v2_table)
    v1 = summarize(v1_rows)
    v2 = summarize(v2_rows)

    backfill = json.loads(args.backfill_audit.read_text(encoding="utf-8"))
    build = json.loads(args.build_report.read_text(encoding="utf-8"))

    accepted_records = backfill.get("records", [])
    accepted_count = int(backfill.get("accepted_count", len(accepted_records)))
    accepted_with_bad_fields = 0
    for rec in accepted_records:
        if not normalize(rec.get("match_strategy", "")) or not normalize(rec.get("confidence", "")):
            accepted_with_bad_fields += 1

    chembl_delta = v2["rows_with_chembl"] - v1["rows_with_chembl"]
    pubchem_delta = v2["rows_with_pubchem"] - v1["rows_with_pubchem"]
    positive_required = bool(args.require_positive_delta)

    build_consistency_pass = (
        True
        if not positive_required
        else (
            int(build.get("metrics", {}).get("chembl_rows_delta_abs", chembl_delta)) == chembl_delta
            and int(build.get("metrics", {}).get("pubchem_rows_delta_abs", pubchem_delta)) == pubchem_delta
        )
    )

    gates = [
        Gate(
            gate_id="inchikey_unique_100pct",
            passed=(v2["inchikey_duplicates"] == 0 and v2["inchikey_bad_format"] == 0),
            detail={
                "inchikey_duplicates": v2["inchikey_duplicates"],
                "inchikey_bad_format": v2["inchikey_bad_format"],
            },
        ),
        Gate(
            gate_id="chembl_coverage_delta_positive",
            passed=(chembl_delta > 0) if positive_required else True,
            detail={
                "v1_rows_with_chembl": v1["rows_with_chembl"],
                "v2_rows_with_chembl": v2["rows_with_chembl"],
                "delta_abs": chembl_delta,
                "v1_rate": v1["chembl_row_rate"],
                "v2_rate": v2["chembl_row_rate"],
                "delta_rate_abs": v2["chembl_row_rate"] - v1["chembl_row_rate"],
                "required_positive": positive_required,
            },
        ),
        Gate(
            gate_id="pubchem_coverage_delta_positive",
            passed=(pubchem_delta > 0) if positive_required else True,
            detail={
                "v1_rows_with_pubchem": v1["rows_with_pubchem"],
                "v2_rows_with_pubchem": v2["rows_with_pubchem"],
                "delta_abs": pubchem_delta,
                "v1_rate": v1["pubchem_row_rate"],
                "v2_rate": v2["pubchem_row_rate"],
                "delta_rate_abs": v2["pubchem_row_rate"] - v1["pubchem_row_rate"],
                "required_positive": positive_required,
            },
        ),
        Gate(
            gate_id="backfill_records_have_strategy_and_confidence",
            passed=(accepted_with_bad_fields == 0 and accepted_count == len(accepted_records)),
            detail={
                "accepted_count": accepted_count,
                "records_length": len(accepted_records),
                "records_missing_fields": accepted_with_bad_fields,
            },
        ),
        Gate(
            gate_id="backfill_rows_in_v2_have_confidence",
            passed=(v2["backfill_rows_missing_conf"] == 0),
            detail={
                "backfill_rows": v2["backfill_rows"],
                "backfill_rows_missing_conf": v2["backfill_rows_missing_conf"],
            },
        ),
        Gate(
            gate_id="chembl_token_format_valid",
            passed=(v2["chembl_bad_token"] == 0),
            detail={"chembl_bad_token": v2["chembl_bad_token"]},
        ),
        Gate(
            gate_id="build_report_consistency",
            passed=build_consistency_pass,
            detail={
                "build_chembl_delta": int(build.get("metrics", {}).get("chembl_rows_delta_abs", chembl_delta)),
                "qa_chembl_delta": chembl_delta,
                "build_pubchem_delta": int(build.get("metrics", {}).get("pubchem_rows_delta_abs", pubchem_delta)),
                "qa_pubchem_delta": pubchem_delta,
                "skipped_in_smoke_mode": (not positive_required),
            },
        ),
    ]

    report = {
        "name": "molecule_xref_core_v2.qa",
        "created_at": utc_now(),
        "inputs": {
            "v1_table": str(args.v1_table),
            "v2_table": str(args.v2_table),
            "backfill_audit": str(args.backfill_audit),
            "build_report": str(args.build_report),
        },
        "require_positive_delta": positive_required,
        "metrics": {
            "v1": v1,
            "v2": v2,
            "chembl_delta_abs": chembl_delta,
            "pubchem_delta_abs": pubchem_delta,
            "chembl_delta_rate_abs": v2["chembl_row_rate"] - v1["chembl_row_rate"],
            "pubchem_delta_rate_abs": v2["pubchem_row_rate"] - v1["pubchem_row_rate"],
        },
        "gates": [{"id": g.gate_id, "passed": g.passed, "detail": g.detail} for g in gates],
        "passed": all(g.passed for g in gates),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if report["passed"] else "FAIL"
    print(f"[{status}] molecule_xref_core_v2 QA -> {args.out}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
