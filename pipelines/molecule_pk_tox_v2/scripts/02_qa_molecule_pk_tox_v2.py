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
FLOAT_RE = re.compile(r"^-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?$")
FIELDS = [
    "bioavailability",
    "clearance",
    "half_life",
    "volume_of_distribution",
    "ld50",
    "mutagenicity",
]
EVIDENCE_ALLOWED = {"structured/high", "text_mined/medium", "inferred/low", "none"}
MUT_LABEL_ALLOWED = {"positive", "negative", "mixed", "unknown"}
CLEARANCE_UNIT_ALLOWED = {"", "mL/min/kg", "mL/min", "L/h"}
HALF_LIFE_UNIT_ALLOWED = {"", "h"}
VOLUME_UNIT_ALLOWED = {"", "L/kg", "L", "L/1.73m2"}
LD50_UNIT_ALLOWED = {"", "mg/kg"}


@dataclass
class Gate:
    gate_id: str
    passed: bool
    detail: Dict[str, Any]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def norm(x: Any) -> str:
    return str(x or "").strip()


def load_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        return list(r)


def parse_float(x: str) -> Tuple[bool, float]:
    t = norm(x)
    if not t or not FLOAT_RE.match(t):
        return False, 0.0
    try:
        return True, float(t)
    except Exception:
        return False, 0.0


def legal_numeric_unit(value: str, unit: str, allowed_units: Set[str], min_positive: bool = True) -> bool:
    v = norm(value)
    u = norm(unit)
    if not v and not u:
        return True
    if (v and not u) or (u and not v):
        return False
    if u not in allowed_units:
        return False
    ok, fv = parse_float(v)
    if not ok:
        return False
    if min_positive and fv <= 0:
        return False
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xref", type=Path, required=True)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--build-report", type=Path, required=True)
    ap.add_argument("--validation-report", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--baseline-rows-with-any", type=int, default=2370)
    ap.add_argument("--min-normalized-legal-rate", type=float, default=0.98)
    args = ap.parse_args()

    for p in [args.xref, args.table, args.build_report, args.validation_report]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    xref_rows = load_tsv(args.xref)
    table_rows = load_tsv(args.table)
    build = json.loads(args.build_report.read_text(encoding="utf-8"))
    validation = json.loads(args.validation_report.read_text(encoding="utf-8"))

    xref_ik = {norm(r.get("inchikey", "")).upper() for r in xref_rows if norm(r.get("inchikey", ""))}

    bad_ik = 0
    dup_ik = 0
    seen: Set[str] = set()
    missing_backlink = 0
    evidence_violations = 0
    anchor_violations = 0
    field_non_empty = {f: 0 for f in FIELDS}
    rows_with_any = 0

    source_non_empty = 0
    source_version_non_empty = 0
    fetch_date_non_empty = 0

    normalized_checks = 0
    normalized_legal = 0
    normalized_by_field = {
        "bioavailability_pct": {"checked": 0, "legal": 0},
        "clearance": {"checked": 0, "legal": 0},
        "half_life": {"checked": 0, "legal": 0},
        "volume_distribution": {"checked": 0, "legal": 0},
        "ld50": {"checked": 0, "legal": 0},
        "mutagenicity": {"checked": 0, "legal": 0},
    }

    for r in table_rows:
        ik = norm(r.get("inchikey", "")).upper()
        if not INCHIKEY_RE.match(ik):
            bad_ik += 1
        if ik in seen:
            dup_ik += 1
        seen.add(ik)
        if ik not in xref_ik:
            missing_backlink += 1

        if norm(r.get("source", "")):
            source_non_empty += 1
        if norm(r.get("source_version", "")):
            source_version_non_empty += 1
        if norm(r.get("fetch_date", "")):
            fetch_date_non_empty += 1

        any_field = False
        for f in FIELDS:
            v = norm(r.get(f, ""))
            ev = norm(r.get(f"{f}_evidence", "none"))
            src = norm(r.get(f"{f}_source", "none"))
            if ev not in EVIDENCE_ALLOWED:
                evidence_violations += 1
            if v:
                any_field = True
                field_non_empty[f] += 1
                if ev == "none" or src in {"", "none"}:
                    anchor_violations += 1

        if any_field:
            rows_with_any += 1

        # normalized legality checks
        normalized_checks += 1
        normalized_by_field["bioavailability_pct"]["checked"] += 1
        b = norm(r.get("bioavailability_pct", ""))
        ok_b = True
        if b:
            ok, fv = parse_float(b)
            ok_b = ok and (0.0 <= fv <= 1000.0)
        if ok_b:
            normalized_legal += 1
            normalized_by_field["bioavailability_pct"]["legal"] += 1

        normalized_checks += 1
        normalized_by_field["clearance"]["checked"] += 1
        ok_c = legal_numeric_unit(
            r.get("clearance_value", ""),
            r.get("clearance_unit_std", ""),
            CLEARANCE_UNIT_ALLOWED,
            min_positive=True,
        )
        if ok_c:
            normalized_legal += 1
            normalized_by_field["clearance"]["legal"] += 1

        normalized_checks += 1
        normalized_by_field["half_life"]["checked"] += 1
        ok_h = legal_numeric_unit(
            r.get("half_life_value", ""),
            r.get("half_life_unit_std", ""),
            HALF_LIFE_UNIT_ALLOWED,
            min_positive=True,
        )
        if ok_h:
            normalized_legal += 1
            normalized_by_field["half_life"]["legal"] += 1

        normalized_checks += 1
        normalized_by_field["volume_distribution"]["checked"] += 1
        ok_v = legal_numeric_unit(
            r.get("volume_distribution_value", ""),
            r.get("volume_distribution_unit_std", ""),
            VOLUME_UNIT_ALLOWED,
            min_positive=True,
        )
        if ok_v:
            normalized_legal += 1
            normalized_by_field["volume_distribution"]["legal"] += 1

        normalized_checks += 1
        normalized_by_field["ld50"]["checked"] += 1
        ok_l = legal_numeric_unit(
            r.get("ld50_value", ""),
            r.get("ld50_unit_std", ""),
            LD50_UNIT_ALLOWED,
            min_positive=True,
        )
        if ok_l:
            normalized_legal += 1
            normalized_by_field["ld50"]["legal"] += 1

        normalized_checks += 1
        normalized_by_field["mutagenicity"]["checked"] += 1
        ok_m = norm(r.get("mutagenicity_label_std", "unknown")) in MUT_LABEL_ALLOWED
        if ok_m:
            normalized_legal += 1
            normalized_by_field["mutagenicity"]["legal"] += 1

    rows = len(table_rows)
    backlink_rate = ((rows - missing_backlink) / rows) if rows else 0.0
    source_rate = source_non_empty / rows if rows else 0.0
    source_version_rate = source_version_non_empty / rows if rows else 0.0
    fetch_rate = fetch_date_non_empty / rows if rows else 0.0

    normalized_legal_rate = (normalized_legal / normalized_checks) if normalized_checks else 0.0
    normalized_legal_rate_by_field = {
        k: (v["legal"] / v["checked"] if v["checked"] else 0.0) for k, v in normalized_by_field.items()
    }

    build_rows_with_any = int(build.get("metrics", {}).get("rows_with_any_pk_tox", rows_with_any))
    second_source_contrib = int(build.get("metrics", {}).get("rows_with_second_source_contrib", 0))
    validation_passed = bool(validation.get("passed"))

    gates = [
        Gate(
            gate_id="inchikey_unique_and_valid",
            passed=(bad_ik == 0 and dup_ik == 0),
            detail={"bad_inchikey": bad_ik, "duplicate_inchikey": dup_ik},
        ),
        Gate(
            gate_id="xref_backlink_rate",
            passed=(backlink_rate >= 0.999),
            detail={
                "backlink_rate": backlink_rate,
                "min_required": 0.999,
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
            gate_id="validation_pass",
            passed=validation_passed,
            detail={"validation_passed": validation_passed},
        ),
        Gate(
            gate_id="v2_non_degradation_rows_with_any",
            passed=(rows_with_any > args.baseline_rows_with_any),
            detail={
                "rows_with_any_pk_tox": rows_with_any,
                "baseline_rows_with_any_pk_tox": args.baseline_rows_with_any,
                "operator": ">",
                "build_report_rows_with_any": build_rows_with_any,
            },
        ),
        Gate(
            gate_id="second_source_contribution_positive",
            passed=(second_source_contrib > 0),
            detail={"rows_with_second_source_contrib": second_source_contrib, "operator": ">"},
        ),
        Gate(
            gate_id="normalized_fields_legal_rate",
            passed=(normalized_legal_rate >= args.min_normalized_legal_rate),
            detail={
                "normalized_legal_rate": normalized_legal_rate,
                "min_required": args.min_normalized_legal_rate,
                "by_field": normalized_legal_rate_by_field,
            },
        ),
        Gate(
            gate_id="source_triplet_complete",
            passed=(source_rate == 1.0 and source_version_rate == 1.0 and fetch_rate == 1.0),
            detail={
                "source_non_empty_rate": source_rate,
                "source_version_non_empty_rate": source_version_rate,
                "fetch_date_non_empty_rate": fetch_rate,
                "required": 1.0,
            },
        ),
    ]

    out = {
        "name": "molecule_pk_tox_v2.qa",
        "created_at": utc_now(),
        "inputs": {
            "xref": str(args.xref),
            "table": str(args.table),
            "build_report": str(args.build_report),
            "validation_report": str(args.validation_report),
            "baseline_rows_with_any": args.baseline_rows_with_any,
            "min_normalized_legal_rate": args.min_normalized_legal_rate,
        },
        "metrics": {
            "rows": rows,
            "xref_rows": len(xref_rows),
            "rows_with_any_pk_tox": rows_with_any,
            "field_non_empty_counts": field_non_empty,
            "backlink_rate": backlink_rate,
            "source_non_empty_rate": source_rate,
            "source_version_non_empty_rate": source_version_rate,
            "fetch_date_non_empty_rate": fetch_rate,
            "normalized_legal_rate": normalized_legal_rate,
            "normalized_legal_rate_by_field": normalized_legal_rate_by_field,
            "rows_with_second_source_contrib": second_source_contrib,
            "anchor_violations": anchor_violations,
            "evidence_violations": evidence_violations,
        },
        "gates": [{"id": g.gate_id, "passed": g.passed, "detail": g.detail} for g in gates],
        "passed": all(g.passed for g in gates),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if out["passed"] else "FAIL"
    print(f"[{status}] molecule_pk_tox_v2 QA -> {args.out}")
    for g in gates:
        print(f"  - {'PASS' if g.passed else 'FAIL'} {g.gate_id}")
    return 0 if out["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
