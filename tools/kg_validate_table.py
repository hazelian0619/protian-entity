#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple


@dataclass
class RuleResult:
    rule_id: str
    rule_type: str
    passed: bool
    checked: int
    passed_count: int
    rate: Optional[float]
    detail: Dict[str, Any]


def _load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _where_matches(where: Optional[Dict[str, Any]], row: Dict[str, str]) -> bool:
    if not where:
        return True
    col = where.get("column")
    if col is None:
        return True
    if "equals" in where:
        return row.get(col, "") == str(where["equals"])
    raise ValueError(f"Unsupported where: {where}")


def _is_non_empty(value: str) -> bool:
    v = (value or "").strip()
    return v not in ("", "NA", "N/A", "None", "null")


def _iter_rows(tsv: Path) -> Tuple[Iterable[Dict[str, str]], list[str]]:
    f = tsv.open("r", encoding="utf-8", newline="")
    reader = csv.DictReader(f, delimiter="\t")
    if reader.fieldnames is None:
        raise SystemExit("Missing header")
    fieldnames = list(reader.fieldnames)

    def gen():
        try:
            for row in reader:
                yield row
        finally:
            f.close()

    return gen(), fieldnames


def validate(tsv: Path, contract: Dict[str, Any]) -> Dict[str, Any]:
    rows, fieldnames = _iter_rows(tsv)

    required_columns = contract.get("required_columns", [])
    missing = [c for c in required_columns if c not in fieldnames]
    schema_ok = len(missing) == 0

    # Prepare rules
    rules = contract.get("rules", [])
    compiled: list[Dict[str, Any]] = []
    for rule in rules:
        r = dict(rule)
        if r.get("type") == "regex_rate":
            r["_re"] = re.compile(r["pattern"])
        if r.get("type") == "sequence_charset_rate":
            r["_allowed"] = set(r.get("allowed_chars", []))
        compiled.append(r)

    # State for rules
    stats: Dict[str, Dict[str, Any]] = {}
    seen_for_unique: Dict[str, set[str]] = {}

    for rule in compiled:
        rid = rule["id"]
        rtype = rule["type"]
        stats[rid] = {
            "rule_id": rid,
            "rule_type": rtype,
            "checked": 0,
            "passed": 0,
            "detail": {},
        }
        if rtype == "unique":
            seen_for_unique[rid] = set()
            stats[rid]["detail"]["duplicates"] = 0

    total_rows = 0

    for row in rows:
        total_rows += 1
        for rule in compiled:
            rid = rule["id"]
            rtype = rule["type"]

            if not _where_matches(rule.get("where"), row):
                continue

            col = rule.get("column")
            if col and col not in row:
                continue

            stats[rid]["checked"] += 1

            value = row.get(col, "") if col else ""

            ok = True
            if rtype == "non_empty_rate":
                ok = _is_non_empty(value)
            elif rtype == "equals_value_rate":
                ok = str(value) == str(rule["value"])
            elif rtype == "allowed_values_rate":
                ok = str(value) in set(str(x) for x in rule.get("allowed", []))
            elif rtype == "regex_rate":
                ok = bool(rule["_re"].match(str(value)))
            elif rtype == "sequence_charset_rate":
                seq = (value or "").strip().upper()
                ok = len(seq) > 0 and set(seq) <= rule["_allowed"]
            elif rtype == "unique":
                v = str(value)
                if v in seen_for_unique[rid]:
                    ok = False
                    stats[rid]["detail"]["duplicates"] += 1
                else:
                    seen_for_unique[rid].add(v)
                    ok = True
            else:
                raise ValueError(f"Unsupported rule type: {rtype}")

            if ok:
                stats[rid]["passed"] += 1

    # Finalize
    results: list[RuleResult] = []
    passed_all = schema_ok

    for rule in compiled:
        rid = rule["id"]
        checked = int(stats[rid]["checked"])
        passed_cnt = int(stats[rid]["passed"])
        rate: Optional[float]
        if rule["type"] == "unique":
            rate = 1.0 if stats[rid]["detail"].get("duplicates", 0) == 0 else (passed_cnt / checked if checked else None)
        else:
            rate = (passed_cnt / checked) if checked else None

        min_rate = rule.get("min_rate")
        ok = True
        if min_rate is not None:
            ok = rate is not None and rate >= float(min_rate)
        if rule["type"] == "unique":
            ok = stats[rid]["detail"].get("duplicates", 0) == 0

        passed_all = passed_all and ok
        results.append(
            RuleResult(
                rule_id=rid,
                rule_type=rule["type"],
                passed=ok,
                checked=checked,
                passed_count=passed_cnt,
                rate=rate,
                detail=stats[rid]["detail"],
            )
        )

    out = {
        "table": str(tsv),
        "contract": contract.get("name"),
        "row_count": total_rows,
        "columns": fieldnames,
        "schema": {
            "required_columns": required_columns,
            "missing_required": missing,
            "passed": schema_ok,
        },
        "rules": [
            {
                "id": r.rule_id,
                "type": r.rule_type,
                "passed": r.passed,
                "checked": r.checked,
                "passed_count": r.passed_count,
                "rate": r.rate,
                **({"detail": r.detail} if r.detail else {}),
            }
            for r in results
        ],
        "passed": passed_all,
    }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contract", type=Path, required=True)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    contract = _load_json(args.contract)
    report = validate(args.table, contract)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if report["passed"] else "FAIL"
    print(f"[{status}] {contract.get('name')} -> {args.table} (rows={report['row_count']})")
    for r in report["rules"]:
        if r["type"] == "unique":
            dup = r.get("detail", {}).get("duplicates", 0)
            if dup:
                print(f"  [FAIL] {r['id']}: duplicates={dup}")
            else:
                print(f"  [OK]   {r['id']}: unique")
        else:
            rate = r.get("rate")
            min_rate = None
            for rr in contract.get("rules", []):
                if rr.get("id") == r["id"]:
                    min_rate = rr.get("min_rate")
                    break
            if rate is None:
                print(f"  [WARN] {r['id']}: checked=0")
            else:
                ok = r["passed"]
                tag = "OK" if ok else "FAIL"
                exp = f">={min_rate}" if min_rate is not None else ""
                print(f"  [{tag}] {r['id']}: rate={rate:.4f} {exp}")

    return 0 if report["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
