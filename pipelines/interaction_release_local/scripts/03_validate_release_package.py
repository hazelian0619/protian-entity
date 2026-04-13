#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List


def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def check(cond: bool, cid: str, detail: str) -> Dict[str, Any]:
    return {"id": cid, "passed": bool(cond), "detail": detail}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Validate interaction release package report")
    p.add_argument("--contract", type=Path, required=True)
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    c = load_json(args.contract)
    r = load_json(args.report)

    checks: List[Dict[str, Any]] = []

    for k in c.get("required_top_level_keys", []):
        checks.append(check(k in r, f"top_level.{k}", f"required top-level key: {k}"))

    req_checks = set(c.get("required_checks", []))
    got_checks = set((r.get("checks") or {}).keys()) if isinstance(r.get("checks"), dict) else set()
    checks.append(check(req_checks.issubset(got_checks), "checks.keys", f"required checks present: {sorted(req_checks)}"))

    req_tables = set(c.get("required_tables", []))
    got_tables = set((r.get("table_stats") or {}).keys()) if isinstance(r.get("table_stats"), dict) else set()
    checks.append(check(req_tables.issubset(got_tables), "table_stats.keys", f"required tables present: {sorted(req_tables)}"))

    if isinstance(r.get("table_stats"), dict):
        for t, obj in r["table_stats"].items():
            if isinstance(obj, dict):
                checks.append(check("path" in obj, f"{t}.path", "table stats contains path"))
                checks.append(check("rows" in obj, f"{t}.rows", "table stats contains rows"))

    passed = all(x["passed"] for x in checks)
    out = {
        "name": "interaction_release_local_package_validation",
        "passed": passed,
        "contract": str(args.contract),
        "report": str(args.report),
        "checks": checks,
    }
    write_json(args.out, out)
    print(f"[{'PASS' if passed else 'FAIL'}] release package validation -> {args.out}")
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
