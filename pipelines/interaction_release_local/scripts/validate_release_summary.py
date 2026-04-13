#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Validate interaction release summary json")
    p.add_argument("--contract", type=Path, required=True)
    p.add_argument("--summary", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    return p.parse_args()


def write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    c = json.loads(args.contract.read_text(encoding="utf-8"))
    s = json.loads(args.summary.read_text(encoding="utf-8"))

    checks: List[Dict[str, Any]] = []

    for k in c.get("required_top_level_keys", []):
        checks.append({"id": f"top.{k}", "passed": k in s, "detail": f"required key: {k}"})

    status_ok = s.get("status") in set(c.get("allowed_status", []))
    checks.append({"id": "status.allowed", "passed": status_ok, "detail": f"status={s.get('status')}"})

    required_checks = c.get("required_checks", [])
    got_checks = s.get("checks", {}) if isinstance(s.get("checks"), dict) else {}
    for k in required_checks:
        checks.append({"id": f"checks.{k}.exists", "passed": k in got_checks, "detail": f"check key: {k}"})

    tables_ok = isinstance(s.get("tables"), dict) and len(s.get("tables")) > 0
    qa_ok = isinstance(s.get("qa_reports"), dict) and len(s.get("qa_reports")) > 0
    checks.append({"id": "tables.non_empty", "passed": tables_ok, "detail": "tables must be non-empty dict"})
    checks.append({"id": "qa_reports.non_empty", "passed": qa_ok, "detail": "qa_reports must be non-empty dict"})

    passed = all(x["passed"] for x in checks)
    out = {
        "name": "interaction_release_local_summary_validation",
        "passed": passed,
        "contract": str(args.contract),
        "summary": str(args.summary),
        "checks": checks,
    }
    write_json(args.out, out)

    print(f"[{'PASS' if passed else 'FAIL'}] release summary validation -> {args.out}")
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
