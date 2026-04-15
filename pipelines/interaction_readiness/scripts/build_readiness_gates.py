#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build compact gates summary for interaction readiness")
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    return p.parse_args()


def write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    report = json.loads(args.report.read_text(encoding="utf-8"))
    gates = report.get("preupload_gates", [])
    passed = all(bool(g.get("passed")) for g in gates) if gates else False

    payload = {
        "pipeline": "interaction_readiness",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if passed and report.get("overall_status") == "PASS" else "FAIL",
        "overall_status": report.get("overall_status"),
        "gate_pass_count": sum(1 for g in gates if g.get("passed")),
        "gate_total": len(gates),
        "blocking_gap_count": report.get("summary", {}).get("blocking_gap_count", None),
        "gates": gates,
    }
    write_json(args.out, payload)
    print(f"[OK] readiness gates -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
