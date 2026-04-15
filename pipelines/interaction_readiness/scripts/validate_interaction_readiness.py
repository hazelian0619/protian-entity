#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Validate interaction readiness JSON report")
    p.add_argument("--contract", type=Path, required=True)
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    return p.parse_args()


def write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def check(cond: bool, cid: str, detail: str) -> Dict[str, Any]:
    return {"id": cid, "passed": bool(cond), "detail": detail}


def main() -> int:
    args = parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    report = json.loads(args.report.read_text(encoding="utf-8"))

    checks: List[Dict[str, Any]] = []

    for k in contract.get("required_top_level_keys", []):
        checks.append(check(k in report, f"top_level.{k}", f"required key present: {k}"))

    required_interactions = set(contract.get("required_interactions", []))
    interaction_map = {x.get("interaction"): x for x in report.get("interactions", []) if isinstance(x, dict)}
    checks.append(
        check(
            required_interactions.issubset(set(interaction_map.keys())),
            "interaction.coverage",
            f"required interactions present: {sorted(required_interactions)}",
        )
    )

    required_dimensions = contract.get("required_dimensions", [])
    allowed_status = set(contract.get("allowed_status_values", []))

    for i_name in sorted(required_interactions):
        item = interaction_map.get(i_name, {})
        dim = item.get("dimension_status", {}) if isinstance(item, dict) else {}
        for d in required_dimensions:
            checks.append(check(d in dim, f"{i_name}.dimension.{d}", f"dimension key exists: {d}"))
            if d in dim:
                checks.append(
                    check(
                        str(dim.get(d)) in allowed_status,
                        f"{i_name}.dimension.{d}.status_value",
                        f"status in allowed set: {dim.get(d)}",
                    )
                )

    gates = report.get("preupload_gates", [])
    required_gate_keys = set(contract.get("required_gate_keys", []))
    checks.append(check(isinstance(gates, list) and len(gates) > 0, "gates.non_empty", "preupload_gates must be non-empty list"))
    if isinstance(gates, list):
        for idx, gate in enumerate(gates):
            if not isinstance(gate, dict):
                checks.append(check(False, f"gates.{idx}.dict", "gate item must be object"))
                continue
            miss = [k for k in required_gate_keys if k not in gate]
            checks.append(check(len(miss) == 0, f"gates.{idx}.keys", f"missing gate keys: {miss}"))

    top10 = report.get("top10_gaps", [])
    checks.append(check(isinstance(top10, list), "top10.type", "top10_gaps must be list"))
    if isinstance(top10, list):
        checks.append(check(len(top10) <= 10, "top10.max_len", f"top10_gaps length <=10, got {len(top10)}"))

    passed = all(c["passed"] for c in checks)
    output = {
        "name": "interaction_l2_readiness_validation",
        "passed": passed,
        "report": str(args.report),
        "contract": str(args.contract),
        "checks": checks,
    }
    write_json(args.out, output)

    print(f"[{'PASS' if passed else 'FAIL'}] interaction readiness validation -> {args.out}")
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
