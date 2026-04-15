#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Aggregate metrics + contract validations into final gate report")
    p.add_argument("--metrics", type=Path, required=True)
    p.add_argument("--validations", type=Path, nargs="+", required=True)
    p.add_argument("--out", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()

    metrics = load_json(args.metrics)
    validations: List[Dict] = [load_json(p) for p in args.validations]

    contract_checks = {str(p): bool(v.get("passed", False)) for p, v in zip(args.validations, validations)}

    metric_checks = metrics.get("gates", {}).get("checks", {})
    checks = {
        "all_three_interaction_types_present": bool(metric_checks.get("all_three_interaction_types_present", False)),
        "aggregate_score_not_all_zero_or_one": bool(metric_checks.get("aggregate_score_not_all_zero_or_one", False)),
        "aggregate_score_not_constant": bool(metric_checks.get("aggregate_score_not_constant", False)),
        "cross_table_edge_join_rate_ge_0_99": bool(metric_checks.get("cross_table_edge_join_rate_ge_0_99", False)),
        "contract_validation_pass": all(contract_checks.values()),
    }

    payload = {
        "pipeline": "interaction_cross_validation_v2",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "row_count": metrics.get("row_count", {}),
        "aggregate_score_distribution": metrics.get("distribution", {}).get("aggregate_score", {}),
        "contract_reports": contract_checks,
        "metrics_report": str(args.metrics),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"[{payload['status']}] gates -> {args.out}")
    return 0 if payload["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
