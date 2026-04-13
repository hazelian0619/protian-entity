#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Set


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize(x: str) -> str:
    return (x or "").strip()


def quantile(values: List[float], q: float) -> float:
    if not values:
        return 0.0
    s = sorted(values)
    idx = int(round((len(s) - 1) * q))
    idx = max(0, min(idx, len(s) - 1))
    return float(s[idx])


def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cross", type=Path, required=True)
    ap.add_argument("--aggregate", type=Path, required=True)
    ap.add_argument("--build-report", type=Path, required=True)
    ap.add_argument("--cross-validation", type=Path, required=True)
    ap.add_argument("--aggregate-validation", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-spread", type=float, default=0.05)
    args = ap.parse_args()

    for p in [
        args.cross,
        args.aggregate,
        args.build_report,
        args.cross_validation,
        args.aggregate_validation,
    ]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    build = load_json(args.build_report)
    val_cross = load_json(args.cross_validation)
    val_agg = load_json(args.aggregate_validation)

    types_in_cross: Set[str] = set()
    cross_edge_ids: Set[str] = set()
    agg_edge_ids: Set[str] = set()
    agg_scores: List[float] = []
    agg_bucket = {"high": 0, "medium": 0, "low": 0}

    with args.cross.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        need = {"interaction_type", "edge_id", "conflict_flag", "consistent_across_n"}
        miss = need - set(r.fieldnames or [])
        if miss:
            raise SystemExit(f"[ERROR] cross table missing columns: {sorted(miss)}")
        for row in r:
            types_in_cross.add(normalize(row.get("interaction_type", "")))
            eid = normalize(row.get("edge_id", ""))
            if eid:
                cross_edge_ids.add(eid)

    with args.aggregate.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        need = {"interaction_type", "edge_id", "aggregate_score", "score_bucket"}
        miss = need - set(r.fieldnames or [])
        if miss:
            raise SystemExit(f"[ERROR] aggregate table missing columns: {sorted(miss)}")
        for row in r:
            eid = normalize(row.get("edge_id", ""))
            if eid:
                agg_edge_ids.add(eid)
            try:
                agg_scores.append(float(normalize(row.get("aggregate_score", "0")) or 0.0))
            except ValueError:
                agg_scores.append(0.0)
            b = normalize(row.get("score_bucket", ""))
            if b in agg_bucket:
                agg_bucket[b] += 1

    score_min = min(agg_scores) if agg_scores else 0.0
    score_max = max(agg_scores) if agg_scores else 0.0
    score_p10 = quantile(agg_scores, 0.10) if agg_scores else 0.0
    score_p90 = quantile(agg_scores, 0.90) if agg_scores else 0.0

    join_rate = (len(cross_edge_ids & agg_edge_ids) / len(cross_edge_ids)) if cross_edge_ids else 0.0

    checks = {
        "all_three_interaction_types_present": {"PPI", "PSI", "RPI"}.issubset(types_in_cross),
        "aggregate_score_not_all_zero_or_one": (score_min < 0.999) and (score_max > 0.001),
        "aggregate_score_not_constant": (score_p90 - score_p10) >= args.min_spread,
        "cross_aggregate_join_rate_ge_0_99": join_rate >= 0.99,
        "contract_cross_pass": bool(val_cross.get("passed", False)),
        "contract_aggregate_pass": bool(val_agg.get("passed", False)),
    }

    report = {
        "name": "interaction_cross_validation_v2.qa",
        "created_at": utc_now(),
        "inputs": {
            "cross": str(args.cross),
            "aggregate": str(args.aggregate),
            "build_report": str(args.build_report),
            "cross_validation": str(args.cross_validation),
            "aggregate_validation": str(args.aggregate_validation),
        },
        "metrics": {
            "cross_rows": len(cross_edge_ids),
            "aggregate_rows": len(agg_edge_ids),
            "types_in_cross": sorted(types_in_cross),
            "cross_aggregate_join_rate": join_rate,
            "aggregate_score": {
                "min": score_min,
                "p10": score_p10,
                "p90": score_p90,
                "max": score_max,
            },
            "score_bucket": agg_bucket,
            "build_summary": build.get("distribution", {}).get("aggregate_score", {}),
            "psi_b_update_integration": build.get("psi_b_update_integration", {}),
        },
        "gates": [{"id": k, "passed": bool(v)} for k, v in checks.items()],
        "passed": all(checks.values()),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if report["passed"] else "FAIL"
    print(f"[{status}] interaction_cross_validation_v2 QA -> {args.out}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
