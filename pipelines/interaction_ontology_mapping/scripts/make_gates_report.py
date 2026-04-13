#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Aggregate metrics + contract validation into final gate report")
    p.add_argument("--metrics", type=Path, required=True)
    p.add_argument("--validation", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    m = load_json(args.metrics)
    v = load_json(args.validation)

    checks = {
        "mappable_record_ratio_ge_0_85": bool(m.get("gates", {}).get("checks", {}).get("mappable_record_ratio_ge_0_85", False)),
        "ontology_uri_valid_rate_ge_0_99": bool(m.get("gates", {}).get("checks", {}).get("ontology_uri_valid_rate_ge_0_99", False)),
        "three_interaction_types_present": bool(m.get("gates", {}).get("checks", {}).get("three_interaction_types_present", False)),
        "contract_validation_pass": bool(v.get("passed", False)),
    }

    payload = {
        "pipeline": "interaction_ontology_mapping_v2",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "coverage": m.get("coverage", {}),
        "row_count": m.get("row_count", {}),
        "metrics_report": str(args.metrics),
        "validation_report": str(args.validation),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"[{payload['status']}] gates -> {args.out}")
    return 0 if payload["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
