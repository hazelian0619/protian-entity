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
    p.add_argument(
        "--metrics",
        type=Path,
        required=True,
        help="metrics json from build_rpi_site_domain_enrichment.py",
    )
    p.add_argument(
        "--validations",
        type=Path,
        nargs="+",
        required=True,
        help="validation report json paths",
    )
    p.add_argument("--out", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()

    metrics = load_json(args.metrics)
    validations: List[Dict] = [load_json(p) for p in args.validations]

    contract_checks = {
        str(p): bool(v.get("passed", False)) for p, v in zip(args.validations, validations)
    }

    coverage = metrics.get("coverage", {})
    checks = {
        "edge_id_join_rate_ge_0_99": bool(metrics.get("gates", {}).get("checks", {}).get("edge_id_join_rate_ge_0_99", False)),
        "method_field_coverage_ge_0_90": bool(metrics.get("gates", {}).get("checks", {}).get("method_field_coverage_ge_0_90", False)),
        "site_or_domain_coverage_ge_0_60": bool(metrics.get("gates", {}).get("checks", {}).get("site_or_domain_coverage_ge_0_60", False)),
        "contract_validation_pass": all(contract_checks.values()),
    }

    payload = {
        "pipeline": "rpi_site_domain_enrichment_v2",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "coverage": {
            "edge_id_join_rate": coverage.get("edge_id_join_rate", 0.0),
            "method_field_coverage": coverage.get("method_field_coverage", 0.0),
            "site_edge_coverage": coverage.get("site_edge_coverage", 0.0),
            "domain_edge_coverage": coverage.get("domain_edge_coverage", 0.0),
            "site_or_domain_coverage": coverage.get("site_or_domain_coverage", 0.0),
        },
        "contract_reports": contract_checks,
        "metrics_report": str(args.metrics),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"[{payload['status']}] gates -> {args.out}")
    return 0 if payload["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
