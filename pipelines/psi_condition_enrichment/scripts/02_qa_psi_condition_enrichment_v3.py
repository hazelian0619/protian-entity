#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
from pathlib import Path
from typing import Dict

ALLOWED_REL = {"=", "<", ">", "<=", ">=", "~"}

HARD_THRESHOLDS = {
    "condition_context": 0.55,
    "condition_pH": 0.22,
    "condition_temperature_c": 0.15,
}

TARGET_THRESHOLDS = {
    "condition_context": 0.65,
    "condition_pH": 0.30,
    "condition_temperature_c": 0.22,
}


def utc_now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def non_empty(v: object) -> bool:
    s = str(v or "").strip()
    return s not in {"", "NA", "N/A", "None", "null"}


def safe_rate(num: int, den: int) -> float:
    if den <= 0:
        return 0.0
    return num / den


def scan_table(path: Path) -> Dict[str, float]:
    stats = {
        "rows": 0,
        "edge_id_non_empty": 0,
        "standard_type_non_empty": 0,
        "standard_relation_allowed": 0,
        "assay_type_non_empty": 0,
        "condition_context_non_empty": 0,
        "condition_pH_non_empty": 0,
        "condition_temperature_non_empty": 0,
    }
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        need = {
            "edge_id",
            "standard_type",
            "standard_relation",
            "assay_type",
            "condition_context",
            "condition_pH",
            "condition_temperature_c",
        }
        miss = need - set(r.fieldnames or [])
        if miss:
            raise SystemExit(f"[ERROR] missing columns in {path}: {sorted(miss)}")

        for row in r:
            stats["rows"] += 1
            if non_empty(row.get("edge_id")):
                stats["edge_id_non_empty"] += 1
            if non_empty(row.get("standard_type")):
                stats["standard_type_non_empty"] += 1
            if str(row.get("standard_relation") or "").strip() in ALLOWED_REL:
                stats["standard_relation_allowed"] += 1
            if non_empty(row.get("assay_type")):
                stats["assay_type_non_empty"] += 1
            if non_empty(row.get("condition_context")):
                stats["condition_context_non_empty"] += 1
            if non_empty(row.get("condition_pH")):
                stats["condition_pH_non_empty"] += 1
            if non_empty(row.get("condition_temperature_c")):
                stats["condition_temperature_non_empty"] += 1

    rows = stats["rows"]
    return {
        "rows": rows,
        "edge_id": safe_rate(stats["edge_id_non_empty"], rows),
        "standard_type": safe_rate(stats["standard_type_non_empty"], rows),
        "standard_relation_allowed": safe_rate(stats["standard_relation_allowed"], rows),
        "assay_type": safe_rate(stats["assay_type_non_empty"], rows),
        "condition_context": safe_rate(stats["condition_context_non_empty"], rows),
        "condition_pH": safe_rate(stats["condition_pH_non_empty"], rows),
        "condition_temperature_c": safe_rate(stats["condition_temperature_non_empty"], rows),
    }


def gate(gate_id: str, passed: bool, detail: Dict) -> Dict:
    return {"id": gate_id, "passed": bool(passed), "detail": detail}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--v2", type=Path, required=True)
    ap.add_argument("--v3", type=Path, required=True)
    ap.add_argument("--v3-validation", type=Path, required=True)
    ap.add_argument("--audit-validation", type=Path, required=True)
    ap.add_argument("--build-report", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    for p in [args.v2, args.v3, args.v3_validation, args.audit_validation, args.build_report]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    m2 = scan_table(args.v2)
    m3 = scan_table(args.v3)

    v3_validation = json.loads(args.v3_validation.read_text(encoding="utf-8"))
    audit_validation = json.loads(args.audit_validation.read_text(encoding="utf-8"))
    build = json.loads(args.build_report.read_text(encoding="utf-8"))

    hard_gates = [
        gate(
            "edge_id_coverage_keep_1_0",
            (m3["edge_id"] >= 1.0) and (m3["edge_id"] >= m2["edge_id"]),
            {"v2": m2["edge_id"], "v3": m3["edge_id"], "required": 1.0},
        ),
        gate(
            "standard_type_no_regression",
            m3["standard_type"] >= m2["standard_type"],
            {"v2": m2["standard_type"], "v3": m3["standard_type"]},
        ),
        gate(
            "standard_relation_no_regression",
            m3["standard_relation_allowed"] >= m2["standard_relation_allowed"],
            {"v2": m2["standard_relation_allowed"], "v3": m3["standard_relation_allowed"]},
        ),
        gate(
            "assay_type_no_regression",
            m3["assay_type"] >= m2["assay_type"],
            {"v2": m2["assay_type"], "v3": m3["assay_type"]},
        ),
        gate(
            "condition_context_min_0_55",
            m3["condition_context"] >= HARD_THRESHOLDS["condition_context"],
            {"v3": m3["condition_context"], "threshold": HARD_THRESHOLDS["condition_context"]},
        ),
        gate(
            "condition_pH_min_0_22",
            m3["condition_pH"] >= HARD_THRESHOLDS["condition_pH"],
            {"v3": m3["condition_pH"], "threshold": HARD_THRESHOLDS["condition_pH"]},
        ),
        gate(
            "condition_temperature_min_0_15",
            m3["condition_temperature_c"] >= HARD_THRESHOLDS["condition_temperature_c"],
            {"v3": m3["condition_temperature_c"], "threshold": HARD_THRESHOLDS["condition_temperature_c"]},
        ),
        gate(
            "contracts_validation_pass",
            bool(v3_validation.get("passed")) and bool(audit_validation.get("passed")),
            {
                "psi_activity_context_v3_contract_pass": bool(v3_validation.get("passed")),
                "psi_condition_parse_audit_v3_contract_pass": bool(audit_validation.get("passed")),
            },
        ),
    ]

    target_gates = [
        gate(
            "target_condition_context_0_65",
            m3["condition_context"] >= TARGET_THRESHOLDS["condition_context"],
            {"v3": m3["condition_context"], "target": TARGET_THRESHOLDS["condition_context"]},
        ),
        gate(
            "target_condition_pH_0_30",
            m3["condition_pH"] >= TARGET_THRESHOLDS["condition_pH"],
            {"v3": m3["condition_pH"], "target": TARGET_THRESHOLDS["condition_pH"]},
        ),
        gate(
            "target_condition_temperature_0_22",
            m3["condition_temperature_c"] >= TARGET_THRESHOLDS["condition_temperature_c"],
            {"v3": m3["condition_temperature_c"], "target": TARGET_THRESHOLDS["condition_temperature_c"]},
        ),
    ]

    passed = all(g["passed"] for g in hard_gates)

    report = {
        "name": "psi_condition_enrichment_v3.qa",
        "created_at": utc_now_iso(),
        "inputs": {
            "v2": str(args.v2),
            "v3": str(args.v3),
            "v3_validation": str(args.v3_validation),
            "audit_validation": str(args.audit_validation),
            "build_report": str(args.build_report),
        },
        "metrics": {
            "v2": m2,
            "v3": m3,
            "delta": {k: (m3.get(k, 0.0) - m2.get(k, 0.0)) for k in [
                "edge_id",
                "standard_type",
                "standard_relation_allowed",
                "assay_type",
                "condition_context",
                "condition_pH",
                "condition_temperature_c",
            ]},
            "conflict_rows": build.get("counts", {}).get("conflict_rows", 0),
            "rule_hits": build.get("rule_hits", {}),
        },
        "hard_gates": hard_gates,
        "target_gates": target_gates,
        "passed": passed,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if passed else "FAIL"
    print(
        f"[{status}] psi_condition_enrichment_v3 QA "
        f"ctx={m3['condition_context']:.4f} pH={m3['condition_pH']:.4f} temp={m3['condition_temperature_c']:.4f}"
    )
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
