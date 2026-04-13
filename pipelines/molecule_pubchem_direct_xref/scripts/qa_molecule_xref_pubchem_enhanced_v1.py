#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize(x: str) -> str:
    return (x or "").strip()


def split_multi(x: str) -> List[str]:
    text = normalize(x)
    if not text:
        return []
    return [t.strip() for t in text.split(";") if t.strip()]


def load_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        return list(reader)


def coverage(rows: List[Dict[str, str]]) -> Dict[str, Any]:
    n = len(rows)
    c = sum(1 for row in rows if split_multi(row.get("pubchem_cid", "")))
    return {"rows": n, "covered_rows": c, "rate": (c / n) if n else 0.0}


@dataclass
class Gate:
    gate_id: str
    passed: bool
    detail: Dict[str, Any]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--core-tsv", type=Path, required=True)
    ap.add_argument("--enhanced-tsv", type=Path, required=True)
    ap.add_argument("--build-report", type=Path, required=True)
    ap.add_argument("--conflict-audit", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--baseline-min-rate", type=float, default=0.244)
    ap.add_argument("--allow-row-subset", type=int, default=0)
    args = ap.parse_args()

    for p in [args.core_tsv, args.enhanced_tsv, args.build_report, args.conflict_audit]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    core_rows = load_tsv(args.core_tsv)
    enhanced_rows = load_tsv(args.enhanced_tsv)
    allow_row_subset = bool(args.allow_row_subset)

    core_by_ik = {normalize(r.get("inchikey", "")).upper(): r for r in core_rows}
    enhanced_by_ik = {normalize(r.get("inchikey", "")).upper(): r for r in enhanced_rows}
    comparable_core_rows = [core_by_ik[ik] for ik in enhanced_by_ik.keys() if ik in core_by_ik]

    if allow_row_subset:
        core = coverage(comparable_core_rows)
    else:
        core = coverage(core_rows)
    enh = coverage(enhanced_rows)

    newly_filled = []
    high_conf_overwrite = 0
    missing_required_on_added = 0
    for ik, src in core_by_ik.items():
        dst = enhanced_by_ik.get(ik)
        if not dst:
            continue

        src_pub = split_multi(src.get("pubchem_cid", ""))
        dst_pub = split_multi(dst.get("pubchem_cid", ""))

        if src_pub and normalize(src.get("confidence", "")).lower() == "high":
            if ";".join(src_pub) != ";".join(dst_pub):
                high_conf_overwrite += 1

        if (not src_pub) and dst_pub:
            newly_filled.append(ik)
            if not normalize(dst.get("match_strategy", "")) or not normalize(dst.get("confidence", "")) or not normalize(dst.get("source_version", "")):
                missing_required_on_added += 1

    build = json.loads(args.build_report.read_text(encoding="utf-8"))
    conflict = json.loads(args.conflict_audit.read_text(encoding="utf-8"))

    gates = [
        Gate(
            gate_id="row_count_preserved",
            passed=((len(core_rows) == len(enhanced_rows)) if not allow_row_subset else (len(enhanced_rows) <= len(core_rows))),
            detail={
                "core_rows": len(core_rows),
                "enhanced_rows": len(enhanced_rows),
                "allow_row_subset": allow_row_subset,
            },
        ),
        Gate(
            gate_id="pubchem_coverage_improved_vs_core",
            passed=(enh["rate"] > core["rate"]),
            detail={"core_rate": core["rate"], "enhanced_rate": enh["rate"], "delta": enh["rate"] - core["rate"]},
        ),
        Gate(
            gate_id="pubchem_coverage_above_baseline_24_4pct",
            passed=(enh["rate"] > args.baseline_min_rate),
            detail={"enhanced_rate": enh["rate"], "baseline_min_rate": args.baseline_min_rate},
        ),
        Gate(
            gate_id="added_rows_have_strategy_confidence_source_version",
            passed=(missing_required_on_added == 0),
            detail={"newly_filled_rows": len(newly_filled), "missing_required_fields_count": missing_required_on_added},
        ),
        Gate(
            gate_id="high_confidence_existing_mapping_not_overwritten",
            passed=(high_conf_overwrite == 0),
            detail={"high_conf_overwrite_count": high_conf_overwrite},
        ),
        Gate(
            gate_id="conflict_audit_exists_and_structured",
            passed=isinstance(conflict.get("records", []), list) and "tie_break_counts" in conflict,
            detail={"conflict_rows": len(conflict.get("records", [])), "tie_break_keys": sorted(list(conflict.get("tie_break_counts", {}).keys()))},
        ),
        Gate(
            gate_id="build_report_coverage_consistency",
            passed=abs(float(build.get("metrics", {}).get("after", {}).get("coverage_rate", enh["rate"])) - enh["rate"]) < 1e-9,
            detail={
                "qa_enhanced_rate": enh["rate"],
                "build_enhanced_rate": float(build.get("metrics", {}).get("after", {}).get("coverage_rate", enh["rate"])),
            },
        ),
    ]

    report = {
        "name": "molecule_xref_pubchem_enhanced_v1.qa",
        "created_at": utc_now(),
        "inputs": {
            "core_tsv": str(args.core_tsv),
            "enhanced_tsv": str(args.enhanced_tsv),
            "build_report": str(args.build_report),
            "conflict_audit": str(args.conflict_audit),
        },
        "metrics": {
            "core": core,
            "enhanced": enh,
            "comparable_core_rows": len(comparable_core_rows),
            "newly_filled_rows": len(newly_filled),
            "high_conf_overwrite_count": high_conf_overwrite,
        },
        "gates": [{"id": g.gate_id, "passed": g.passed, "detail": g.detail} for g in gates],
        "passed": all(g.passed for g in gates),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    status = "PASS" if report["passed"] else "FAIL"
    print(f"[{status}] molecule_xref_pubchem_enhanced_v1 QA -> {args.out}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
