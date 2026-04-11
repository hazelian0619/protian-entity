#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def parse_pdb_ids(raw: str) -> List[str]:
    if not raw:
        return []
    out: List[str] = []
    seen = set()
    for token in raw.split(";"):
        t = token.strip().upper()
        if not t:
            continue
        if t in seen:
            continue
        seen.add(t)
        out.append(t)
    return out


def load_master_pairs(path: Path) -> Tuple[Set[Tuple[str, str]], int]:
    pairs: Set[Tuple[str, str]] = set()
    rows = 0
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {"uniprot_id", "pdb_ids"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] master schema missing columns: {sorted(missing)}")

        for row in r:
            rows += 1
            uid = (row.get("uniprot_id") or "").strip()
            if not uid:
                continue
            for pdb_id in parse_pdb_ids((row.get("pdb_ids") or "").strip()):
                pairs.add((uid, pdb_id))
    return pairs, rows


@dataclass
class Gate:
    gate_id: str
    passed: bool
    detail: Dict[str, Any]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", type=Path, required=True)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--build-report", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--resolution-min", type=float, default=0.1)
    ap.add_argument("--resolution-max-non-em", type=float, default=25.0)
    ap.add_argument("--resolution-max-em", type=float, default=60.0)
    args = ap.parse_args()

    if not args.master.exists():
        raise SystemExit(f"[ERROR] missing master: {args.master}")
    if not args.table.exists():
        raise SystemExit(f"[ERROR] missing table: {args.table}")
    if not args.build_report.exists():
        raise SystemExit(f"[ERROR] missing build report: {args.build_report}")

    master_pairs, master_rows = load_master_pairs(args.master)
    build_report = json.loads(args.build_report.read_text(encoding="utf-8"))

    out_pairs: Set[Tuple[str, str]] = set()
    duplicate_pairs = 0
    orphan_pairs = 0
    orphan_examples: List[Dict[str, str]] = []

    resolution_non_empty = 0
    resolution_bad_format = 0
    resolution_bad_range = 0
    resolution_bad_range_non_em = 0
    resolution_bad_range_em = 0
    resolution_bad_examples: List[Dict[str, Any]] = []
    resolution_min_seen = None
    resolution_max_seen = None

    rows = 0

    with args.table.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {
            "pdb_id",
            "uniprot_id",
            "experimental_method",
            "resolution",
            "release_date",
            "ligand_count",
            "source",
            "fetch_date",
        }
        missing = required - set(r.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] output schema missing columns: {sorted(missing)}")

        for row in r:
            rows += 1
            uid = (row.get("uniprot_id") or "").strip()
            pdb_id = (row.get("pdb_id") or "").strip().upper()
            pair = (uid, pdb_id)
            if pair in out_pairs:
                duplicate_pairs += 1
            else:
                out_pairs.add(pair)

            if pair not in master_pairs:
                orphan_pairs += 1
                if len(orphan_examples) < 20:
                    orphan_examples.append({"uniprot_id": uid, "pdb_id": pdb_id})

            raw_res = (row.get("resolution") or "").strip()
            if raw_res:
                resolution_non_empty += 1
                try:
                    value = float(raw_res)
                    method = (row.get("experimental_method") or "").upper()
                    allowed_max = args.resolution_max_em if "ELECTRON MICROSCOPY" in method else args.resolution_max_non_em
                    if value < args.resolution_min or value > allowed_max:
                        resolution_bad_range += 1
                        if "ELECTRON MICROSCOPY" in method:
                            resolution_bad_range_em += 1
                        else:
                            resolution_bad_range_non_em += 1
                        if len(resolution_bad_examples) < 20:
                            resolution_bad_examples.append(
                                {
                                    "pdb_id": pdb_id,
                                    "uniprot_id": uid,
                                    "experimental_method": row.get("experimental_method") or "",
                                    "resolution": value,
                                    "allowed_min": args.resolution_min,
                                    "allowed_max": allowed_max,
                                }
                            )
                    if resolution_min_seen is None or value < resolution_min_seen:
                        resolution_min_seen = value
                    if resolution_max_seen is None or value > resolution_max_seen:
                        resolution_max_seen = value
                except Exception:
                    resolution_bad_format += 1

    missing_pairs = len(master_pairs - out_pairs)
    expected_pairs = int(build_report.get("selected_pairs") or len(master_pairs))
    sample_mode = bool(build_report.get("sample_mode"))

    api_failure_rate = float(build_report.get("api_failure_rate") or 0.0)
    api_fail_threshold = float(build_report.get("fail_threshold") or 0.05)

    gates: List[Gate] = [
        Gate(
            gate_id="row_count_matches_expected_pairs",
            passed=(rows == expected_pairs),
            detail={
                "output_rows": rows,
                "expected_pairs": expected_pairs,
                "sample_mode": sample_mode,
                "master_pairs": len(master_pairs),
                "missing_pairs_vs_master": missing_pairs,
            },
        ),
        Gate(
            gate_id="pair_traceability_to_master",
            passed=(orphan_pairs == 0),
            detail={"orphans": orphan_pairs, "orphan_examples": orphan_examples},
        ),
        Gate(
            gate_id="pair_unique",
            passed=(duplicate_pairs == 0),
            detail={"duplicate_pairs": duplicate_pairs},
        ),
        Gate(
            gate_id="resolution_legal_range",
            passed=(resolution_bad_format == 0 and resolution_bad_range == 0),
            detail={
                "non_empty_resolution_rows": resolution_non_empty,
                "bad_format": resolution_bad_format,
                "bad_range": resolution_bad_range,
                "bad_range_non_em": resolution_bad_range_non_em,
                "bad_range_em": resolution_bad_range_em,
                "allowed_min": args.resolution_min,
                "allowed_max_non_em": args.resolution_max_non_em,
                "allowed_max_em": args.resolution_max_em,
                "observed_min": resolution_min_seen,
                "observed_max": resolution_max_seen,
                "examples": resolution_bad_examples,
            },
        ),
        Gate(
            gate_id="api_failure_rate_below_threshold",
            passed=(api_failure_rate <= api_fail_threshold),
            detail={"api_failure_rate": api_failure_rate, "threshold": api_fail_threshold},
        ),
    ]

    passed = all(g.passed for g in gates)
    report: Dict[str, Any] = {
        "name": "pdb_structures_v1_qa",
        "created_at": utc_now_iso(),
        "inputs": {
            "master": str(args.master),
            "table": str(args.table),
            "build_report": str(args.build_report),
        },
        "metrics": {
            "master_rows": master_rows,
            "master_pairs": len(master_pairs),
            "expected_pairs": expected_pairs,
            "sample_mode": sample_mode,
            "output_rows": rows,
            "output_unique_pairs": len(out_pairs),
            "missing_pairs": missing_pairs,
            "orphan_pairs": orphan_pairs,
            "duplicate_pairs": duplicate_pairs,
            "resolution_non_empty_rows": resolution_non_empty,
        },
        "gates": [{"id": g.gate_id, "passed": g.passed, "detail": g.detail} for g in gates],
        "passed": passed,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if passed else "FAIL"
    print(f"[{status}] pdb_structures_v1 QA -> {args.out}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
