#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

PDB_ID_RE = re.compile(r"^[0-9A-Z]{4}$")


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def parse_iso_date(s: str):
    try:
        return datetime.strptime(s, "%Y-%m-%d").date()
    except Exception:
        return None


def load_master_ids(path: Path) -> Tuple[Set[str], int]:
    ids: Set[str] = set()
    rows = 0
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {"rna_id"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] master schema missing columns: {sorted(missing)}")

        for row in r:
            rows += 1
            rid = (row.get("rna_id") or "").strip()
            if rid:
                ids.add(rid)
    return ids, rows


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
    ap.add_argument("--release-date-min", default="1970-01-01")
    args = ap.parse_args()

    if not args.master.exists():
        raise SystemExit(f"[ERROR] missing master: {args.master}")
    if not args.table.exists():
        raise SystemExit(f"[ERROR] missing table: {args.table}")
    if not args.build_report.exists():
        raise SystemExit(f"[ERROR] missing build report: {args.build_report}")

    master_ids, master_rows = load_master_ids(args.master)
    build_report = json.loads(args.build_report.read_text(encoding="utf-8"))

    out_pair_keys: Set[Tuple[str, str, str]] = set()
    duplicate_pairs = 0

    orphan_rna = 0
    orphan_examples: List[Dict[str, str]] = []

    pdb_format_bad = 0
    pdb_format_bad_examples: List[Dict[str, str]] = []

    release_non_empty = 0
    release_bad_format = 0
    release_bad_range = 0
    release_min_seen = None
    release_max_seen = None
    release_bad_examples: List[Dict[str, Any]] = []

    resolution_non_empty = 0
    resolution_bad_format = 0
    resolution_bad_range = 0
    resolution_bad_range_non_em = 0
    resolution_bad_range_em = 0
    resolution_min_seen = None
    resolution_max_seen = None
    resolution_bad_examples: List[Dict[str, Any]] = []

    unique_rna_with_rows: Set[str] = set()
    rows = 0

    release_date_floor = parse_iso_date(args.release_date_min)
    if release_date_floor is None:
        raise SystemExit(f"[ERROR] invalid --release-date-min: {args.release_date_min}")
    release_date_ceiling = datetime.now(timezone.utc).date()

    with args.table.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {
            "rna_id",
            "rna_type",
            "urs_id",
            "pdb_id",
            "pdb_entity_id",
            "mapping_strategy",
            "experimental_method",
            "resolution",
            "release_date",
            "source",
            "source_version",
            "fetch_date",
        }
        missing = required - set(r.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] output schema missing columns: {sorted(missing)}")

        for row in r:
            rows += 1
            rna_id = (row.get("rna_id") or "").strip()
            pdb_id = (row.get("pdb_id") or "").strip().upper()
            pdb_entity_id = (row.get("pdb_entity_id") or "").strip().upper()

            if rna_id:
                unique_rna_with_rows.add(rna_id)
            key = (rna_id, pdb_id, pdb_entity_id)
            if key in out_pair_keys:
                duplicate_pairs += 1
            else:
                out_pair_keys.add(key)

            if rna_id not in master_ids:
                orphan_rna += 1
                if len(orphan_examples) < 20:
                    orphan_examples.append({"rna_id": rna_id, "pdb_id": pdb_id})

            if not PDB_ID_RE.match(pdb_id):
                pdb_format_bad += 1
                if len(pdb_format_bad_examples) < 20:
                    pdb_format_bad_examples.append({"rna_id": rna_id, "pdb_id": pdb_id})

            release_raw = (row.get("release_date") or "").strip()
            if release_raw:
                release_non_empty += 1
                d = parse_iso_date(release_raw)
                if d is None:
                    release_bad_format += 1
                    if len(release_bad_examples) < 20:
                        release_bad_examples.append(
                            {"rna_id": rna_id, "pdb_id": pdb_id, "release_date": release_raw, "error": "bad_format"}
                        )
                else:
                    if d < release_date_floor or d > release_date_ceiling:
                        release_bad_range += 1
                        if len(release_bad_examples) < 20:
                            release_bad_examples.append(
                                {
                                    "rna_id": rna_id,
                                    "pdb_id": pdb_id,
                                    "release_date": release_raw,
                                    "allowed_min": str(release_date_floor),
                                    "allowed_max": str(release_date_ceiling),
                                }
                            )
                    if release_min_seen is None or d < release_min_seen:
                        release_min_seen = d
                    if release_max_seen is None or d > release_max_seen:
                        release_max_seen = d

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
                                    "rna_id": rna_id,
                                    "pdb_id": pdb_id,
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

    expected_pairs = int(build_report.get("selected_pairs") or rows)
    sample_mode = bool(build_report.get("sample_mode"))

    api_failure_rate = float(build_report.get("api_failure_rate") or 0.0)
    api_fail_threshold = float(build_report.get("fail_threshold") or 0.05)

    rna_backlink_coverage = (len(unique_rna_with_rows) / master_rows) if master_rows else 0.0

    gates: List[Gate] = [
        Gate(
            gate_id="row_count_matches_expected_pairs",
            passed=(rows == expected_pairs),
            detail={
                "output_rows": rows,
                "expected_pairs": expected_pairs,
                "sample_mode": sample_mode,
            },
        ),
        Gate(
            gate_id="rna_id_traceable_to_master",
            passed=(orphan_rna == 0),
            detail={"orphans": orphan_rna, "orphan_examples": orphan_examples},
        ),
        Gate(
            gate_id="pair_unique",
            passed=(duplicate_pairs == 0),
            detail={"duplicate_pairs": duplicate_pairs},
        ),
        Gate(
            gate_id="pdb_id_format",
            passed=(pdb_format_bad == 0),
            detail={"bad_rows": pdb_format_bad, "examples": pdb_format_bad_examples},
        ),
        Gate(
            gate_id="release_date_format_and_range",
            passed=(release_bad_format == 0 and release_bad_range == 0),
            detail={
                "non_empty_release_rows": release_non_empty,
                "bad_format": release_bad_format,
                "bad_range": release_bad_range,
                "observed_min": str(release_min_seen) if release_min_seen else None,
                "observed_max": str(release_max_seen) if release_max_seen else None,
                "allowed_min": str(release_date_floor),
                "allowed_max": str(release_date_ceiling),
                "examples": release_bad_examples,
            },
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
        "name": "rna_pdb_structures_v1_qa",
        "created_at": utc_now_iso(),
        "inputs": {
            "master": str(args.master),
            "table": str(args.table),
            "build_report": str(args.build_report),
        },
        "metrics": {
            "master_rows": master_rows,
            "output_rows": rows,
            "output_unique_pairs": len(out_pair_keys),
            "output_unique_rna_ids": len(unique_rna_with_rows),
            "rna_backlink_coverage_rate": rna_backlink_coverage,
            "orphan_rna_rows": orphan_rna,
            "duplicate_pairs": duplicate_pairs,
            "pdb_format_bad_rows": pdb_format_bad,
            "release_date_non_empty_rows": release_non_empty,
            "release_date_bad_format": release_bad_format,
            "release_date_bad_range": release_bad_range,
            "resolution_non_empty_rows": resolution_non_empty,
            "resolution_bad_format": resolution_bad_format,
            "resolution_bad_range": resolution_bad_range,
            "api_failure_rate": api_failure_rate,
            "api_failure_threshold": api_fail_threshold,
            "expected_pairs": expected_pairs,
            "sample_mode": sample_mode,
        },
        "gates": [{"id": g.gate_id, "passed": g.passed, "detail": g.detail} for g in gates],
        "passed": passed,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if passed else "FAIL"
    print(f"[{status}] rna_pdb_structures_v1 QA -> {args.out}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
