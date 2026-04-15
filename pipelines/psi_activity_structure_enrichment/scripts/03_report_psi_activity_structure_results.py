#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import random
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple


def non_empty(v: object) -> bool:
    return str(v or "").strip() not in {"", "NA", "N/A", "None", "null"}


def pick_value(row: Dict[str, str], names: List[str]) -> str:
    for n in names:
        if n in row:
            return str(row.get(n) or "").strip()
    return ""


def reservoir_sample(samples: List[Dict[str, str]], row: Dict[str, str], n: int, idx: int, rng: random.Random) -> None:
    if len(samples) < n:
        samples.append(dict(row))
        return
    j = rng.randint(1, idx)
    if j <= n:
        samples[j - 1] = dict(row)


def write_sample(path: Path, fieldnames: List[str], rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def scan_activity(path: Path, sample_n: int, seed: int) -> Tuple[Dict, List[str], List[Dict[str, str]]]:
    counts = {
        "rows": 0,
        "activity_type_non_empty": 0,
        "relation_non_empty": 0,
        "assay_type_non_empty": 0,
        "assay_id_non_empty": 0,
        "conditions_non_empty": 0,
    }
    dist = Counter()
    samples: List[Dict[str, str]] = []
    rng = random.Random(seed)
    fieldnames: List[str] = []

    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        fieldnames = list(r.fieldnames or [])
        for i, row in enumerate(r, start=1):
            counts["rows"] += 1
            activity_type = pick_value(row, ["activity_type", "standard_type"]).upper()
            relation = pick_value(row, ["relation", "standard_relation"])
            assay_type = pick_value(row, ["assay_type"])
            assay_id = pick_value(row, ["assay_id"])
            conditions = pick_value(row, ["conditions", "condition_context"])

            if non_empty(activity_type):
                counts["activity_type_non_empty"] += 1
            if non_empty(relation):
                counts["relation_non_empty"] += 1
            if non_empty(assay_type):
                counts["assay_type_non_empty"] += 1
            if non_empty(assay_id):
                counts["assay_id_non_empty"] += 1
            if non_empty(conditions):
                counts["conditions_non_empty"] += 1

            if activity_type in {"IC50", "KI", "KD", "EC50"}:
                dist[activity_type] += 1
            elif non_empty(activity_type):
                dist["OTHER"] += 1

            reservoir_sample(samples, row, sample_n, i, rng)

    rows = counts["rows"] or 1
    coverage = {
        "activity_type": counts["activity_type_non_empty"] / rows,
        "relation": counts["relation_non_empty"] / rows,
        "assay_type": counts["assay_type_non_empty"] / rows,
        "assay_id": counts["assay_id_non_empty"] / rows,
        "conditions": counts["conditions_non_empty"] / rows,
    }

    return {
        "rows": counts["rows"],
        "coverage": coverage,
        "activity_type_distribution": {
            "IC50": dist.get("IC50", 0),
            "Ki": dist.get("KI", 0),
            "Kd": dist.get("KD", 0),
            "EC50": dist.get("EC50", 0),
            "OTHER": dist.get("OTHER", 0),
        },
    }, fieldnames, samples


def scan_structure(path: Path, sample_n: int, seed: int) -> Tuple[Dict, List[str], List[Dict[str, str]]]:
    counts = {
        "rows": 0,
        "pdb_id_non_empty": 0,
        "affinity_non_empty": 0,
    }
    samples: List[Dict[str, str]] = []
    rng = random.Random(seed)
    fieldnames: List[str] = []

    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        fieldnames = list(r.fieldnames or [])
        for i, row in enumerate(r, start=1):
            counts["rows"] += 1
            pdb_id = pick_value(row, ["pdb_id"])
            affinity = pick_value(row, ["affinity", "structure_affinity_score"])

            if non_empty(pdb_id):
                counts["pdb_id_non_empty"] += 1
            if non_empty(affinity):
                counts["affinity_non_empty"] += 1

            reservoir_sample(samples, row, sample_n, i, rng)

    rows = counts["rows"] or 1
    coverage = {
        "pdb_id": counts["pdb_id_non_empty"] / rows,
        "affinity": counts["affinity_non_empty"] / rows,
    }

    return {
        "rows": counts["rows"],
        "coverage": coverage,
    }, fieldnames, samples


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--activity", type=Path, required=True)
    ap.add_argument("--structure", type=Path, required=True)
    ap.add_argument("--qa-report", type=Path, required=True)
    ap.add_argument("--out-report", type=Path, required=True)
    ap.add_argument("--activity-sample", type=Path, required=True)
    ap.add_argument("--structure-sample", type=Path, required=True)
    ap.add_argument("--sample-size", type=int, default=20)
    ap.add_argument("--sample-seed", type=int, default=20260412)
    args = ap.parse_args()

    for p in [args.activity, args.structure, args.qa_report]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    qa = json.loads(args.qa_report.read_text(encoding="utf-8"))

    activity_metrics, activity_fields, activity_samples = scan_activity(args.activity, args.sample_size, args.sample_seed)
    structure_metrics, structure_fields, structure_samples = scan_structure(
        args.structure, args.sample_size, args.sample_seed + 17
    )

    write_sample(args.activity_sample, activity_fields, activity_samples)
    write_sample(args.structure_sample, structure_fields, structure_samples)

    report = {
        "name": "psi_activity_structure_enrichment_v2.result_evidence",
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "tables": {
            "psi_activity_context_v2": {
                "path": str(args.activity),
                **activity_metrics,
                "sample20_path": str(args.activity_sample),
            },
            "psi_structure_evidence_v2": {
                "path": str(args.structure),
                **structure_metrics,
                "sample20_path": str(args.structure_sample),
            },
        },
        "acceptance_metrics": {
            "edge_id_join_rate_activity": qa.get("metrics", {}).get("edge_id_join_rate_activity", 0.0),
            "edge_id_join_rate_structure": qa.get("metrics", {}).get("edge_id_join_rate_structure", 0.0),
            "activity_type_coverage": qa.get("metrics", {}).get("activity_type_coverage", 0.0),
            "assay_field_coverage": qa.get("metrics", {}).get("assay_field_coverage", 0.0),
            "structure_edge_coverage": qa.get("metrics", {}).get("structure_edge_coverage", 0.0),
            "contracts_pass": qa.get("metrics", {}).get("contracts_pass", False),
            "gates_pass": qa.get("passed", False),
        },
    }

    args.out_report.parent.mkdir(parents=True, exist_ok=True)
    args.out_report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(
        "[OK] result evidence report generated "
        f"activity_rows={activity_metrics['rows']:,} structure_rows={structure_metrics['rows']:,}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
