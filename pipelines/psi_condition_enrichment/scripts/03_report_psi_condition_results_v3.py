#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import random
from pathlib import Path
from typing import Dict, List


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


def count_rows(path: Path) -> int:
    n = 0
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for _ in r:
            n += 1
    return n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--v3", type=Path, required=True)
    ap.add_argument("--audit", type=Path, required=True)
    ap.add_argument("--build", type=Path, required=True)
    ap.add_argument("--qa", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--sample", type=Path, required=True)
    ap.add_argument("--sample-size", type=int, default=20)
    ap.add_argument("--sample-seed", type=int, default=20260414)
    args = ap.parse_args()

    for p in [args.v3, args.audit, args.build, args.qa]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    rng = random.Random(args.sample_seed)
    samples: List[Dict[str, str]] = []
    fields: List[str] = []

    with args.v3.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        fields = list(r.fieldnames or [])
        for i, row in enumerate(r, start=1):
            reservoir_sample(samples, row, args.sample_size, i, rng)

    write_sample(args.sample, fields, samples)

    qa = json.loads(args.qa.read_text(encoding="utf-8"))
    build = json.loads(args.build.read_text(encoding="utf-8"))

    report = {
        "name": "psi_condition_enrichment_v3.result_evidence",
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "tables": {
            "psi_activity_context_v3": {
                "path": str(args.v3),
                "rows": count_rows(args.v3),
                "sample20_path": str(args.sample),
            },
            "psi_condition_parse_audit_v3": {
                "path": str(args.audit),
                "rows": count_rows(args.audit),
            },
        },
        "v2_to_v3_metrics": qa.get("metrics", {}),
        "extraction_rule_hits": build.get("rule_hits", {}),
        "source_field_hits": build.get("source_field_hits", {}),
        "qa_passed": bool(qa.get("passed", False)),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(
        "[OK] psi condition result evidence generated "
        f"sample={args.sample} rows={report['tables']['psi_activity_context_v3']['rows']:,}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
