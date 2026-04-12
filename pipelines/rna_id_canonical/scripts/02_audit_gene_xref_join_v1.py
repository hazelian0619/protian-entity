#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def write_json(path: Path, obj: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--canonical-map", type=Path, required=True)
    ap.add_argument("--gene-xref-rna", type=Path, default=Path("data/output/gene_xref_rna_v1.tsv"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--max-rows", type=int, default=None)
    args = ap.parse_args()

    if not args.canonical_map.exists():
        raise SystemExit(f"[ERROR] missing canonical map: {args.canonical_map}")
    if not args.gene_xref_rna.exists():
        raise SystemExit(f"[ERROR] missing gene_xref_rna: {args.gene_xref_rna}")

    legacy_to_canonical: Dict[str, str] = {}
    legacy_to_type: Dict[str, str] = {}
    canonical_set = set()
    canonical_rows = 0

    with args.canonical_map.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"legacy_rna_id", "canonical_rna_id", "rna_type"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] missing required columns in canonical map: {sorted(missing)}")

        for row in reader:
            legacy = (row.get("legacy_rna_id") or "").strip()
            canonical = (row.get("canonical_rna_id") or "").strip()
            rna_type = (row.get("rna_type") or "").strip().lower()
            if legacy == "":
                continue
            legacy_to_canonical[legacy] = canonical
            legacy_to_type[legacy] = rna_type
            if canonical != "":
                canonical_set.add(canonical)
            canonical_rows += 1

    total_gene_rows = 0
    matched_by_legacy = 0
    matched_by_canonical_direct = 0
    changed_after_canonicalization = 0
    unchanged_after_canonicalization = 0
    per_type = defaultdict(lambda: Counter())
    unmatched_examples = []

    with args.gene_xref_rna.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "rna_id" not in (reader.fieldnames or []):
            raise SystemExit("[ERROR] gene_xref_rna missing column: rna_id")

        for row in reader:
            rid = (row.get("rna_id") or "").strip()
            total_gene_rows += 1

            if rid in legacy_to_canonical:
                matched_by_legacy += 1
                canonical = legacy_to_canonical[rid]
                rt = legacy_to_type.get(rid, "unknown")
                per_type[rt]["total"] += 1
                per_type[rt]["matched_legacy"] += 1
                if canonical == rid:
                    unchanged_after_canonicalization += 1
                    per_type[rt]["unchanged"] += 1
                else:
                    changed_after_canonicalization += 1
                    per_type[rt]["changed"] += 1
            else:
                if len(unmatched_examples) < 20:
                    unmatched_examples.append(rid)

            if rid in canonical_set:
                matched_by_canonical_direct += 1

            if args.max_rows is not None and total_gene_rows >= args.max_rows:
                break

    by_type = {}
    for rt, cnt in sorted(per_type.items()):
        n = cnt.get("total", 0)
        by_type[rt] = {
            "rows": n,
            "legacy_join_success_rate": (cnt.get("matched_legacy", 0) / n) if n else 0.0,
            "changed_after_canonicalization": cnt.get("changed", 0),
            "changed_rate": (cnt.get("changed", 0) / n) if n else 0.0,
        }

    report = {
        "name": "rna_id_canonical_map_v1.gene_xref_join_audit",
        "generated_at": utc_now_iso(),
        "inputs": {
            "canonical_map": str(args.canonical_map),
            "gene_xref_rna": str(args.gene_xref_rna),
            "sample_mode": args.max_rows is not None,
            "max_rows": args.max_rows,
        },
        "counts": {
            "canonical_rows": canonical_rows,
            "gene_xref_rows_scanned": total_gene_rows,
            "matched_by_legacy": matched_by_legacy,
            "matched_by_canonical_direct": matched_by_canonical_direct,
            "changed_after_canonicalization": changed_after_canonicalization,
            "unchanged_after_canonicalization": unchanged_after_canonicalization,
        },
        "rates": {
            "legacy_join_success_rate": (matched_by_legacy / total_gene_rows) if total_gene_rows else 0.0,
            "canonical_direct_join_rate": (matched_by_canonical_direct / total_gene_rows) if total_gene_rows else 0.0,
            "changed_after_canonicalization_rate": (changed_after_canonicalization / matched_by_legacy) if matched_by_legacy else 0.0,
        },
        "by_rna_type": by_type,
        "unmatched_examples": unmatched_examples,
        "gates": {
            "legacy_join_success_100": matched_by_legacy == total_gene_rows,
        },
    }
    write_json(args.out, report)
    print(f"[OK] wrote join audit: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
