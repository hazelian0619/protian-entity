#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def write_tsv(path: Path, header: List[str], rows: Iterable[List[str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    n = 0
    with tmp.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for row in rows:
            w.writerow(row)
            n += 1
    tmp.replace(path)
    return n


def write_json(path: Path, obj: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def load_xref(xref_path: Path) -> Tuple[Dict[str, Dict[str, str]], Dict[str, object]]:
    xref: Dict[str, Dict[str, str]] = {}
    dup_conflicts: List[Dict[str, str]] = []
    match_strategy_counts = Counter()
    xref_rows = 0

    with xref_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"rna_id", "xref_id", "match_strategy", "source_version", "fetch_date"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] missing required columns in xref: {sorted(missing)}")

        for row in reader:
            xref_rows += 1
            rid = (row.get("rna_id") or "").strip()
            xref_id = (row.get("xref_id") or "").strip()
            match_strategy = (row.get("match_strategy") or "").strip()
            source_version = (row.get("source_version") or "").strip()
            fetch_date = (row.get("fetch_date") or "").strip()
            if rid == "" or xref_id == "":
                continue

            match_strategy_counts[match_strategy] += 1

            prev = xref.get(rid)
            if prev is not None and prev["xref_id"] != xref_id:
                dup_conflicts.append(
                    {
                        "rna_id": rid,
                        "xref_id_existing": prev["xref_id"],
                        "xref_id_new": xref_id,
                    }
                )
                continue

            xref[rid] = {
                "xref_id": xref_id,
                "match_strategy": match_strategy,
                "source_version": source_version,
                "fetch_date": fetch_date,
            }

    stats = {
        "xref_rows": xref_rows,
        "xref_unique_rna_ids": len(xref),
        "xref_conflict_count": len(dup_conflicts),
        "xref_conflicts": dup_conflicts,
        "xref_match_strategy_counts": dict(match_strategy_counts),
    }
    return xref, stats


def canonicalize_row(
    legacy_rna_id: str,
    rna_type: str,
    taxon_id: str,
    row_source_version: str,
    xref_map: Dict[str, Dict[str, str]],
    default_source: str,
    fetch_date: str,
) -> Tuple[List[str], str]:
    # strategy labels for report
    if legacy_rna_id in xref_map:
        x = xref_map[legacy_rna_id]
        canonical_rna_id = x["xref_id"]
        canonical_system = "RNAcentral"
        confidence = "high"
        strategy = "xref_enst_to_urs_v2"
        evidence = f"xref:{x.get('match_strategy', '')}"
        sv = row_source_version.strip()
        xv = x.get("source_version", "").strip()
        source_version = f"{sv};xref:{xv}" if sv and xv else (sv or xv or "unknown")
    elif legacy_rna_id.startswith("URS"):
        canonical_rna_id = legacy_rna_id
        canonical_system = "RNAcentral"
        confidence = "high"
        strategy = "legacy_is_urs"
        evidence = "legacy_id_is_rnacentral_urs"
        source_version = row_source_version.strip() or "unknown"
    elif legacy_rna_id.startswith("ENST"):
        canonical_rna_id = legacy_rna_id
        canonical_system = "EnsemblTranscript"
        confidence = "medium"
        strategy = "fallback_legacy_enst"
        evidence = "no_xref_match_keep_legacy_enst"
        source_version = row_source_version.strip() or "unknown"
    else:
        canonical_rna_id = legacy_rna_id
        canonical_system = "LegacyID"
        confidence = "low"
        strategy = "fallback_legacy_other"
        evidence = "unsupported_legacy_prefix_keep_legacy_id"
        source_version = row_source_version.strip() or "unknown"

    out = [
        legacy_rna_id,
        canonical_rna_id,
        canonical_system,
        rna_type,
        taxon_id,
        confidence,
        evidence,
        default_source,
        fetch_date,
        source_version,
    ]
    return out, strategy


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", type=Path, default=Path("data/output/rna_master_v1.tsv"))
    ap.add_argument("--xref", type=Path, default=Path("data/output/rna_xref_mrna_enst_urs_v2.tsv"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--fetch-date", default=date.today().isoformat())
    ap.add_argument("--source", default="rna_master_v1;rna_xref_mrna_enst_urs_v2")
    ap.add_argument("--max-rows", type=int, default=None)
    args = ap.parse_args()

    if not args.master.exists():
        raise SystemExit(f"[ERROR] missing master: {args.master}")
    if not args.xref.exists():
        raise SystemExit(
            "[ERROR] missing xref (assistant A output): "
            f"{args.xref}\n"
            "Please run: bash pipelines/rna_xref_enst_urs/run.sh"
        )

    xref_map, xref_stats = load_xref(args.xref)

    required_master_cols = {"rna_id", "rna_type", "taxon_id", "source_version"}
    header = [
        "legacy_rna_id",
        "canonical_rna_id",
        "canonical_system",
        "rna_type",
        "taxon_id",
        "confidence",
        "evidence",
        "source",
        "fetch_date",
        "source_version",
    ]

    rows_out: List[List[str]] = []
    strategy_counts = Counter()
    canonical_system_counts = Counter()
    confidence_counts = Counter()
    rna_type_counts = Counter()
    canonical_non_empty_by_type = defaultdict(int)
    changed_by_type = defaultdict(int)
    total = 0
    duplicate_legacy = 0
    seen_legacy = set()

    with args.master.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = required_master_cols - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] missing required columns in master: {sorted(missing)}")

        for row in reader:
            legacy_rna_id = (row.get("rna_id") or "").strip()
            rna_type = (row.get("rna_type") or "").strip().lower()
            taxon_id = (row.get("taxon_id") or "").strip()
            row_source_version = (row.get("source_version") or "").strip()

            if legacy_rna_id in seen_legacy:
                duplicate_legacy += 1
                continue
            seen_legacy.add(legacy_rna_id)

            out, strategy = canonicalize_row(
                legacy_rna_id=legacy_rna_id,
                rna_type=rna_type,
                taxon_id=taxon_id,
                row_source_version=row_source_version,
                xref_map=xref_map,
                default_source=args.source,
                fetch_date=args.fetch_date,
            )
            rows_out.append(out)
            total += 1

            canonical_rna_id = out[1]
            canonical_system = out[2]
            confidence = out[5]

            strategy_counts[strategy] += 1
            canonical_system_counts[canonical_system] += 1
            confidence_counts[confidence] += 1
            rna_type_counts[rna_type] += 1
            if canonical_rna_id != "":
                canonical_non_empty_by_type[rna_type] += 1
            if canonical_rna_id != legacy_rna_id:
                changed_by_type[rna_type] += 1

            if args.max_rows is not None and total >= args.max_rows:
                break

    rows_out.sort(key=lambda x: x[0])
    out_rows = write_tsv(args.out, header=header, rows=rows_out)

    by_type = {}
    for rt, n in sorted(rna_type_counts.items()):
        non_empty = canonical_non_empty_by_type.get(rt, 0)
        changed = changed_by_type.get(rt, 0)
        by_type[rt] = {
            "rows": n,
            "canonical_non_empty": non_empty,
            "canonical_non_empty_rate": (non_empty / n) if n else 0.0,
            "canonical_changed": changed,
            "canonical_changed_rate": (changed / n) if n else 0.0,
        }

    report = {
        "name": "rna_id_canonical_map_v1.build",
        "generated_at": utc_now_iso(),
        "inputs": {
            "master": str(args.master),
            "xref": str(args.xref),
            "sample_mode": args.max_rows is not None,
            "max_rows": args.max_rows,
        },
        "counts": {
            "rows_output": out_rows,
            "legacy_coverage_total": total,
            "legacy_unique_rows": len(seen_legacy),
            "legacy_duplicate_rows_skipped": duplicate_legacy,
            "canonical_non_empty_total": sum(canonical_non_empty_by_type.values()),
            "canonical_non_empty_rate": (sum(canonical_non_empty_by_type.values()) / total) if total else 0.0,
        },
        "distribution": {
            "strategy_counts": dict(strategy_counts),
            "canonical_system_counts": dict(canonical_system_counts),
            "confidence_counts": dict(confidence_counts),
            "by_rna_type": by_type,
        },
        "xref_stats": xref_stats,
        "gates": {
            "legacy_coverage_100": total == len(seen_legacy),
            "canonical_non_empty_100": sum(canonical_non_empty_by_type.values()) == total,
        },
    }
    write_json(args.report, report)

    print(f"[OK] wrote canonical map: {args.out} rows={out_rows}")
    print(f"[OK] wrote build report: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
