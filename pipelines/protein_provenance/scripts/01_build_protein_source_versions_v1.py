#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


@dataclass(frozen=True)
class DatasetSpec:
    dataset_name: str
    rel_path: str
    scan_columns: List[str]


DATASET_SPECS: List[DatasetSpec] = [
    DatasetSpec("protein_master_v6_clean.tsv", "data/processed/protein_master_v6_clean.tsv", ["source", "fetch_date", "date_modified"]),
    DatasetSpec("protein_edges.tsv", "data/processed/protein_edges.tsv", ["source", "fetch_date"]),
    DatasetSpec("ptm_sites.tsv", "data/processed/ptm_sites.tsv", ["source", "fetch_date"]),
    DatasetSpec("pathway_members.tsv", "data/processed/pathway_members.tsv", ["source", "fetch_date"]),
    DatasetSpec("alphafold_quality.tsv", "data/processed/alphafold_quality.tsv", ["alphafold_version"]),
    DatasetSpec("hgnc_core.tsv", "data/processed/hgnc_core.tsv", ["date_modified"]),
]


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _non_empty(v: str) -> bool:
    return (v or "").strip() not in ("", "NA", "N/A", "None", "null")


def _scan_tsv(path: Path, scan_columns: Iterable[str], max_rows: Optional[int]) -> Dict[str, Any]:
    scan_columns = list(scan_columns)
    counters: Dict[str, Counter[str]] = {c: Counter() for c in scan_columns}

    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        missing_cols = [c for c in scan_columns if c not in (r.fieldnames or [])]
        if missing_cols:
            raise SystemExit(f"[ERROR] {path} missing required columns for provenance scan: {missing_cols}")

        rows_scanned = 0
        for row in r:
            rows_scanned += 1
            for col in scan_columns:
                v = (row.get(col) or "").strip()
                if _non_empty(v):
                    counters[col][v] += 1
            if max_rows is not None and rows_scanned >= max_rows:
                break

    return {
        "path": str(path),
        "rows_scanned": rows_scanned,
        "sampled": max_rows is not None,
        "max_rows": max_rows,
        "columns": {
            c: {
                "non_empty": int(sum(counters[c].values())),
                "unique_non_empty": int(len(counters[c])),
                "top_values": [{"value": k, "count": v} for k, v in counters[c].most_common(10)],
                "max_value": (max(counters[c]) if counters[c] else None),
                "min_value": (min(counters[c]) if counters[c] else None),
            }
            for c in scan_columns
        },
    }


def _single_value(profile: Dict[str, Any], col: str) -> Optional[str]:
    top = profile["columns"][col]["top_values"]
    if len(top) == 1:
        return str(top[0]["value"])
    return None


def _max_date_like(profile: Dict[str, Any], col: str) -> Optional[str]:
    values = [x["value"] for x in profile["columns"][col]["top_values"]]
    max_value = profile["columns"][col]["max_value"]
    if max_value and DATE_RE.match(str(max_value)):
        return str(max_value)
    for v in sorted(values, reverse=True):
        if DATE_RE.match(str(v)):
            return str(v)
    return None


def _safe_version(value: Optional[str], fallback: str = "unknown") -> str:
    v = (value or "").strip()
    return v if v else fallback


def _build_rows(scans: Dict[str, Dict[str, Any]]) -> List[Dict[str, str]]:
    master = scans["protein_master_v6_clean.tsv"]
    edges = scans["protein_edges.tsv"]
    ptm = scans["ptm_sites.tsv"]
    pathway = scans["pathway_members.tsv"]
    alphafold = scans["alphafold_quality.tsv"]
    hgnc = scans["hgnc_core.tsv"]

    master_fetch_date = _single_value(master, "fetch_date") or "unknown"
    master_uniprot_max_date = _max_date_like(master, "date_modified") or "unknown"

    edges_source_version = _single_value(edges, "source")
    edges_fetch_date = _single_value(edges, "fetch_date") or master_fetch_date

    ptm_fetch_date = _single_value(ptm, "fetch_date") or master_fetch_date
    pathway_fetch_date = _single_value(pathway, "fetch_date") or master_fetch_date

    alphafold_version = _single_value(alphafold, "alphafold_version")

    hgnc_max_date = _max_date_like(hgnc, "date_modified")
    hgnc_month = hgnc_max_date[:7] if hgnc_max_date and DATE_RE.match(hgnc_max_date) else "unknown"

    master_source_version = (
        "v6_clean"
        f"|UniProt:max_date_modified={_safe_version(master_uniprot_max_date)}"
        f"|AlphaFold:{_safe_version(alphafold_version)}"
        f"|STRING:{_safe_version(edges_source_version)}"
        f"|HGNC:snapshot_{hgnc_month}"
    )

    rows: List[Dict[str, str]] = [
        {
            "dataset": "protein_master_v6_clean.tsv",
            "primary_source": "UniProtKB (Swiss-Prot-centric, integrated)",
            "source_version": master_source_version,
            "fetch_date": master_fetch_date,
            "evidence_field": "source,date_modified,fetch_date + alphafold_quality.alphafold_version + protein_edges.source + hgnc_core.date_modified",
            "notes": "Composite anchor rule: one row aggregates component versions; keep source_version in sidecar table to avoid changing existing protein_master contracts/joins.",
        },
        {
            "dataset": "protein_edges.tsv",
            "primary_source": "STRING",
            "source_version": _safe_version(edges_source_version),
            "fetch_date": edges_fetch_date,
            "evidence_field": "source,fetch_date",
            "notes": "Raw source already carries explicit version in source (e.g., STRING_v12.0).",
        },
        {
            "dataset": "ptm_sites.tsv",
            "primary_source": "PhosphoSitePlus",
            "source_version": f"snapshot@{_safe_version(ptm_fetch_date)}",
            "fetch_date": ptm_fetch_date,
            "evidence_field": "source,fetch_date",
            "notes": "No explicit upstream release code in table; lock by fetch_date snapshot.",
        },
        {
            "dataset": "pathway_members.tsv",
            "primary_source": "Reactome",
            "source_version": f"snapshot@{_safe_version(pathway_fetch_date)}",
            "fetch_date": pathway_fetch_date,
            "evidence_field": "source,fetch_date",
            "notes": "No explicit release code in table; lock by fetch_date snapshot.",
        },
        {
            "dataset": "alphafold_quality.tsv",
            "primary_source": "AlphaFold DB",
            "source_version": _safe_version(alphafold_version),
            "fetch_date": master_fetch_date,
            "evidence_field": "alphafold_version (+protein_master.fetch_date fallback)",
            "notes": "alphafold_quality has explicit version but no fetch_date column; inherit fetch_date from protein_master_v6_clean.tsv.",
        },
        {
            "dataset": "hgnc_core.tsv",
            "primary_source": "HGNC",
            "source_version": f"snapshot_{hgnc_month}(max_date_modified={_safe_version(hgnc_max_date)})",
            "fetch_date": master_fetch_date,
            "evidence_field": "date_modified (+protein_master.fetch_date fallback)",
            "notes": "HGNC file has per-gene date_modified, not a single release field; use latest observed date as snapshot anchor and inherit fetch_date.",
        },
    ]

    return rows


def _write_tsv(path: Path, rows: List[Dict[str, str]]) -> None:
    fieldnames = ["dataset", "primary_source", "source_version", "fetch_date", "evidence_field", "notes"]
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fieldnames})
    tmp.replace(path)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--out-build-report", type=Path, required=True)
    ap.add_argument("--max-rows-per-input", type=int, default=None, help="Optional smoke mode: scan only first N rows per input TSV")
    args = ap.parse_args()

    repo_root = args.repo_root.resolve()

    scans: Dict[str, Dict[str, Any]] = {}
    for spec in DATASET_SPECS:
        path = repo_root / spec.rel_path
        if not path.exists():
            raise SystemExit(f"[ERROR] missing input: {path}")
        scans[spec.dataset_name] = _scan_tsv(path, spec.scan_columns, max_rows=args.max_rows_per_input)

    rows = _build_rows(scans)
    _write_tsv(args.out, rows)

    build_report = {
        "name": "protein_source_versions_v1",
        "created_at": _utc_now(),
        "output": str(args.out),
        "rows_written": len(rows),
        "sampled_mode": args.max_rows_per_input is not None,
        "max_rows_per_input": args.max_rows_per_input,
        "scan_profiles": scans,
        "preview_rows": rows,
    }

    args.out_build_report.parent.mkdir(parents=True, exist_ok=True)
    args.out_build_report.write_text(json.dumps(build_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    mode = "SMOKE" if args.max_rows_per_input is not None else "FULL"
    print(f"[OK] {mode} build -> {args.out} (rows={len(rows)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
