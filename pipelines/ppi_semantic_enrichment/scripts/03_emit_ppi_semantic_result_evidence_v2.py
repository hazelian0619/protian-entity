#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


SOURCE_CANDIDATE_URLS: Dict[str, List[str]] = {
    "string_detailed": [
        "https://stringdb-static.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz",
    ],
    "intact_mitab": [
        "https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/human.zip",
        "https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip",
    ],
    "biogrid_tab3": [
        "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip",
        "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ALL-LATEST.tab3.zip",
    ],
}

SOURCE_VERSION_HINTS: Dict[str, str] = {
    "string_detailed": "STRING v12.0",
    "intact_mitab": "IntAct current human MITAB",
    "biogrid_tab3": "BioGRID Latest TAB3",
}


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def safe_rate(num: int, den: int) -> float:
    if den <= 0:
        return 1.0
    return num / den


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(list(header))
        for row in rows:
            w.writerow(list(row))
            n += 1
    return n


def reservoir_update(reservoir: List[Dict[str, str]], row: Dict[str, str], seen: int, k: int, rng: random.Random) -> None:
    if len(reservoir) < k:
        reservoir.append(row)
        return
    j = rng.randint(1, seen)
    if j <= k:
        reservoir[j - 1] = row


def load_qa_metrics(path: Path) -> Dict[str, object]:
    obj = json.loads(path.read_text(encoding="utf-8"))
    return {
        "passed": bool(obj.get("passed")),
        "metrics": obj.get("metrics", {}),
        "gates": obj.get("gates", []),
    }


def analyze_method_table(path: Path, sample_size: int, seed: int) -> Dict[str, object]:
    coverage_fields = [
        "method",
        "throughput",
        "pmid",
        "doi",
        "text_mining_score",
        "experimental_score",
    ]
    counts = Counter()
    combination_counter = Counter()
    exploded_counter = Counter()
    source_missing = 0
    sample_rows: List[Dict[str, str]] = []
    rng = random.Random(seed)

    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = ["edge_id", *coverage_fields, "source_databases"]
        missing = sorted(set(required) - set(r.fieldnames or []))
        if missing:
            raise SystemExit(f"[ERROR] method table missing columns: {missing}")

        for i, row in enumerate(r, start=1):
            counts["rows"] += 1
            for field in coverage_fields:
                if (row.get(field) or "").strip():
                    counts[f"{field}_non_empty"] += 1

            raw_sources = (row.get("source_databases") or "").strip()
            if raw_sources:
                parts = sorted({x.strip() for x in raw_sources.split(";") if x.strip()})
                if parts:
                    combination_counter["+".join(parts)] += 1
                    for part in parts:
                        exploded_counter[part] += 1
                else:
                    source_missing += 1
            else:
                source_missing += 1

            reservoir_update(sample_rows, row, i, sample_size, rng)

    rows = counts["rows"]
    coverage = {f: safe_rate(counts[f"{f}_non_empty"], rows) for f in coverage_fields}
    return {
        "rows": rows,
        "coverage": coverage,
        "source_distribution_combination": dict(combination_counter.most_common()),
        "source_distribution_exploded": dict(exploded_counter.most_common()),
        "source_missing_rows": source_missing,
        "sample_rows": sample_rows,
    }


def analyze_function_table(path: Path, sample_size: int, seed: int) -> Dict[str, object]:
    sample_rows: List[Dict[str, str]] = []
    source_counter = Counter()
    counts = Counter()
    rng = random.Random(seed)

    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = ["edge_id", "source"]
        missing = sorted(set(required) - set(r.fieldnames or []))
        if missing:
            raise SystemExit(f"[ERROR] function table missing columns: {missing}")

        for i, row in enumerate(r, start=1):
            counts["rows"] += 1
            source_counter[(row.get("source") or "").strip() or ""] += 1
            reservoir_update(sample_rows, row, i, sample_size, rng)

    return {
        "rows": counts["rows"],
        "source_distribution": dict(source_counter.most_common()),
        "sample_rows": sample_rows,
    }


def build_blocked_inputs_report(build_report_path: Path) -> Dict[str, object]:
    build = json.loads(build_report_path.read_text(encoding="utf-8"))
    blocked_items: List[Dict[str, object]] = []

    for item in build.get("external_sources", []):
        name = str(item.get("name") or "")
        status = str(item.get("status") or "")
        path = Path(str(item.get("path") or ""))
        exists = path.exists() and path.is_file() and path.stat().st_size > 0 if path else False

        if status == "missing" or not exists:
            blocked_items.append(
                {
                    "name": name,
                    "status": status or "missing",
                    "path": str(path),
                    "exists": exists,
                    "version_hint": SOURCE_VERSION_HINTS.get(name, ""),
                    "download_urls": SOURCE_CANDIDATE_URLS.get(name, []),
                    "sha256": sha256_file(path) if exists else None,
                    "size_bytes": path.stat().st_size if exists else 0,
                }
            )

    return {
        "name": "ppi_semantic_enrichment_v2.blocked_missing_inputs",
        "created_at": utc_now_iso(),
        "blocked_count": len(blocked_items),
        "blocked_items": blocked_items,
        "status": "blocked" if blocked_items else "clear",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--method-table", type=Path, required=True)
    ap.add_argument("--function-table", type=Path, required=True)
    ap.add_argument("--qa-report", type=Path, required=True)
    ap.add_argument("--build-report", type=Path, required=True)
    ap.add_argument("--sample-size", type=int, default=20)
    ap.add_argument("--sample-seed", type=int, default=20260412)
    ap.add_argument("--out-evidence-report", type=Path, required=True)
    ap.add_argument("--out-blocked-report", type=Path, required=True)
    ap.add_argument("--out-method-sample", type=Path, required=True)
    ap.add_argument("--out-function-sample", type=Path, required=True)
    args = ap.parse_args()

    for p in [args.method_table, args.function_table, args.qa_report, args.build_report]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    method_info = analyze_method_table(args.method_table, sample_size=args.sample_size, seed=args.sample_seed)
    function_info = analyze_function_table(args.function_table, sample_size=args.sample_size, seed=args.sample_seed + 1)
    qa_info = load_qa_metrics(args.qa_report)
    blocked = build_blocked_inputs_report(args.build_report)

    method_sample_header = list((method_info["sample_rows"][0].keys() if method_info["sample_rows"] else []))
    function_sample_header = list((function_info["sample_rows"][0].keys() if function_info["sample_rows"] else []))

    method_sample_rows_written = write_tsv(
        args.out_method_sample,
        method_sample_header,
        ([row.get(c, "") for c in method_sample_header] for row in method_info["sample_rows"]),
    )
    function_sample_rows_written = write_tsv(
        args.out_function_sample,
        function_sample_header,
        ([row.get(c, "") for c in function_sample_header] for row in function_info["sample_rows"]),
    )

    evidence_report = {
        "name": "ppi_semantic_enrichment_v2.result_evidence",
        "created_at": utc_now_iso(),
        "outputs": {
            "method_table": str(args.method_table),
            "function_table": str(args.function_table),
            "method_sample20": str(args.out_method_sample),
            "function_sample20": str(args.out_function_sample),
        },
        "row_counts": {
            "ppi_method_context_v2": method_info["rows"],
            "ppi_function_context_v2": function_info["rows"],
            "method_sample_rows": method_sample_rows_written,
            "function_sample_rows": function_sample_rows_written,
        },
        "coverage_rates": {
            "method": method_info["coverage"],
        },
        "source_distribution": {
            "method_source_databases_combination": method_info["source_distribution_combination"],
            "method_source_databases_exploded": method_info["source_distribution_exploded"],
            "function_source": function_info["source_distribution"],
        },
        "acceptance": {
            "qa_passed": qa_info["passed"],
            "qa_metrics": qa_info["metrics"],
            "qa_gates": qa_info["gates"],
            "blocked_missing_inputs_status": blocked["status"],
            "blocked_missing_inputs_count": blocked["blocked_count"],
        },
    }

    args.out_evidence_report.parent.mkdir(parents=True, exist_ok=True)
    args.out_evidence_report.write_text(json.dumps(evidence_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    args.out_blocked_report.parent.mkdir(parents=True, exist_ok=True)
    args.out_blocked_report.write_text(json.dumps(blocked, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"[OK] result evidence report -> {args.out_evidence_report}")
    print(f"[OK] blocked missing inputs -> {args.out_blocked_report} ({blocked['status']})")
    print(f"[OK] sample rows -> {args.out_method_sample} ({method_sample_rows_written}), {args.out_function_sample} ({function_sample_rows_written})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
