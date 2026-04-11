#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import json
import re
import threading
import time
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

try:
    import requests
except ImportError as exc:  # pragma: no cover
    raise SystemExit("Missing dependency: requests. Install with `python3 -m pip install requests`.") from exc


OUTPUT_COLUMNS = [
    "uniprot_id",
    "interpro_id",
    "pfam_id",
    "entry_name",
    "start",
    "end",
    "source_version",
    "fetch_date",
]

PFAM_RE = re.compile(r"^PF\d{5}$")
INTERPRO_RE = re.compile(r"^IPR\d{6}$")
RETRYABLE_STATUSES = {408, 425, 429, 500, 502, 503, 504}
MANUAL_DOWNLOAD_FAILURE_RATE_THRESHOLD = 0.05

_thread_local = threading.local()


@dataclass
class FetchResult:
    uniprot_id: str
    status: int
    rows: List[Dict[str, str]]
    source_version: Optional[str]
    error: Optional[str] = None


def _session() -> requests.Session:
    s = getattr(_thread_local, "session", None)
    if s is None:
        s = requests.Session()
        s.headers.update(
            {
                "User-Agent": "protian-entity/protein_domains_pipeline (contact: local)",
                "Accept": "application/json",
            }
        )
        _thread_local.session = s
    return s


def _iter_input_rows(input_path: Path) -> Iterable[Dict[str, str]]:
    with input_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"Missing header: {input_path}")
        required = {"uniprot_id", "domains"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise SystemExit(f"Input missing required columns: {sorted(missing)}")
        for row in reader:
            yield row


def load_uniprot_ids(input_path: Path, limit: Optional[int]) -> tuple[List[str], Dict[str, str]]:
    ids: List[str] = []
    domains_lookup: Dict[str, str] = {}
    seen: set[str] = set()

    for row in _iter_input_rows(input_path):
        uid = (row.get("uniprot_id") or "").strip().upper()
        if not uid:
            continue
        if uid in seen:
            continue
        seen.add(uid)
        ids.append(uid)
        domains_lookup[uid] = (row.get("domains") or "")
        if limit is not None and len(ids) >= limit:
            break

    return ids, domains_lookup


def _normalize_interpro_id(raw: Any) -> str:
    val = str(raw or "").strip()
    return val if INTERPRO_RE.match(val) else ""


def _fetch_one(
    uniprot_id: str,
    timeout: float,
    max_retries: int,
    sleep_base: float,
) -> FetchResult:
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{uniprot_id}/"

    last_error: Optional[str] = None
    for attempt in range(max_retries + 1):
        try:
            resp = _session().get(url, timeout=timeout)
            status = int(resp.status_code)
            version = resp.headers.get("InterPro-Version")

            if status == 204:
                return FetchResult(uniprot_id=uniprot_id, status=204, rows=[], source_version=version)

            if status == 200:
                payload = resp.json()
                results = payload.get("results", []) or []
                rows: List[Dict[str, str]] = []

                for item in results:
                    metadata = item.get("metadata", {}) or {}
                    pfam_id = str(metadata.get("accession") or "").strip()
                    if not PFAM_RE.match(pfam_id):
                        continue
                    interpro_id = _normalize_interpro_id(metadata.get("integrated"))
                    entry_name = str(metadata.get("name") or "").strip()

                    proteins = item.get("proteins", []) or []
                    for protein in proteins:
                        locations = protein.get("entry_protein_locations", []) or []
                        for loc in locations:
                            fragments = loc.get("fragments", []) or []
                            for frag in fragments:
                                start = frag.get("start")
                                end = frag.get("end")
                                if start is None or end is None:
                                    continue
                                try:
                                    start_i = int(start)
                                    end_i = int(end)
                                except Exception:
                                    continue
                                if start_i <= 0 or end_i <= 0 or end_i < start_i:
                                    continue

                                rows.append(
                                    {
                                        # Keep queried UniProt ID for stable join-back to master table.
                                        "uniprot_id": uniprot_id,
                                        "interpro_id": interpro_id,
                                        "pfam_id": pfam_id,
                                        "entry_name": entry_name,
                                        "start": str(start_i),
                                        "end": str(end_i),
                                    }
                                )

                return FetchResult(uniprot_id=uniprot_id, status=200, rows=rows, source_version=version)

            if status in RETRYABLE_STATUSES and attempt < max_retries:
                time.sleep(sleep_base * (attempt + 1))
                continue

            return FetchResult(
                uniprot_id=uniprot_id,
                status=status,
                rows=[],
                source_version=version,
                error=f"HTTP {status}",
            )

        except Exception as exc:
            last_error = str(exc)
            if attempt < max_retries:
                time.sleep(sleep_base * (attempt + 1))
                continue

    return FetchResult(
        uniprot_id=uniprot_id,
        status=-1,
        rows=[],
        source_version=None,
        error=last_error or "unknown error",
    )


def write_tsv(rows: List[Dict[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in OUTPUT_COLUMNS})


def _safe_rate(numerator: int, denominator: int) -> float:
    return (numerator / denominator) if denominator else 0.0


def build_manual_download_plan() -> Dict[str, Any]:
    return {
        "target_local_dir": "data/raw/protein_domains/",
        "files": [
            {
                "name": "protein2ipr.dat.gz",
                "url": "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz",
                "expected_size": "16G",
                "purpose": "Protein->InterPro/Pfam mapping with residue start/end positions",
            },
            {
                "name": "entry.list",
                "url": "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list",
                "expected_size": "2.7M",
                "purpose": "InterPro entry metadata (accession/name/type)",
            },
            {
                "name": "ParentChildTreeFile.txt",
                "url": "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/ParentChildTreeFile.txt",
                "expected_size": "632K",
                "purpose": "InterPro hierarchy helper file for downstream QA (optional)",
            },
        ],
        "sha256_check_commands": [
            "mkdir -p data/raw/protein_domains",
            "cd data/raw/protein_domains",
            "sha256sum protein2ipr.dat.gz entry.list ParentChildTreeFile.txt > SHA256SUMS.txt",
            "sha256sum -c SHA256SUMS.txt",
        ],
        "note": "File sizes are from InterPro current_release index at pipeline build time and may vary by release.",
    }


def build_metrics(
    input_ids: List[str],
    domains_lookup: Dict[str, str],
    output_rows: List[Dict[str, str]],
    status_counter: Counter,
    failed: Dict[str, str],
    source_versions: List[str],
    elapsed_seconds: float,
    sample_mode: bool,
) -> Dict[str, Any]:
    input_set = set(input_ids)

    output_ids = [r["uniprot_id"] for r in output_rows if r.get("uniprot_id")]
    output_set = set(output_ids)

    # join quality
    joinable_ids = {uid for uid in output_set if uid in input_set}
    non_joinable_ids = sorted(output_set - input_set)
    join_pass_rate = _safe_rate(len(joinable_ids), len(output_set)) if output_set else 1.0

    # coverage quality
    master_coverage_rate = _safe_rate(len(joinable_ids), len(input_set))

    with_domain_text_ids = {
        uid
        for uid in input_ids
        if (domains_lookup.get(uid) or "").strip() not in ("", "NA", "N/A", "None", "null")
    }
    covered_with_domain_text = with_domain_text_ids & output_set
    with_domain_text_coverage_rate = _safe_rate(len(covered_with_domain_text), len(with_domain_text_ids))

    missing_ids = sorted(input_set - output_set)

    # field format checks
    pfam_values = [r.get("pfam_id", "") for r in output_rows]
    interpro_values = [r.get("interpro_id", "") for r in output_rows]

    pfam_valid = sum(1 for v in pfam_values if PFAM_RE.match(v or ""))
    interpro_present = [v for v in interpro_values if v]
    interpro_valid = sum(1 for v in interpro_present if INTERPRO_RE.match(v or ""))
    start_end_valid = sum(
        1
        for r in output_rows
        if str(r.get("start", "")).isdigit()
        and str(r.get("end", "")).isdigit()
        and int(r["start"]) > 0
        and int(r["end"]) >= int(r["start"])
    )

    rps = _safe_rate(len(input_ids), int(elapsed_seconds)) if elapsed_seconds >= 1 else float(len(input_ids))

    acceptance = {
        "join_pass_rate_gte_0_99": join_pass_rate >= 0.99,
        "pfam_format_pass": (pfam_valid == len(pfam_values)) if pfam_values else True,
        "interpro_format_pass": (interpro_valid == len(interpro_present)) if interpro_present else True,
    }

    return {
        "pipeline": "protein_domains_interpro_v1",
        "sample_mode": sample_mode,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "row_count": {
            "input_uniprot_ids": len(input_ids),
            "output_rows": len(output_rows),
            "output_uniprot_ids": len(output_set),
        },
        "api": {
            "status_counts": dict(status_counter),
            "failed_count": len(failed),
            "failed_uniprot_ids_top50": [
                {"uniprot_id": uid, "error": failed[uid]} for uid in sorted(failed)[:50]
            ],
            "source_versions": sorted(set(v for v in source_versions if v)),
            "elapsed_seconds": elapsed_seconds,
            "throughput_uniprot_per_sec": rps,
        },
        "coverage": {
            "join_pass_rate": join_pass_rate,
            "master_coverage_rate": master_coverage_rate,
            "domain_text_input_count": len(with_domain_text_ids),
            "domain_text_covered_count": len(covered_with_domain_text),
            "domain_text_coverage_rate": with_domain_text_coverage_rate,
        },
        "format_validation": {
            "pfam_valid_rows": pfam_valid,
            "pfam_total_rows": len(pfam_values),
            "pfam_valid_rate": _safe_rate(pfam_valid, len(pfam_values)),
            "interpro_present_rows": len(interpro_present),
            "interpro_valid_rows": interpro_valid,
            "interpro_valid_rate": _safe_rate(interpro_valid, len(interpro_present)),
            "start_end_valid_rows": start_end_valid,
            "start_end_valid_rate": _safe_rate(start_end_valid, len(output_rows)),
        },
        "missing_analysis": {
            "missing_uniprot_count": len(missing_ids),
            "missing_uniprot_ids_top50": missing_ids[:50],
            "non_joinable_uniprot_count": len(non_joinable_ids),
            "non_joinable_uniprot_ids_top50": non_joinable_ids[:50],
        },
        "acceptance": acceptance,
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build Protein InterPro/Pfam domain mapping table.")
    p.add_argument(
        "--input",
        type=Path,
        default=Path("data/processed/protein_master_v6_clean.tsv"),
    )
    p.add_argument(
        "--output",
        type=Path,
        default=Path("data/output/protein/protein_domains_interpro_v1.tsv"),
    )
    p.add_argument(
        "--report",
        type=Path,
        default=Path("pipelines/protein_domains/reports/protein_domains_interpro_v1.metrics.json"),
    )
    p.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional: only process first N unique UniProt IDs (sample run).",
    )
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--timeout", type=float, default=25.0)
    p.add_argument("--max-retries", type=int, default=3)
    p.add_argument("--sleep-base", type=float, default=0.6)
    p.add_argument(
        "--audit",
        type=Path,
        default=Path("pipelines/protein_domains/reports/protein_domains_interpro_v1.audit.json"),
    )
    p.add_argument(
        "--min-throughput",
        type=float,
        default=5.0,
        help="If full-run throughput (uniprot/sec) is lower than this threshold, mark API mode as too slow.",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if not args.input.exists():
        raise SystemExit(f"Input not found: {args.input}")

    input_ids, domains_lookup = load_uniprot_ids(args.input, args.limit)
    if not input_ids:
        raise SystemExit("No UniProt IDs found in input.")

    started = time.time()

    status_counter: Counter = Counter()
    failed: Dict[str, str] = {}
    source_versions: List[str] = []
    rows_out: List[Dict[str, str]] = []

    total = len(input_ids)
    print(f"[INFO] Fetching InterPro/Pfam for {total} UniProt IDs with workers={args.workers}...")

    done = 0
    with cf.ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futures = {
            ex.submit(_fetch_one, uid, args.timeout, args.max_retries, args.sleep_base): uid
            for uid in input_ids
        }

        for fut in cf.as_completed(futures):
            uid = futures[fut]
            done += 1
            try:
                res = fut.result()
            except Exception as exc:  # pragma: no cover
                status_counter[-1] += 1
                failed[uid] = f"executor exception: {exc}"
                continue

            status_counter[res.status] += 1
            if res.source_version:
                source_versions.append(res.source_version)
            if res.error:
                failed[uid] = res.error

            for row in res.rows:
                row["source_version"] = (res.source_version or "")
                rows_out.append(row)

            if done % 500 == 0 or done == total:
                print(f"[INFO] progress: {done}/{total}")

    fetch_date = datetime.now(timezone.utc).date().isoformat()
    for row in rows_out:
        row["fetch_date"] = fetch_date

    rows_out.sort(
        key=lambda r: (
            r.get("uniprot_id", ""),
            int(r.get("start", "0")),
            int(r.get("end", "0")),
            r.get("pfam_id", ""),
            r.get("interpro_id", ""),
        )
    )

    write_tsv(rows_out, args.output)

    elapsed = time.time() - started
    metrics = build_metrics(
        input_ids=input_ids,
        domains_lookup=domains_lookup,
        output_rows=rows_out,
        status_counter=status_counter,
        failed=failed,
        source_versions=source_versions,
        elapsed_seconds=elapsed,
        sample_mode=args.limit is not None,
    )

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(metrics, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    input_set = set(input_ids)
    output_set = {r["uniprot_id"] for r in rows_out if r.get("uniprot_id")}
    missing_uniprot_ids = sorted(input_set - output_set)
    non_joinable_uniprot_ids = sorted(output_set - input_set)
    failed_rate = _safe_rate(len(failed), len(input_ids))
    throughput = float(metrics["api"]["throughput_uniprot_per_sec"])
    manual_download_required = (
        args.limit is None
        and (failed_rate > MANUAL_DOWNLOAD_FAILURE_RATE_THRESHOLD or throughput < args.min_throughput)
    )
    manual_download_plan = build_manual_download_plan() if manual_download_required else None

    audit = {
        "pipeline": "protein_domains_interpro_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_uniprot_count": len(input_ids),
        "output_uniprot_count": len(output_set),
        "status_counts": dict(status_counter),
        "failed_count": len(failed),
        "failed_uniprot_ids": [{"uniprot_id": uid, "error": failed[uid]} for uid in sorted(failed)],
        "missing_uniprot_count": len(missing_uniprot_ids),
        "missing_uniprot_ids": missing_uniprot_ids,
        "non_joinable_uniprot_count": len(non_joinable_uniprot_ids),
        "non_joinable_uniprot_ids": non_joinable_uniprot_ids,
        "manual_download_required": manual_download_required,
        "manual_download_reason": {
            "failed_rate": failed_rate,
            "failed_rate_threshold": MANUAL_DOWNLOAD_FAILURE_RATE_THRESHOLD,
            "throughput_uniprot_per_sec": throughput,
            "min_throughput": args.min_throughput,
        },
        "manual_download_plan": manual_download_plan,
    }
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    args.audit.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    ok = (
        metrics["acceptance"]["join_pass_rate_gte_0_99"]
        and metrics["acceptance"]["pfam_format_pass"]
        and metrics["acceptance"]["interpro_format_pass"]
    )

    status = "PASS" if ok else "FAIL"
    print(
        f"[{status}] rows={metrics['row_count']['output_rows']} "
        f"join={metrics['coverage']['join_pass_rate']:.4f} "
        f"master_coverage={metrics['coverage']['master_coverage_rate']:.4f} "
        f"pfam_valid={metrics['format_validation']['pfam_valid_rate']:.4f} "
        f"interpro_valid={metrics['format_validation']['interpro_valid_rate']:.4f}"
    )

    if failed:
        print(f"[WARN] {len(failed)} UniProt IDs failed after retry; see report for details.")

    if manual_download_required:
        print("[INTERRUPT] API mode is too slow or unstable. Please switch to bulk files:")
        if manual_download_plan:
            print(json.dumps(manual_download_plan, ensure_ascii=False, indent=2))
        return 3

    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
