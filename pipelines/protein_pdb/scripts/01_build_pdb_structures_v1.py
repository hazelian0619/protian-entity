#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"
RCSB_CORE_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
USER_AGENT = "kg-protein-pdb/1.0"
SOURCE = "RCSB_DATA_API"

GRAPHQL_QUERY = """
query($ids:[String!]!) {
  entries(entry_ids:$ids) {
    rcsb_id
    exptl { method }
    rcsb_entry_info {
      resolution_combined
      nonpolymer_entity_count
    }
    rcsb_accession_info {
      initial_release_date
    }
  }
}
""".strip()


@dataclass
class PairRow:
    uniprot_id: str
    pdb_id: str


@dataclass
class EntryMeta:
    experimental_method: str
    resolution: str
    release_date: str
    ligand_count: str


@dataclass
class FetchStats:
    total_unique_pdb: int
    batches_total: int = 0
    batches_failed: int = 0
    graphql_requests: int = 0
    graphql_failures: int = 0
    core_fallback_requests: int = 0
    core_fallback_failures: int = 0


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def utc_today() -> str:
    return datetime.now(timezone.utc).date().isoformat()


def _ssl_context() -> ssl.SSLContext:
    try:
        import certifi  # type: ignore

        return ssl.create_default_context(cafile=certifi.where())
    except Exception:
        return ssl.create_default_context()


def _json_post(url: str, payload: Dict[str, Any], timeout: int, ctx: ssl.SSLContext) -> Dict[str, Any]:
    req = urllib.request.Request(
        url,
        data=json.dumps(payload, ensure_ascii=False).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        },
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _json_get(url: str, timeout: int, ctx: ssl.SSLContext) -> Dict[str, Any]:
    req = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        },
        method="GET",
    )
    with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
        return json.loads(resp.read().decode("utf-8"))


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


def load_master_pairs(path: Path) -> Tuple[List[PairRow], List[str], Dict[str, int]]:
    pairs: List[PairRow] = []
    unique_pdb_order: List[str] = []
    unique_seen = set()
    pair_seen = set()

    metrics = {
        "master_rows": 0,
        "rows_with_pdb": 0,
        "pairs_raw": 0,
        "pairs_deduped": 0,
        "pairs_dropped_duplicate": 0,
    }

    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"uniprot_id", "pdb_ids"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] input schema missing columns: {sorted(missing)}")

        for row in reader:
            metrics["master_rows"] += 1
            uniprot_id = (row.get("uniprot_id") or "").strip()
            if not uniprot_id:
                continue
            pdb_ids = parse_pdb_ids((row.get("pdb_ids") or "").strip())
            if not pdb_ids:
                continue
            metrics["rows_with_pdb"] += 1
            for pdb_id in pdb_ids:
                metrics["pairs_raw"] += 1
                pair_key = (uniprot_id, pdb_id)
                if pair_key in pair_seen:
                    metrics["pairs_dropped_duplicate"] += 1
                    continue
                pair_seen.add(pair_key)
                pairs.append(PairRow(uniprot_id=uniprot_id, pdb_id=pdb_id))
                metrics["pairs_deduped"] += 1
                if pdb_id not in unique_seen:
                    unique_seen.add(pdb_id)
                    unique_pdb_order.append(pdb_id)

    return pairs, unique_pdb_order, metrics


def normalize_entry(entry: Dict[str, Any]) -> EntryMeta:
    methods: List[str] = []
    for item in entry.get("exptl") or []:
        m = (item or {}).get("method")
        if m:
            mm = str(m).strip()
            if mm and mm not in methods:
                methods.append(mm)
    experimental_method = ";".join(methods)

    resolution = ""
    resolution_raw = (entry.get("rcsb_entry_info") or {}).get("resolution_combined")
    vals: List[float] = []
    if isinstance(resolution_raw, list):
        for x in resolution_raw:
            if x is None:
                continue
            try:
                vals.append(float(x))
            except Exception:
                continue
    elif resolution_raw is not None:
        try:
            vals.append(float(resolution_raw))
        except Exception:
            pass
    if vals:
        best = min(vals)
        resolution = f"{best:.3f}".rstrip("0").rstrip(".")

    release_date = ""
    release_raw = (entry.get("rcsb_accession_info") or {}).get("initial_release_date")
    if release_raw:
        release_date = str(release_raw).strip()[:10]

    ligand_count = ""
    ligand_raw = (entry.get("rcsb_entry_info") or {}).get("nonpolymer_entity_count")
    if ligand_raw is not None:
        try:
            ligand_count = str(int(ligand_raw))
        except Exception:
            ligand_count = str(ligand_raw).strip()

    return EntryMeta(
        experimental_method=experimental_method,
        resolution=resolution,
        release_date=release_date,
        ligand_count=ligand_count,
    )


def batched(items: Sequence[str], size: int) -> Iterable[List[str]]:
    for i in range(0, len(items), size):
        yield list(items[i : i + size])


def fetch_graphql_batch(
    pdb_ids: List[str],
    timeout: int,
    retries: int,
    retry_sleep: float,
    ctx: ssl.SSLContext,
    stats: FetchStats,
) -> Tuple[Optional[Dict[str, EntryMeta]], Optional[str]]:
    payload = {"query": GRAPHQL_QUERY, "variables": {"ids": pdb_ids}}
    last_error: Optional[str] = None

    for attempt in range(1, retries + 1):
        stats.graphql_requests += 1
        try:
            resp = _json_post(RCSB_GRAPHQL_URL, payload=payload, timeout=timeout, ctx=ctx)
            if resp.get("errors"):
                last_error = f"graphql_errors={resp.get('errors')}"
                raise RuntimeError(last_error)

            out: Dict[str, EntryMeta] = {}
            for entry in (resp.get("data") or {}).get("entries") or []:
                rid = str((entry or {}).get("rcsb_id") or "").upper().strip()
                if not rid:
                    continue
                out[rid] = normalize_entry(entry or {})
            return out, None
        except Exception as e:
            last_error = f"{type(e).__name__}: {e}"
            if attempt < retries:
                time.sleep(retry_sleep * attempt)

    stats.graphql_failures += 1
    return None, (last_error or "unknown_graphql_error")


def fetch_with_partition(
    pdb_ids: List[str],
    timeout: int,
    retries: int,
    retry_sleep: float,
    ctx: ssl.SSLContext,
    stats: FetchStats,
) -> Tuple[Dict[str, EntryMeta], Dict[str, str], Dict[str, str]]:
    """
    Robust fetch strategy:
    1) Try GraphQL on whole batch.
    2) If batch fails, recursively split the batch.
    3) For single ID failure, fallback to core/entry endpoint.
    Returns:
      - found details map
      - missing map {pdb_id: reason}
      - api_error map {pdb_id: reason}
    """

    found: Dict[str, EntryMeta] = {}
    missing: Dict[str, str] = {}
    api_error: Dict[str, str] = {}

    batch_map, batch_err = fetch_graphql_batch(
        pdb_ids=pdb_ids,
        timeout=timeout,
        retries=retries,
        retry_sleep=retry_sleep,
        ctx=ctx,
        stats=stats,
    )

    if batch_map is not None:
        found.update(batch_map)
        returned = set(batch_map.keys())
        for pdb_id in pdb_ids:
            if pdb_id not in returned:
                missing[pdb_id] = "not_found_in_graphql"
        return found, missing, api_error

    # GraphQL failed for this batch
    if len(pdb_ids) == 1:
        pdb_id = pdb_ids[0]
        meta, status = fetch_core_entry(
            pdb_id=pdb_id,
            timeout=timeout,
            retries=max(1, retries - 1),
            retry_sleep=retry_sleep,
            ctx=ctx,
            stats=stats,
        )
        if status == "ok" and meta is not None:
            found[pdb_id] = meta
        elif status == "missing":
            missing[pdb_id] = "not_found_in_core"
        else:
            api_error[pdb_id] = f"graphql_error={batch_err}; core_error={status}"
        return found, missing, api_error

    mid = len(pdb_ids) // 2
    left = pdb_ids[:mid]
    right = pdb_ids[mid:]

    left_found, left_missing, left_api_error = fetch_with_partition(
        pdb_ids=left,
        timeout=timeout,
        retries=retries,
        retry_sleep=retry_sleep,
        ctx=ctx,
        stats=stats,
    )
    right_found, right_missing, right_api_error = fetch_with_partition(
        pdb_ids=right,
        timeout=timeout,
        retries=retries,
        retry_sleep=retry_sleep,
        ctx=ctx,
        stats=stats,
    )

    found.update(left_found)
    found.update(right_found)
    missing.update(left_missing)
    missing.update(right_missing)
    api_error.update(left_api_error)
    api_error.update(right_api_error)
    return found, missing, api_error


def fetch_core_entry(
    pdb_id: str,
    timeout: int,
    retries: int,
    retry_sleep: float,
    ctx: ssl.SSLContext,
    stats: FetchStats,
) -> Tuple[Optional[EntryMeta], str]:
    url = RCSB_CORE_ENTRY_URL.format(pdb_id=urllib.parse.quote(pdb_id))
    last_error = ""
    for attempt in range(1, retries + 1):
        stats.core_fallback_requests += 1
        try:
            entry = _json_get(url=url, timeout=timeout, ctx=ctx)
            return normalize_entry(entry), "ok"
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None, "missing"
            last_error = f"HTTPError:{e.code}"
        except Exception as e:  # noqa: BLE001
            last_error = f"{type(e).__name__}: {e}"

        if attempt < retries:
            time.sleep(retry_sleep * attempt)

    stats.core_fallback_failures += 1
    return None, (last_error or "core_fallback_failed")


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


def write_json(path: Path, obj: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def manual_download_plan() -> Dict[str, Any]:
    return {
        "reason": "RCSB API batch failure rate exceeded threshold; switch to bulk download mode.",
        "download_checklist": [
            {
                "name": "RCSB entries index",
                "url": "https://files.rcsb.org/pub/pdb/derived_data/index/entries.idx",
                "expected_size_hint": "~70MB (varies by release)",
            },
            {
                "name": "RCSB resolution index",
                "url": "https://files.rcsb.org/pub/pdb/derived_data/index/resolu.idx",
                "expected_size_hint": "~10MB (varies by release)",
            },
            {
                "name": "RCSB compound index (for ligand counting fallback)",
                "url": "https://files.rcsb.org/pub/pdb/derived_data/index/compound.idx",
                "expected_size_hint": "~100MB (varies by release)",
            },
        ],
        "place_under": "data/raw/protein_pdb/rcsb_bulk/",
        "checksum": {
            "method": "sha256",
            "commands": [
                "mkdir -p data/raw/protein_pdb/rcsb_bulk",
                "cd data/raw/protein_pdb/rcsb_bulk",
                "shasum -a 256 entries.idx resolu.idx compound.idx > SHA256SUMS.txt",
                "shasum -a 256 -c SHA256SUMS.txt",
            ],
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--report-build", type=Path, required=True)
    ap.add_argument("--report-audit", type=Path, required=True)
    ap.add_argument("--max-unique-pdb", type=int, default=None)
    ap.add_argument("--batch-size", type=int, default=200)
    ap.add_argument("--timeout", type=int, default=30)
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--retry-sleep", type=float, default=1.0)
    ap.add_argument("--sleep-between-batches", type=float, default=0.05)
    ap.add_argument("--fail-threshold", type=float, default=0.05)
    ap.add_argument("--fetch-date", default=utc_today())
    args = ap.parse_args()

    if not args.input.exists():
        raise SystemExit(f"[ERROR] missing input: {args.input}")
    if args.batch_size <= 0:
        raise SystemExit("[ERROR] --batch-size must be > 0")

    created_at = utc_now_iso()
    pairs, unique_pdb_ids, source_metrics = load_master_pairs(args.input)

    selected_unique_ids = unique_pdb_ids
    sample_mode = False
    if args.max_unique_pdb is not None and args.max_unique_pdb >= 0 and args.max_unique_pdb < len(unique_pdb_ids):
        sample_mode = True
        selected_unique_ids = unique_pdb_ids[: args.max_unique_pdb]

    selected_set = set(selected_unique_ids)
    selected_pairs = [p for p in pairs if p.pdb_id in selected_set]

    stats = FetchStats(total_unique_pdb=len(selected_unique_ids))
    ctx = _ssl_context()

    details: Dict[str, EntryMeta] = {}
    missing_ids: Dict[str, str] = {}
    api_error_ids: Dict[str, str] = {}

    batches = list(batched(selected_unique_ids, args.batch_size))
    stats.batches_total = len(batches)

    for idx, batch in enumerate(batches, start=1):
        graphql_failures_before = stats.graphql_failures
        batch_map, batch_missing, batch_api_error = fetch_with_partition(
            pdb_ids=batch,
            timeout=args.timeout,
            retries=args.retries,
            retry_sleep=args.retry_sleep,
            ctx=ctx,
            stats=stats,
        )

        details.update(batch_map)
        missing_ids.update(batch_missing)
        api_error_ids.update(batch_api_error)

        # GraphQL fail counter increments when at least one GraphQL attempt fully failed.
        if stats.graphql_failures > graphql_failures_before:
            stats.batches_failed += 1

        if idx % 10 == 0 or idx == len(batches):
            print(f"[INFO] fetched {idx}/{len(batches)} batches | found={len(details)} missing={len(missing_ids)} api_error={len(api_error_ids)}")

        if args.sleep_between_batches > 0 and idx < len(batches):
            time.sleep(args.sleep_between_batches)

    api_failure_rate = (len(api_error_ids) / len(selected_unique_ids)) if selected_unique_ids else 0.0

    audit_report: Dict[str, Any] = {
        "name": "pdb_structures_v1_api_audit",
        "created_at": created_at,
        "input": str(args.input),
        "output": str(args.output),
        "sample_mode": sample_mode,
        "max_unique_pdb": args.max_unique_pdb,
        "selected_unique_pdb": len(selected_unique_ids),
        "selected_pairs": len(selected_pairs),
        "found_unique_pdb": len(details),
        "missing_unique_pdb": len(missing_ids),
        "api_error_unique_pdb": len(api_error_ids),
        "api_failure_rate": api_failure_rate,
        "fail_threshold": args.fail_threshold,
        "stats": {
            "batches_total": stats.batches_total,
            "batches_failed": stats.batches_failed,
            "graphql_requests": stats.graphql_requests,
            "graphql_failures": stats.graphql_failures,
            "core_fallback_requests": stats.core_fallback_requests,
            "core_fallback_failures": stats.core_fallback_failures,
        },
        "missing_ids": sorted(missing_ids.keys()),
        "missing_details": [{"pdb_id": k, "reason": v} for k, v in sorted(missing_ids.items())],
        "api_error_ids": [{"pdb_id": k, "error": v} for k, v in sorted(api_error_ids.items())],
    }

    if api_failure_rate > args.fail_threshold:
        audit_report["manual_download_required"] = True
        audit_report["manual_download_plan"] = manual_download_plan()

        build_report = {
            "name": "pdb_structures_v1_build",
            "created_at": created_at,
            "status": "aborted",
            "reason": "api_failure_rate_exceeded_threshold",
            "api_failure_rate": api_failure_rate,
            "fail_threshold": args.fail_threshold,
            "manual_download_plan": audit_report["manual_download_plan"],
            "source_metrics": source_metrics,
        }
        write_json(args.report_audit, audit_report)
        write_json(args.report_build, build_report)
        print(
            "[ABORT] API failure rate exceeded threshold "
            f"({api_failure_rate:.4f} > {args.fail_threshold:.4f}). "
            "See report for manual download checklist."
        )
        return 3

    header = [
        "pdb_id",
        "uniprot_id",
        "experimental_method",
        "resolution",
        "release_date",
        "ligand_count",
        "source",
        "fetch_date",
    ]

    def iter_rows() -> Iterable[List[str]]:
        for pair in selected_pairs:
            meta = details.get(pair.pdb_id)
            if meta is None:
                yield [pair.pdb_id, pair.uniprot_id, "", "", "", "", SOURCE, args.fetch_date]
            else:
                yield [
                    pair.pdb_id,
                    pair.uniprot_id,
                    meta.experimental_method,
                    meta.resolution,
                    meta.release_date,
                    meta.ligand_count,
                    SOURCE,
                    args.fetch_date,
                ]

    out_rows = write_tsv(args.output, header=header, rows=iter_rows())

    build_report = {
        "name": "pdb_structures_v1_build",
        "created_at": created_at,
        "status": "ok",
        "input": str(args.input),
        "output": str(args.output),
        "sample_mode": sample_mode,
        "max_unique_pdb": args.max_unique_pdb,
        "source_metrics": source_metrics,
        "selected_pairs": len(selected_pairs),
        "selected_unique_pdb": len(selected_unique_ids),
        "output_rows": out_rows,
        "found_unique_pdb": len(details),
        "missing_unique_pdb": len(missing_ids),
        "api_error_unique_pdb": len(api_error_ids),
        "api_failure_rate": api_failure_rate,
        "fail_threshold": args.fail_threshold,
        "stats": {
            "batches_total": stats.batches_total,
            "batches_failed": stats.batches_failed,
            "graphql_requests": stats.graphql_requests,
            "graphql_failures": stats.graphql_failures,
            "core_fallback_requests": stats.core_fallback_requests,
            "core_fallback_failures": stats.core_fallback_failures,
        },
    }

    audit_report["manual_download_required"] = False

    write_json(args.report_build, build_report)
    write_json(args.report_audit, audit_report)

    print(f"[OK] output -> {args.output} (rows={out_rows})")
    print(
        "[OK] unique_pdb "
        f"found={len(details)} missing={len(missing_ids)} api_error={len(api_error_ids)} "
        f"failure_rate={api_failure_rate:.4f}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
