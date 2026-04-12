#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"
RCSB_CORE_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
USER_AGENT = "kg-rna-pdb/1.0"
SOURCE = "RCSB_DATA_API"
PDB_ID_RE = re.compile(r"^[0-9A-Z]{4}$")

GRAPHQL_QUERY = """
query($ids:[String!]!) {
  entries(entry_ids:$ids) {
    rcsb_id
    exptl { method }
    rcsb_entry_info {
      resolution_combined
    }
    rcsb_accession_info {
      initial_release_date
    }
  }
}
""".strip()


@dataclass
class MasterRow:
    rna_id: str
    rna_type: str
    taxon_id: str
    enst_base: str
    direct_urs: str


@dataclass
class PairRow:
    rna_id: str
    rna_type: str
    urs_id: str
    pdb_id: str
    pdb_entity_id: str
    mapping_strategy: str


@dataclass
class EntryMeta:
    experimental_method: str
    resolution: str
    release_date: str


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


def open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def resolve_tsv_or_gz(path: Path) -> Optional[Path]:
    if path.exists():
        return path
    if path.suffix == ".gz":
        alt = Path(str(path)[:-3])
        if alt.exists():
            return alt
    else:
        alt = Path(str(path) + ".gz")
        if alt.exists():
            return alt
    return None


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


def parse_pdb_external_id(raw: str) -> Optional[Tuple[str, str]]:
    token = (raw or "").strip()
    if not token:
        return None

    head, sep, tail = token.partition("_")
    pdb_id = head.strip().upper()
    if not PDB_ID_RE.match(pdb_id):
        return None
    pdb_entity_id = tail.strip().upper() if sep else ""
    return pdb_id, pdb_entity_id


def read_master_rows(path: Path, taxon_id: str, max_rows: Optional[int]) -> Tuple[List[MasterRow], Dict[str, int], Set[str], Set[str]]:
    rows: List[MasterRow] = []
    enst_bases: Set[str] = set()
    master_ids: Set[str] = set()

    metrics = {
        "master_rows_scanned": 0,
        "master_rows_target_taxon": 0,
        "master_rows_with_urs_id": 0,
        "master_rows_with_enst_id": 0,
    }

    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"rna_id", "rna_type", "taxon_id"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] rna_master schema missing columns: {sorted(missing)}")

        for row in reader:
            metrics["master_rows_scanned"] += 1
            tx = (row.get("taxon_id") or "").strip()
            if tx != taxon_id:
                continue
            metrics["master_rows_target_taxon"] += 1

            rna_id = (row.get("rna_id") or "").strip()
            if not rna_id:
                continue
            master_ids.add(rna_id)

            rna_type = (row.get("rna_type") or "").strip().lower()
            direct_urs = ""
            enst_base = ""

            if rna_id.startswith("URS") and rna_id.endswith(f"_{taxon_id}"):
                direct_urs = rna_id
                metrics["master_rows_with_urs_id"] += 1
            elif rna_id.startswith("ENST") and rna_id.endswith(f"_{taxon_id}"):
                enst_base = rna_id[: -(len(taxon_id) + 1)]
                enst_bases.add(enst_base)
                metrics["master_rows_with_enst_id"] += 1

            rows.append(
                MasterRow(
                    rna_id=rna_id,
                    rna_type=rna_type,
                    taxon_id=tx,
                    enst_base=enst_base,
                    direct_urs=direct_urs,
                )
            )

            if max_rows is not None and len(rows) >= max_rows:
                break

    return rows, metrics, enst_bases, master_ids


def load_xref_map(path: Optional[Path], master_ids: Set[str], taxon_id: str) -> Tuple[Dict[str, str], Dict[str, Any]]:
    mapping: Dict[str, str] = {}
    metrics: Dict[str, Any] = {
        "xref_used": False,
        "xref_path": str(path) if path else None,
        "xref_rows_scanned": 0,
        "xref_rows_target_taxon": 0,
        "xref_rows_for_master": 0,
        "xref_rows_with_urs": 0,
        "xref_conflicts": 0,
    }

    if path is None:
        return mapping, metrics
    if not path.exists():
        raise SystemExit(f"[ERROR] xref file not found: {path}")

    metrics["xref_used"] = True
    vote: Dict[str, Counter] = defaultdict(Counter)

    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"rna_id", "xref_id"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] xref schema missing columns: {sorted(missing)}")

        for row in reader:
            metrics["xref_rows_scanned"] += 1
            rna_id = (row.get("rna_id") or "").strip()
            if not rna_id or rna_id not in master_ids:
                continue
            metrics["xref_rows_for_master"] += 1

            tx = (row.get("taxon_id") or "").strip()
            if tx and tx != taxon_id:
                continue
            metrics["xref_rows_target_taxon"] += 1

            urs = (row.get("xref_id") or "").strip()
            if not (urs.startswith("URS") and urs.endswith(f"_{taxon_id}")):
                continue
            metrics["xref_rows_with_urs"] += 1
            vote[rna_id][urs] += 1

    for rna_id, counter in vote.items():
        ranked = sorted(counter.items(), key=lambda kv: (-kv[1], kv[0]))
        mapping[rna_id] = ranked[0][0]
        if len(ranked) > 1:
            metrics["xref_conflicts"] += 1

    metrics["xref_rows_resolved"] = len(mapping)
    return mapping, metrics


def scan_id_mapping(
    path: Path,
    taxon_id: str,
    enst_bases: Set[str],
    max_lines: Optional[int],
) -> Tuple[Dict[str, Counter], Dict[str, Set[Tuple[str, str]]], Dict[str, Any]]:
    allowed_dbs = {"ENSEMBL", "ENSEMBL_GENCODE", "GENCODE"}

    enst_to_urs: Dict[str, Counter] = defaultdict(Counter)
    urs_to_pdb: Dict[str, Set[Tuple[str, str]]] = defaultdict(set)

    metrics: Dict[str, Any] = {
        "idmap_lines_scanned": 0,
        "idmap_lines_target_taxon": 0,
        "idmap_ensembl_rows": 0,
        "idmap_ensembl_rows_hit_master_enst": 0,
        "idmap_pdb_rows": 0,
        "idmap_pdb_rows_malformed": 0,
        "idmap_pdb_rows_valid": 0,
    }
    malformed_examples: List[Dict[str, str]] = []

    with open_maybe_gzip(path) as f:
        for line in f:
            metrics["idmap_lines_scanned"] += 1
            if max_lines is not None and metrics["idmap_lines_scanned"] > max_lines:
                break
            if line.strip() == "":
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue

            urs_raw = parts[0].strip()
            db = parts[1].strip().upper()
            ext_id = parts[2].strip()
            tx = parts[3].strip()

            if tx != taxon_id:
                continue
            metrics["idmap_lines_target_taxon"] += 1

            urs = urs_raw if urs_raw.endswith(f"_{taxon_id}") else f"{urs_raw}_{taxon_id}"

            if db in allowed_dbs:
                metrics["idmap_ensembl_rows"] += 1
                if ext_id.startswith("ENST"):
                    enst_base = ext_id.split(".")[0]
                    if enst_base in enst_bases:
                        enst_to_urs[enst_base][urs] += 1
                        metrics["idmap_ensembl_rows_hit_master_enst"] += 1
                continue

            if db == "PDB":
                metrics["idmap_pdb_rows"] += 1
                parsed = parse_pdb_external_id(ext_id)
                if parsed is None:
                    metrics["idmap_pdb_rows_malformed"] += 1
                    if len(malformed_examples) < 20:
                        malformed_examples.append(
                            {
                                "urs_id": urs,
                                "external_id": ext_id,
                            }
                        )
                    continue
                metrics["idmap_pdb_rows_valid"] += 1
                urs_to_pdb[urs].add(parsed)

    metrics["idmap_enst_keys_resolved"] = len(enst_to_urs)
    metrics["idmap_urs_with_pdb"] = len(urs_to_pdb)
    metrics["idmap_pdb_malformed_examples"] = malformed_examples
    return enst_to_urs, urs_to_pdb, metrics


def resolve_canonical_urs(
    master_rows: List[MasterRow],
    xref_map: Dict[str, str],
    enst_to_urs: Dict[str, Counter],
) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, Any]]:
    rna_to_urs: Dict[str, str] = {}
    strategy_by_rna: Dict[str, str] = {}

    strategy_counts: Counter = Counter()
    unresolved = 0
    unresolved_examples: List[str] = []
    ambiguous_enst = 0

    for row in master_rows:
        rid = row.rna_id

        if rid in xref_map:
            rna_to_urs[rid] = xref_map[rid]
            strategy_by_rna[rid] = "xref_file"
            strategy_counts["xref_file"] += 1
            continue

        if row.direct_urs:
            rna_to_urs[rid] = row.direct_urs
            strategy_by_rna[rid] = "rna_id_direct_urs"
            strategy_counts["rna_id_direct_urs"] += 1
            continue

        if row.enst_base and row.enst_base in enst_to_urs:
            ranked = sorted(enst_to_urs[row.enst_base].items(), key=lambda kv: (-kv[1], kv[0]))
            rna_to_urs[rid] = ranked[0][0]
            strategy_by_rna[rid] = "id_mapping_ensembl_vote"
            strategy_counts["id_mapping_ensembl_vote"] += 1
            if len(ranked) > 1:
                ambiguous_enst += 1
            continue

        unresolved += 1
        if len(unresolved_examples) < 20:
            unresolved_examples.append(rid)

    metrics: Dict[str, Any] = {
        "resolved_rna_rows": len(rna_to_urs),
        "unresolved_rna_rows": unresolved,
        "strategy_counts": dict(strategy_counts),
        "ambiguous_enst_rows": ambiguous_enst,
        "unresolved_examples": unresolved_examples,
    }
    return rna_to_urs, strategy_by_rna, metrics


def build_pair_rows(
    master_rows: List[MasterRow],
    rna_to_urs: Dict[str, str],
    strategy_by_rna: Dict[str, str],
    urs_to_pdb: Dict[str, Set[Tuple[str, str]]],
) -> Tuple[List[PairRow], List[str], Dict[str, Any]]:
    pairs: List[PairRow] = []
    pair_seen: Set[Tuple[str, str, str, str]] = set()
    rna_rows_with_pdb: Set[str] = set()
    unique_pdb_order: List[str] = []
    unique_pdb_seen: Set[str] = set()

    metrics: Dict[str, Any] = {
        "pairs_raw": 0,
        "pairs_deduped": 0,
        "pairs_dropped_duplicate": 0,
        "rna_rows_with_pdb": 0,
        "rna_rows_without_pdb": 0,
    }

    for row in master_rows:
        rid = row.rna_id
        urs = rna_to_urs.get(rid)
        if not urs:
            continue

        pdb_items = urs_to_pdb.get(urs) or set()
        if not pdb_items:
            metrics["rna_rows_without_pdb"] += 1
            continue

        rna_rows_with_pdb.add(rid)
        for pdb_id, pdb_entity_id in sorted(pdb_items, key=lambda x: (x[0], x[1])):
            metrics["pairs_raw"] += 1
            key = (rid, urs, pdb_id, pdb_entity_id)
            if key in pair_seen:
                metrics["pairs_dropped_duplicate"] += 1
                continue
            pair_seen.add(key)
            pairs.append(
                PairRow(
                    rna_id=rid,
                    rna_type=row.rna_type,
                    urs_id=urs,
                    pdb_id=pdb_id,
                    pdb_entity_id=pdb_entity_id,
                    mapping_strategy=strategy_by_rna.get(rid, "unknown"),
                )
            )
            metrics["pairs_deduped"] += 1
            if pdb_id not in unique_pdb_seen:
                unique_pdb_seen.add(pdb_id)
                unique_pdb_order.append(pdb_id)

    metrics["rna_rows_with_pdb"] = len(rna_rows_with_pdb)
    metrics["unique_pdb_total"] = len(unique_pdb_seen)
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

    return EntryMeta(
        experimental_method=experimental_method,
        resolution=resolution,
        release_date=release_date,
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
        except Exception as e:  # noqa: BLE001
            last_error = f"{type(e).__name__}: {e}"
            if attempt < retries:
                time.sleep(retry_sleep * attempt)

    stats.graphql_failures += 1
    return None, (last_error or "unknown_graphql_error")


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


def fetch_with_partition(
    pdb_ids: List[str],
    timeout: int,
    retries: int,
    retry_sleep: float,
    ctx: ssl.SSLContext,
    stats: FetchStats,
) -> Tuple[Dict[str, EntryMeta], Dict[str, str], Dict[str, str]]:
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
        "reason": "RCSB API batch failure rate exceeded threshold; switch to manual bulk acquisition.",
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
                "name": "RCSB release notes/checksums",
                "url": "https://files.rcsb.org/pub/pdb/derived_data/index/",
                "expected_size_hint": "index directory",
            },
        ],
        "place_under": "data/raw/rna/pdb_bulk/",
        "checksum": {
            "method": "sha256",
            "commands": [
                "mkdir -p data/raw/rna/pdb_bulk",
                "cd data/raw/rna/pdb_bulk",
                "shasum -a 256 entries.idx resolu.idx > SHA256SUMS.txt",
                "shasum -a 256 -c SHA256SUMS.txt",
            ],
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna-master", type=Path, required=True)
    ap.add_argument("--id-mapping", type=Path, required=True)
    ap.add_argument("--xref", type=Path, default=None)
    ap.add_argument("--taxon-id", default="9606")
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--report-build", type=Path, required=True)
    ap.add_argument("--report-audit", type=Path, required=True)
    ap.add_argument("--source-version", default="RNAcentral:25;RCSB:live")
    ap.add_argument("--fetch-date", default=utc_today())

    ap.add_argument("--max-master-rows", type=int, default=None)
    ap.add_argument("--max-idmap-lines", type=int, default=None)
    ap.add_argument("--max-unique-pdb", type=int, default=None)

    ap.add_argument("--batch-size", type=int, default=500)
    ap.add_argument("--timeout", type=int, default=40)
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--retry-sleep", type=float, default=1.0)
    ap.add_argument("--sleep-between-batches", type=float, default=0.05)
    ap.add_argument("--fail-threshold", type=float, default=0.05)
    args = ap.parse_args()

    created_at = utc_now_iso()

    if not args.rna_master.exists():
        raise SystemExit(f"[ERROR] missing rna_master: {args.rna_master}")

    id_mapping_path = resolve_tsv_or_gz(args.id_mapping)
    if id_mapping_path is None:
        raise SystemExit(f"[ERROR] missing id_mapping (.tsv/.gz): {args.id_mapping}")

    if args.batch_size <= 0:
        raise SystemExit("[ERROR] --batch-size must be > 0")

    master_rows, master_metrics, enst_bases, master_ids = read_master_rows(
        path=args.rna_master,
        taxon_id=args.taxon_id,
        max_rows=args.max_master_rows,
    )
    if not master_rows:
        raise SystemExit("[ERROR] no master rows after taxon filter")

    xref_map, xref_metrics = load_xref_map(
        path=args.xref,
        master_ids=master_ids,
        taxon_id=args.taxon_id,
    )

    enst_to_urs, urs_to_pdb, idmap_metrics = scan_id_mapping(
        path=id_mapping_path,
        taxon_id=args.taxon_id,
        enst_bases=enst_bases,
        max_lines=args.max_idmap_lines,
    )

    rna_to_urs, strategy_by_rna, resolve_metrics = resolve_canonical_urs(
        master_rows=master_rows,
        xref_map=xref_map,
        enst_to_urs=enst_to_urs,
    )

    pairs, unique_pdb_ids, pair_metrics = build_pair_rows(
        master_rows=master_rows,
        rna_to_urs=rna_to_urs,
        strategy_by_rna=strategy_by_rna,
        urs_to_pdb=urs_to_pdb,
    )

    selected_unique_ids = unique_pdb_ids
    sample_mode = False
    if args.max_unique_pdb is not None and args.max_unique_pdb >= 0 and args.max_unique_pdb < len(unique_pdb_ids):
        sample_mode = True
        selected_unique_ids = unique_pdb_ids[: args.max_unique_pdb]

    selected_set = set(selected_unique_ids)
    selected_pairs = [p for p in pairs if p.pdb_id in selected_set]
    selected_rna_with_pdb = len({p.rna_id for p in selected_pairs})

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

        if stats.graphql_failures > graphql_failures_before:
            stats.batches_failed += 1

        if idx % 10 == 0 or idx == len(batches):
            print(
                f"[INFO] fetched {idx}/{len(batches)} batches | "
                f"found={len(details)} missing={len(missing_ids)} api_error={len(api_error_ids)}"
            )

        if args.sleep_between_batches > 0 and idx < len(batches):
            time.sleep(args.sleep_between_batches)

    api_failure_rate = (len(api_error_ids) / len(selected_unique_ids)) if selected_unique_ids else 0.0

    master_total = len(master_rows)
    canonical_coverage_rate = (len(rna_to_urs) / master_total) if master_total else 0.0
    pdb_coverage_rate = (selected_rna_with_pdb / master_total) if master_total else 0.0

    audit_report: Dict[str, Any] = {
        "name": "rna_pdb_structures_v1_api_audit",
        "created_at": created_at,
        "inputs": {
            "rna_master": str(args.rna_master),
            "id_mapping": str(id_mapping_path),
            "xref": str(args.xref) if args.xref else None,
            "taxon_id": args.taxon_id,
        },
        "sample_mode": sample_mode,
        "max_unique_pdb": args.max_unique_pdb,
        "selected_unique_pdb": len(selected_unique_ids),
        "selected_pairs": len(selected_pairs),
        "found_unique_pdb": len(details),
        "missing_unique_pdb": len(missing_ids),
        "api_error_unique_pdb": len(api_error_ids),
        "api_failure_rate": api_failure_rate,
        "fail_threshold": args.fail_threshold,
        "coverage": {
            "master_rows": master_total,
            "resolved_urs_rows": len(rna_to_urs),
            "resolved_urs_rate": canonical_coverage_rate,
            "rna_rows_with_pdb": selected_rna_with_pdb,
            "rna_rows_with_pdb_rate": pdb_coverage_rate,
        },
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
            "name": "rna_pdb_structures_v1.build",
            "created_at": created_at,
            "status": "aborted",
            "reason": "api_failure_rate_exceeded_threshold",
            "api_failure_rate": api_failure_rate,
            "fail_threshold": args.fail_threshold,
            "manual_download_plan": audit_report["manual_download_plan"],
            "metrics": {
                "master": master_metrics,
                "xref": xref_metrics,
                "id_mapping": idmap_metrics,
                "resolve": resolve_metrics,
                "pairs": pair_metrics,
                "coverage": {
                    "resolved_urs_rate": canonical_coverage_rate,
                    "rna_with_pdb_rate": pdb_coverage_rate,
                },
            },
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
    ]

    def iter_rows() -> Iterable[List[str]]:
        for pair in selected_pairs:
            meta = details.get(pair.pdb_id)
            if meta is None:
                yield [
                    pair.rna_id,
                    pair.rna_type,
                    pair.urs_id,
                    pair.pdb_id,
                    pair.pdb_entity_id,
                    pair.mapping_strategy,
                    "",
                    "",
                    "",
                    SOURCE,
                    args.source_version,
                    args.fetch_date,
                ]
            else:
                yield [
                    pair.rna_id,
                    pair.rna_type,
                    pair.urs_id,
                    pair.pdb_id,
                    pair.pdb_entity_id,
                    pair.mapping_strategy,
                    meta.experimental_method,
                    meta.resolution,
                    meta.release_date,
                    SOURCE,
                    args.source_version,
                    args.fetch_date,
                ]

    out_rows = write_tsv(args.output, header=header, rows=iter_rows())

    build_report = {
        "name": "rna_pdb_structures_v1_build",
        "created_at": created_at,
        "status": "ok",
        "input": {
            "rna_master": str(args.rna_master),
            "id_mapping": str(id_mapping_path),
            "xref": str(args.xref) if args.xref else None,
            "taxon_id": args.taxon_id,
        },
        "output": str(args.output),
        "sample_mode": sample_mode,
        "max_unique_pdb": args.max_unique_pdb,
        "max_master_rows": args.max_master_rows,
        "max_idmap_lines": args.max_idmap_lines,
        "metrics": {
            "master": master_metrics,
            "xref": xref_metrics,
            "id_mapping": idmap_metrics,
            "resolve": resolve_metrics,
            "pairs": pair_metrics,
            "coverage": {
                "master_rows": master_total,
                "resolved_urs_rows": len(rna_to_urs),
                "resolved_urs_rate": canonical_coverage_rate,
                "rna_rows_with_pdb": selected_rna_with_pdb,
                "rna_rows_with_pdb_rate": pdb_coverage_rate,
            },
        },
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
        "[OK] coverage "
        f"resolved_urs={len(rna_to_urs)}/{master_total} ({canonical_coverage_rate:.4%}), "
        f"rna_with_pdb={selected_rna_with_pdb}/{master_total} ({pdb_coverage_rate:.4%})"
    )
    print(
        "[OK] unique_pdb "
        f"found={len(details)} missing={len(missing_ids)} api_error={len(api_error_ids)} "
        f"failure_rate={api_failure_rate:.4f}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
