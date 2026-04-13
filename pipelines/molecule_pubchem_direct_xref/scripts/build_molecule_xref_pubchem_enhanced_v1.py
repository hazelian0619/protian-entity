#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sqlite3
import threading
import time
import warnings
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import requests
from requests.exceptions import RequestException, RequestsDependencyWarning

INCHIKEY_PREFIX = "IK::"
SMILES_PREFIX = "SMILES::"
PUBCHEM_API_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
CONF_RANK = {"": 0, "low": 1, "medium": 2, "high": 3}

warnings.filterwarnings("ignore", category=RequestsDependencyWarning)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def utc_today() -> str:
    return datetime.now(timezone.utc).date().isoformat()


def normalize(x: str) -> str:
    return (x or "").strip()


def split_multi(x: str) -> List[str]:
    text = normalize(x)
    if not text:
        return []
    return [token.strip() for token in text.split(";") if token.strip()]


def join_multi(values: Iterable[str]) -> str:
    return ";".join(sorted({normalize(v) for v in values if normalize(v)}))


def merge_confidence(current: str, incoming: str) -> str:
    c = normalize(current).lower()
    i = normalize(incoming).lower()
    return c if CONF_RANK.get(c, 0) >= CONF_RANK.get(i, 0) else i


def parse_pubchem_cids(payload: Optional[Dict[str, Any]]) -> List[str]:
    if not payload:
        return []
    raw = payload.get("IdentifierList", {}).get("CID", [])
    out: List[str] = []
    for cid in raw:
        token = normalize(str(cid))
        if token.isdigit():
            out.append(token)
    return sorted(set(out), key=lambda x: int(x))


def resolve_pubchem_candidates(
    inchikey_cids: Sequence[str],
    smiles_cids: Sequence[str],
    has_smiles: bool,
) -> Dict[str, Any]:
    ik = sorted(set([normalize(x) for x in inchikey_cids if normalize(x)]), key=lambda x: int(x))
    sm = sorted(set([normalize(x) for x in smiles_cids if normalize(x)]), key=lambda x: int(x))

    selected: List[str] = []
    strategy = ""
    confidence = ""
    tie_break = "none"

    if len(ik) == 1:
        selected = ik
        strategy = "pubchem_direct_inchikey_exact_unique"
        confidence = "high"
        tie_break = "not_needed"
    elif len(ik) > 1:
        if sm:
            inter = sorted(set(ik).intersection(sm), key=lambda x: int(x))
            if len(inter) == 1:
                selected = inter
                strategy = "pubchem_direct_inchikey_multi_tiebreak_smiles"
                confidence = "high"
                tie_break = "resolved_unique_by_smiles_intersection"
            elif len(inter) > 1:
                selected = inter
                strategy = "pubchem_direct_inchikey_multi_intersection"
                confidence = "medium"
                tie_break = "narrowed_multi_by_smiles_intersection"
            else:
                selected = ik
                strategy = "pubchem_direct_inchikey_multi_cid"
                confidence = "low"
                tie_break = "no_overlap_keep_inchikey_candidates"
        else:
            selected = ik
            strategy = "pubchem_direct_inchikey_multi_cid"
            confidence = "low"
            tie_break = "smiles_missing_keep_inchikey_candidates"
    else:
        if len(sm) == 1:
            selected = sm
            strategy = "pubchem_direct_smiles_fallback_unique"
            confidence = "medium"
            tie_break = "fallback_smiles_unique"
        elif len(sm) > 1:
            selected = sm
            strategy = "pubchem_direct_smiles_fallback_multi_cid"
            confidence = "low"
            tie_break = "fallback_smiles_multi_keep_all"

    return {
        "selected_cids": selected,
        "match_strategy": strategy,
        "confidence": confidence,
        "tie_break": tie_break,
        "inchikey_candidates": ik,
        "smiles_candidates": sm,
        "smiles_used": bool(has_smiles),
    }


def apply_pubchem_enhancement(
    row: Dict[str, str],
    resolved: Dict[str, Any],
    source_tag: str,
    source_version_tag: str,
    today: Optional[str] = None,
) -> Tuple[Dict[str, str], bool, str]:
    out = dict(row)
    selected = resolved.get("selected_cids", [])
    if not selected:
        return out, False, "no_new_cid"

    existing_cids = split_multi(out.get("pubchem_cid", ""))
    existing_conf = normalize(out.get("confidence", "")).lower()
    if existing_cids and existing_conf == "high":
        return out, False, "locked_existing_high_confidence"
    if existing_cids:
        return out, False, "existing_pubchem_present_skip"

    out["pubchem_cid"] = join_multi(selected)
    out["match_strategy"] = join_multi(split_multi(out.get("match_strategy", "")) + [resolved.get("match_strategy", "")])
    out["confidence"] = merge_confidence(out.get("confidence", ""), resolved.get("confidence", ""))
    out["xref_source"] = join_multi(split_multi(out.get("xref_source", "")) + [source_tag])
    out["source_version"] = join_multi(split_multi(out.get("source_version", "")) + [source_version_tag])
    if today is not None:
        out["fetch_date"] = today
    return out, True, "updated"


def _read_tsv(path: Path, max_rows: Optional[int] = None) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        rows: List[Dict[str, str]] = []
        for idx, row in enumerate(reader, start=1):
            rows.append(row)
            if max_rows is not None and idx >= max_rows:
                break
    return list(reader.fieldnames), rows


def _write_tsv(path: Path, header: Sequence[str], rows: Iterable[Dict[str, str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    written = 0
    with tmp.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=list(header), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in header})
            written += 1
    tmp.replace(path)
    return written


def _table_exists(conn: sqlite3.Connection, name: str) -> bool:
    q = "SELECT 1 FROM sqlite_master WHERE type IN ('table','view') AND name=?"
    return conn.execute(q, (name,)).fetchone() is not None


def load_structure_map(m1_db: Path, inchikeys: Sequence[str]) -> Dict[str, Dict[str, str]]:
    if not m1_db.exists() or not inchikeys:
        return {}

    candidates = [
        "molecule_entity_strict_chembl36",
        "molecule_entity_chembl36",
        "molecule_raw_chembl36",
    ]
    out: Dict[str, Dict[str, str]] = {}
    with sqlite3.connect(m1_db) as conn:
        tables = [t for t in candidates if _table_exists(conn, t)]
        if not tables:
            return {}

        uniq = sorted(set([normalize(x).upper() for x in inchikeys if normalize(x)]))
        if not uniq:
            return {}

        unresolved = set(uniq)
        batch_size = 500
        for table in tables:
            if not unresolved:
                break
            pending = sorted(unresolved)
            for i in range(0, len(pending), batch_size):
                batch = pending[i : i + batch_size]
                placeholders = ",".join(["?"] * len(batch))
                q = f"""
                    SELECT inchikey,
                           COALESCE(NULLIF(TRIM(smiles), ''), '') AS smiles,
                           COALESCE(NULLIF(TRIM(inchi), ''), '') AS inchi
                    FROM {table}
                    WHERE inchikey IN ({placeholders})
                """
                found_in_batch: List[str] = []
                for inchikey, smiles, inchi in conn.execute(q, batch):
                    ik = normalize(inchikey).upper()
                    if not ik:
                        continue
                    current = out.get(ik, {})
                    s = normalize(smiles) or current.get("smiles", "")
                    ic = normalize(inchi) or current.get("inchi", "")
                    out[ik] = {"smiles": s, "inchi": ic}
                    found_in_batch.append(ik)
                if found_in_batch:
                    unresolved.difference_update(found_in_batch)

        # Tiny fallback only for rare format-noise leftovers.
        if unresolved and len(unresolved) <= 200:
            conn.execute("DROP TABLE IF EXISTS temp.tmp_ik")
            conn.execute("CREATE TEMP TABLE tmp_ik (inchikey TEXT PRIMARY KEY)")
            conn.executemany("INSERT OR IGNORE INTO tmp_ik(inchikey) VALUES (?)", [(x,) for x in sorted(unresolved)])
            for table in tables:
                q = f"""
                    SELECT t.inchikey,
                           COALESCE(MAX(NULLIF(TRIM(m.smiles), '')), '') AS smiles,
                           COALESCE(MAX(NULLIF(TRIM(m.inchi), '')), '') AS inchi
                    FROM tmp_ik t
                    LEFT JOIN {table} m ON UPPER(TRIM(m.inchikey)) = t.inchikey
                    GROUP BY t.inchikey
                """
                for inchikey, smiles, inchi in conn.execute(q):
                    ik = normalize(inchikey).upper()
                    current = out.get(ik, {})
                    s = normalize(smiles) or current.get("smiles", "")
                    ic = normalize(inchi) or current.get("inchi", "")
                    out[ik] = {"smiles": s, "inchi": ic}
    return out


class JsonlCache:
    def __init__(self, path: Path):
        self.path = path
        self.lock = threading.Lock()
        self.store: Dict[str, List[str]] = {}
        self.loaded = 0
        self._load()

    def _load(self) -> None:
        if not self.path.exists():
            return
        with self.path.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    obj = json.loads(line)
                    key = str(obj.get("key", ""))
                    cids = [normalize(str(x)) for x in obj.get("cids", []) if normalize(str(x))]
                    if key:
                        self.store[key] = sorted(set(cids), key=lambda x: int(x))
                        self.loaded += 1
                except Exception:
                    continue

    def get(self, key: str) -> Optional[List[str]]:
        with self.lock:
            value = self.store.get(key)
            return list(value) if value is not None else None

    def set(self, key: str, cids: Sequence[str], meta: Optional[Dict[str, Any]] = None) -> None:
        norm = sorted(set([normalize(x) for x in cids if normalize(x)]), key=lambda x: int(x))
        payload = {
            "key": key,
            "cids": norm,
            "updated_at": utc_now(),
        }
        if meta:
            payload["meta"] = meta
        with self.lock:
            self.store[key] = norm
            self.path.parent.mkdir(parents=True, exist_ok=True)
            with self.path.open("a", encoding="utf-8") as f:
                f.write(json.dumps(payload, ensure_ascii=False) + "\n")


class PubChemClient:
    def __init__(
        self,
        cache: JsonlCache,
        timeout: int,
        retries: int,
        qps: float,
        user_agent: str = "protian-entity-molecule-pubchem-direct-xref/1.0",
    ):
        self.cache = cache
        self.timeout = timeout
        self.retries = retries
        self.min_interval = (1.0 / qps) if qps > 0 else 0.0
        self.user_agent = user_agent

        self._throttle_lock = threading.Lock()
        self._local = threading.local()
        self._next_allowed = 0.0
        self.stats = Counter()

    def _session(self) -> requests.Session:
        sess = getattr(self._local, "session", None)
        if sess is None:
            sess = requests.Session()
            sess.headers.update({"Accept": "application/json", "User-Agent": self.user_agent})
            self._local.session = sess
        return sess

    def _sleep_for_rate_limit(self) -> None:
        if self.min_interval <= 0:
            return
        with self._throttle_lock:
            now = time.monotonic()
            if now < self._next_allowed:
                time.sleep(self._next_allowed - now)
                now = time.monotonic()
            self._next_allowed = now + self.min_interval

    def _fetch_json(self, url: str, data: Optional[Dict[str, str]] = None) -> Tuple[Optional[Dict[str, Any]], str]:
        method = "POST" if data is not None else "GET"
        session = self._session()

        for attempt in range(self.retries + 1):
            self._sleep_for_rate_limit()
            try:
                if method == "GET":
                    resp = session.get(url, timeout=self.timeout)
                else:
                    resp = session.post(url, data=data, timeout=self.timeout)
                code = int(resp.status_code)
                if code == 200:
                    self.stats["http_success"] += 1
                    return resp.json(), "ok"
                if code in {400, 404}:
                    self.stats[f"http_{code}_empty"] += 1
                    return None, "not_found"
                self.stats[f"http_{code}_retryable"] += 1
                if attempt >= self.retries:
                    self.stats["http_failed_final"] += 1
                    return None, "failed"
            except RequestException:
                self.stats["http_error_retryable"] += 1
                if attempt >= self.retries:
                    self.stats["http_failed_final"] += 1
                    return None, "failed"
            time.sleep(min(4.0, 0.4 * (2 ** attempt)))
        return None, "failed"

    def _cache_key_smiles(self, smiles: str) -> str:
        digest = hashlib.sha1(smiles.encode("utf-8")).hexdigest()
        return f"{SMILES_PREFIX}{digest}"

    def query_inchikey(self, inchikey: str) -> List[str]:
        key = f"{INCHIKEY_PREFIX}{inchikey}"
        cached = self.cache.get(key)
        if cached is not None:
            self.stats["cache_hit_inchikey"] += 1
            return cached

        self.stats["cache_miss_inchikey"] += 1
        self.stats["request_inchikey"] += 1
        url = f"{PUBCHEM_API_BASE}/inchikey/{inchikey}/cids/JSON"
        payload, outcome = self._fetch_json(url=url)
        cids = parse_pubchem_cids(payload)
        if outcome in {"ok", "not_found"}:
            self.cache.set(key=key, cids=cids, meta={"type": "inchikey", "inchikey": inchikey, "outcome": outcome})
        else:
            self.stats["cache_skip_failed_inchikey"] += 1
        return cids

    def query_smiles(self, smiles: str) -> List[str]:
        key = self._cache_key_smiles(smiles)
        cached = self.cache.get(key)
        if cached is not None:
            self.stats["cache_hit_smiles"] += 1
            return cached

        self.stats["cache_miss_smiles"] += 1
        self.stats["request_smiles"] += 1
        url = f"{PUBCHEM_API_BASE}/smiles/cids/JSON"
        payload, outcome = self._fetch_json(url=url, data={"smiles": smiles})
        cids = parse_pubchem_cids(payload)
        if outcome in {"ok", "not_found"}:
            self.cache.set(key=key, cids=cids, meta={"type": "smiles", "outcome": outcome})
        else:
            self.stats["cache_skip_failed_smiles"] += 1
        return cids


def coverage_stats(rows: Sequence[Dict[str, str]]) -> Dict[str, Any]:
    total = len(rows)
    covered = sum(1 for r in rows if split_multi(r.get("pubchem_cid", "")))
    return {
        "rows": total,
        "rows_with_pubchem_cid": covered,
        "coverage_rate": (covered / total) if total else 0.0,
    }


def _process_single_lookup(
    inchikey: str,
    smiles: str,
    client: PubChemClient,
) -> Dict[str, Any]:
    ik_cids = client.query_inchikey(inchikey) if normalize(inchikey) else []
    smiles_cids: List[str] = []

    if smiles and (len(ik_cids) == 0 or len(ik_cids) > 1):
        smiles_cids = client.query_smiles(smiles)

    resolved = resolve_pubchem_candidates(
        inchikey_cids=ik_cids,
        smiles_cids=smiles_cids,
        has_smiles=bool(smiles),
    )

    return {
        "inchikey": inchikey,
        "smiles": smiles,
        "inchikey_cids": ik_cids,
        "smiles_cids": smiles_cids,
        "resolved": resolved,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--core-tsv", type=Path, required=True)
    ap.add_argument("--m1-db", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--build-report", type=Path, required=True)
    ap.add_argument("--conflict-audit", type=Path, required=True)
    ap.add_argument("--max-rows", type=int, default=None)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--qps", type=float, default=5.0)
    ap.add_argument("--timeout", type=int, default=20)
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument(
        "--cache-path",
        type=Path,
        default=Path("pipelines/molecule_pubchem_direct_xref/.cache/pubchem_lookup_cache_v1.jsonl"),
    )
    ap.add_argument("--progress-every", type=int, default=200)
    args = ap.parse_args()

    if not args.core_tsv.exists():
        raise SystemExit(f"[ERROR] missing input: {args.core_tsv}")
    if not args.m1_db.exists():
        raise SystemExit(f"[ERROR] missing input: {args.m1_db}")

    created_at = utc_now()
    today = utc_today()
    source_tag = "pubchem_pug_rest"
    source_version_tag = f"PubChem:PUGREST:{today}"

    header, rows = _read_tsv(args.core_tsv, max_rows=args.max_rows)
    required = [
        "inchikey",
        "chembl_id",
        "drugbank_id",
        "pubchem_cid",
        "match_strategy",
        "confidence",
        "xref_source",
        "fetch_date",
        "source_version",
    ]
    for col in required:
        if col not in header:
            header.append(col)

    before = coverage_stats(rows)
    rows_missing_pubchem = [r for r in rows if not split_multi(r.get("pubchem_cid", ""))]
    target_inchikeys = sorted({normalize(r.get("inchikey", "")).upper() for r in rows_missing_pubchem if normalize(r.get("inchikey", ""))})
    structure_map = load_structure_map(args.m1_db, target_inchikeys)

    cache = JsonlCache(args.cache_path)
    client = PubChemClient(
        cache=cache,
        timeout=args.timeout,
        retries=args.retries,
        qps=args.qps,
    )

    lookups: Dict[str, Dict[str, Any]] = {}
    print(f"[INFO] target missing-pubchem rows: {len(rows_missing_pubchem)}")
    print(f"[INFO] structure map loaded: {len(structure_map)} inchikey records")

    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        future_map = {}
        for row in rows_missing_pubchem:
            ik = normalize(row.get("inchikey", "")).upper()
            smiles = normalize(structure_map.get(ik, {}).get("smiles", ""))
            future = pool.submit(_process_single_lookup, ik, smiles, client)
            future_map[future] = ik

        for idx, future in enumerate(as_completed(future_map), start=1):
            result = future.result()
            lookups[result["inchikey"]] = result
            if args.progress_every > 0 and idx % args.progress_every == 0:
                print(f"[INFO] lookup progress {idx}/{len(future_map)}")

    out_rows: List[Dict[str, str]] = []
    update_reasons = Counter()
    strategy_counts = Counter()
    conflict_records: List[Dict[str, Any]] = []
    cids_added_total = 0

    for row in rows:
        ik = normalize(row.get("inchikey", "")).upper()
        lookup = lookups.get(ik)
        if not lookup:
            out_rows.append(dict(row))
            continue

        resolved = lookup["resolved"]
        updated_row, changed, reason = apply_pubchem_enhancement(
            row=row,
            resolved=resolved,
            source_tag=source_tag,
            source_version_tag=source_version_tag,
            today=today,
        )
        update_reasons[reason] += 1

        if changed:
            strategy_counts[resolved.get("match_strategy", "")] += 1
            cids_added_total += len(split_multi(updated_row.get("pubchem_cid", "")))
        if len(resolved.get("inchikey_candidates", [])) > 1 or len(resolved.get("smiles_candidates", [])) > 1:
            conflict_records.append(
                {
                    "inchikey": ik,
                    "match_strategy": resolved.get("match_strategy", ""),
                    "confidence": resolved.get("confidence", ""),
                    "tie_break": resolved.get("tie_break", ""),
                    "selected_cids": resolved.get("selected_cids", []),
                    "inchikey_candidates": resolved.get("inchikey_candidates", []),
                    "smiles_candidates": resolved.get("smiles_candidates", []),
                    "smiles_used": resolved.get("smiles_used", False),
                    "smiles_present_in_m1": bool(lookup.get("smiles")),
                }
            )

        out_rows.append(updated_row)

    rows_written = _write_tsv(args.out, header=header, rows=out_rows)
    after = coverage_stats(out_rows)

    conflict_counter = Counter()
    for rec in conflict_records:
        tie_break = normalize(rec.get("tie_break", "none")) or "none"
        conflict_counter[tie_break] += 1

    build_report = {
        "name": "molecule_xref_pubchem_enhanced_v1.build",
        "created_at": created_at,
        "inputs": {
            "core_tsv": str(args.core_tsv),
            "m1_db": str(args.m1_db),
            "cache_path": str(args.cache_path),
            "cache_preloaded_entries": cache.loaded,
        },
        "output": str(args.out),
        "metrics": {
            "rows_written": rows_written,
            "before": before,
            "after": after,
            "coverage_delta_abs": after["coverage_rate"] - before["coverage_rate"],
            "coverage_delta_rows": after["rows_with_pubchem_cid"] - before["rows_with_pubchem_cid"],
            "target_rows_missing_pubchem_before": len(rows_missing_pubchem),
            "rows_updated": int(update_reasons.get("updated", 0)),
            "cids_added_total": cids_added_total,
            "strategy_counts": dict(strategy_counts),
            "update_reason_counts": dict(update_reasons),
            "lookup_stats": dict(client.stats),
            "conflict_rows": len(conflict_records),
        },
        "source_version_tag": source_version_tag,
        "sample_mode": args.max_rows is not None,
        "max_rows": args.max_rows,
    }

    conflict_audit = {
        "name": "molecule_xref_pubchem_enhanced_v1.conflict_audit",
        "created_at": created_at,
        "total_conflicts": len(conflict_records),
        "tie_break_counts": dict(conflict_counter),
        "records": conflict_records,
    }

    args.build_report.parent.mkdir(parents=True, exist_ok=True)
    args.build_report.write_text(json.dumps(build_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    args.conflict_audit.parent.mkdir(parents=True, exist_ok=True)
    args.conflict_audit.write_text(json.dumps(conflict_audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(
        "[OK] enhanced rows={} coverage {:.4f} -> {:.4f} (delta rows={})".format(
            rows_written,
            before["coverage_rate"],
            after["coverage_rate"],
            after["rows_with_pubchem_cid"] - before["rows_with_pubchem_cid"],
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
