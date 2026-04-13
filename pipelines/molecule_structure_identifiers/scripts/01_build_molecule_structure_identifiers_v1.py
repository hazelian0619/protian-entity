#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
import threading
import time
import warnings
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import requests
from requests.exceptions import RequestException, RequestsDependencyWarning

warnings.filterwarnings("ignore", category=RequestsDependencyWarning)

INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def utc_today() -> str:
    return datetime.now(timezone.utc).date().isoformat()


def normalize(x: str) -> str:
    return (x or "").strip()


def split_multi(x: str) -> List[str]:
    t = normalize(x)
    if not t:
        return []
    return [p.strip() for p in t.split(";") if p.strip()]


def join_multi(values: Iterable[str]) -> str:
    return ";".join(sorted({normalize(v) for v in values if normalize(v)}))


def valid_inchikey(x: str) -> bool:
    return bool(INCHIKEY_RE.match(normalize(x).upper()))


def read_tsv(path: Path, max_rows: Optional[int] = None) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        rows: List[Dict[str, str]] = []
        for i, row in enumerate(reader, start=1):
            rows.append(row)
            if max_rows is not None and i >= max_rows:
                break
    return list(reader.fieldnames), rows


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Dict[str, str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    written = 0
    with tmp.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=list(header), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in header})
            written += 1
    tmp.replace(path)
    return written


def load_chembl_structures(chembl_db: Path, inchikeys: Sequence[str]) -> Dict[str, Dict[str, str]]:
    if not chembl_db.exists() or not inchikeys:
        return {}

    unique_iks = sorted({normalize(x).upper() for x in inchikeys if valid_inchikey(x)})
    if not unique_iks:
        return {}

    out: Dict[str, Dict[str, str]] = {}
    batch_size = 1000

    with sqlite3.connect(chembl_db) as conn:
        for i in range(0, len(unique_iks), batch_size):
            batch = unique_iks[i : i + batch_size]
            placeholders = ",".join(["?"] * len(batch))
            query = f"""
                SELECT
                    standard_inchi_key,
                    COALESCE(NULLIF(TRIM(standard_inchi), ''), '') AS standard_inchi,
                    COALESCE(NULLIF(TRIM(canonical_smiles), ''), '') AS canonical_smiles
                FROM compound_structures
                WHERE standard_inchi_key IN ({placeholders})
            """
            for ik, inchi, smiles in conn.execute(query, batch):
                key = normalize(ik).upper()
                if not valid_inchikey(key):
                    continue

                candidate = {
                    "inchi": normalize(inchi),
                    "canonical_smiles": normalize(smiles),
                }
                prev = out.get(key)
                if prev is None:
                    out[key] = candidate
                    continue

                prev_score = int(bool(prev.get("inchi"))) + int(bool(prev.get("canonical_smiles")))
                cand_score = int(bool(candidate.get("inchi"))) + int(bool(candidate.get("canonical_smiles")))
                if cand_score > prev_score:
                    out[key] = candidate

    return out


class JsonlCache:
    def __init__(self, path: Path):
        self.path = path
        self.lock = threading.Lock()
        self.store: Dict[str, Dict[str, str]] = {}
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
                except Exception:
                    continue
                key = normalize(str(obj.get("key", "")))
                if not key:
                    continue
                self.store[key] = {
                    "inchi": normalize(str(obj.get("inchi", ""))),
                    "canonical_smiles": normalize(str(obj.get("canonical_smiles", ""))),
                    "status": normalize(str(obj.get("status", ""))),
                }
                self.loaded += 1

    def get(self, key: str) -> Optional[Dict[str, str]]:
        with self.lock:
            v = self.store.get(key)
            return dict(v) if v is not None else None

    def set(self, key: str, value: Dict[str, str]) -> None:
        payload = {
            "key": key,
            "inchi": normalize(value.get("inchi", "")),
            "canonical_smiles": normalize(value.get("canonical_smiles", "")),
            "status": normalize(value.get("status", "")),
            "updated_at": utc_now(),
        }
        with self.lock:
            self.store[key] = {
                "inchi": payload["inchi"],
                "canonical_smiles": payload["canonical_smiles"],
                "status": payload["status"],
            }
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
        user_agent: str = "kg-molecule-structure-identifiers/1.0",
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
        session = getattr(self._local, "session", None)
        if session is None:
            session = requests.Session()
            session.headers.update({"Accept": "application/json", "User-Agent": self.user_agent})
            self._local.session = session
        return session

    def _throttle(self) -> None:
        if self.min_interval <= 0:
            return
        with self._throttle_lock:
            now = time.monotonic()
            if now < self._next_allowed:
                time.sleep(self._next_allowed - now)
                now = time.monotonic()
            self._next_allowed = now + self.min_interval

    def _fetch_json(self, url: str) -> Tuple[Optional[Dict[str, Any]], str]:
        session = self._session()
        for attempt in range(self.retries + 1):
            self._throttle()
            try:
                resp = session.get(url, timeout=self.timeout)
                code = int(resp.status_code)
                if code == 200:
                    self.stats["http_success"] += 1
                    return resp.json(), "ok"
                if code in {400, 404}:
                    self.stats[f"http_{code}_not_found"] += 1
                    return None, "not_found"
                self.stats[f"http_{code}_retryable"] += 1
            except RequestException:
                self.stats["http_error_retryable"] += 1

            if attempt < self.retries:
                time.sleep(min(4.0, 0.4 * (2 ** attempt)))
            else:
                self.stats["http_failed_final"] += 1
                return None, "failed"
        return None, "failed"

    @staticmethod
    def _extract_fields(payload: Optional[Dict[str, Any]]) -> Dict[str, str]:
        if not payload:
            return {"inchi": "", "canonical_smiles": ""}

        props = payload.get("PropertyTable", {}).get("Properties", [])
        if not isinstance(props, list) or not props:
            return {"inchi": "", "canonical_smiles": ""}

        row = props[0] if isinstance(props[0], dict) else {}
        inchi = normalize(str(row.get("InChI", "")))
        smiles = normalize(
            str(
                row.get("ConnectivitySMILES", "")
                or row.get("CanonicalSMILES", "")
                or row.get("SMILES", "")
            )
        )
        return {"inchi": inchi, "canonical_smiles": smiles}

    def query_by_inchikey(self, inchikey: str) -> Dict[str, str]:
        ik = normalize(inchikey).upper()
        cache_key = f"IK::{ik}"
        cached = self.cache.get(cache_key)
        if cached is not None:
            self.stats["cache_hit"] += 1
            return cached

        self.stats["cache_miss"] += 1
        self.stats["request_inchikey"] += 1
        url = f"{PUBCHEM_BASE}/inchikey/{ik}/property/InChI,CanonicalSMILES/JSON"
        payload, outcome = self._fetch_json(url)
        values = self._extract_fields(payload)
        result = {
            "inchi": values["inchi"],
            "canonical_smiles": values["canonical_smiles"],
            "status": outcome,
        }

        if outcome in {"ok", "not_found"}:
            self.cache.set(cache_key, result)
        else:
            self.stats["cache_skip_failed"] += 1
        return result


def non_empty_rate(rows: Sequence[Dict[str, str]], column: str) -> float:
    if not rows:
        return 0.0
    filled = sum(1 for r in rows if normalize(r.get(column, "")))
    return filled / len(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--core-tsv", type=Path, required=True)
    ap.add_argument("--chembl-db", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--build-report", type=Path, required=True)
    ap.add_argument("--max-rows", type=int, default=None)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--qps", type=float, default=6.0)
    ap.add_argument("--timeout", type=int, default=20)
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--progress-every", type=int, default=300)
    ap.add_argument(
        "--cache-path",
        type=Path,
        default=Path("pipelines/molecule_structure_identifiers/.cache/pubchem_structure_lookup_v1.jsonl"),
    )
    args = ap.parse_args()

    for p in [args.core_tsv, args.chembl_db]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    created_at = utc_now()
    fetch_date = utc_today()
    chembl_version = f"ChEMBL:{args.chembl_db.name}"
    pubchem_version = f"PubChem:PUGREST:{fetch_date}"

    header_in, input_rows = read_tsv(args.core_tsv, max_rows=args.max_rows)
    if "inchikey" not in header_in:
        raise SystemExit(f"[ERROR] inchikey not found in input: {args.core_tsv}")

    input_rows_norm: List[Dict[str, str]] = []
    invalid_inchikey_rows = 0
    for row in input_rows:
        record = dict(row)
        ik = normalize(record.get("inchikey", "")).upper()
        record["inchikey"] = ik
        if not valid_inchikey(ik):
            invalid_inchikey_rows += 1
        input_rows_norm.append(record)

    unique_iks = sorted({r["inchikey"] for r in input_rows_norm if valid_inchikey(r.get("inchikey", ""))})
    chembl_map = load_chembl_structures(args.chembl_db, unique_iks)

    need_pubchem: List[str] = []
    for ik in unique_iks:
        chem = chembl_map.get(ik, {})
        if not normalize(chem.get("inchi", "")) or not normalize(chem.get("canonical_smiles", "")):
            need_pubchem.append(ik)

    cache = JsonlCache(args.cache_path)
    client = PubChemClient(
        cache=cache,
        timeout=args.timeout,
        retries=args.retries,
        qps=args.qps,
    )

    pubchem_map: Dict[str, Dict[str, str]] = {}
    total_lookup = len(need_pubchem)
    if total_lookup > 0:
        print(f"[INFO] pubchem lookups needed: {total_lookup}")
        with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
            futures = {pool.submit(client.query_by_inchikey, ik): ik for ik in need_pubchem}
            for idx, future in enumerate(as_completed(futures), start=1):
                ik = futures[future]
                pubchem_map[ik] = future.result()
                if args.progress_every > 0 and idx % args.progress_every == 0:
                    print(f"[INFO] pubchem progress {idx}/{total_lookup}")

    need_pubchem_set = set(need_pubchem)
    output_rows: List[Dict[str, str]] = []

    source_counter = Counter()
    inchi_source_counter = Counter()
    smiles_source_counter = Counter()

    for row in input_rows_norm:
        ik = normalize(row.get("inchikey", "")).upper()
        chembl = chembl_map.get(ik, {})
        pubchem = pubchem_map.get(ik, {})

        chembl_inchi = normalize(chembl.get("inchi", ""))
        chembl_smiles = normalize(chembl.get("canonical_smiles", ""))
        pubchem_inchi = normalize(pubchem.get("inchi", ""))
        pubchem_smiles = normalize(pubchem.get("canonical_smiles", ""))

        inchi = chembl_inchi or pubchem_inchi
        canonical_smiles = chembl_smiles or pubchem_smiles

        if chembl_inchi:
            inchi_source = "chembl_compound_structures.standard_inchi"
        elif pubchem_inchi:
            inchi_source = "pubchem_pug_rest.inchi"
        else:
            inchi_source = ""

        if chembl_smiles:
            smiles_source = "chembl_compound_structures.canonical_smiles"
        elif pubchem_smiles:
            smiles_source = "pubchem_pug_rest.canonical_smiles"
        else:
            smiles_source = ""

        selected_sources = [s for s in [inchi_source, smiles_source] if s]
        selected_versions = []
        if any(s.startswith("chembl_") for s in selected_sources):
            selected_versions.append(chembl_version)
        if any(s.startswith("pubchem_") for s in selected_sources):
            selected_versions.append(pubchem_version)

        if selected_sources:
            row_source = join_multi(selected_sources)
            row_source_version = join_multi(selected_versions)
        else:
            attempts = ["chembl_compound_structures"]
            versions = [chembl_version]
            if ik in need_pubchem_set:
                attempts.append("pubchem_pug_rest")
                versions.append(pubchem_version)
            row_source = join_multi(attempts)
            row_source_version = join_multi(versions)

        out_row = {
            "inchikey": ik,
            "chembl_id": normalize(row.get("chembl_id", "")),
            "pubchem_cid": normalize(row.get("pubchem_cid", "")),
            "inchi": inchi,
            "canonical_smiles": canonical_smiles,
            "inchi_source": inchi_source,
            "canonical_smiles_source": smiles_source,
            "source": row_source,
            "source_version": row_source_version,
            "fetch_date": fetch_date,
        }
        output_rows.append(out_row)

        source_counter[row_source] += 1
        inchi_source_counter[inchi_source or "missing"] += 1
        smiles_source_counter[smiles_source or "missing"] += 1

    header_out = [
        "inchikey",
        "chembl_id",
        "pubchem_cid",
        "inchi",
        "canonical_smiles",
        "inchi_source",
        "canonical_smiles_source",
        "source",
        "source_version",
        "fetch_date",
    ]

    rows_written = write_tsv(args.out, header=header_out, rows=output_rows)

    inchi_rate = non_empty_rate(output_rows, "inchi")
    smiles_rate = non_empty_rate(output_rows, "canonical_smiles")
    source_rate = non_empty_rate(output_rows, "source")
    source_version_rate = non_empty_rate(output_rows, "source_version")
    fetch_date_rate = non_empty_rate(output_rows, "fetch_date")

    mappable_rows = [r for r in output_rows if normalize(r.get("chembl_id", "")) or normalize(r.get("pubchem_cid", ""))]
    inchi_mappable_rate = non_empty_rate(mappable_rows, "inchi") if mappable_rows else 0.0

    build_report = {
        "name": "molecule_structure_identifiers_v1.build",
        "created_at": created_at,
        "inputs": {
            "core_tsv": str(args.core_tsv),
            "chembl_db": str(args.chembl_db),
            "cache_path": str(args.cache_path),
            "cache_preloaded_entries": cache.loaded,
        },
        "output": str(args.out),
        "metrics": {
            "rows_input": len(input_rows_norm),
            "rows_written": rows_written,
            "invalid_inchikey_rows": invalid_inchikey_rows,
            "unique_valid_inchikeys": len(unique_iks),
            "chembl_matches": len(chembl_map),
            "pubchem_lookup_rows": total_lookup,
            "pubchem_hit_any_field": sum(
                1
                for rec in pubchem_map.values()
                if normalize(rec.get("inchi", "")) or normalize(rec.get("canonical_smiles", ""))
            ),
            "inchi_non_empty_rate": inchi_rate,
            "inchi_non_empty_rate_mappable_subset": inchi_mappable_rate,
            "canonical_smiles_non_empty_rate": smiles_rate,
            "source_non_empty_rate": source_rate,
            "source_version_non_empty_rate": source_version_rate,
            "fetch_date_non_empty_rate": fetch_date_rate,
            "inchi_source_counts": dict(inchi_source_counter),
            "canonical_smiles_source_counts": dict(smiles_source_counter),
            "top_row_source_counts": dict(source_counter.most_common(8)),
            "pubchem_client_stats": dict(client.stats),
        },
        "source_version_tags": {
            "chembl": chembl_version,
            "pubchem": pubchem_version,
        },
        "sample_mode": args.max_rows is not None,
        "max_rows": args.max_rows,
    }

    args.build_report.parent.mkdir(parents=True, exist_ok=True)
    args.build_report.write_text(json.dumps(build_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(
        "[OK] molecule_structure_identifiers_v1 rows={} inchi_rate={:.4f} smiles_rate={:.4f}".format(
            rows_written,
            inchi_rate,
            smiles_rate,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
