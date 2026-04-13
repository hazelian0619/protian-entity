#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import requests
from requests.exceptions import RequestException, RequestsDependencyWarning
import warnings

warnings.filterwarnings("ignore", category=RequestsDependencyWarning)

INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
CHEMBL_RE = re.compile(r"^CHEMBL\d+$")
DIGIT_RE = re.compile(r"\d+")
CHEMBL_VERSION_RE = re.compile(r"chembl[_-]?(\d+)", re.IGNORECASE)

DESCRIPTOR_FIELDS = [
    "molecular_weight",
    "logp",
    "tpsa",
    "hbd",
    "hba",
    "rotatable_bonds",
]

UNIT_STRATEGY = {
    "molecular_weight": "Da",
    "logp": "unitless",
    "tpsa": "Å²",
    "hbd": "count",
    "hba": "count",
    "rotatable_bonds": "count",
}
UNIT_STRATEGY_TEXT = "MW:Da|LogP:unitless|TPSA:Å²|HBD:count|HBA:count|RotatableBonds:count"

SOURCE_CHEMBL = "chembl_compound_properties"
SOURCE_PUBCHEM = "pubchem_pugrest_property"
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
PUBCHEM_PROPERTIES = [
    "MolecularWeight",
    "XLogP",
    "TPSA",
    "HBondDonorCount",
    "HBondAcceptorCount",
    "RotatableBondCount",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def utc_today() -> str:
    return datetime.now(timezone.utc).date().isoformat()


def normalize(x: object) -> str:
    return str(x or "").strip()


def split_multi(x: object) -> List[str]:
    text = normalize(x)
    if not text:
        return []
    return [t.strip() for t in text.split(";") if t.strip()]


def normalize_inchikey(x: object) -> str:
    return normalize(x).upper()


def valid_inchikey(x: str) -> bool:
    return bool(INCHIKEY_RE.match(normalize_inchikey(x)))


def parse_pubchem_cid(token: object) -> Optional[str]:
    text = normalize(token)
    if not text:
        return None
    if text.isdigit():
        return text
    m = DIGIT_RE.search(text)
    if m:
        return m.group(0)
    return None


def fmt_float(value: float, digits: int = 6) -> str:
    s = f"{value:.{digits}f}".rstrip("0").rstrip(".")
    return s if s else "0"


def parse_float(value: object) -> Optional[float]:
    text = normalize(value)
    if not text:
        return None
    try:
        return float(text)
    except Exception:
        return None


def parse_int_like(value: object) -> Optional[int]:
    text = normalize(value)
    if not text:
        return None
    try:
        f = float(text)
    except Exception:
        return None
    if f != f:
        return None
    return int(round(f))


def empty_descriptor_record() -> Dict[str, str]:
    return {k: "" for k in DESCRIPTOR_FIELDS}


def non_empty_count(record: Dict[str, str]) -> int:
    return sum(1 for k in DESCRIPTOR_FIELDS if normalize(record.get(k, "")))


def choose_best_record(records: Sequence[Dict[str, str]]) -> Dict[str, str]:
    best: Dict[str, str] = empty_descriptor_record()
    best_score = -1
    for rec in records:
        score = non_empty_count(rec)
        if score > best_score:
            best = {k: normalize(rec.get(k, "")) for k in DESCRIPTOR_FIELDS}
            best_score = score
    return best


def merge_descriptors(primary: Dict[str, str], fallback: Dict[str, str]) -> Dict[str, str]:
    merged = empty_descriptor_record()
    for field in DESCRIPTOR_FIELDS:
        p = normalize(primary.get(field, ""))
        f = normalize(fallback.get(field, ""))
        merged[field] = p if p else f
    return merged


def parse_pubchem_property_entry(entry: Dict[str, Any]) -> Dict[str, str]:
    out = empty_descriptor_record()

    mw = parse_float(entry.get("MolecularWeight"))
    if mw is not None:
        out["molecular_weight"] = fmt_float(mw)

    logp = parse_float(entry.get("XLogP"))
    if logp is not None:
        out["logp"] = fmt_float(logp)

    tpsa = parse_float(entry.get("TPSA"))
    if tpsa is not None:
        out["tpsa"] = fmt_float(tpsa)

    hbd = parse_int_like(entry.get("HBondDonorCount"))
    if hbd is not None:
        out["hbd"] = str(hbd)

    hba = parse_int_like(entry.get("HBondAcceptorCount"))
    if hba is not None:
        out["hba"] = str(hba)

    rb = parse_int_like(entry.get("RotatableBondCount"))
    if rb is not None:
        out["rotatable_bonds"] = str(rb)

    return out


def parse_chembl_property_row(row: Tuple[Any, ...]) -> Dict[str, str]:
    out = empty_descriptor_record()
    mw = parse_float(row[0])
    if mw is not None:
        out["molecular_weight"] = fmt_float(mw)

    logp = parse_float(row[1])
    if logp is not None:
        out["logp"] = fmt_float(logp)

    tpsa = parse_float(row[2])
    if tpsa is not None:
        out["tpsa"] = fmt_float(tpsa)

    hbd = parse_int_like(row[3])
    if hbd is not None:
        out["hbd"] = str(hbd)

    hba = parse_int_like(row[4])
    if hba is not None:
        out["hba"] = str(hba)

    rb = parse_int_like(row[5])
    if rb is not None:
        out["rotatable_bonds"] = str(rb)
    return out


def detect_chembl_source_version(chembl_db: Path) -> str:
    tokens = [chembl_db.name, str(chembl_db.parent.name), str(chembl_db)]
    for token in tokens:
        m = CHEMBL_VERSION_RE.search(token)
        if m:
            return f"ChEMBL:{m.group(1)}:compound_properties"
    return "ChEMBL:unknown:compound_properties"


def read_xref(path: Path, max_rows: Optional[int] = None) -> List[Dict[str, Any]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        required = {"inchikey", "chembl_id", "pubchem_cid"}
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise SystemExit(f"[ERROR] xref missing columns: {missing}")

        rows: List[Dict[str, Any]] = []
        for idx, row in enumerate(reader, start=1):
            inchikey = normalize_inchikey(row.get("inchikey", ""))
            chembl_ids = [c.upper() for c in split_multi(row.get("chembl_id", "")) if CHEMBL_RE.match(c.upper())]
            pubchem_ids = []
            for token in split_multi(row.get("pubchem_cid", "")):
                cid = parse_pubchem_cid(token)
                if cid:
                    pubchem_ids.append(cid)

            rows.append(
                {
                    "inchikey": inchikey,
                    "chembl_ids": chembl_ids,
                    "pubchem_ids": list(dict.fromkeys(pubchem_ids)),
                }
            )
            if max_rows is not None and idx >= max_rows:
                break
    return rows


def load_chembl_properties(chembl_db: Path, chembl_ids: Sequence[str]) -> Dict[str, Dict[str, str]]:
    if not chembl_db.exists():
        raise SystemExit(f"[ERROR] missing ChEMBL db: {chembl_db}")
    if not chembl_ids:
        return {}

    uniq = sorted(set([normalize(x).upper() for x in chembl_ids if CHEMBL_RE.match(normalize(x).upper())]))
    if not uniq:
        return {}

    out: Dict[str, Dict[str, str]] = {}
    with sqlite3.connect(chembl_db) as conn:
        conn.execute("DROP TABLE IF EXISTS temp.tmp_chembl_id")
        conn.execute("CREATE TEMP TABLE tmp_chembl_id (chembl_id TEXT PRIMARY KEY)")
        conn.executemany("INSERT OR IGNORE INTO tmp_chembl_id(chembl_id) VALUES (?)", [(x,) for x in uniq])

        q = """
            SELECT UPPER(TRIM(md.chembl_id)) AS chembl_id,
                   cp.full_mwt,
                   cp.alogp,
                   cp.psa,
                   cp.hbd,
                   cp.hba,
                   cp.rtb
            FROM tmp_chembl_id t
            JOIN molecule_dictionary md ON UPPER(TRIM(md.chembl_id)) = t.chembl_id
            JOIN compound_properties cp ON cp.molregno = md.molregno
        """
        for chembl_id, mw, logp, tpsa, hbd, hba, rtb in conn.execute(q):
            rec = parse_chembl_property_row((mw, logp, tpsa, hbd, hba, rtb))
            prev = out.get(chembl_id)
            if prev is None or non_empty_count(rec) > non_empty_count(prev):
                out[chembl_id] = rec
    return out


class PubChemCache:
    def __init__(self, path: Path):
        self.path = path
        self.data: Dict[str, Dict[str, str]] = {}
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
                cid = normalize(obj.get("cid", ""))
                props = obj.get("props", {})
                if cid and isinstance(props, dict):
                    self.data[cid] = {k: normalize(v) for k, v in props.items()}

    def get(self, cid: str) -> Optional[Dict[str, str]]:
        rec = self.data.get(cid)
        return dict(rec) if rec is not None else None

    def set(self, cid: str, props: Dict[str, str]) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        cleaned = {k: normalize(v) for k, v in props.items()}
        payload = {
            "cid": cid,
            "props": cleaned,
            "updated_at": utc_now(),
        }
        with self.path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(payload, ensure_ascii=False) + "\n")
        self.data[cid] = cleaned


class InChIKeyCidCache:
    def __init__(self, path: Path):
        self.path = path
        self.data: Dict[str, List[str]] = {}
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
                inchikey = normalize_inchikey(obj.get("inchikey", ""))
                cids_raw = obj.get("cids", [])
                cids = [normalize(x) for x in cids_raw if normalize(x).isdigit()]
                if inchikey:
                    self.data[inchikey] = list(dict.fromkeys(cids))

    def get(self, inchikey: str) -> Optional[List[str]]:
        rec = self.data.get(normalize_inchikey(inchikey))
        return list(rec) if rec is not None else None

    def set(self, inchikey: str, cids: Sequence[str]) -> None:
        ik = normalize_inchikey(inchikey)
        uniq = list(dict.fromkeys([normalize(c) for c in cids if normalize(c).isdigit()]))
        self.path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "inchikey": ik,
            "cids": uniq,
            "updated_at": utc_now(),
        }
        with self.path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(payload, ensure_ascii=False) + "\n")
        self.data[ik] = uniq


def chunked(items: Sequence[str], size: int) -> Iterable[List[str]]:
    for i in range(0, len(items), size):
        yield list(items[i : i + size])


def fetch_pubchem_properties(
    cids: Sequence[str],
    cache: PubChemCache,
    timeout: int,
    retries: int,
    batch_size: int,
    qps: float,
) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Any]]:
    out: Dict[str, Dict[str, str]] = {}
    uniq = sorted(set([normalize(c) for c in cids if normalize(c)]), key=lambda x: int(x))

    cache_hits = 0
    for cid in uniq:
        cached = cache.get(cid)
        if cached is not None:
            out[cid] = cached
            cache_hits += 1

    pending = [cid for cid in uniq if cid not in out]

    batch_size = max(1, int(batch_size))
    sleep_s = 0.0 if qps <= 0 else (1.0 / qps)
    api_batches = 0
    api_errors: List[str] = []
    api_success_cids = 0

    session = requests.Session()
    props = ",".join(PUBCHEM_PROPERTIES)

    for batch in chunked(pending, batch_size):
        api_batches += 1
        url = f"{PUBCHEM_BASE}/cid/{','.join(batch)}/property/{props}/JSON"
        payload: Optional[Dict[str, Any]] = None

        for attempt in range(1, retries + 2):
            try:
                resp = session.get(url, timeout=timeout)
                if resp.status_code == 404:
                    payload = {"PropertyTable": {"Properties": []}}
                    break
                resp.raise_for_status()
                payload = resp.json()
                break
            except (RequestException, ValueError) as exc:
                if attempt >= retries + 1:
                    api_errors.append(f"batch[{batch[0]}..{batch[-1]}]: {exc}")
                    payload = {"PropertyTable": {"Properties": []}}
                else:
                    time.sleep(min(2.0, 0.3 * attempt))

        props_list = payload.get("PropertyTable", {}).get("Properties", []) if payload else []
        got: Dict[str, Dict[str, str]] = {}
        for entry in props_list:
            cid = normalize(entry.get("CID", ""))
            if not cid:
                continue
            parsed = parse_pubchem_property_entry(entry)
            got[cid] = parsed

        for cid in batch:
            parsed = got.get(cid, empty_descriptor_record())
            cache.set(cid, parsed)
            out[cid] = parsed
            if non_empty_count(parsed) > 0:
                api_success_cids += 1

        if sleep_s > 0:
            time.sleep(sleep_s)

    metrics = {
        "requested_unique_cids": len(uniq),
        "cache_hits": cache_hits,
        "fetched_from_api": len(pending),
        "api_batches": api_batches,
        "api_success_cids_with_any_value": api_success_cids,
        "api_error_count": len(api_errors),
        "api_error_samples": api_errors[:20],
    }
    return out, metrics


def resolve_pubchem_cids_from_inchikey(
    inchikeys: Sequence[str],
    cache: InChIKeyCidCache,
    timeout: int,
    retries: int,
    qps: float,
) -> Tuple[Dict[str, List[str]], Dict[str, Any]]:
    uniq = [ik for ik in dict.fromkeys([normalize_inchikey(x) for x in inchikeys if valid_inchikey(x)])]
    out: Dict[str, List[str]] = {}
    cache_hits = 0
    api_calls = 0
    api_error_samples: List[str] = []
    resolved_non_empty = 0
    sleep_s = 0.0 if qps <= 0 else (1.0 / qps)
    session = requests.Session()

    for ik in uniq:
        cached = cache.get(ik)
        if cached is not None:
            out[ik] = list(cached)
            cache_hits += 1
            if cached:
                resolved_non_empty += 1
            continue

        api_calls += 1
        url = f"{PUBCHEM_BASE}/inchikey/{ik}/cids/JSON"
        cids: List[str] = []
        for attempt in range(1, retries + 2):
            try:
                resp = session.get(url, timeout=timeout)
                if resp.status_code == 404:
                    cids = []
                    break
                resp.raise_for_status()
                payload = resp.json()
                raw = payload.get("IdentifierList", {}).get("CID", [])
                cids = [normalize(x) for x in raw if normalize(x).isdigit()]
                cids = list(dict.fromkeys(cids))
                break
            except (RequestException, ValueError) as exc:
                if attempt >= retries + 1:
                    api_error_samples.append(f"{ik}: {exc}")
                    cids = []
                else:
                    time.sleep(min(3.0, 0.4 * attempt))

        cache.set(ik, cids)
        out[ik] = cids
        if cids:
            resolved_non_empty += 1
        if sleep_s > 0:
            time.sleep(sleep_s)

    metrics = {
        "requested_inchikeys": len(uniq),
        "cache_hits": cache_hits,
        "api_calls": api_calls,
        "resolved_non_empty_cid_lists": resolved_non_empty,
        "api_error_count": len(api_error_samples),
        "api_error_samples": api_error_samples[:20],
    }
    return out, metrics


def write_tsv(path: Path, rows: Iterable[Dict[str, str]], header: Sequence[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    n = 0
    with tmp.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=list(header), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in header})
            n += 1
    tmp.replace(path)
    return n


def build_table(
    xref_rows: Sequence[Dict[str, Any]],
    chembl_map: Dict[str, Dict[str, str]],
    pubchem_map: Dict[str, Dict[str, str]],
    pubchem_inchikey_map: Dict[str, Dict[str, str]],
    chembl_source_version: str,
    pubchem_source_version: str,
    fetch_date: str,
) -> Tuple[List[Dict[str, str]], Dict[str, Any]]:
    rows_out: List[Dict[str, str]] = []

    descriptor_non_empty = Counter()
    mappable_non_empty = Counter()
    source_usage = Counter()

    mappable_rows = 0
    rows_with_any = 0
    rows_with_core_triple = 0

    for rec in xref_rows:
        chembl_candidates = [chembl_map[cid] for cid in rec["chembl_ids"] if cid in chembl_map]
        pubchem_candidates = [pubchem_map[cid] for cid in rec["pubchem_ids"] if cid in pubchem_map]
        pubchem_ik_best = pubchem_inchikey_map.get(rec["inchikey"], empty_descriptor_record())

        chembl_best = choose_best_record(chembl_candidates)
        pubchem_best = choose_best_record(pubchem_candidates)
        pubchem_combined = merge_descriptors(pubchem_best, pubchem_ik_best)
        merged = merge_descriptors(chembl_best, pubchem_combined)

        has_mapping = bool(rec["chembl_ids"] or rec["pubchem_ids"])
        if has_mapping:
            mappable_rows += 1

        used_chembl = False
        used_pubchem = False
        for field in DESCRIPTOR_FIELDS:
            value = normalize(merged.get(field, ""))
            if value:
                descriptor_non_empty[field] += 1
                if has_mapping:
                    mappable_non_empty[field] += 1

                if normalize(chembl_best.get(field, "")):
                    used_chembl = True
                elif normalize(pubchem_combined.get(field, "")):
                    used_pubchem = True

        if non_empty_count(merged) > 0:
            rows_with_any += 1

        if all(normalize(merged.get(f, "")) for f in ["molecular_weight", "logp", "tpsa"]):
            rows_with_core_triple += 1

        source_tokens: List[str] = []
        source_version_tokens: List[str] = []
        if used_chembl:
            source_tokens.append(SOURCE_CHEMBL)
            source_version_tokens.append(chembl_source_version)
        if used_pubchem:
            source_tokens.append(SOURCE_PUBCHEM)
            source_version_tokens.append(pubchem_source_version)

        if used_chembl and used_pubchem:
            source_usage["hybrid"] += 1
        elif used_chembl:
            source_usage["chembl_only"] += 1
        elif used_pubchem:
            source_usage["pubchem_only"] += 1
        else:
            source_usage["none"] += 1

        rows_out.append(
            {
                "inchikey": rec["inchikey"],
                **merged,
                "unit_strategy": UNIT_STRATEGY_TEXT,
                "source": ";".join(source_tokens) if source_tokens else "none",
                "source_version": ";".join(source_version_tokens) if source_version_tokens else "none",
                "fetch_date": fetch_date,
            }
        )

    n = len(rows_out)
    metrics = {
        "rows_written": n,
        "mappable_subset_rows": mappable_rows,
        "rows_with_any_descriptor": rows_with_any,
        "rows_with_mw_logp_tpsa": rows_with_core_triple,
        "core_triple_coverage_rate": (rows_with_core_triple / n) if n else 0.0,
        "descriptor_non_empty_counts": {f: int(descriptor_non_empty[f]) for f in DESCRIPTOR_FIELDS},
        "descriptor_non_empty_rates_all_rows": {
            f: ((descriptor_non_empty[f] / n) if n else 0.0) for f in DESCRIPTOR_FIELDS
        },
        "descriptor_non_empty_rates_mappable_subset": {
            f: ((mappable_non_empty[f] / mappable_rows) if mappable_rows else 0.0) for f in DESCRIPTOR_FIELDS
        },
        "source_usage": dict(source_usage),
    }
    return rows_out, metrics


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Build molecule physicochemical descriptor layer v1")
    ap.add_argument("--xref", type=Path, default=Path("data/output/molecules/molecule_xref_core_v2.tsv"))
    ap.add_argument("--chembl-db", type=Path, default=Path("data/raw/molecules/chembl_36/chembl_36.db"))
    ap.add_argument("--out", type=Path, default=Path("data/output/molecules/molecule_physchem_descriptors_v1.tsv"))
    ap.add_argument(
        "--report",
        type=Path,
        default=Path("pipelines/molecule_physchem_descriptors/reports/molecule_physchem_descriptors_v1.build.json"),
    )
    ap.add_argument(
        "--pubchem-cache",
        type=Path,
        default=Path("pipelines/molecule_physchem_descriptors/.cache/pubchem_properties_v1.jsonl"),
    )
    ap.add_argument(
        "--pubchem-inchikey-cache",
        type=Path,
        default=Path("pipelines/molecule_physchem_descriptors/.cache/pubchem_inchikey_to_cid_v1.jsonl"),
    )
    ap.add_argument("--max-rows", type=int, default=None)
    ap.add_argument("--disable-pubchem", action="store_true")
    ap.add_argument("--disable-pubchem-inchikey-fallback", action="store_true")
    ap.add_argument("--pubchem-timeout", type=int, default=30)
    ap.add_argument("--pubchem-retries", type=int, default=2)
    ap.add_argument("--pubchem-batch-size", type=int, default=100)
    ap.add_argument("--pubchem-qps", type=float, default=5.0)
    ap.add_argument("--pubchem-inchikey-qps", type=float, default=8.0)
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    if not args.xref.exists():
        raise SystemExit(f"[ERROR] missing xref: {args.xref}")

    xref_rows = read_xref(args.xref, max_rows=args.max_rows)
    chembl_ids: List[str] = []
    pubchem_ids: List[str] = []
    bad_inchikey = 0
    for row in xref_rows:
        if not valid_inchikey(row["inchikey"]):
            bad_inchikey += 1
        chembl_ids.extend(row["chembl_ids"])
        pubchem_ids.extend(row["pubchem_ids"])

    chembl_map = load_chembl_properties(args.chembl_db, chembl_ids)
    chembl_source_version = detect_chembl_source_version(args.chembl_db)

    pubchem_map: Dict[str, Dict[str, str]] = {}
    pubchem_inchikey_map: Dict[str, Dict[str, str]] = {}
    pubchem_fallback_metrics: Dict[str, Any] = {
        "candidate_inchikeys": 0,
        "inchikey_to_cid": {
            "requested_inchikeys": 0,
            "cache_hits": 0,
            "api_calls": 0,
            "resolved_non_empty_cid_lists": 0,
            "api_error_count": 0,
            "api_error_samples": [],
        },
        "resolved_unique_cids": 0,
        "property_fetch_for_resolved_cids": {
            "requested_unique_cids": 0,
            "cache_hits": 0,
            "fetched_from_api": 0,
            "api_batches": 0,
            "api_success_cids_with_any_value": 0,
            "api_error_count": 0,
            "api_error_samples": [],
        },
        "rows_with_any_inchikey_fallback_value": 0,
    }
    pubchem_metrics = {
        "requested_unique_cids": 0,
        "cache_hits": 0,
        "fetched_from_api": 0,
        "api_batches": 0,
        "api_success_cids_with_any_value": 0,
        "api_error_count": 0,
        "api_error_samples": [],
    }

    if not args.disable_pubchem:
        cache = PubChemCache(args.pubchem_cache)
        pubchem_map, pubchem_metrics = fetch_pubchem_properties(
            cids=pubchem_ids,
            cache=cache,
            timeout=args.pubchem_timeout,
            retries=max(0, args.pubchem_retries),
            batch_size=max(1, args.pubchem_batch_size),
            qps=max(0.0, args.pubchem_qps),
        )

        if not args.disable_pubchem_inchikey_fallback:
            fallback_candidates: List[str] = []
            for rec in xref_rows:
                if not rec["chembl_ids"] or not valid_inchikey(rec["inchikey"]):
                    continue
                chembl_best = choose_best_record([chembl_map[cid] for cid in rec["chembl_ids"] if cid in chembl_map])
                pubchem_best = choose_best_record([pubchem_map[cid] for cid in rec["pubchem_ids"] if cid in pubchem_map])
                merged_pre = merge_descriptors(chembl_best, pubchem_best)
                if any(not normalize(merged_pre.get(field, "")) for field in DESCRIPTOR_FIELDS):
                    fallback_candidates.append(rec["inchikey"])

            pubchem_fallback_metrics["candidate_inchikeys"] = len(fallback_candidates)
            ik_cache = InChIKeyCidCache(args.pubchem_inchikey_cache)
            ik_to_cids, ik_metrics = resolve_pubchem_cids_from_inchikey(
                inchikeys=fallback_candidates,
                cache=ik_cache,
                timeout=args.pubchem_timeout,
                retries=max(0, args.pubchem_retries),
                qps=max(0.0, args.pubchem_inchikey_qps),
            )
            pubchem_fallback_metrics["inchikey_to_cid"] = ik_metrics

            resolved_cids = sorted(
                set([cid for cids in ik_to_cids.values() for cid in cids if cid]),
                key=lambda x: int(x),
            )
            pubchem_fallback_metrics["resolved_unique_cids"] = len(resolved_cids)

            if resolved_cids:
                resolved_props, resolved_prop_metrics = fetch_pubchem_properties(
                    cids=resolved_cids,
                    cache=cache,
                    timeout=args.pubchem_timeout,
                    retries=max(0, args.pubchem_retries),
                    batch_size=max(1, args.pubchem_batch_size),
                    qps=max(0.0, args.pubchem_qps),
                )
                pubchem_map.update(resolved_props)
                pubchem_fallback_metrics["property_fetch_for_resolved_cids"] = resolved_prop_metrics

            rows_with_fallback_any = 0
            for ik, cids in ik_to_cids.items():
                rec = choose_best_record([pubchem_map[cid] for cid in cids if cid in pubchem_map])
                if non_empty_count(rec) > 0:
                    pubchem_inchikey_map[ik] = rec
                    rows_with_fallback_any += 1
            pubchem_fallback_metrics["rows_with_any_inchikey_fallback_value"] = rows_with_fallback_any

    fetch_date = utc_today()
    pubchem_source_version = f"PubChem:PUGREST:compound/property:{fetch_date}"

    rows_out, table_metrics = build_table(
        xref_rows=xref_rows,
        chembl_map=chembl_map,
        pubchem_map=pubchem_map,
        pubchem_inchikey_map=pubchem_inchikey_map,
        chembl_source_version=chembl_source_version,
        pubchem_source_version=pubchem_source_version,
        fetch_date=fetch_date,
    )

    header = [
        "inchikey",
        *DESCRIPTOR_FIELDS,
        "unit_strategy",
        "source",
        "source_version",
        "fetch_date",
    ]
    written = write_tsv(args.out, rows_out, header)

    report = {
        "name": "molecule_physchem_descriptors_v1.build",
        "created_at": utc_now(),
        "inputs": {
            "xref": str(args.xref),
            "chembl_db": str(args.chembl_db),
            "pubchem_enabled": (not args.disable_pubchem),
            "pubchem_cache": str(args.pubchem_cache),
            "pubchem_inchikey_cache": str(args.pubchem_inchikey_cache),
            "pubchem_inchikey_fallback_enabled": (not args.disable_pubchem_inchikey_fallback),
        },
        "output": str(args.out),
        "unit_strategy": UNIT_STRATEGY,
        "sources": {
            "chembl": SOURCE_CHEMBL,
            "pubchem": SOURCE_PUBCHEM,
            "chembl_source_version": chembl_source_version,
            "pubchem_source_version": pubchem_source_version,
        },
        "metrics": {
            "input_rows": len(xref_rows),
            "rows_written": written,
            "bad_inchikey_rows": bad_inchikey,
            "unique_chembl_ids_in_xref": len(set([x for x in chembl_ids if x])),
            "unique_pubchem_cids_in_xref": len(set([x for x in pubchem_ids if x])),
            "chembl_ids_with_properties": len(chembl_map),
            **table_metrics,
            "pubchem_fetch": pubchem_metrics,
            "pubchem_inchikey_fallback": pubchem_fallback_metrics,
        },
        "acceptance_preview": {
            "backlink_expected_rows_equal_input": (written == len(xref_rows)),
            "core_triple_coverage_rate": table_metrics["core_triple_coverage_rate"],
            "target_core_triple_coverage_min": 0.70,
            "target_mappable_non_empty_rate_min": 0.90,
        },
    }

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(
        "[OK] molecule_physchem_descriptors_v1 build "
        f"rows={written} core_coverage={table_metrics['core_triple_coverage_rate']:.4f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
