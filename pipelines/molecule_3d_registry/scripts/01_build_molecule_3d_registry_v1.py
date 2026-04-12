#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import requests

INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
CHEMBL_RE = re.compile(r"^CHEMBL\d+$")
CID_RE = re.compile(r"\d+")
ZINC_RE = re.compile(r"ZINC\s*([0-9]{1,15})", flags=re.IGNORECASE)


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
    return [s.strip() for s in t.split(";") if s.strip()]


def valid_inchikey(x: str) -> bool:
    return bool(INCHIKEY_RE.match(normalize(x).upper()))


def valid_chembl(x: str) -> bool:
    return bool(CHEMBL_RE.match(normalize(x).upper()))


def parse_cid(x: str) -> Optional[str]:
    s = normalize(x)
    if not s:
        return None
    if s.isdigit():
        return s
    m = CID_RE.search(s)
    return m.group(0) if m else None


def parse_zinc_id(x: str) -> Optional[str]:
    text = normalize(x).upper()
    if not text:
        return None
    m = ZINC_RE.search(text)
    if not m:
        return None
    digits = m.group(1)
    return f"ZINC{digits.zfill(12)}"


def join_values(values: Iterable[str]) -> str:
    return ";".join(sorted({normalize(v) for v in values if normalize(v)}))


def write_tsv(path: Path, rows: Iterable[Dict[str, str]], header: Sequence[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    n = 0
    with tmp.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(header), lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in header})
            n += 1
    tmp.replace(path)
    return n


def read_tsv(path: Path, max_rows: Optional[int] = None) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        rows: List[Dict[str, str]] = []
        for i, row in enumerate(r, start=1):
            rows.append(row)
            if max_rows is not None and i >= max_rows:
                break
    return list(r.fieldnames), rows


@dataclass
class XrefRow:
    inchikey: str
    chembl_ids: Set[str]
    pubchem_cids: Set[str]


def load_xref_rows(path: Path, max_rows: Optional[int]) -> Tuple[List[XrefRow], Dict[str, Any]]:
    cols, raw_rows = read_tsv(path, max_rows=max_rows)
    if "inchikey" not in cols:
        raise SystemExit(f"[ERROR] xref table missing inchikey: {path}")

    rows: List[XrefRow] = []
    skipped_bad_ik = 0

    for row in raw_rows:
        ik = normalize(row.get("inchikey", "")).upper()
        if not valid_inchikey(ik):
            skipped_bad_ik += 1
            continue

        chembl_ids = {c.upper() for c in split_multi(row.get("chembl_id", "")) if valid_chembl(c.upper())}
        pubchem_cids = {cid for cid in (parse_cid(x) for x in split_multi(row.get("pubchem_cid", ""))) if cid}

        rows.append(XrefRow(inchikey=ik, chembl_ids=chembl_ids, pubchem_cids=pubchem_cids))

    stats = {
        "input_rows": len(raw_rows),
        "valid_inchikey_rows": len(rows),
        "skipped_bad_inchikey_rows": skipped_bad_ik,
        "rows_with_pubchem_cid": sum(1 for r in rows if r.pubchem_cids),
        "rows_with_chembl_id": sum(1 for r in rows if r.chembl_ids),
    }
    return rows, stats


def _fetch_pubchem_chunk(session: requests.Session, cids: Sequence[str], timeout: int) -> Dict[str, int]:
    cid_csv = ",".join(cids)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_csv}/property/ConformerCount3D/CSV"
    resp = session.get(url, timeout=timeout)
    if resp.status_code != 200:
        raise RuntimeError(f"HTTP {resp.status_code}")
    data: Dict[str, int] = {}
    reader = csv.DictReader(StringIO(resp.text))
    for row in reader:
        cid = normalize(row.get("CID", ""))
        val = normalize(row.get("ConformerCount3D", ""))
        if not cid:
            continue
        try:
            data[cid] = int(float(val)) if val else 0
        except Exception:
            data[cid] = 0
    return data


def _fetch_pubchem_registry_chunk(session: requests.Session, cids: Sequence[str], timeout: int) -> Dict[str, Set[str]]:
    cid_csv = ",".join(cids)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_csv}/xrefs/RegistryID/JSON"
    resp = session.get(url, timeout=timeout)
    if resp.status_code != 200:
        raise RuntimeError(f"HTTP {resp.status_code}")

    out: Dict[str, Set[str]] = {}
    payload = resp.json()
    infos = payload.get("InformationList", {}).get("Information", []) or []
    for info in infos:
        cid = normalize(str(info.get("CID", "")))
        regs = info.get("RegistryID", []) or []
        zids: Set[str] = set()
        for reg in regs:
            zid = parse_zinc_id(str(reg))
            if zid:
                zids.add(zid)
        if cid:
            out[cid] = zids
    return out


def fetch_pubchem_conformer_counts(
    cids: Sequence[str],
    *,
    chunk_size: int = 80,
    timeout: int = 25,
    sleep_seconds: float = 0.12,
    retries: int = 2,
) -> Tuple[Dict[str, int], Dict[str, Any]]:
    session = requests.Session()
    session.headers.update({"User-Agent": "kg-molecule-3d-registry/1.0"})

    cid_list = sorted({c for c in cids if c})
    out: Dict[str, int] = {}
    failed_singletons: Set[str] = set()

    metrics = Counter()

    def fetch_recursive(chunk: List[str]) -> None:
        if not chunk:
            return

        for attempt in range(retries + 1):
            metrics["pubchem_http_calls"] += 1
            try:
                got = _fetch_pubchem_chunk(session, chunk, timeout=timeout)
                out.update(got)
                metrics["pubchem_http_success_calls"] += 1
                return
            except Exception:
                metrics["pubchem_http_failed_calls"] += 1
                if attempt < retries:
                    time.sleep(0.35 * (attempt + 1))
                else:
                    break

        # failed after retries: split for fault isolation
        if len(chunk) == 1:
            failed_singletons.add(chunk[0])
            return

        mid = len(chunk) // 2
        fetch_recursive(chunk[:mid])
        fetch_recursive(chunk[mid:])

    for i in range(0, len(cid_list), chunk_size):
        chunk = cid_list[i : i + chunk_size]
        fetch_recursive(chunk)
        if sleep_seconds > 0:
            time.sleep(sleep_seconds)

    report = {
        "requested_cids": len(cid_list),
        "resolved_cids": len(out),
        "failed_cids": len(failed_singletons),
        "failed_cid_sample": sorted(list(failed_singletons))[:50],
        **{k: int(v) for k, v in metrics.items()},
    }
    return out, report


def fetch_pubchem_zinc_ids(
    cids: Sequence[str],
    *,
    chunk_size: int = 40,
    timeout: int = 30,
    sleep_seconds: float = 0.10,
    retries: int = 2,
) -> Tuple[Dict[str, Set[str]], Dict[str, Any]]:
    session = requests.Session()
    session.headers.update({"User-Agent": "kg-molecule-3d-registry/1.0"})

    cid_list = sorted({c for c in cids if c})
    out: Dict[str, Set[str]] = {}
    failed_singletons: Set[str] = set()
    metrics = Counter()

    def fetch_recursive(chunk: List[str]) -> None:
        if not chunk:
            return

        for attempt in range(retries + 1):
            metrics["pubchem_registry_http_calls"] += 1
            try:
                got = _fetch_pubchem_registry_chunk(session, chunk, timeout=timeout)
                out.update(got)
                metrics["pubchem_registry_http_success_calls"] += 1
                return
            except Exception:
                metrics["pubchem_registry_http_failed_calls"] += 1
                if attempt < retries:
                    time.sleep(0.35 * (attempt + 1))
                else:
                    break

        if len(chunk) == 1:
            failed_singletons.add(chunk[0])
            return

        mid = len(chunk) // 2
        fetch_recursive(chunk[:mid])
        fetch_recursive(chunk[mid:])

    for i in range(0, len(cid_list), chunk_size):
        chunk = cid_list[i : i + chunk_size]
        fetch_recursive(chunk)
        if sleep_seconds > 0:
            time.sleep(sleep_seconds)

    cid_with_zinc = sum(1 for cid in cid_list if out.get(cid))
    total_zinc = int(sum(len(v) for v in out.values()))
    report = {
        "requested_cids": len(cid_list),
        "resolved_cids": len(out),
        "cids_with_zinc_id": cid_with_zinc,
        "total_zinc_ids": total_zinc,
        "failed_cids": len(failed_singletons),
        "failed_cid_sample": sorted(list(failed_singletons))[:50],
        **{k: int(v) for k, v in metrics.items()},
    }
    return out, report


def load_zinc_map_from_chembl(chembl_db: Path) -> Tuple[Dict[str, Set[str]], Dict[str, Any]]:
    if not chembl_db.exists():
        raise FileNotFoundError(str(chembl_db))

    chembl_to_zinc: Dict[str, Set[str]] = defaultdict(set)
    metrics = Counter()

    with sqlite3.connect(chembl_db) as conn:
        # release anchor
        rel = conn.execute(
            "SELECT chembl_release FROM chembl_release ORDER BY chembl_release_id DESC LIMIT 1"
        ).fetchone()
        chembl_release = rel[0] if rel else "CHEMBL_UNKNOWN"

        q = """
            SELECT UPPER(TRIM(md.chembl_id)) AS chembl_id, ms.synonyms
            FROM molecule_synonyms ms
            JOIN molecule_dictionary md ON md.molregno = ms.molregno
            WHERE md.chembl_id IS NOT NULL
              AND ms.synonyms IS NOT NULL
              AND UPPER(ms.synonyms) LIKE 'ZINC%'
        """
        for chembl_id, syn in conn.execute(q):
            metrics["zinc_synonym_rows"] += 1
            c = normalize(chembl_id).upper()
            if not valid_chembl(c):
                metrics["zinc_bad_chembl"] += 1
                continue
            zid = parse_zinc_id(syn)
            if not zid:
                metrics["zinc_bad_synonym"] += 1
                continue
            chembl_to_zinc[c].add(zid)

    metrics["zinc_mapped_chembl_ids"] = len(chembl_to_zinc)
    metrics["zinc_total_ids"] = int(sum(len(v) for v in chembl_to_zinc.values()))

    return chembl_to_zinc, {
        "chembl_release": chembl_release,
        **{k: int(v) for k, v in metrics.items()},
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xref", type=Path, required=True)
    ap.add_argument("--chembl-db", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--coverage-report", type=Path, required=True)
    ap.add_argument("--max-rows", type=int, default=None)
    ap.add_argument("--pubchem-chunk-size", type=int, default=80)
    ap.add_argument("--pubchem-timeout", type=int, default=25)
    args = ap.parse_args()

    if not args.xref.exists():
        raise SystemExit(f"[ERROR] missing xref input: {args.xref}")

    created_at = utc_now()
    fetch_date = utc_today()

    xref_rows, xref_stats = load_xref_rows(args.xref, args.max_rows)

    # PubChem 3D metadata (real-time lookup)
    all_pubchem_cids: Set[str] = set()
    for r in xref_rows:
        all_pubchem_cids.update(r.pubchem_cids)

    cid_to_conformer_count, pubchem_report = fetch_pubchem_conformer_counts(
        sorted(all_pubchem_cids),
        chunk_size=max(1, args.pubchem_chunk_size),
        timeout=max(5, args.pubchem_timeout),
    )

    # Choose per-row primary CID (prefer one with 3D conformer)
    primary_cid_by_inchikey: Dict[str, str] = {}
    for r in xref_rows:
        cids = sorted(r.pubchem_cids)
        if not cids:
            primary_cid_by_inchikey[r.inchikey] = ""
            continue
        chosen = ""
        for cid in cids:
            if cid_to_conformer_count.get(cid, 0) > 0:
                chosen = cid
                break
        if not chosen:
            chosen = cids[0]
        primary_cid_by_inchikey[r.inchikey] = chosen

    # PubChem RegistryID -> ZINC IDs (on primary CIDs only)
    primary_cids = sorted({cid for cid in primary_cid_by_inchikey.values() if cid})
    cid_to_pubchem_zinc, pubchem_zinc_report = fetch_pubchem_zinc_ids(
        primary_cids,
        chunk_size=min(40, max(10, args.pubchem_chunk_size // 2 or 40)),
        timeout=max(10, args.pubchem_timeout),
    )

    # ZINC metadata (local mapping from ChEMBL synonyms)
    try:
        chembl_to_zinc, zinc_report = load_zinc_map_from_chembl(args.chembl_db)
        zinc_map_source = "chembl_molecule_synonyms"
    except FileNotFoundError:
        chembl_to_zinc = {}
        zinc_report = {"chembl_release": "CHEMBL_UNKNOWN", "warning": f"missing chembl db: {args.chembl_db}"}
        zinc_map_source = "missing"

    header = [
        "inchikey",
        "xref_version",
        "pubchem_cid",
        "pubchem_3d_available",
        "pubchem_conformer_count",
        "pubchem_formats",
        "pubchem_download_url",
        "pubchem_source",
        "pubchem_source_type",
        "zinc_id",
        "zinc_3d_available",
        "zinc_conformer_count",
        "zinc_formats",
        "zinc_download_url",
        "zinc_source",
        "zinc_source_type",
        "source_providers",
        "structure_formats",
        "availability_status",
        "source_types",
        "fetch_date",
        "source_version",
    ]

    out_rows: List[Dict[str, str]] = []
    counts = Counter()
    source_type_counter = Counter()

    for row in xref_rows:
        inchikey = row.inchikey
        cids = sorted(row.pubchem_cids)
        chembl_ids = sorted(row.chembl_ids)

        # PubChem aggregation
        pubchem_counts = [cid_to_conformer_count.get(cid, 0) for cid in cids]
        pubchem_has_3d = any(c > 0 for c in pubchem_counts)
        pubchem_count = max(pubchem_counts) if pubchem_counts else 0

        first_pubchem_cid = primary_cid_by_inchikey.get(inchikey, "")

        pubchem_formats = "SDF;JSON" if pubchem_has_3d else "NA"
        pubchem_url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{first_pubchem_cid}/record/SDF/?record_type=3d"
            if first_pubchem_cid
            else ""
        )

        # ZINC aggregation (from chembl->zinc)
        zinc_ids: Set[str] = set()
        for c in chembl_ids:
            zinc_ids.update(chembl_to_zinc.get(c, set()))
        if first_pubchem_cid:
            zinc_ids.update(cid_to_pubchem_zinc.get(first_pubchem_cid, set()))

        zinc_has_3d = bool(zinc_ids)
        zinc_id_str = join_values(zinc_ids)
        first_zinc = sorted(zinc_ids)[0] if zinc_ids else ""
        zinc_formats = "MOL2;SDF;PDBQT" if zinc_has_3d else "NA"
        zinc_url = (
            f"https://zinc20.docking.org/substances/{first_zinc}/"
            if first_zinc
            else f"https://zinc20.docking.org/substances/search/?q={inchikey}"
        )

        providers: List[str] = []
        formats: Set[str] = set()
        row_types: Set[str] = set()

        if pubchem_has_3d:
            providers.append("PubChem")
            formats.update(["SDF", "JSON"])
            row_types.add("computational")
            counts["rows_with_pubchem_3d"] += 1

        if zinc_has_3d:
            providers.append("ZINC")
            formats.update(["MOL2", "SDF", "PDBQT"])
            row_types.add("computational")
            counts["rows_with_zinc_3d"] += 1

        if providers:
            availability = "available"
            counts["rows_with_any_3d"] += 1
        else:
            availability = "unavailable"
            providers = ["none"]
            formats = {"none"}
            row_types = {"unknown"}
            counts["rows_without_3d"] += 1

        source_types = join_values(row_types)
        source_type_counter[source_types] += 1

        out_rows.append(
            {
                "inchikey": inchikey,
                "xref_version": "v2" if "_v2" in args.xref.name else "v1",
                "pubchem_cid": join_values(cids),
                "pubchem_3d_available": "1" if pubchem_has_3d else "0",
                "pubchem_conformer_count": str(pubchem_count),
                "pubchem_formats": pubchem_formats,
                "pubchem_download_url": pubchem_url,
                "pubchem_source": "PubChem",
                "pubchem_source_type": "computational",
                "zinc_id": zinc_id_str,
                "zinc_3d_available": "1" if zinc_has_3d else "0",
                "zinc_conformer_count": "1" if zinc_has_3d else "0",
                "zinc_formats": zinc_formats,
                "zinc_download_url": zinc_url,
                "zinc_source": "ZINC20",
                "zinc_source_type": "computational",
                "source_providers": join_values(providers),
                "structure_formats": join_values(formats),
                "availability_status": availability,
                "source_types": source_types,
                "fetch_date": fetch_date,
                "source_version": (
                    f"PubChem:PUGREST_ConformerCount3D;"
                    f"PubChem:PUGREST_RegistryID;"
                    f"ZINC:{zinc_map_source};"
                    f"ChEMBL:{zinc_report.get('chembl_release', 'CHEMBL_UNKNOWN')}"
                ),
            }
        )

    rows_written = write_tsv(args.out, out_rows, header)

    # coverage report
    non_empty_rates = {}
    for col in ["source_providers", "structure_formats", "availability_status", "source_types"]:
        non_empty = sum(1 for r in out_rows if normalize(r.get(col, "")) != "")
        non_empty_rates[col] = {
            "non_empty": non_empty,
            "total": rows_written,
            "rate": (non_empty / rows_written) if rows_written else 0.0,
        }

    coverage_report = {
        "name": "molecule_3d_registry_v1.coverage",
        "created_at": created_at,
        "xref_input": str(args.xref),
        "rows": rows_written,
        "row_metrics": {
            "rows_with_pubchem_3d": int(counts["rows_with_pubchem_3d"]),
            "rows_with_zinc_3d": int(counts["rows_with_zinc_3d"]),
            "rows_with_any_3d": int(counts["rows_with_any_3d"]),
            "rows_without_3d": int(counts["rows_without_3d"]),
            "pubchem_3d_rate": (counts["rows_with_pubchem_3d"] / rows_written) if rows_written else 0.0,
            "zinc_3d_rate": (counts["rows_with_zinc_3d"] / rows_written) if rows_written else 0.0,
            "any_3d_rate": (counts["rows_with_any_3d"] / rows_written) if rows_written else 0.0,
        },
        "non_empty_rates": non_empty_rates,
        "source_type_distribution": dict(source_type_counter),
    }

    args.coverage_report.parent.mkdir(parents=True, exist_ok=True)
    args.coverage_report.write_text(json.dumps(coverage_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    build_report: Dict[str, Any] = {
        "name": "molecule_3d_registry_v1.build",
        "created_at": created_at,
        "sample_mode": args.max_rows is not None,
        "max_rows": args.max_rows,
        "inputs": {
            "xref": str(args.xref),
            "chembl_db": str(args.chembl_db),
            "chembl_db_exists": args.chembl_db.exists(),
        },
        "output": str(args.out),
        "metrics": {
            "rows_written": rows_written,
            **xref_stats,
            **coverage_report["row_metrics"],
        },
        "pubchem_lookup": pubchem_report,
        "pubchem_zinc_lookup": pubchem_zinc_report,
        "zinc_mapping": zinc_report,
        "source_type_distribution": dict(source_type_counter),
        "coverage_report": str(args.coverage_report),
    }

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(build_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"[OK] build -> {args.out} (rows={rows_written})")
    print(
        "[OK] 3D coverage: pubchem={} zinc={} any={}".format(
            counts["rows_with_pubchem_3d"],
            counts["rows_with_zinc_3d"],
            counts["rows_with_any_3d"],
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
