#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
CHEMBL_RE = re.compile(r"^CHEMBL\d+$")
DIGIT_RE = re.compile(r"\d+")


# ----------------------------- common utils -----------------------------

def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def utc_today() -> str:
    return datetime.now(timezone.utc).date().isoformat()


def normalize(x: str) -> str:
    return (x or "").strip()


def canonical_term(x: str) -> str:
    return " ".join(normalize(x).split()).lower()


def valid_inchikey(x: str) -> bool:
    return bool(INCHIKEY_RE.match(normalize(x).upper()))


def valid_chembl(x: str) -> bool:
    return bool(CHEMBL_RE.match(normalize(x).upper()))


def parse_pubchem_cid(raw: str) -> Optional[str]:
    text = normalize(raw)
    if not text:
        return None
    if text.isdigit():
        return text
    m = DIGIT_RE.search(text)
    if m:
        return m.group(0)
    return None


def split_multi(x: str) -> List[str]:
    text = normalize(x)
    if not text:
        return []
    out: List[str] = []
    for token in text.split(";"):
        t = token.strip()
        if t:
            out.append(t)
    return out


def split_synonyms(x: str) -> List[str]:
    return split_multi(x)


def join_values(values: Iterable[str]) -> str:
    return ";".join(sorted({normalize(v) for v in values if normalize(v)}))


def _read_tsv(path: Path, max_rows: Optional[int] = None) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        out: List[Dict[str, str]] = []
        for i, row in enumerate(r, start=1):
            out.append(row)
            if max_rows is not None and i >= max_rows:
                break
    return list(r.fieldnames), out


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


# ----------------------------- data model -----------------------------

@dataclass
class XrefRec:
    inchikey: str
    drugbank_ids: Set[str] = field(default_factory=set)
    chembl_ids: Set[str] = field(default_factory=set)
    pubchem_cids: Set[str] = field(default_factory=set)
    sources: Set[str] = field(default_factory=set)
    match_strategies: Set[str] = field(default_factory=set)


@dataclass
class BackfillRecord:
    drugbank_id: str
    strategy: str
    confidence: str
    matched_term: str
    chembl_id: str
    inchikey: str


@dataclass
class BuildCtx:
    records: Dict[str, XrefRec] = field(default_factory=dict)
    source_versions: Set[str] = field(default_factory=set)
    warnings: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    metrics: Counter = field(default_factory=Counter)

    def ensure_record(self, inchikey: str) -> XrefRec:
        rec = self.records.get(inchikey)
        if rec is None:
            rec = XrefRec(inchikey=inchikey)
            self.records[inchikey] = rec
        return rec


# ----------------------------- sqlite helpers -----------------------------

def _table_exists(conn: sqlite3.Connection, name: str) -> bool:
    cur = conn.execute("SELECT 1 FROM sqlite_master WHERE type IN ('table','view') AND name=?", (name,))
    return cur.fetchone() is not None


def _columns(conn: sqlite3.Connection, table: str) -> Set[str]:
    return {str(r[1]) for r in conn.execute(f"PRAGMA table_info({table})")}


def extract_chembl_from_m1(m1_db: Path, ctx: BuildCtx) -> Dict[str, Set[str]]:
    if not m1_db.exists():
        ctx.warnings.append(f"optional input missing: {m1_db} (m1 chembl enrichment skipped)")
        return {}

    candidates: List[Tuple[str, str]] = [
        ("molecule_entity_strict_chembl36", "rep_chembl_id"),
        ("molecule_entity_chembl36", "rep_chembl_id"),
        ("molecule_idmap_chembl36", "chembl_id"),
        ("molecule_raw_chembl36", "chembl_id"),
    ]

    out: Dict[str, Set[str]] = defaultdict(set)
    with sqlite3.connect(m1_db) as conn:
        chosen: Optional[Tuple[str, str]] = None
        for table, chembl_col in candidates:
            if not _table_exists(conn, table):
                continue
            cols = _columns(conn, table)
            if "inchikey" in cols and chembl_col in cols:
                chosen = (table, chembl_col)
                break

        if chosen is None:
            ctx.warnings.append(f"m1 sqlite has no usable inchikey/chembl table: {m1_db}")
            return {}

        table, chembl_col = chosen
        q = (
            f"SELECT UPPER(TRIM(inchikey)) AS inchikey, UPPER(TRIM({chembl_col})) AS chembl_id "
            f"FROM {table} WHERE inchikey IS NOT NULL AND {chembl_col} IS NOT NULL"
        )
        for ik, chembl in conn.execute(q):
            ik_s = normalize(ik).upper()
            c_s = normalize(chembl).upper()
            if valid_inchikey(ik_s) and valid_chembl(c_s):
                out[ik_s].add(c_s)

    ctx.metrics["m1_chembl_inchikey_count"] = len(out)
    ctx.metrics["m1_chembl_pairs"] = int(sum(len(v) for v in out.values()))
    ctx.notes.append("chembl enrichment from molecules_m1 sqlite")
    ctx.source_versions.add("ChEMBL:molecules_m1.sqlite")
    return out


def extract_chembl_by_inchikey_from_chembl_db(chembl_db: Path, inchikeys: Sequence[str], ctx: BuildCtx) -> Dict[str, Set[str]]:
    if not chembl_db.exists():
        ctx.warnings.append(f"optional input missing: {chembl_db} (chembl-by-inchikey enrichment skipped)")
        return {}
    if not inchikeys:
        return {}

    out: Dict[str, Set[str]] = defaultdict(set)
    with sqlite3.connect(chembl_db) as conn:
        if not _table_exists(conn, "compound_structures") or not _table_exists(conn, "molecule_dictionary"):
            ctx.warnings.append("chembl db missing required tables (compound_structures/molecule_dictionary)")
            return {}

        conn.execute("DROP TABLE IF EXISTS temp.tmp_ik")
        conn.execute("CREATE TEMP TABLE tmp_ik (inchikey TEXT PRIMARY KEY)")
        conn.executemany("INSERT OR IGNORE INTO tmp_ik(inchikey) VALUES (?)", [(x,) for x in set(inchikeys)])

        q = """
            SELECT t.inchikey, UPPER(TRIM(md.chembl_id)) AS chembl_id
            FROM tmp_ik t
            JOIN compound_structures cs ON cs.standard_inchi_key = t.inchikey
            JOIN molecule_dictionary md ON md.molregno = cs.molregno
            WHERE md.chembl_id IS NOT NULL
        """
        for ik, chembl in conn.execute(q):
            ik_s = normalize(ik).upper()
            c_s = normalize(chembl).upper()
            if valid_inchikey(ik_s) and valid_chembl(c_s):
                out[ik_s].add(c_s)

    ctx.metrics["chembl_db_inchikey_match_count"] = len(out)
    ctx.notes.append("chembl enrichment from chembl_36.db compound_structures")
    ctx.source_versions.add("ChEMBL:chembl_36.db")
    return out


def extract_pubchem_for_chembl_ids(chembl_db: Path, chembl_ids: Sequence[str], ctx: BuildCtx) -> Dict[str, Set[str]]:
    if not chembl_db.exists() or not chembl_ids:
        return {}

    out: Dict[str, Set[str]] = defaultdict(set)
    with sqlite3.connect(chembl_db) as conn:
        required = {"molecule_dictionary", "compound_records", "source"}
        if any(not _table_exists(conn, t) for t in required):
            ctx.warnings.append("chembl db missing tables for pubchem bridge")
            return {}

        conn.execute("DROP TABLE IF EXISTS temp.tmp_chembl")
        conn.execute("CREATE TEMP TABLE tmp_chembl (chembl_id TEXT PRIMARY KEY)")
        conn.executemany("INSERT OR IGNORE INTO tmp_chembl(chembl_id) VALUES (?)", [(x.upper(),) for x in set(chembl_ids)])

        source_cols = _columns(conn, "source")
        desc_expr = "LOWER(COALESCE(s.src_description,''))" if "src_description" in source_cols else "''"
        short_expr = "LOWER(COALESCE(s.src_short_name,''))" if "src_short_name" in source_cols else "''"
        url_expr = "LOWER(COALESCE(s.src_url,''))" if "src_url" in source_cols else "''"

        q = f"""
            SELECT UPPER(TRIM(md.chembl_id)) AS chembl_id, TRIM(cr.src_compound_id) AS src_compound_id
            FROM tmp_chembl t
            JOIN molecule_dictionary md ON UPPER(md.chembl_id) = t.chembl_id
            JOIN compound_records cr ON md.molregno = cr.molregno
            JOIN source s ON cr.src_id = s.src_id
            WHERE cr.src_compound_id IS NOT NULL
              AND (
                    {desc_expr} LIKE '%pubchem%'
                 OR {short_expr} LIKE '%pubchem%'
                 OR {url_expr} LIKE '%pubchem%'
              )
        """
        for chembl, raw_cid in conn.execute(q):
            c = normalize(chembl).upper()
            cid = parse_pubchem_cid(raw_cid)
            if valid_chembl(c) and cid:
                out[c].add(cid)

    if out:
        ctx.source_versions.add("PubChemCID:from_chembl_36_compound_records")
    ctx.metrics["pubchem_chembl_map_count"] = len(out)
    return out


def _query_name_hits_pref_name(chembl_db: Path, terms: Sequence[str], ctx: BuildCtx) -> Dict[str, Set[Tuple[str, str]]]:
    out: Dict[str, Set[Tuple[str, str]]] = defaultdict(set)
    if not chembl_db.exists() or not terms:
        return out

    with sqlite3.connect(chembl_db) as conn:
        conn.execute("DROP TABLE IF EXISTS temp.tmp_terms_name")
        conn.execute("CREATE TEMP TABLE tmp_terms_name (term TEXT PRIMARY KEY)")
        conn.executemany("INSERT OR IGNORE INTO tmp_terms_name(term) VALUES (?)", [(t,) for t in set(terms)])

        q = """
            SELECT t.term, UPPER(TRIM(md.chembl_id)) AS chembl_id, UPPER(TRIM(cs.standard_inchi_key)) AS inchikey
            FROM tmp_terms_name t
            JOIN molecule_dictionary md ON md.pref_name = t.term COLLATE NOCASE
            JOIN compound_structures cs ON cs.molregno = md.molregno
            WHERE md.chembl_id IS NOT NULL AND cs.standard_inchi_key IS NOT NULL
        """
        for term, chembl_id, inchikey in conn.execute(q):
            if valid_chembl(chembl_id) and valid_inchikey(inchikey):
                out[canonical_term(term)].add((chembl_id, inchikey))

    return out


def _query_name_hits_synonym(chembl_db: Path, terms: Sequence[str], ctx: BuildCtx) -> Dict[str, Set[Tuple[str, str]]]:
    out: Dict[str, Set[Tuple[str, str]]] = defaultdict(set)
    if not chembl_db.exists() or not terms:
        return out

    with sqlite3.connect(chembl_db) as conn:
        conn.execute("DROP TABLE IF EXISTS temp.tmp_terms_syn")
        conn.execute("CREATE TEMP TABLE tmp_terms_syn (term TEXT PRIMARY KEY)")
        conn.executemany("INSERT OR IGNORE INTO tmp_terms_syn(term) VALUES (?)", [(t,) for t in set(terms)])

        q = """
            SELECT t.term, UPPER(TRIM(md.chembl_id)) AS chembl_id, UPPER(TRIM(cs.standard_inchi_key)) AS inchikey
            FROM tmp_terms_syn t
            JOIN molecule_synonyms ms ON ms.synonyms = t.term COLLATE NOCASE
            JOIN molecule_dictionary md ON md.molregno = ms.molregno
            JOIN compound_structures cs ON cs.molregno = md.molregno
            WHERE md.chembl_id IS NOT NULL AND cs.standard_inchi_key IS NOT NULL
        """
        for term, chembl_id, inchikey in conn.execute(q):
            if valid_chembl(chembl_id) and valid_inchikey(inchikey):
                out[canonical_term(term)].add((chembl_id, inchikey))

    return out


# ----------------------------- business logic -----------------------------

def ingest_v1_rows(rows: List[Dict[str, str]], ctx: BuildCtx) -> Dict[str, Any]:
    v1_metrics: Dict[str, Any] = {}
    v1_rows = 0
    v1_chembl_rows = 0
    v1_pubchem_rows = 0
    v1_dbid_all: Set[str] = set()
    v1_dbid_chembl: Set[str] = set()
    v1_dbid_pubchem: Set[str] = set()

    for row in rows:
        v1_rows += 1
        ik = normalize(row.get("inchikey", "")).upper()
        if not valid_inchikey(ik):
            ctx.metrics["v1_bad_inchikey_rows"] += 1
            continue

        rec = ctx.ensure_record(ik)
        for dbid in split_multi(row.get("drugbank_id", "")):
            rec.drugbank_ids.add(dbid)
        for chembl in split_multi(row.get("chembl_id", "")):
            c = chembl.upper()
            if valid_chembl(c):
                rec.chembl_ids.add(c)
        for cid in split_multi(row.get("pubchem_cid", "")):
            pc = parse_pubchem_cid(cid)
            if pc:
                rec.pubchem_cids.add(pc)
        for src in split_multi(row.get("xref_source", "")):
            rec.sources.add(src)
        rec.match_strategies.add("seed_v1")

        src_ver = normalize(row.get("source_version", ""))
        if src_ver:
            ctx.source_versions.add(f"v1:{src_ver}")

    for rec in ctx.records.values():
        if not rec.drugbank_ids:
            continue
        v1_dbid_all.update(rec.drugbank_ids)
        if rec.chembl_ids:
            v1_chembl_rows += 1
            v1_dbid_chembl.update(rec.drugbank_ids)
        if rec.pubchem_cids:
            v1_pubchem_rows += 1
            v1_dbid_pubchem.update(rec.drugbank_ids)

    v1_metrics.update(
        {
            "rows": v1_rows,
            "rows_with_chembl": v1_chembl_rows,
            "rows_with_pubchem": v1_pubchem_rows,
            "chembl_row_rate": (v1_chembl_rows / v1_rows) if v1_rows else 0.0,
            "pubchem_row_rate": (v1_pubchem_rows / v1_rows) if v1_rows else 0.0,
            "drugbank_ids": len(v1_dbid_all),
            "drugbank_ids_with_chembl": len(v1_dbid_chembl),
            "drugbank_ids_with_pubchem": len(v1_dbid_pubchem),
            "drugbank_chembl_rate": (len(v1_dbid_chembl) / len(v1_dbid_all)) if v1_dbid_all else 0.0,
            "drugbank_pubchem_rate": (len(v1_dbid_pubchem) / len(v1_dbid_all)) if v1_dbid_all else 0.0,
        }
    )
    return v1_metrics


def build_drugbank_maps(
    drug_master_rows: List[Dict[str, str]],
    drug_xref_rows: List[Dict[str, str]],
    ctx: BuildCtx,
) -> Tuple[Dict[str, str], Dict[str, Dict[str, Any]], Set[str]]:
    dbid_to_inchikey: Dict[str, str] = {}
    master_info: Dict[str, Dict[str, Any]] = {}
    all_dbids: Set[str] = set()

    for row in drug_master_rows:
        dbid = normalize(row.get("drugbank_id", ""))
        if not dbid:
            continue
        all_dbids.add(dbid)
        ik = normalize(row.get("inchikey", "")).upper()
        name = normalize(row.get("name", ""))
        synonyms = split_synonyms(row.get("synonyms", ""))
        master_info[dbid] = {
            "name": name,
            "synonyms": synonyms,
            "source_version": normalize(row.get("source_version", "")),
        }
        if valid_inchikey(ik):
            dbid_to_inchikey[dbid] = ik

        sv = normalize(row.get("source_version", ""))
        if sv:
            ctx.source_versions.add(f"DrugBank:{sv}")

    for row in drug_xref_rows:
        dbid = normalize(row.get("drugbank_id", ""))
        if not dbid:
            continue
        all_dbids.add(dbid)
        ik = normalize(row.get("inchikey", "")).upper()
        if valid_inchikey(ik):
            dbid_to_inchikey[dbid] = ik

    return dbid_to_inchikey, master_info, all_dbids


def add_seed_drugbank_rows(dbid_to_inchikey: Dict[str, str], ctx: BuildCtx) -> None:
    for dbid, ik in dbid_to_inchikey.items():
        rec = ctx.ensure_record(ik)
        rec.drugbank_ids.add(dbid)
        rec.sources.add("drug_master_v1")
        rec.sources.add("drug_xref_molecules_v1")
        rec.match_strategies.add("drugbank_inchikey_seed")


def apply_backfill(
    dbids_missing_inchikey: Set[str],
    master_info: Dict[str, Dict[str, Any]],
    chembl_db: Path,
    ctx: BuildCtx,
) -> Tuple[List[BackfillRecord], List[Dict[str, Any]]]:
    accepted: List[BackfillRecord] = []
    unresolved: List[Dict[str, Any]] = []

    if not chembl_db.exists():
        for dbid in sorted(dbids_missing_inchikey):
            info = master_info.get(dbid, {})
            unresolved.append(
                {
                    "drugbank_id": dbid,
                    "reason": "chembl_db_missing",
                    "name": info.get("name", ""),
                    "candidate_count": 0,
                }
            )
        return accepted, unresolved

    # Stage 1: exact (case-insensitive) DrugBank name -> ChEMBL pref_name
    stage1_terms = sorted({normalize(master_info.get(dbid, {}).get("name", "")) for dbid in dbids_missing_inchikey if normalize(master_info.get(dbid, {}).get("name", ""))})
    pref_hits = _query_name_hits_pref_name(chembl_db, stage1_terms, ctx)

    stage1_resolved: Set[str] = set()
    stage2_pool: Set[str] = set()

    for dbid in sorted(dbids_missing_inchikey):
        info = master_info.get(dbid, {})
        name = normalize(info.get("name", ""))
        hits = pref_hits.get(canonical_term(name), set()) if name else set()
        if len(hits) == 1:
            chembl_id, ik = next(iter(hits))
            accepted.append(
                BackfillRecord(
                    drugbank_id=dbid,
                    strategy="backfill_stage1_name_exact_unique",
                    confidence="high",
                    matched_term=name,
                    chembl_id=chembl_id,
                    inchikey=ik,
                )
            )
            stage1_resolved.add(dbid)
        elif len(hits) > 1:
            unresolved.append(
                {
                    "drugbank_id": dbid,
                    "reason": "stage1_name_ambiguous",
                    "name": name,
                    "candidate_count": len(hits),
                    "candidate_sample": [{"chembl_id": c, "inchikey": i} for c, i in sorted(hits)[:5]],
                }
            )
        else:
            stage2_pool.add(dbid)

    # Stage 2: exact (case-insensitive) DrugBank synonym -> ChEMBL synonym
    stage2_terms: Set[str] = set()
    dbid_to_syn_terms: Dict[str, List[str]] = {}
    for dbid in sorted(stage2_pool):
        syns = [normalize(x) for x in master_info.get(dbid, {}).get("synonyms", []) if normalize(x)]
        if syns:
            dbid_to_syn_terms[dbid] = syns
            stage2_terms.update(syns)

    syn_hits = _query_name_hits_synonym(chembl_db, sorted(stage2_terms), ctx)

    for dbid in sorted(stage2_pool):
        syns = dbid_to_syn_terms.get(dbid, [])
        hit_set: Set[Tuple[str, str]] = set()
        matched_terms: Set[str] = set()
        for syn in syns:
            h = syn_hits.get(canonical_term(syn), set())
            if h:
                hit_set.update(h)
                matched_terms.add(syn)

        if len(hit_set) == 1:
            chembl_id, ik = next(iter(hit_set))
            matched_term = sorted(matched_terms)[0] if matched_terms else ""
            accepted.append(
                BackfillRecord(
                    drugbank_id=dbid,
                    strategy="backfill_stage2_synonym_exact_unique",
                    confidence="medium",
                    matched_term=matched_term,
                    chembl_id=chembl_id,
                    inchikey=ik,
                )
            )
        elif len(hit_set) > 1:
            unresolved.append(
                {
                    "drugbank_id": dbid,
                    "reason": "stage2_synonym_ambiguous",
                    "name": master_info.get(dbid, {}).get("name", ""),
                    "candidate_count": len(hit_set),
                    "candidate_sample": [{"chembl_id": c, "inchikey": i} for c, i in sorted(hit_set)[:5]],
                    "matched_terms_sample": sorted(matched_terms)[:5],
                }
            )
        else:
            unresolved.append(
                {
                    "drugbank_id": dbid,
                    "reason": "no_match_after_two_stage_backfill",
                    "name": master_info.get(dbid, {}).get("name", ""),
                    "candidate_count": 0,
                }
            )

    return accepted, unresolved


def infer_confidence(rec: XrefRec) -> str:
    if "backfill_stage2_synonym_exact_unique" in rec.match_strategies:
        return "medium"
    if rec.drugbank_ids and rec.chembl_ids:
        return "high"
    if rec.drugbank_ids:
        return "medium"
    if rec.chembl_ids:
        return "medium"
    return "low"


def summarize_v2(records: Dict[str, XrefRec]) -> Dict[str, Any]:
    rows = 0
    rows_chembl = 0
    rows_pubchem = 0
    dbid_all: Set[str] = set()
    dbid_chembl: Set[str] = set()
    dbid_pubchem: Set[str] = set()

    for rec in records.values():
        if not rec.drugbank_ids:
            continue
        rows += 1
        dbid_all.update(rec.drugbank_ids)
        if rec.chembl_ids:
            rows_chembl += 1
            dbid_chembl.update(rec.drugbank_ids)
        if rec.pubchem_cids:
            rows_pubchem += 1
            dbid_pubchem.update(rec.drugbank_ids)

    return {
        "rows": rows,
        "rows_with_chembl": rows_chembl,
        "rows_with_pubchem": rows_pubchem,
        "chembl_row_rate": (rows_chembl / rows) if rows else 0.0,
        "pubchem_row_rate": (rows_pubchem / rows) if rows else 0.0,
        "drugbank_ids": len(dbid_all),
        "drugbank_ids_with_chembl": len(dbid_chembl),
        "drugbank_ids_with_pubchem": len(dbid_pubchem),
        "drugbank_chembl_rate": (len(dbid_chembl) / len(dbid_all)) if dbid_all else 0.0,
        "drugbank_pubchem_rate": (len(dbid_pubchem) / len(dbid_all)) if dbid_all else 0.0,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--v1-core", type=Path, required=True)
    ap.add_argument("--drug-master", type=Path, required=True)
    ap.add_argument("--drug-xref", type=Path, required=True)
    ap.add_argument("--m1-db", type=Path, required=True)
    ap.add_argument("--chembl-db", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--backfill-audit", type=Path, required=True)
    ap.add_argument("--missing-audit", type=Path, required=True)
    ap.add_argument("--max-rows", type=int, default=None)
    ap.add_argument("--missing-sample-n", type=int, default=200)
    args = ap.parse_args()

    for p in [args.v1_core, args.drug_master, args.drug_xref]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing required input: {p}")

    created_at = utc_now()
    fetch_date = utc_today()

    ctx = BuildCtx()

    _, v1_rows = _read_tsv(args.v1_core, args.max_rows)
    _, master_rows = _read_tsv(args.drug_master, args.max_rows)
    _, drug_xref_rows = _read_tsv(args.drug_xref, args.max_rows)

    v1_metrics = ingest_v1_rows(v1_rows, ctx)

    dbid_to_inchikey, master_info, all_dbids = build_drugbank_maps(master_rows, drug_xref_rows, ctx)
    add_seed_drugbank_rows(dbid_to_inchikey, ctx)

    # Chembl enrichment (stage A): from m1 sqlite
    m1_map = extract_chembl_from_m1(args.m1_db, ctx)
    chembl_added_from_m1 = 0
    for ik, rec in ctx.records.items():
        before = set(rec.chembl_ids)
        rec.chembl_ids.update(m1_map.get(ik, set()))
        if rec.chembl_ids != before:
            rec.sources.add("molecules_m1.sqlite")
            rec.match_strategies.add("m1_inchikey_to_chembl")
            chembl_added_from_m1 += len(rec.chembl_ids - before)

    # Chembl enrichment (stage B): direct inchikey -> chembl from chembl DB
    target_iks = sorted([ik for ik, rec in ctx.records.items() if rec.drugbank_ids and not rec.chembl_ids])
    chembl_by_ik_db = extract_chembl_by_inchikey_from_chembl_db(args.chembl_db, target_iks, ctx)
    chembl_added_from_chembl_db = 0
    for ik, chembls in chembl_by_ik_db.items():
        rec = ctx.records.get(ik)
        if rec is None:
            continue
        before = set(rec.chembl_ids)
        rec.chembl_ids.update(chembls)
        if rec.chembl_ids != before:
            rec.sources.add("chembl_36.db")
            rec.match_strategies.add("chembl_db_inchikey_exact")
            chembl_added_from_chembl_db += len(rec.chembl_ids - before)

    # Two-stage backfill for DrugBank entries missing InChIKey
    missing_dbids = {dbid for dbid in all_dbids if dbid not in dbid_to_inchikey}
    backfill_accepted, backfill_unresolved = apply_backfill(missing_dbids, master_info, args.chembl_db, ctx)

    backfill_strategy_counter: Counter[str] = Counter()
    for rec_audit in backfill_accepted:
        rec = ctx.ensure_record(rec_audit.inchikey)
        rec.drugbank_ids.add(rec_audit.drugbank_id)
        rec.chembl_ids.add(rec_audit.chembl_id)
        rec.sources.update({"drug_master_v1", "chembl_36.db"})
        rec.match_strategies.add(rec_audit.strategy)
        backfill_strategy_counter[rec_audit.strategy] += 1

    # PubChem bridge from all chembl IDs after backfill
    all_chembl_ids: Set[str] = set()
    for rec in ctx.records.values():
        all_chembl_ids.update(rec.chembl_ids)
    pubchem_map = extract_pubchem_for_chembl_ids(args.chembl_db, sorted(all_chembl_ids), ctx)

    pubchem_added = 0
    for rec in ctx.records.values():
        if not rec.chembl_ids:
            continue
        before = set(rec.pubchem_cids)
        for chembl_id in rec.chembl_ids:
            rec.pubchem_cids.update(pubchem_map.get(chembl_id, set()))
        if rec.pubchem_cids != before:
            rec.sources.add("chembl_36.db")
            rec.match_strategies.add("chembl_to_pubchem_bridge")
            pubchem_added += len(rec.pubchem_cids - before)

    # Output rows (DrugBank-anchored)
    source_version = join_values(ctx.source_versions) or "molecule_xref_enrichment_v2|source_unknown"
    out_rows: List[Dict[str, str]] = []
    duplicate_inchikey = 0
    seen_ik: Set[str] = set()

    for ik in sorted(ctx.records.keys()):
        rec = ctx.records[ik]
        if not rec.drugbank_ids:
            continue
        if ik in seen_ik:
            duplicate_inchikey += 1
        seen_ik.add(ik)

        out_rows.append(
            {
                "inchikey": ik,
                "chembl_id": join_values(rec.chembl_ids),
                "drugbank_id": join_values(rec.drugbank_ids),
                "pubchem_cid": join_values(rec.pubchem_cids),
                "match_strategy": join_values(rec.match_strategies),
                "confidence": infer_confidence(rec),
                "xref_source": join_values(rec.sources),
                "fetch_date": fetch_date,
                "source_version": source_version,
            }
        )

    header = [
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
    rows_written = write_tsv(args.out, out_rows, header)

    v2_metrics = summarize_v2(ctx.records)

    # Audit outputs
    backfill_audit = {
        "name": "molecule_xref_core_v2.backfill_audit",
        "created_at": created_at,
        "accepted_count": len(backfill_accepted),
        "strategy_counts": dict(backfill_strategy_counter),
        "records": [
            {
                "drugbank_id": x.drugbank_id,
                "match_strategy": x.strategy,
                "confidence": x.confidence,
                "matched_term": x.matched_term,
                "chembl_id": x.chembl_id,
                "inchikey": x.inchikey,
            }
            for x in sorted(backfill_accepted, key=lambda z: z.drugbank_id)
        ],
    }

    unresolved_counter = Counter([x.get("reason", "unknown") for x in backfill_unresolved])
    missing_audit = {
        "name": "molecule_xref_core_v2.missing_audit",
        "created_at": created_at,
        "total": len(backfill_unresolved),
        "reason_counts": dict(unresolved_counter),
        "sample_top_n": backfill_unresolved[: args.missing_sample_n],
    }

    args.backfill_audit.parent.mkdir(parents=True, exist_ok=True)
    args.backfill_audit.write_text(json.dumps(backfill_audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    args.missing_audit.parent.mkdir(parents=True, exist_ok=True)
    args.missing_audit.write_text(json.dumps(missing_audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    report: Dict[str, Any] = {
        "name": "molecule_xref_core_v2.build",
        "created_at": created_at,
        "sample_mode": args.max_rows is not None,
        "max_rows": args.max_rows,
        "inputs": {
            "v1_core": str(args.v1_core),
            "drug_master": str(args.drug_master),
            "drug_xref": str(args.drug_xref),
            "m1_db": str(args.m1_db),
            "m1_db_exists": args.m1_db.exists(),
            "chembl_db": str(args.chembl_db),
            "chembl_db_exists": args.chembl_db.exists(),
        },
        "output": str(args.out),
        "reports": {
            "backfill_audit": str(args.backfill_audit),
            "missing_audit": str(args.missing_audit),
        },
        "metrics": {
            "rows_written": rows_written,
            "duplicate_inchikey": duplicate_inchikey,
            "inchikey_unique_rate": 1.0 if duplicate_inchikey == 0 else max(0.0, 1.0 - duplicate_inchikey / max(rows_written, 1)),
            "v1_rows": v1_metrics["rows"],
            "v2_rows": v2_metrics["rows"],
            "chembl_rows_v1": v1_metrics["rows_with_chembl"],
            "chembl_rows_v2": v2_metrics["rows_with_chembl"],
            "chembl_rows_delta_abs": v2_metrics["rows_with_chembl"] - v1_metrics["rows_with_chembl"],
            "chembl_rate_v1": v1_metrics["chembl_row_rate"],
            "chembl_rate_v2": v2_metrics["chembl_row_rate"],
            "chembl_rate_delta_abs": v2_metrics["chembl_row_rate"] - v1_metrics["chembl_row_rate"],
            "pubchem_rows_v1": v1_metrics["rows_with_pubchem"],
            "pubchem_rows_v2": v2_metrics["rows_with_pubchem"],
            "pubchem_rows_delta_abs": v2_metrics["rows_with_pubchem"] - v1_metrics["rows_with_pubchem"],
            "pubchem_rate_v1": v1_metrics["pubchem_row_rate"],
            "pubchem_rate_v2": v2_metrics["pubchem_row_rate"],
            "pubchem_rate_delta_abs": v2_metrics["pubchem_row_rate"] - v1_metrics["pubchem_row_rate"],
            "drugbank_ids_v1": v1_metrics["drugbank_ids"],
            "drugbank_ids_v2": v2_metrics["drugbank_ids"],
            "drugbank_missing_inchikey_input": len(missing_dbids),
            "backfill_accepted": len(backfill_accepted),
            "backfill_stage1_accepted": int(backfill_strategy_counter.get("backfill_stage1_name_exact_unique", 0)),
            "backfill_stage2_accepted": int(backfill_strategy_counter.get("backfill_stage2_synonym_exact_unique", 0)),
            "backfill_unresolved": len(backfill_unresolved),
            "chembl_ids_added_from_m1": chembl_added_from_m1,
            "chembl_ids_added_from_chembl_db_inchikey": chembl_added_from_chembl_db,
            "pubchem_ids_added_from_bridge": pubchem_added,
            **{k: int(v) for k, v in ctx.metrics.items()},
        },
        "warnings": ctx.warnings,
        "notes": ctx.notes,
        "source_version": source_version,
    }

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"[OK] build -> {args.out} (rows={rows_written})")
    print(
        "[OK] chembl rows v1={} v2={} delta={} | pubchem rows v1={} v2={} delta={}".format(
            v1_metrics["rows_with_chembl"],
            v2_metrics["rows_with_chembl"],
            v2_metrics["rows_with_chembl"] - v1_metrics["rows_with_chembl"],
            v1_metrics["rows_with_pubchem"],
            v2_metrics["rows_with_pubchem"],
            v2_metrics["rows_with_pubchem"] - v1_metrics["rows_with_pubchem"],
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
