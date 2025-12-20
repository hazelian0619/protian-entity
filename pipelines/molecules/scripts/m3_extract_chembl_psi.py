#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""M3: Extract ChEMBL PSI (protein–small molecule) activity evidence edges.

Design goals:
- Offline + standard library only (sqlite3).
- Read input ChEMBL SQLite (user-downloaded) in read-only mode.
- Write outputs to a NEW SQLite DB (avoid locking `molecules_m1.sqlite`).
- Produce both evidence-granularity rows (activity-level) and aggregated edges.

Default paths follow the repo conventions:
- Input:  data/raw/molecules/chembl_36/chembl_36.db
- Output: data/output/molecules/chembl_m3.sqlite

This script is intentionally conservative: it focuses on a minimal, traceable v1.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import sqlite3
import sys
import time
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


def utc_now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def display_path(path: str) -> str:
    kg_root = os.environ.get("KG_ROOT")
    if not kg_root:
        return path
    try:
        return os.path.relpath(os.path.abspath(path), kg_root)
    except Exception:
        return path


def ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)


def connect_sqlite_ro(db_path: str) -> sqlite3.Connection:
    # Use SQLite URI for read-only.
    uri = f"file:{os.path.abspath(db_path)}?mode=ro"
    conn = sqlite3.connect(uri, uri=True)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA query_only=ON;")
    return conn


def connect_sqlite_rw(db_path: str) -> sqlite3.Connection:
    ensure_parent_dir(db_path)
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.execute("PRAGMA temp_store=MEMORY;")
    conn.execute("PRAGMA foreign_keys=ON;")
    conn.execute("PRAGMA cache_size=-200000;")

    def sha1_hex(s: Optional[str]) -> Optional[str]:
        if s is None:
            return None
        return hashlib.sha1(s.encode("utf-8")).hexdigest()

    conn.create_function("sha1_hex", 1, sha1_hex)
    return conn


def sql_scalar(conn: sqlite3.Connection, sql: str, params: Sequence[Any] = ()) -> Any:
    cur = conn.execute(sql, params)
    row = cur.fetchone()
    return None if row is None else row[0]


def table_exists(conn: sqlite3.Connection, table: str) -> bool:
    return (
        sql_scalar(
            conn,
            "SELECT 1 FROM sqlite_master WHERE type IN ('table','view') AND name = ? LIMIT 1",
            (table,),
        )
        == 1
    )


def compute_file_sha256_head(path: str, nbytes: int = 10 * 1024 * 1024) -> Optional[str]:
    try:
        with open(path, "rb") as f:
            data = f.read(nbytes)
    except OSError:
        return None
    h = hashlib.sha256()
    h.update(data)
    return h.hexdigest()


def parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="M3: Extract ChEMBL PSI evidence + aggregated edges (offline)")

    p.add_argument(
        "--chembl-db",
        default="data/raw/molecules/chembl_36/chembl_36.db",
        help="Path to ChEMBL SQLite DB (default: data/raw/molecules/chembl_36/chembl_36.db)",
    )
    p.add_argument(
        "--out-db",
        default="data/output/molecules/chembl_m3.sqlite",
        help="Output SQLite path (default: data/output/molecules/chembl_m3.sqlite)",
    )
    p.add_argument(
        "--outdir",
        default="pipelines/molecules/reports/m3",
        help="Directory for QC outputs (default: pipelines/molecules/reports/m3)",
    )

    p.add_argument(
        "--mode",
        choices=["all", "aggregate"],
        default="all",
        help="Run mode: all=extract evidence + aggregate edges; aggregate=aggregate/QC only from existing evidence",
    )
    p.add_argument(
        "--qc-prefix",
        default="qc_m3_chembl",
        help="QC output file prefix written to --outdir (default: qc_m3_chembl)",
    )

    p.add_argument(
        "--types",
        default="IC50,Ki,Kd,EC50",
        help="Comma-separated standard_type allowlist (default: IC50,Ki,Kd,EC50)",
    )
    p.add_argument(
        "--min-confidence",
        type=int,
        default=0,
        help="Minimum assay confidence_score (default: 0)",
    )
    p.add_argument(
        "--require-inchikey",
        action="store_true",
        help="Drop evidence rows missing standard_inchi_key",
    )
    p.add_argument(
        "--require-uniprot",
        action="store_true",
        help="Drop evidence rows missing Uniprot accession",
    )

    p.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Limit evidence rows extracted (0 = no limit). Useful for smoke tests.",
    )
    p.add_argument("--batch-size", type=int, default=5000, help="Insert batch size")
    p.add_argument("--progress-every", type=int, default=100000, help="Progress print frequency")

    p.add_argument(
        "--start-activity-id",
        type=int,
        default=0,
        help="Optional lower bound (inclusive) for activity_id (default: 0)",
    )
    p.add_argument(
        "--end-activity-id",
        type=int,
        default=0,
        help="Optional upper bound (inclusive) for activity_id (default: 0 = no bound)",
    )

    return p.parse_args(argv)


def create_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS meta_run (
            run_id TEXT PRIMARY KEY,
            started_at_utc TEXT NOT NULL,
            finished_at_utc TEXT,
            input_path TEXT NOT NULL,
            input_sha256_head TEXT,
            script_version TEXT NOT NULL,
            rules_version TEXT NOT NULL,
            notes TEXT
        );

        CREATE TABLE IF NOT EXISTS psi_evidence_v1 (
            activity_id INTEGER PRIMARY KEY,
            assay_id INTEGER,
            doc_id INTEGER,

            compound_chembl_id TEXT,
            molregno INTEGER,
            compound_inchikey TEXT,
            compound_canonical_smiles TEXT,

            target_tid INTEGER,
            target_chembl_id TEXT,
            target_uniprot_accession TEXT,

            assay_type TEXT,
            assay_confidence_score INTEGER,
            assay_description TEXT,
            bao_format TEXT,

            standard_type TEXT,
            standard_relation TEXT,
            standard_value REAL,
            standard_units TEXT,
            standard_value_nM REAL,

            pchembl_value REAL,
            pchembl_value_eff REAL,

            data_validity_comment TEXT,
            activity_comment TEXT,

            doi TEXT,
            pubmed_id INTEGER,
            source TEXT NOT NULL,

            evidence_score_v1 REAL,
            extracted_at_utc TEXT NOT NULL
        );

        CREATE INDEX IF NOT EXISTS idx_psi_ev_compound_inchikey ON psi_evidence_v1(compound_inchikey);
        CREATE INDEX IF NOT EXISTS idx_psi_ev_target_uniprot ON psi_evidence_v1(target_uniprot_accession);
        CREATE INDEX IF NOT EXISTS idx_psi_ev_type ON psi_evidence_v1(standard_type);
        CREATE INDEX IF NOT EXISTS idx_psi_ev_conf ON psi_evidence_v1(assay_confidence_score);

        CREATE TABLE IF NOT EXISTS psi_edges_v1 (
            edge_id TEXT PRIMARY KEY,

            compound_inchikey TEXT NOT NULL,
            compound_chembl_id TEXT,
            target_uniprot_accession TEXT,
            target_chembl_id TEXT NOT NULL,
            standard_type TEXT NOT NULL,

            n_evidence INTEGER NOT NULL,
            n_docs INTEGER NOT NULL,

            pchembl_max REAL,
            pchembl_mean REAL,
            assay_confidence_score_max INTEGER,
            evidence_score_max REAL,

            best_activity_id INTEGER,
            best_assay_id INTEGER,
            best_doc_id INTEGER,
            best_doi TEXT,
            best_pubmed_id INTEGER,

            aggregated_at_utc TEXT NOT NULL
        );

        CREATE INDEX IF NOT EXISTS idx_psi_edges_compound ON psi_edges_v1(compound_inchikey);
        CREATE INDEX IF NOT EXISTS idx_psi_edges_target ON psi_edges_v1(target_chembl_id);
        CREATE INDEX IF NOT EXISTS idx_psi_edges_uniprot ON psi_edges_v1(target_uniprot_accession);
        CREATE INDEX IF NOT EXISTS idx_psi_edges_type ON psi_edges_v1(standard_type);
        CREATE INDEX IF NOT EXISTS idx_psi_edges_best_activity ON psi_edges_v1(best_activity_id);

        CREATE TABLE IF NOT EXISTS target_uniprot_map_v1 (
            target_tid INTEGER NOT NULL,
            target_chembl_id TEXT,
            component_id INTEGER NOT NULL,
            target_uniprot_accession TEXT NOT NULL,
            PRIMARY KEY (target_tid, component_id)
        );

        CREATE INDEX IF NOT EXISTS idx_tu_target_chembl_id ON target_uniprot_map_v1(target_chembl_id);
        CREATE INDEX IF NOT EXISTS idx_tu_accession ON target_uniprot_map_v1(target_uniprot_accession);
        """
    )


def populate_target_uniprot_map(out_conn: sqlite3.Connection, in_conn: sqlite3.Connection) -> None:
    """Populate target→UniProt mapping table from ChEMBL tables.

    This is independent of the evidence extraction filters; it enables QA and
    downstream target normalization.
    """

    out_conn.execute("DELETE FROM target_uniprot_map_v1")
    out_conn.commit()

    sql = """
    SELECT
      td.tid AS target_tid,
      td.chembl_id AS target_chembl_id,
      tc.component_id AS component_id,
      csq.accession AS target_uniprot_accession
    FROM target_dictionary td
    JOIN target_components tc
      ON td.tid = tc.tid
    JOIN component_sequences csq
      ON tc.component_id = csq.component_id
    WHERE csq.accession IS NOT NULL
    """

    cur = in_conn.execute(sql)
    batch: list[tuple[Any, ...]] = []

    out_conn.execute("BEGIN")
    while True:
        rows = cur.fetchmany(5000)
        if not rows:
            break
        for r in rows:
            batch.append((r["target_tid"], r["target_chembl_id"], r["component_id"], r["target_uniprot_accession"]))
        out_conn.executemany(
            """
            INSERT OR IGNORE INTO target_uniprot_map_v1(
                target_tid, target_chembl_id, component_id, target_uniprot_accession
            ) VALUES (?, ?, ?, ?)
            """,
            batch,
        )
        batch.clear()
    out_conn.execute("COMMIT")



def unit_to_nM(value: Optional[float], units: Optional[str]) -> Optional[float]:
    if value is None or units is None:
        return None
    u = units.strip().lower()
    if u in {"nm", "nanomolar"}:
        return float(value)
    if u in {"pm", "picomolar"}:
        return float(value) / 1000.0
    if u in {"um", "µm", "micromolar"}:
        return float(value) * 1000.0
    if u in {"mm", "millimolar"}:
        return float(value) * 1_000_000.0
    if u in {"m", "mol/l", "molar"}:
        return float(value) * 1_000_000_000.0
    return None


def pchembl_from_nM(value_nM: Optional[float]) -> Optional[float]:
    if value_nM is None:
        return None
    if value_nM <= 0:
        return None
    # pChEMBL for concentration-like measures: p = -log10(M)
    # M = nM * 1e-9  => p = 9 - log10(nM)
    return 9.0 - math.log10(float(value_nM))


def evidence_score_v1(
    assay_confidence_score: Optional[int],
    assay_type: Optional[str],
    standard_relation: Optional[str],
    standard_units: Optional[str],
) -> float:
    conf = float(assay_confidence_score or 0)
    conf_term = max(0.0, min(1.0, conf / 9.0))

    at = (assay_type or "").strip().upper()
    assay_type_term = {
        "B": 1.0,
        "F": 0.9,
        "A": 0.6,
        "T": 0.4,
        "P": 0.2,
    }.get(at, 0.5)

    rel = (standard_relation or "").strip()
    relation_term = {
        "=": 1.0,
        "<=": 0.9,
        ">=": 0.9,
        "<": 0.7,
        ">": 0.7,
        "~": 0.5,
    }.get(rel, 0.6)

    u = (standard_units or "").strip().lower()
    units_term = {
        "pm": 1.0,
        "nm": 1.0,
        "um": 0.8,
        "µm": 0.8,
        "mm": 0.6,
        "m": 0.5,
    }.get(u, 0.5)

    score = 0.45 * conf_term + 0.25 * assay_type_term + 0.2 * relation_term + 0.1 * units_term
    return float(max(0.0, min(1.0, score)))


def required_input_tables() -> List[str]:
    return [
        "activities",
        "assays",
        "docs",
        "molecule_dictionary",
        "compound_structures",
        "target_dictionary",
        "target_components",
        "component_sequences",
    ]


def build_evidence_query(standard_types: List[str], require_inchikey: bool, require_uniprot: bool) -> str:
    placeholders = ",".join(["?"] * len(standard_types))

    where = [
        "a.standard_value IS NOT NULL",
        "a.standard_type IS NOT NULL",
        f"a.standard_type IN ({placeholders})",
    ]

    # ChEMBL stores confidence_score on assays.
    where.append("COALESCE(ass.confidence_score, 0) >= ?")

    if require_inchikey:
        where.append("cs.standard_inchi_key IS NOT NULL")
    if require_uniprot:
        where.append(
            "EXISTS (SELECT 1 FROM target_components tc JOIN component_sequences csq ON tc.component_id=csq.component_id WHERE tc.tid=ass.tid AND csq.accession IS NOT NULL)"
        )

    # Optional activity_id bounds will be appended later.

    where_sql = "\n      AND ".join(where)

    # Notes:
    # - target_uniprot_accession is a deterministic representative (MIN accession) for the target.
    # - Full mapping (target tid -> component -> UniProt) is stored in target_uniprot_map_v1.
    return f"""
    SELECT
        a.activity_id AS activity_id,
        a.assay_id AS assay_id,
        ass.doc_id AS doc_id,

        md.chembl_id AS compound_chembl_id,
        a.molregno AS molregno,
        cs.standard_inchi_key AS compound_inchikey,
        cs.canonical_smiles AS compound_canonical_smiles,

        ass.tid AS target_tid,
        td.chembl_id AS target_chembl_id,
        (
          SELECT COALESCE(
            MIN(CASE WHEN csq.db_source = 'SWISS-PROT' THEN csq.accession END),
            MIN(csq.accession)
          )
          FROM target_components tc
          JOIN component_sequences csq
            ON tc.component_id = csq.component_id
          WHERE tc.tid = ass.tid
            AND csq.accession IS NOT NULL
        ) AS target_uniprot_accession,

        ass.assay_type AS assay_type,
        ass.confidence_score AS assay_confidence_score,
        ass.description AS assay_description,
        ass.bao_format AS bao_format,

        a.standard_type AS standard_type,
        a.standard_relation AS standard_relation,
        a.standard_value AS standard_value,
        a.standard_units AS standard_units,

        a.pchembl_value AS pchembl_value,
        a.data_validity_comment AS data_validity_comment,
        a.activity_comment AS activity_comment,

        d.doi AS doi,
        d.pubmed_id AS pubmed_id
    FROM activities AS a
    JOIN assays AS ass
      ON a.assay_id = ass.assay_id
    JOIN molecule_dictionary AS md
      ON a.molregno = md.molregno
    LEFT JOIN compound_structures AS cs
      ON a.molregno = cs.molregno
    JOIN target_dictionary AS td
      ON ass.tid = td.tid
    LEFT JOIN docs AS d
      ON ass.doc_id = d.doc_id
    WHERE {where_sql}
    ORDER BY a.activity_id
    """


def iter_evidence_rows(
    in_conn: sqlite3.Connection,
    sql: str,
    params: Sequence[Any],
    fetch_size: int = 5000,
) -> Iterable[sqlite3.Row]:
    cur = in_conn.execute(sql, params)
    while True:
        rows = cur.fetchmany(fetch_size)
        if not rows:
            break
        for r in rows:
            yield r


def insert_many(conn: sqlite3.Connection, sql: str, rows: List[Tuple[Any, ...]]) -> None:
    if not rows:
        return
    conn.executemany(sql, rows)


def aggregate_edges(conn: sqlite3.Connection, aggregated_at: str) -> None:
    # v1 aggregation rule:
    # - group by (compound_inchikey, COALESCE(target_uniprot_accession,target_chembl_id), target_chembl_id, standard_type)
    # - pick best evidence by pchembl_value_eff desc, evidence_score desc, confidence desc, activity_id asc
    # This must run in SQL for performance at full scale (millions of rows).
    conn.execute("DELETE FROM psi_edges_v1")
    conn.commit()

    conn.execute(
        """
        CREATE INDEX IF NOT EXISTS idx_psi_ev_group
        ON psi_evidence_v1(
          compound_inchikey,
          target_chembl_id,
          target_uniprot_accession,
          standard_type,
          activity_id
        )
        """
    )

    conn.commit()
    conn.execute("BEGIN")
    conn.execute(
        """
        WITH base AS (
          SELECT
            activity_id,
            assay_id,
            doc_id,
            compound_inchikey,
            compound_chembl_id,
            target_chembl_id,
            target_uniprot_accession,
            COALESCE(target_uniprot_accession, target_chembl_id) AS target_key,
            standard_type,
            pchembl_value_eff,
            assay_confidence_score,
            evidence_score_v1,
            doi,
            pubmed_id
          FROM psi_evidence_v1
          WHERE compound_inchikey IS NOT NULL
            AND target_chembl_id IS NOT NULL
            AND standard_type IS NOT NULL
        ),
        stats AS (
          SELECT
            compound_inchikey,
            target_chembl_id,
            target_uniprot_accession,
            target_key,
            standard_type,
            COUNT(*) AS n_evidence,
            COUNT(DISTINCT doc_id) AS n_docs,
            MAX(pchembl_value_eff) AS pchembl_max,
            AVG(pchembl_value_eff) AS pchembl_mean,
            MAX(assay_confidence_score) AS assay_confidence_score_max,
            MAX(evidence_score_v1) AS evidence_score_max
          FROM base
          GROUP BY compound_inchikey, target_key, target_chembl_id, target_uniprot_accession, standard_type
        ),
        ranked AS (
          SELECT
            *,
            ROW_NUMBER() OVER (
              PARTITION BY compound_inchikey, target_key, target_chembl_id, standard_type
              ORDER BY
                COALESCE(pchembl_value_eff, -1e9) DESC,
                COALESCE(evidence_score_v1, -1e9) DESC,
                COALESCE(assay_confidence_score, -1e9) DESC,
                activity_id ASC
            ) AS rn
          FROM base
        ),
        best AS (
          SELECT * FROM ranked WHERE rn = 1
        )
        INSERT OR REPLACE INTO psi_edges_v1(
          edge_id,
          compound_inchikey, compound_chembl_id,
          target_uniprot_accession, target_chembl_id,
          standard_type,
          n_evidence, n_docs,
          pchembl_max, pchembl_mean,
          assay_confidence_score_max, evidence_score_max,
          best_activity_id, best_assay_id, best_doc_id,
          best_doi, best_pubmed_id,
          aggregated_at_utc
        )
        SELECT
          sha1_hex(stats.compound_inchikey || '|' || stats.target_key || '|' || stats.standard_type) AS edge_id,
          stats.compound_inchikey,
          best.compound_chembl_id,
          stats.target_uniprot_accession,
          stats.target_chembl_id,
          stats.standard_type,
          stats.n_evidence,
          stats.n_docs,
          stats.pchembl_max,
          stats.pchembl_mean,
          stats.assay_confidence_score_max,
          stats.evidence_score_max,
          best.activity_id,
          best.assay_id,
          best.doc_id,
          best.doi,
          best.pubmed_id,
          ?
        FROM stats
        JOIN best
          ON best.compound_inchikey = stats.compound_inchikey
         AND best.target_key = stats.target_key
         AND best.target_chembl_id = stats.target_chembl_id
         AND best.standard_type = stats.standard_type
        """,
        (aggregated_at,),
    )
    conn.execute("COMMIT")


def write_qc(outdir: str, qc_prefix: str, qc: Dict[str, Any]) -> None:
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, f"{qc_prefix}.json"), "w", encoding="utf-8") as f:
        json.dump(qc, f, ensure_ascii=False, indent=2)

    with open(os.path.join(outdir, f"{qc_prefix}.md"), "w", encoding="utf-8") as f:
        f.write("# M3 ChEMBL PSI QC (v1)\n\n")
        f.write(f"- run_id: `{qc.get('run_id')}`\n")
        f.write(f"- started_at_utc: `{qc.get('started_at_utc')}`\n")
        f.write(f"- finished_at_utc: `{qc.get('finished_at_utc')}`\n")
        f.write(f"- input_db: `{qc.get('input', {}).get('chembl_db')}`\n")
        f.write(f"- output_db: `{qc.get('output', {}).get('out_db')}`\n\n")

        f.write("## Counts\n\n")
        f.write("| metric | value |\n|---|---:|\n")
        for k, v in qc.get("counts", {}).items():
            f.write(f"| {k} | {v} |\n")

        f.write("\n## Distributions\n\n")
        dists = qc.get("distributions", {})
        for name, rows in dists.items():
            f.write(f"### {name}\n\n")
            f.write("| key | n |\n|---|---:|\n")
            for key, n in rows:
                f.write(f"| {key} | {n} |\n")
            f.write("\n")


def topk_counts(conn: sqlite3.Connection, sql: str, params: Sequence[Any] = (), k: int = 20) -> List[Tuple[str, int]]:
    cur = conn.execute(sql, params)
    out: List[Tuple[str, int]] = []
    for row in cur.fetchmany(k):
        out.append((str(row[0]), int(row[1])))
    return out


def main(argv: List[str]) -> int:
    args = parse_args(argv)

    standard_types = [t.strip() for t in (args.types or "").split(",") if t.strip()]
    if not standard_types:
        print("ERROR: --types must contain at least one standard_type", file=sys.stderr)
        return 2

    if not os.path.exists(args.chembl_db):
        print(f"ERROR: ChEMBL DB not found: {args.chembl_db}", file=sys.stderr)
        print(
            "Hint: download chembl_36_sqlite.tar.gz and extract chembl_36.db into data/raw/molecules/chembl_36/",
            file=sys.stderr,
        )
        return 2

    os.makedirs(args.outdir, exist_ok=True)

    run_id = dt.datetime.now(dt.timezone.utc).strftime("m3_%Y%m%dT%H%M%SZ")
    started_at = utc_now_iso()
    extracted_at = utc_now_iso()

    in_conn = connect_sqlite_ro(args.chembl_db)
    try:
        missing = [t for t in required_input_tables() if not table_exists(in_conn, t)]
        if missing:
            print("ERROR: input ChEMBL SQLite is missing required tables:", file=sys.stderr)
            for t in missing:
                print(f"  - {t}", file=sys.stderr)
            return 2

        out_conn = connect_sqlite_rw(args.out_db)
        try:
            create_schema(out_conn)

            if args.mode == "all":
                # Make runs idempotent: never mix outputs across different filter settings.
                out_conn.executescript(
                    "DELETE FROM psi_evidence_v1; DELETE FROM psi_edges_v1; DELETE FROM target_uniprot_map_v1;"
                )

            populate_target_uniprot_map(out_conn, in_conn)

            input_fingerprint = compute_file_sha256_head(args.chembl_db)
            out_conn.execute(
                """
                INSERT OR REPLACE INTO meta_run(
                    run_id, started_at_utc, finished_at_utc,
                    input_path, input_sha256_head,
                    script_version, rules_version, notes
                ) VALUES (?, ?, NULL, ?, ?, ?, ?, ?)
                """,
                (
                    run_id,
                    started_at,
                    display_path(args.chembl_db),
                    input_fingerprint,
                    "m3_extract_chembl_psi_v1",
                    "m3_rules_v1",
                    (
                        "PSI evidence edges from ChEMBL activities/assays/targets/docs"
                        f" | mode={args.mode}"
                        f" | types={','.join(standard_types)}"
                        f" | min_conf={int(args.min_confidence)}"
                        f" | require_inchikey={bool(args.require_inchikey)}"
                        f" | require_uniprot={bool(args.require_uniprot)}"
                        f" | limit={int(args.limit)}"
                    ),
                ),
            )
            out_conn.commit()

            t0 = time.time()

            if args.mode == "all":
                q = build_evidence_query(standard_types, args.require_inchikey, args.require_uniprot)
                params: List[Any] = []
                params.extend(standard_types)
                params.append(int(args.min_confidence))

                # Add activity_id bounds as a second-stage filter to keep the main query stable.
                bounds: List[str] = []
                if args.start_activity_id:
                    bounds.append("activity_id >= ?")
                    params.append(int(args.start_activity_id))
                if args.end_activity_id:
                    bounds.append("activity_id <= ?")
                    params.append(int(args.end_activity_id))

                if bounds:
                    q = f"SELECT * FROM ({q}) WHERE {' AND '.join(bounds)} ORDER BY activity_id"

                insert_sql = (
                    "INSERT OR REPLACE INTO psi_evidence_v1("
                    "activity_id, assay_id, doc_id, "
                    "compound_chembl_id, molregno, compound_inchikey, compound_canonical_smiles, "
                    "target_tid, target_chembl_id, target_uniprot_accession, "
                    "assay_type, assay_confidence_score, assay_description, bao_format, "
                    "standard_type, standard_relation, standard_value, standard_units, standard_value_nM, "
                    "pchembl_value, pchembl_value_eff, "
                    "data_validity_comment, activity_comment, "
                    "doi, pubmed_id, source, evidence_score_v1, extracted_at_utc"
                    ") VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
                )

                n_seen = 0
                n_inserted = 0
                buffer: List[Tuple[Any, ...]] = []

                for r in iter_evidence_rows(in_conn, q, params):
                    if args.limit and n_seen >= args.limit:
                        break

                    n_seen += 1
                    std_val = r["standard_value"]
                    std_units = r["standard_units"]
                    std_rel = r["standard_relation"]

                    val_nM = unit_to_nM(std_val, std_units)
                    pchembl_eff = r["pchembl_value"]
                    if pchembl_eff is None:
                        # Compute fallback pChEMBL only when value is convertible.
                        pchembl_eff = pchembl_from_nM(val_nM)

                    score = evidence_score_v1(
                        assay_confidence_score=r["assay_confidence_score"],
                        assay_type=r["assay_type"],
                        standard_relation=std_rel,
                        standard_units=std_units,
                    )

                    buffer.append(
                        (
                            r["activity_id"],
                            r["assay_id"],
                            r["doc_id"],
                            r["compound_chembl_id"],
                            r["molregno"],
                            r["compound_inchikey"],
                            r["compound_canonical_smiles"],
                            r["target_tid"],
                            r["target_chembl_id"],
                            r["target_uniprot_accession"],
                            r["assay_type"],
                            r["assay_confidence_score"],
                            r["assay_description"],
                            r["bao_format"],
                            r["standard_type"],
                            std_rel,
                            std_val,
                            std_units,
                            val_nM,
                            r["pchembl_value"],
                            pchembl_eff,
                            r["data_validity_comment"],
                            r["activity_comment"],
                            r["doi"],
                            r["pubmed_id"],
                            "chembl_36",
                            score,
                            extracted_at,
                        )
                    )

                    if len(buffer) >= args.batch_size:
                        out_conn.execute("BEGIN")
                        insert_many(out_conn, insert_sql, buffer)
                        out_conn.execute("COMMIT")
                        n_inserted += len(buffer)
                        buffer.clear()

                    if args.progress_every and n_seen % args.progress_every == 0:
                        elapsed = time.time() - t0
                        rate = n_seen / elapsed if elapsed > 0 else 0
                        print(f"[{run_id}] seen={n_seen:,} inserted={n_inserted:,} ({rate:,.0f}/s)", file=sys.stderr)

                if buffer:
                    out_conn.execute("BEGIN")
                    insert_many(out_conn, insert_sql, buffer)
                    out_conn.execute("COMMIT")
                    n_inserted += len(buffer)
            else:
                n_seen = 0
                n_inserted = 0

            aggregated_at = utc_now_iso()
            print(f"[{run_id}] aggregating edges...", file=sys.stderr)
            aggregate_edges(out_conn, aggregated_at=aggregated_at)


            finished_at = utc_now_iso()

            out_conn.execute("UPDATE meta_run SET finished_at_utc = ? WHERE run_id = ?", (finished_at, run_id))
            out_conn.commit()

            # QC (computed from output DB to avoid reading input twice)
            qc: Dict[str, Any] = {
                "run_id": run_id,
                "started_at_utc": started_at,
                "finished_at_utc": finished_at,
                "input": {
                    "chembl_db": display_path(args.chembl_db),
                    "sha256_head": input_fingerprint,
                    "types": standard_types,
                    "min_confidence": int(args.min_confidence),
                    "require_inchikey": bool(args.require_inchikey),
                    "require_uniprot": bool(args.require_uniprot),
                    "activity_id_bounds": {
                        "start": int(args.start_activity_id) if args.start_activity_id else None,
                        "end": int(args.end_activity_id) if args.end_activity_id else None,
                    },
                },
                "output": {"out_db": display_path(args.out_db), "outdir": display_path(args.outdir)},
                "counts": {
                    "evidence_rows": int(sql_scalar(out_conn, "SELECT COUNT(*) FROM psi_evidence_v1") or 0),
                    "edges": int(sql_scalar(out_conn, "SELECT COUNT(*) FROM psi_edges_v1") or 0),
                    "target_map_rows": int(sql_scalar(out_conn, "SELECT COUNT(*) FROM target_uniprot_map_v1") or 0),
                    "distinct_compounds_inchikey": int(
                        sql_scalar(out_conn, "SELECT COUNT(DISTINCT compound_inchikey) FROM psi_evidence_v1") or 0
                    ),
                    "distinct_targets_uniprot": int(
                        sql_scalar(out_conn, "SELECT COUNT(DISTINCT target_uniprot_accession) FROM psi_evidence_v1") or 0
                    ),
                    "distinct_targets_chembl": int(
                        sql_scalar(out_conn, "SELECT COUNT(DISTINCT target_chembl_id) FROM psi_evidence_v1") or 0
                    ),
                    "missing_inchikey": int(
                        sql_scalar(out_conn, "SELECT COUNT(*) FROM psi_evidence_v1 WHERE compound_inchikey IS NULL") or 0
                    ),
                    "missing_uniprot": int(
                        sql_scalar(out_conn, "SELECT COUNT(*) FROM psi_evidence_v1 WHERE target_uniprot_accession IS NULL")
                        or 0
                    ),
                    "missing_pchembl": int(
                        sql_scalar(out_conn, "SELECT COUNT(*) FROM psi_evidence_v1 WHERE pchembl_value IS NULL") or 0
                    ),
                    "missing_pchembl_eff": int(
                        sql_scalar(out_conn, "SELECT COUNT(*) FROM psi_evidence_v1 WHERE pchembl_value_eff IS NULL") or 0
                    ),
                },
                "distributions": {
                    "standard_type": topk_counts(
                        out_conn,
                        "SELECT standard_type, COUNT(*) AS n FROM psi_evidence_v1 GROUP BY standard_type ORDER BY n DESC",
                        (),
                        20,
                    ),
                    "standard_units": topk_counts(
                        out_conn,
                        "SELECT COALESCE(standard_units,'(NULL)') AS u, COUNT(*) AS n FROM psi_evidence_v1 GROUP BY u ORDER BY n DESC",
                        (),
                        20,
                    ),
                    "assay_type": topk_counts(
                        out_conn,
                        "SELECT COALESCE(assay_type,'(NULL)') AS t, COUNT(*) AS n FROM psi_evidence_v1 GROUP BY t ORDER BY n DESC",
                        (),
                        20,
                    ),
                    "assay_confidence_score": topk_counts(
                        out_conn,
                        "SELECT COALESCE(assay_confidence_score,-1) AS c, COUNT(*) AS n FROM psi_evidence_v1 GROUP BY c ORDER BY n DESC",
                        (),
                        20,
                    ),
                },
                "runtime": {
                    "seconds": round(time.time() - t0, 2),
                    "seen": int(n_seen),
                    "inserted": int(n_inserted),
                },
            }

            write_qc(args.outdir, args.qc_prefix, qc)
            out_conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
            return 0
        finally:
            out_conn.close()
    finally:
        in_conn.close()


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
