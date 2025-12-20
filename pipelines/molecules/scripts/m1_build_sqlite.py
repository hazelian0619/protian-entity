#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import random
import sqlite3
import sys
import time
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple


def utc_now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat() + "Z"


def display_path(path: str) -> str:
    """Prefer repo-relative paths in QC/metadata for open-source reproducibility."""

    kg_root = os.environ.get("KG_ROOT")
    if not kg_root:
        return path
    try:
        return os.path.relpath(os.path.abspath(path), kg_root)
    except Exception:
        return path


def normalize_cell(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    value = value.strip()
    if value == "" or value == "\\N":
        return None
    return value


@dataclass
class LengthStats:
    n: int = 0
    mean: float = 0.0
    m2: float = 0.0
    min_v: Optional[int] = None
    max_v: Optional[int] = None

    def add(self, x: int) -> None:
        self.n += 1
        if self.min_v is None or x < self.min_v:
            self.min_v = x
        if self.max_v is None or x > self.max_v:
            self.max_v = x
        delta = x - self.mean
        self.mean += delta / self.n
        delta2 = x - self.mean
        self.m2 += delta * delta2

    def variance(self) -> float:
        if self.n < 2:
            return 0.0
        return self.m2 / (self.n - 1)

    def std(self) -> float:
        return self.variance() ** 0.5


class ReservoirSampler:
    def __init__(self, k: int, seed: int = 0) -> None:
        self.k = k
        self.samples: List[int] = []
        self.n_seen = 0
        self.rng = random.Random(seed)

    def add(self, x: int) -> None:
        self.n_seen += 1
        if self.k <= 0:
            return
        if len(self.samples) < self.k:
            self.samples.append(x)
            return
        j = self.rng.randint(1, self.n_seen)
        if j <= self.k:
            self.samples[j - 1] = x

    def percentiles(self, ps: List[float]) -> Dict[str, Optional[float]]:
        if not self.samples:
            return {f"p{int(p*100)}": None for p in ps}
        data = sorted(self.samples)
        out: Dict[str, Optional[float]] = {}
        for p in ps:
            if p <= 0:
                out[f"p{int(p*100)}"] = float(data[0])
                continue
            if p >= 1:
                out[f"p{int(p*100)}"] = float(data[-1])
                continue
            idx = int(round(p * (len(data) - 1)))
            out[f"p{int(p*100)}"] = float(data[idx])
        return out


def ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)


def connect_sqlite(db_path: str) -> sqlite3.Connection:
    ensure_parent_dir(db_path)
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.execute("PRAGMA temp_store=MEMORY;")
    conn.execute("PRAGMA foreign_keys=ON;")
    conn.execute("PRAGMA cache_size=-200000;")
    return conn


def create_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS meta_run (
            run_id TEXT PRIMARY KEY,
            started_at_utc TEXT NOT NULL,
            finished_at_utc TEXT,
            input_path TEXT NOT NULL,
            input_size_bytes INTEGER,
            input_mtime_utc TEXT,
            input_fingerprint_sha256 TEXT,
            script_version TEXT NOT NULL,
            rules_version TEXT NOT NULL,
            notes TEXT
        );

        CREATE TABLE IF NOT EXISTS molecule_raw_chembl36 (
            ingest_rownum INTEGER NOT NULL,
            chembl_id TEXT,
            smiles TEXT,
            inchi TEXT,
            inchikey TEXT,
            smiles_len INTEGER,
            has_dot INTEGER,
            PRIMARY KEY (ingest_rownum)
        );

        CREATE TABLE IF NOT EXISTS molecule_rejects_chembl36 (
            ingest_rownum INTEGER NOT NULL,
            chembl_id TEXT,
            reason TEXT NOT NULL,
            raw_line_prefix TEXT,
            PRIMARY KEY (ingest_rownum)
        );

        DROP VIEW IF EXISTS molecule_idmap_chembl36;
        CREATE VIEW molecule_idmap_chembl36 AS
        SELECT inchikey, chembl_id
        FROM molecule_raw_chembl36
        WHERE inchikey IS NOT NULL AND chembl_id IS NOT NULL;
        """
    )


def parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="M1: build ChEMBL chemreps SQLite + QC (offline, no deps)")
    p.add_argument(
        "--input",
        default="data/raw/molecules/chembl_36_chemreps.txt",
        help="Path to chemreps TSV (default: data/raw/molecules/chembl_36_chemreps.txt)",
    )
    p.add_argument(
        "--db",
        default="data/output/molecules/molecules_m1.sqlite",
        help="Output SQLite path (default: data/output/molecules/molecules_m1.sqlite)",
    )
    p.add_argument(
        "--outdir",
        default="pipelines/molecules/reports/m1",
        help="Output directory for QC reports (default: pipelines/molecules/reports/m1)",
    )
    p.add_argument("--batch-size", type=int, default=20000, help="SQLite insert batch size")
    p.add_argument(
        "--progress-every",
        type=int,
        default=200000,
        help="Print progress every N data lines",
    )
    p.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Only process first N data lines (0 = full file)",
    )
    p.add_argument(
        "--strict-smiles-len",
        type=int,
        default=200,
        help="Default strict cutoff for smiles length (used for strict table & QC)",
    )
    p.add_argument(
        "--strict-no-dot",
        action="store_true",
        default=True,
        help="Strict table excludes multi-component SMILES (has_dot=1)",
    )
    p.add_argument(
        "--fingerprint-bytes",
        type=int,
        default=10 * 1024 * 1024,
        help="Bytes to hash from start of input for fingerprint (default: 10MB)",
    )
    return p.parse_args(argv)


def file_mtime_utc_iso(path: str) -> Optional[str]:
    try:
        ts = os.path.getmtime(path)
    except OSError:
        return None
    return dt.datetime.fromtimestamp(ts, dt.timezone.utc).replace(microsecond=0).isoformat() + "Z"


def compute_head_fingerprint(path: str, nbytes: int) -> Optional[str]:
    try:
        with open(path, "rb") as f:
            data = f.read(nbytes)
    except OSError:
        return None
    h = hashlib.sha256()
    h.update(data)
    return h.hexdigest()


def build_entity_tables(conn: sqlite3.Connection, strict_smiles_len: int, strict_no_dot: bool) -> None:
    conn.executescript(
        """
        DROP TABLE IF EXISTS molecule_entity_chembl36;
        DROP TABLE IF EXISTS molecule_entity_strict_chembl36;
        """
    )

    conn.execute(
        """
        CREATE TABLE molecule_entity_chembl36 AS
        SELECT
            inchikey AS inchikey,
            SUBSTR(inchikey, 1, 14) AS inchikey_connectivity,
            MIN(chembl_id) AS rep_chembl_id,
            MIN(smiles) AS smiles,
            MIN(inchi) AS inchi,
            MIN(smiles_len) AS smiles_len,
            MAX(has_dot) AS has_dot,
            COUNT(*) AS n_rows
        FROM molecule_raw_chembl36
        WHERE inchikey IS NOT NULL AND smiles IS NOT NULL
        GROUP BY inchikey;
        """
    )

    where_clauses = ["inchikey IS NOT NULL", "smiles IS NOT NULL", f"smiles_len <= {int(strict_smiles_len)}"]
    if strict_no_dot:
        where_clauses.append("has_dot = 0")

    strict_where = " AND ".join(where_clauses)
    conn.execute(
        f"""
        CREATE TABLE molecule_entity_strict_chembl36 AS
        SELECT *
        FROM molecule_entity_chembl36
        WHERE {strict_where};
        """
    )

    conn.executescript(
        """
        CREATE INDEX IF NOT EXISTS idx_raw_inchikey ON molecule_raw_chembl36(inchikey);
        CREATE INDEX IF NOT EXISTS idx_raw_chembl_id ON molecule_raw_chembl36(chembl_id);

        CREATE UNIQUE INDEX IF NOT EXISTS idx_entity_inchikey ON molecule_entity_chembl36(inchikey);
        CREATE INDEX IF NOT EXISTS idx_entity_connectivity ON molecule_entity_chembl36(inchikey_connectivity);
        CREATE INDEX IF NOT EXISTS idx_entity_rep_chembl_id ON molecule_entity_chembl36(rep_chembl_id);

        CREATE UNIQUE INDEX IF NOT EXISTS idx_entity_strict_inchikey ON molecule_entity_strict_chembl36(inchikey);
        """
    )


def sql_scalar(conn: sqlite3.Connection, sql: str, params: Tuple[Any, ...] = ()) -> Any:
    cur = conn.execute(sql, params)
    row = cur.fetchone()
    return None if row is None else row[0]


def write_qc_reports(outdir: str, qc: Dict[str, Any]) -> None:
    os.makedirs(outdir, exist_ok=True)

    json_path = os.path.join(outdir, "qc_report.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(qc, f, ensure_ascii=False, indent=2)

    md_path = os.path.join(outdir, "qc_report.md")
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("# M1 QC Report (ChEMBL chemreps)\n\n")
        f.write(f"- run_id: `{qc.get('run_id')}`\n")
        f.write(f"- started_at_utc: `{qc.get('started_at_utc')}`\n")
        f.write(f"- finished_at_utc: `{qc.get('finished_at_utc')}`\n")
        f.write(f"- input_path: `{qc.get('input', {}).get('path')}`\n")
        f.write(f"- output_db: `{qc.get('output', {}).get('db_path')}`\n\n")

        f.write("## Summary\n\n")
        s = qc.get("summary", {})
        f.write("| metric | value |\n|---|---:|\n")
        for k in [
            "n_input_lines",
            "n_raw_rows",
            "n_rejects",
            "n_entity",
            "n_entity_strict",
            "missing_smiles_rate",
            "missing_inchikey_rate",
            "has_dot_rate",
            "duplicate_inchikey_rows",
        ]:
            if k in s:
                f.write(f"| {k} | {s[k]} |\n")

        f.write("\n## SMILES Length (sampled)\n\n")
        lens = qc.get("smiles_len", {})
        f.write("| stat | value |\n|---|---:|\n")
        for k in ["count", "min", "p50", "p95", "p99", "max", "mean", "std"]:
            if k in lens:
                f.write(f"| {k} | {lens[k]} |\n")

        f.write("\n## Longest SMILES examples (prefix)\n\n")
        examples = qc.get("longest_smiles_examples", [])
        if examples:
            f.write("| rank | smiles_len | chembl_id | smiles_prefix |\n|---:|---:|---|---|\n")
            for i, ex in enumerate(examples, start=1):
                f.write(
                    "| {rank} | {smiles_len} | `{chembl_id}` | `{prefix}` |\n".format(
                        rank=i,
                        smiles_len=ex.get("smiles_len"),
                        chembl_id=ex.get("chembl_id") or "",
                        prefix=(ex.get("smiles_prefix") or "").replace("`", "'"),
                    )
                )
        else:
            f.write("(no examples)\n")


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    input_path = args.input
    outdir = args.outdir
    db_path = args.db

    if not os.path.exists(input_path):
        print(f"ERROR: input not found: {input_path}", file=sys.stderr)
        return 2

    run_id = dt.datetime.now(dt.timezone.utc).strftime("m1_%Y%m%dT%H%M%SZ")
    started_at = utc_now_iso()

    input_size = os.path.getsize(input_path)
    input_mtime = file_mtime_utc_iso(input_path)
    fingerprint = compute_head_fingerprint(input_path, args.fingerprint_bytes)

    conn = connect_sqlite(db_path)
    try:
        create_schema(conn)

        script_version = "m1_build_sqlite_v1"
        rules_version = "m1_rules_v1"

        conn.execute(
            """
            INSERT OR REPLACE INTO meta_run(
                run_id, started_at_utc, finished_at_utc,
                input_path, input_size_bytes, input_mtime_utc, input_fingerprint_sha256,
                script_version, rules_version, notes
            ) VALUES (?, ?, NULL, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                run_id,
                started_at,
                display_path(input_path),
                int(input_size),
                input_mtime,
                fingerprint,
                script_version,
                rules_version,
                "M1 offline build from chemreps (no RDKit; flags only)",
            ),
        )
        conn.commit()

        # Fresh build (keep it simple and deterministic)
        conn.executescript(
            """
            DELETE FROM molecule_raw_chembl36;
            DELETE FROM molecule_rejects_chembl36;
            """
        )
        conn.commit()

        total_data_lines = 0
        raw_rows = 0
        rejects = 0
        missing_smiles = 0
        missing_inchikey = 0
        has_dot = 0

        len_stats = LengthStats()
        sampler = ReservoirSampler(k=200_000, seed=42)
        longest: List[Tuple[int, str, Optional[str]]] = []  # (len, smiles, chembl_id)

        insert_raw = []
        insert_rej = []

        t0 = time.time()

        with open(input_path, "r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            expected = {"chembl_id", "canonical_smiles", "standard_inchi", "standard_inchi_key"}
            header = set(reader.fieldnames or [])
            if not expected.issubset(header):
                raise RuntimeError(f"Unexpected header. Missing: {sorted(expected - header)}; got: {sorted(header)}")

            for row in reader:
                total_data_lines += 1
                if args.limit and total_data_lines > args.limit:
                    break

                chembl_id = normalize_cell(row.get("chembl_id"))
                smiles = normalize_cell(row.get("canonical_smiles"))
                inchi = normalize_cell(row.get("standard_inchi"))
                inchikey = normalize_cell(row.get("standard_inchi_key"))

                ingest_rownum = total_data_lines

                reason = None
                if smiles is None:
                    missing_smiles += 1
                    reason = "missing_smiles"
                if inchikey is None:
                    missing_inchikey += 1
                    reason = "missing_inchikey" if reason is None else (reason + "+missing_inchikey")
                elif len(inchikey) != 27:
                    reason = "bad_inchikey_len"

                if reason is not None:
                    rejects += 1
                    raw_prefix = None
                    try:
                        raw_prefix = "\t".join(
                            [
                                row.get("chembl_id", "") or "",
                                row.get("canonical_smiles", "") or "",
                                row.get("standard_inchi", "") or "",
                                row.get("standard_inchi_key", "") or "",
                            ]
                        )[:200]
                    except Exception:
                        raw_prefix = None
                    insert_rej.append((ingest_rownum, chembl_id, reason, raw_prefix))
                else:
                    slen = len(smiles)
                    dot = 1 if "." in smiles else 0
                    if dot:
                        has_dot += 1

                    len_stats.add(slen)
                    sampler.add(slen)

                    if len(longest) < 20:
                        longest.append((slen, smiles, chembl_id))
                        longest.sort(key=lambda x: x[0], reverse=True)
                    else:
                        if slen > longest[-1][0]:
                            longest.append((slen, smiles, chembl_id))
                            longest.sort(key=lambda x: x[0], reverse=True)
                            longest = longest[:20]

                    insert_raw.append((ingest_rownum, chembl_id, smiles, inchi, inchikey, slen, dot))
                    raw_rows += 1

                if (len(insert_raw) + len(insert_rej)) >= args.batch_size:
                    conn.execute("BEGIN")
                    if insert_raw:
                        conn.executemany(
                            """
                            INSERT INTO molecule_raw_chembl36(
                                ingest_rownum, chembl_id, smiles, inchi, inchikey, smiles_len, has_dot
                            ) VALUES (?, ?, ?, ?, ?, ?, ?)
                            """,
                            insert_raw,
                        )
                    if insert_rej:
                        conn.executemany(
                            """
                            INSERT INTO molecule_rejects_chembl36(
                                ingest_rownum, chembl_id, reason, raw_line_prefix
                            ) VALUES (?, ?, ?, ?)
                            """,
                            insert_rej,
                        )
                    conn.execute("COMMIT")
                    insert_raw.clear()
                    insert_rej.clear()

                if args.progress_every and total_data_lines % args.progress_every == 0:
                    elapsed = time.time() - t0
                    rate = total_data_lines / elapsed if elapsed > 0 else 0
                    print(
                        f"[{run_id}] lines={total_data_lines:,} raw={raw_rows:,} rejects={rejects:,} "
                        f"dot_rate={(has_dot/raw_rows if raw_rows else 0):.2%} "
                        f"({rate:,.0f} lines/s)",
                        file=sys.stderr,
                    )

        if insert_raw or insert_rej:
            conn.execute("BEGIN")
            if insert_raw:
                conn.executemany(
                    """
                    INSERT INTO molecule_raw_chembl36(
                        ingest_rownum, chembl_id, smiles, inchi, inchikey, smiles_len, has_dot
                    ) VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    insert_raw,
                )
            if insert_rej:
                conn.executemany(
                    """
                    INSERT INTO molecule_rejects_chembl36(
                        ingest_rownum, chembl_id, reason, raw_line_prefix
                    ) VALUES (?, ?, ?, ?)
                    """,
                    insert_rej,
                )
            conn.execute("COMMIT")

        build_entity_tables(conn, strict_smiles_len=args.strict_smiles_len, strict_no_dot=args.strict_no_dot)
        conn.commit()

        # QC via SQL
        n_raw = sql_scalar(conn, "SELECT COUNT(*) FROM molecule_raw_chembl36")
        n_rej = sql_scalar(conn, "SELECT COUNT(*) FROM molecule_rejects_chembl36")
        n_entity = sql_scalar(conn, "SELECT COUNT(*) FROM molecule_entity_chembl36")
        n_entity_strict = sql_scalar(conn, "SELECT COUNT(*) FROM molecule_entity_strict_chembl36")
        distinct_inchikey = sql_scalar(conn, "SELECT COUNT(DISTINCT inchikey) FROM molecule_raw_chembl36")
        dup_inchikey_rows = (n_raw - distinct_inchikey) if (n_raw is not None and distinct_inchikey is not None) else None

        finished_at = utc_now_iso()
        conn.execute(
            "UPDATE meta_run SET finished_at_utc=? WHERE run_id=?",
            (finished_at, run_id),
        )
        conn.commit()

        pcts = sampler.percentiles([0.50, 0.95, 0.99])
        longest_examples = [
            {
                "smiles_len": l,
                "chembl_id": cid,
                "smiles_prefix": s[:180],
            }
            for (l, s, cid) in longest
        ]

        qc: Dict[str, Any] = {
            "run_id": run_id,
            "started_at_utc": started_at,
            "finished_at_utc": finished_at,
            "input": {
                "path": display_path(input_path),
                "size_bytes": input_size,
                "mtime_utc": input_mtime,
                "head_fingerprint_sha256": fingerprint,
            },
            "output": {
                "db_path": display_path(db_path),
                "outdir": display_path(outdir),
            },
            "parameters": {
                "limit": args.limit,
                "batch_size": args.batch_size,
                "strict_smiles_len": args.strict_smiles_len,
                "strict_no_dot": bool(args.strict_no_dot),
                "fingerprint_bytes": args.fingerprint_bytes,
                "reservoir_k": sampler.k,
                "progress_every": args.progress_every,
            },
            "summary": {
                "n_input_lines": total_data_lines,
                "n_raw_rows": int(n_raw or 0),
                "n_rejects": int(n_rej or 0),
                "n_entity": int(n_entity or 0),
                "n_entity_strict": int(n_entity_strict or 0),
                "missing_smiles_rate": round((missing_smiles / total_data_lines) if total_data_lines else 0.0, 6),
                "missing_inchikey_rate": round((missing_inchikey / total_data_lines) if total_data_lines else 0.0, 6),
                "has_dot_rate": round((has_dot / raw_rows) if raw_rows else 0.0, 6),
                "duplicate_inchikey_rows": int(dup_inchikey_rows or 0),
            },
            "smiles_len": {
                "count": len_stats.n,
                "min": len_stats.min_v,
                "max": len_stats.max_v,
                "mean": round(len_stats.mean, 3),
                "std": round(len_stats.std(), 3),
                "p50": pcts.get("p50"),
                "p95": pcts.get("p95"),
                "p99": pcts.get("p99"),
            },
            "longest_smiles_examples": longest_examples,
        }

        write_qc_reports(outdir, qc)

        elapsed = time.time() - t0
        print(
            f"DONE {run_id}: lines={total_data_lines:,} raw={n_raw:,} rejects={n_rej:,} "
            f"entity={n_entity:,} strict={n_entity_strict:,} elapsed={elapsed:.1f}s",
            file=sys.stderr,
        )
        return 0

    finally:
        conn.close()


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
