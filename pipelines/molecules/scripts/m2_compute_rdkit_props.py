#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import datetime as dt
import json
import os
import sqlite3
import sys
import time
from typing import Any, Dict, Iterable, List, Tuple


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


def connect_sqlite(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.execute("PRAGMA temp_store=MEMORY;")
    conn.execute("PRAGMA foreign_keys=ON;")
    return conn


def create_schema(conn: sqlite3.Connection, table_name: str) -> None:
    conn.executescript(
        f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            inchikey TEXT PRIMARY KEY,
            rdkit_canonical_smiles TEXT,
            mol_wt REAL,
            mol_logp REAL,
            tpsa REAL,
            hbd INTEGER,
            hba INTEGER,
            rot_bonds INTEGER,
            rings INTEGER,
            heavy_atoms INTEGER,
            frac_csp3 REAL,
            rdkit_sanitize_ok INTEGER NOT NULL,
            rdkit_error TEXT,
            computed_at_utc TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS {table_name}_rejects (
            inchikey TEXT,
            smiles TEXT,
            error TEXT,
            computed_at_utc TEXT NOT NULL
        );

        CREATE INDEX IF NOT EXISTS idx_{table_name}_sanitize ON {table_name}(rdkit_sanitize_ok);
        """
    )


def parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="M2: compute RDKit physchem descriptors for molecules_m1.sqlite")
    p.add_argument(
        "--db",
        default="data/output/molecules/molecules_m1.sqlite",
        help="Input SQLite from M1 (default: data/output/molecules/molecules_m1.sqlite)",
    )
    p.add_argument(
        "--source-table",
        default="molecule_entity_strict_chembl36",
        help="Source table with at least (inchikey, smiles)",
    )
    p.add_argument(
        "--outdir",
        default="pipelines/molecules/reports/m2",
        help="Output directory for QC reports (default: pipelines/molecules/reports/m2)",
    )
    p.add_argument(
        "--table",
        default="molecule_props_rdkit_v1",
        help="Output table name (default: molecule_props_rdkit_v1)",
    )
    p.add_argument("--limit", type=int, default=100000, help="Limit rows to compute (0 = all)")
    p.add_argument("--batch-size", type=int, default=5000, help="Insert batch size")
    p.add_argument("--progress-every", type=int, default=20000, help="Progress print frequency")

    p.add_argument(
        "--skip-existing",
        dest="skip_existing",
        action="store_true",
        help="Skip inchikeys already present in output table (default)",
    )
    p.add_argument(
        "--no-skip-existing",
        dest="skip_existing",
        action="store_false",
        help="Recompute even if already present",
    )
    p.set_defaults(skip_existing=True)

    return p.parse_args(argv)


def fetch_rows(conn: sqlite3.Connection, source_table: str, out_table: str, skip_existing: bool) -> Iterable[Tuple[str, str]]:
    if skip_existing:
        sql = f"""
        SELECT s.inchikey, s.smiles
        FROM {source_table} AS s
        LEFT JOIN {out_table} AS o
          ON s.inchikey = o.inchikey
        WHERE s.inchikey IS NOT NULL
          AND s.smiles IS NOT NULL
          AND o.inchikey IS NULL
        """
    else:
        sql = f"SELECT inchikey, smiles FROM {source_table} WHERE inchikey IS NOT NULL AND smiles IS NOT NULL"

    cur = conn.execute(sql)
    while True:
        rows = cur.fetchmany(5000)
        if not rows:
            break
        for inchikey, smiles in rows:
            yield inchikey, smiles


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    try:
        from rdkit import Chem
        from rdkit import RDLogger
        from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
    except Exception as e:
        print("ERROR: RDKit not available. Create a venv and install rdkit first.", file=sys.stderr)
        print(f"Import error: {type(e).__name__}: {e}", file=sys.stderr)
        return 2

    RDLogger.DisableLog("rdApp.*")

    run_id = dt.datetime.now(dt.timezone.utc).strftime("m2_%Y%m%dT%H%M%SZ")
    started_at = utc_now_iso()

    conn = connect_sqlite(args.db)
    try:
        create_schema(conn, args.table)

        existing_before = conn.execute(f"SELECT COUNT(*) FROM {args.table}").fetchone()[0]

        n_selected = 0
        n_ok = 0
        n_fail = 0

        to_insert: List[Tuple[Any, ...]] = []
        to_reject: List[Tuple[Any, ...]] = []

        computed_at = utc_now_iso()
        t0 = time.time()

        def mol_from_smiles(smiles_str: str) -> Tuple[Any, int, str]:
            try:
                mol = Chem.MolFromSmiles(smiles_str, sanitize=True)
                if mol is not None:
                    return mol, 1, ""
            except Exception as e:
                primary_error = f"{type(e).__name__}: {e}"
            else:
                primary_error = "MolFromSmiles returned None"

            mol2 = Chem.MolFromSmiles(smiles_str, sanitize=False)
            if mol2 is None:
                raise ValueError(primary_error)

            try:
                Chem.SanitizeMol(mol2)
                return mol2, 0, primary_error
            except Exception as e:
                try:
                    ops = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
                    Chem.SanitizeMol(mol2, sanitizeOps=ops)
                    return mol2, 0, primary_error + f" | sanitize_no_kekulize: {type(e).__name__}: {e}"
                except Exception as e2:
                    raise ValueError(primary_error + f" | sanitize_failed: {type(e2).__name__}: {e2}")

        for inchikey, smiles in fetch_rows(conn, args.source_table, args.table, args.skip_existing):
            if args.limit and n_selected >= args.limit:
                break
            n_selected += 1

            try:
                mol, sanitize_ok, warn = mol_from_smiles(smiles)

                rdkit_smiles = Chem.MolToSmiles(mol, canonical=True)

                row = (
                    inchikey,
                    rdkit_smiles,
                    float(Descriptors.MolWt(mol)),
                    float(Crippen.MolLogP(mol)),
                    float(rdMolDescriptors.CalcTPSA(mol)),
                    int(Lipinski.NumHDonors(mol)),
                    int(Lipinski.NumHAcceptors(mol)),
                    int(Lipinski.NumRotatableBonds(mol)),
                    int(rdMolDescriptors.CalcNumRings(mol)),
                    int(Lipinski.HeavyAtomCount(mol)),
                    float(rdMolDescriptors.CalcFractionCSP3(mol)),
                    int(sanitize_ok),
                    (warn or None),
                    computed_at,
                )
                to_insert.append(row)
                n_ok += 1
            except Exception as e:
                to_reject.append((inchikey, smiles[:500], f"{type(e).__name__}: {e}", computed_at))
                n_fail += 1

            if len(to_insert) + len(to_reject) >= args.batch_size:
                conn.execute("BEGIN")
                if to_insert:
                    conn.executemany(
                        f"""
                        INSERT OR REPLACE INTO {args.table}(
                            inchikey, rdkit_canonical_smiles,
                            mol_wt, mol_logp, tpsa,
                            hbd, hba, rot_bonds, rings, heavy_atoms, frac_csp3,
                            rdkit_sanitize_ok, rdkit_error, computed_at_utc
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        to_insert,
                    )
                if to_reject:
                    conn.executemany(
                        f"INSERT INTO {args.table}_rejects(inchikey, smiles, error, computed_at_utc) VALUES (?, ?, ?, ?)",
                        to_reject,
                    )
                conn.execute("COMMIT")
                to_insert.clear()
                to_reject.clear()

            if args.progress_every and n_selected % args.progress_every == 0:
                elapsed = time.time() - t0
                rate = n_selected / elapsed if elapsed > 0 else 0
                print(
                    f"[{run_id}] selected={n_selected:,} ok={n_ok:,} fail={n_fail:,} ({rate:,.0f}/s)",
                    file=sys.stderr,
                    flush=True,
                )

        if to_insert or to_reject:
            conn.execute("BEGIN")
            if to_insert:
                conn.executemany(
                    f"""
                    INSERT OR REPLACE INTO {args.table}(
                        inchikey, rdkit_canonical_smiles,
                        mol_wt, mol_logp, tpsa,
                        hbd, hba, rot_bonds, rings, heavy_atoms, frac_csp3,
                        rdkit_sanitize_ok, rdkit_error, computed_at_utc
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    to_insert,
                )
            if to_reject:
                conn.executemany(
                    f"INSERT INTO {args.table}_rejects(inchikey, smiles, error, computed_at_utc) VALUES (?, ?, ?, ?)",
                    to_reject,
                )
            conn.execute("COMMIT")

        conn.execute(
            f"""
            CREATE INDEX IF NOT EXISTS idx_{args.table}_rejects_inchikey
            ON {args.table}_rejects(inchikey);
            """
        )
        conn.execute(
            f"DELETE FROM {args.table}_rejects WHERE inchikey IN (SELECT inchikey FROM {args.table});"
        )
        conn.commit()

        finished_at = utc_now_iso()
        elapsed = time.time() - t0

        qc: Dict[str, Any] = {
            "run_id": run_id,
            "started_at_utc": started_at,
            "finished_at_utc": finished_at,
            "input": {"db": display_path(args.db), "source_table": args.source_table},
            "output": {"db": display_path(args.db), "table": args.table, "outdir": display_path(args.outdir)},
            "parameters": {"limit": args.limit, "batch_size": args.batch_size, "skip_existing": bool(args.skip_existing)},
            "counts": {"existing_before": existing_before, "selected": n_selected, "ok": n_ok, "fail": n_fail},
            "runtime": {"seconds": round(elapsed, 2), "rows_per_second": round((n_selected / elapsed) if elapsed else 0, 2)},
        }

        qc_json_run = os.path.join(args.outdir, f"qc_m2_rdkit_{run_id}.json")
        qc_md_run = os.path.join(args.outdir, f"qc_m2_rdkit_{run_id}.md")
        qc_json_latest = os.path.join(args.outdir, "qc_m2_rdkit_latest.json")
        qc_md_latest = os.path.join(args.outdir, "qc_m2_rdkit_latest.md")

        for path in [qc_json_run, qc_json_latest]:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(qc, f, ensure_ascii=False, indent=2)

        for path in [qc_md_run, qc_md_latest]:
            with open(path, "w", encoding="utf-8") as f:
                f.write("# M2 RDKit QC\n\n")
                f.write(f"- run_id: `{run_id}`\n")
                f.write(f"- started_at_utc: `{started_at}`\n")
                f.write(f"- finished_at_utc: `{finished_at}`\n")
                f.write(f"- db: `{display_path(args.db)}`\n")
                f.write(f"- source_table: `{args.source_table}`\n")
                f.write(f"- output_table: `{args.table}`\n\n")
                f.write("## Counts\n\n")
                f.write("| metric | value |\n|---|---:|\n")
                for k, v in qc["counts"].items():
                    f.write(f"| {k} | {v} |\n")
                f.write("\n## Runtime\n\n")
                f.write("| metric | value |\n|---|---:|\n")
                for k, v in qc["runtime"].items():
                    f.write(f"| {k} | {v} |\n")

        print(
            f"DONE {run_id}: existing_before={existing_before:,} selected={n_selected:,} ok={n_ok:,} fail={n_fail:,} elapsed={elapsed:.1f}s",
            file=sys.stderr,
            flush=True,
        )
        return 0

    finally:
        conn.close()


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
