#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Sanity-check a locally downloaded ChEMBL SQLite dump.

This is a lightweight pre-flight check before running M3 ETL.
It does NOT modify the input DB.

Usage:
  python3 pipelines/molecules/scripts/m3_check_chembl_db.py --chembl-db data/raw/molecules/chembl_36/chembl_36.db
"""

import argparse
import datetime as dt
import hashlib
import os
import sqlite3
import sys
from typing import List


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


def sha256_head(path: str, nbytes: int = 10 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read(nbytes))
    return h.hexdigest()


def parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Check ChEMBL SQLite presence and required tables")
    p.add_argument(
        "--chembl-db",
        default="data/raw/molecules/chembl_36/chembl_36.db",
        help="Path to ChEMBL SQLite DB (default: data/raw/molecules/chembl_36/chembl_36.db)",
    )
    return p.parse_args(argv)


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    path = args.chembl_db

    if not os.path.exists(path):
        print(f"ERROR: not found: {path}", file=sys.stderr)
        print("Expected a locally downloaded ChEMBL SQLite dump (.db).", file=sys.stderr)
        return 2

    size = os.path.getsize(path)
    mtime = dt.datetime.fromtimestamp(os.path.getmtime(path), dt.timezone.utc).replace(microsecond=0)
    fp = sha256_head(path)

    uri = f"file:{os.path.abspath(path)}?mode=ro"
    conn = sqlite3.connect(uri, uri=True)
    try:
        conn.execute("PRAGMA query_only=ON;")
        tables = [r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name").fetchall()]

        required = [
            "activities",
            "assays",
            "docs",
            "molecule_dictionary",
            "compound_structures",
            "target_dictionary",
            "target_components",
            "component_sequences",
        ]
        missing = [t for t in required if t not in set(tables)]

        print("ChEMBL SQLite check")
        print(f"- checked_at_utc: {utc_now_iso()}")
        print(f"- path: {display_path(path)}")
        print(f"- size_bytes: {size}")
        print(f"- mtime_utc: {mtime.isoformat().replace('+00:00','Z')}")
        print(f"- sha256_head_10mb: {fp}")
        print(f"- n_tables: {len(tables)}")

        if missing:
            print("- status: FAIL")
            print("- missing_tables:")
            for t in missing:
                print(f"  - {t}")
            return 3

        print("- status: OK")
        return 0

    finally:
        conn.close()


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
