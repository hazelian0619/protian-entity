#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


@dataclass(frozen=True)
class GateResult:
    ok: bool
    details: Dict[str, Any]


def _connect_ro(db_path: Path) -> sqlite3.Connection:
    uri = f"file:{db_path.resolve()}?mode=ro"
    conn = sqlite3.connect(uri, uri=True)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA query_only=ON;")
    return conn


def _table_exists(conn: sqlite3.Connection, name: str, schema: str = "main") -> bool:
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type IN ('table','view') AND name=? LIMIT 1", (name,)
    ).fetchone()
    if row is not None:
        return True
    if schema != "main":
        row = conn.execute(
            f"SELECT 1 FROM {schema}.sqlite_master WHERE type IN ('table','view') AND name=? LIMIT 1",
            (name,),
        ).fetchone()
        return row is not None
    return False


def _scalar(conn: sqlite3.Connection, sql: str, args: Tuple[Any, ...] = ()) -> Any:
    row = conn.execute(sql, args).fetchone()
    return row[0] if row else None


def _missing_rate(conn: sqlite3.Connection, table: str, col: str) -> float:
    # Missing = NULL or empty string.
    return float(
        _scalar(
            conn,
            f"""
            SELECT AVG(CASE WHEN {col} IS NULL OR {col}='' THEN 1.0 ELSE 0.0 END)
            FROM {table}
            """.strip(),
        )
        or 0.0
    )


def _dup_count(conn: sqlite3.Connection, table: str, col: str) -> int:
    v = _scalar(conn, f"SELECT COUNT(*) - COUNT(DISTINCT {col}) FROM {table}")
    return int(v or 0)


def _disallowed_values(conn: sqlite3.Connection, table: str, col: str, allowed: Iterable[str]) -> List[Tuple[str, int]]:
    allowed_list = [a for a in allowed if a]
    if not allowed_list:
        return []
    placeholders = ",".join(["?"] * len(allowed_list))
    rows = conn.execute(
        f"SELECT {col} AS v, COUNT(*) AS n FROM {table} WHERE {col} NOT IN ({placeholders}) GROUP BY {col} ORDER BY n DESC",
        tuple(allowed_list),
    ).fetchall()
    return [(str(r["v"]) if r["v"] is not None else "(NULL)", int(r["n"])) for r in rows]


def run_gates(
    *,
    m3_db: Path,
    m1_db: Optional[Path],
    min_inchikey_rate: float,
    min_uniprot_rate: float,
    standard_types: List[str],
    check_join_coverage: bool,
    refs_sample: int,
    deep_refs: bool,
) -> GateResult:
    details: Dict[str, Any] = {"m3_db": str(m3_db)}

    if not m3_db.exists():
        return GateResult(False, {**details, "error": "m3_db_not_found"})

    conn = _connect_ro(m3_db)
    try:
        required_tables = ["psi_evidence_v1", "psi_edges_v1", "target_uniprot_map_v1"]
        missing_tables = [t for t in required_tables if not _table_exists(conn, t)]
        details["tables"] = {"required": required_tables, "missing": missing_tables}
        if missing_tables:
            return GateResult(False, details)

        evidence_rows = int(_scalar(conn, "SELECT COUNT(*) FROM psi_evidence_v1") or 0)
        edges = int(_scalar(conn, "SELECT COUNT(*) FROM psi_edges_v1") or 0)
        details["counts"] = {"evidence_rows": evidence_rows, "edges": edges}
        if evidence_rows <= 0 or edges <= 0:
            return GateResult(False, details)

        missing_inchikey_rate = _missing_rate(conn, "psi_evidence_v1", "compound_inchikey")
        missing_uniprot_rate = _missing_rate(conn, "psi_evidence_v1", "target_uniprot_accession")
        details["rates"] = {
            "inchikey_present_rate": 1.0 - missing_inchikey_rate,
            "uniprot_present_rate": 1.0 - missing_uniprot_rate,
        }
        if (1.0 - missing_inchikey_rate) < min_inchikey_rate:
            details["rates"]["inchikey_gate"] = {
                "min_required": min_inchikey_rate,
                "ok": False,
            }
            return GateResult(False, details)
        if (1.0 - missing_uniprot_rate) < min_uniprot_rate:
            details["rates"]["uniprot_gate"] = {"min_required": min_uniprot_rate, "ok": False}
            return GateResult(False, details)

        evidence_dups = _dup_count(conn, "psi_evidence_v1", "activity_id")
        edges_dups = _dup_count(conn, "psi_edges_v1", "edge_id")
        details["uniqueness"] = {
            "activity_id_duplicates": evidence_dups,
            "edge_id_duplicates": edges_dups,
        }
        if evidence_dups != 0 or edges_dups != 0:
            return GateResult(False, details)

        disallowed = _disallowed_values(conn, "psi_evidence_v1", "standard_type", standard_types)
        details["standard_type"] = {"allowed": standard_types, "disallowed": disallowed}
        if disallowed:
            return GateResult(False, details)

        # Referential integrity: edges.best_activity_id -> evidence.activity_id
        # Full verification can be expensive on multi-million edge tables, so the
        # default is a range check + sample join.
        evidence_min = _scalar(conn, "SELECT MIN(activity_id) FROM psi_evidence_v1")
        evidence_max = _scalar(conn, "SELECT MAX(activity_id) FROM psi_evidence_v1")
        edge_min = _scalar(conn, "SELECT MIN(best_activity_id) FROM psi_edges_v1")
        edge_max = _scalar(conn, "SELECT MAX(best_activity_id) FROM psi_edges_v1")
        refs: Dict[str, Any] = {
            "evidence_activity_id_min": int(evidence_min) if evidence_min is not None else None,
            "evidence_activity_id_max": int(evidence_max) if evidence_max is not None else None,
            "edges_best_activity_id_min": int(edge_min) if edge_min is not None else None,
            "edges_best_activity_id_max": int(edge_max) if edge_max is not None else None,
        }
        if (
            evidence_min is not None
            and edge_min is not None
            and int(edge_min) < int(evidence_min)
            or evidence_max is not None
            and edge_max is not None
            and int(edge_max) > int(evidence_max)
        ):
            refs["range_check_ok"] = False
            details["refs"] = refs
            return GateResult(False, details)
        refs["range_check_ok"] = True

        sample_n = max(0, int(refs_sample or 0))
        if sample_n > 0:
            missing_sample = int(
                _scalar(
                    conn,
                    f"""
                    SELECT COUNT(*)
                    FROM (SELECT best_activity_id FROM psi_edges_v1 LIMIT {sample_n}) s
                    LEFT JOIN psi_evidence_v1 ev ON ev.activity_id = s.best_activity_id
                    WHERE ev.activity_id IS NULL
                    """.strip(),
                )
                or 0
            )
            refs["sample_n"] = sample_n
            refs["missing_best_activity_id_in_sample"] = missing_sample
            if missing_sample != 0:
                details["refs"] = refs
                return GateResult(False, details)

        if deep_refs:
            missing_best = int(
                _scalar(
                    conn,
                    """
                    SELECT COUNT(*)
                    FROM psi_edges_v1 e
                    LEFT JOIN psi_evidence_v1 ev ON ev.activity_id = e.best_activity_id
                    WHERE ev.activity_id IS NULL
                    """.strip(),
                )
                or 0
            )
            refs["missing_best_activity_id"] = missing_best
            if missing_best != 0:
                details["refs"] = refs
                return GateResult(False, details)

        details["refs"] = refs

        if check_join_coverage and m1_db is not None:
            join_details: Dict[str, Any] = {"m1_db": str(m1_db)}
            if not m1_db.exists():
                join_details["skipped"] = "m1_db_not_found"
            else:
                conn.execute(f"ATTACH DATABASE 'file:{m1_db.resolve()}?mode=ro' AS m1")
                if not _table_exists(conn, "molecule_props_rdkit_v1", schema="m1"):
                    join_details["skipped"] = "m1_table_not_found"
                else:
                    matched = int(
                        _scalar(
                            conn,
                            """
                            SELECT COUNT(*)
                            FROM psi_evidence_v1 e
                            JOIN m1.molecule_props_rdkit_v1 p ON p.inchikey = e.compound_inchikey
                            """.strip(),
                        )
                        or 0
                    )
                    join_details["join_props_rate"] = matched / float(evidence_rows)
            details["join"] = join_details

        return GateResult(True, details)
    finally:
        conn.close()


def main() -> int:
    ap = argparse.ArgumentParser(description="Validate M3 (ChEMBL PSI) SQLite against quality gates")
    ap.add_argument(
        "--m3-db",
        type=Path,
        default=Path("data/output/molecules/chembl_m3.sqlite"),
        help="M3 output sqlite (default: data/output/molecules/chembl_m3.sqlite)",
    )
    ap.add_argument(
        "--m1-db",
        type=Path,
        default=Path("data/output/molecules/molecules_m1.sqlite"),
        help="Optional M1 sqlite for join coverage (default: data/output/molecules/molecules_m1.sqlite)",
    )
    ap.add_argument(
        "--skip-join-coverage",
        action="store_true",
        help="Skip read-only join coverage check against M1",
    )
    ap.add_argument("--min-inchikey-rate", type=float, default=0.95)
    ap.add_argument("--min-uniprot-rate", type=float, default=0.95)
    ap.add_argument(
        "--types",
        type=str,
        default="IC50,Ki,Kd,EC50",
        help="Comma-separated allowlist of standard_type values",
    )
    ap.add_argument(
        "--refs-sample",
        type=int,
        default=1000,
        help="Sample size for best_activity_id reference check (0 disables sampling)",
    )
    ap.add_argument(
        "--deep-refs",
        action="store_true",
        help="Run a full edge->evidence reference scan (can be slow)",
    )
    ap.add_argument("--json-out", type=Path, default=None, help="Optional JSON report output path")
    args = ap.parse_args()

    standard_types = [t.strip() for t in (args.types or "").split(",") if t.strip()]
    result = run_gates(
        m3_db=args.m3_db,
        m1_db=None if args.skip_join_coverage else args.m1_db,
        min_inchikey_rate=args.min_inchikey_rate,
        min_uniprot_rate=args.min_uniprot_rate,
        standard_types=standard_types,
        check_join_coverage=(not args.skip_join_coverage),
        refs_sample=args.refs_sample,
        deep_refs=args.deep_refs,
    )

    payload: Dict[str, Any] = {"ok": result.ok, **result.details}
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    stream = sys.stdout if result.ok else sys.stderr
    stream.write(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    return 0 if result.ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
