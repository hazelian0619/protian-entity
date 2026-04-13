#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List


THRESH_EDGE_JOIN = 0.99
THRESH_ACTIVITY_TYPE = 0.95
THRESH_ASSAY = 0.90
THRESH_STRUCTURE = 0.20


def utc_now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def non_empty(v: object) -> bool:
    return str(v or "").strip() not in {"", "NA", "N/A", "None", "null"}


def safe_rate(num: int, den: int) -> float:
    if den <= 0:
        return 1.0
    return num / den


@dataclass
class Gate:
    gate_id: str
    passed: bool
    detail: Dict[str, Any]


def load_rows_into_temp(
    conn: sqlite3.Connection,
    tsv_path: Path,
    table_name: str,
    edge_col: str = "edge_id",
) -> Dict[str, int]:
    conn.execute(f"DROP TABLE IF EXISTS {table_name}")
    conn.execute(f"CREATE TEMP TABLE {table_name}(edge_id TEXT)")
    conn.execute(f"CREATE INDEX idx_{table_name}_edge_id ON {table_name}(edge_id)")

    rows_total = 0
    rows_with_edge = 0
    batch = []

    with tsv_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or edge_col not in reader.fieldnames:
            raise SystemExit(f"[ERROR] {tsv_path} missing column: {edge_col}")

        for row in reader:
            rows_total += 1
            eid = (row.get(edge_col) or "").strip()
            if eid:
                rows_with_edge += 1
                batch.append((eid,))
            if len(batch) >= 20000:
                conn.executemany(f"INSERT INTO {table_name}(edge_id) VALUES (?)", batch)
                batch.clear()

    if batch:
        conn.executemany(f"INSERT INTO {table_name}(edge_id) VALUES (?)", batch)

    return {"rows_total": rows_total, "rows_with_edge": rows_with_edge}


def scan_activity(tsv_path: Path) -> Dict[str, int]:
    out = {
        "rows_total": 0,
        "standard_type_non_empty": 0,
        "assay_field_non_empty": 0,
        "relation_allowed": 0,
    }
    allowed_rel = {"=", "<", ">", "<=", ">=", "~"}

    with tsv_path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {"standard_type", "assay_type", "assay_description", "bao_format", "standard_relation"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] activity context missing columns: {sorted(missing)}")

        for row in r:
            out["rows_total"] += 1
            if non_empty(row.get("standard_type")):
                out["standard_type_non_empty"] += 1
            if (
                non_empty(row.get("assay_type"))
                or non_empty(row.get("assay_description"))
                or non_empty(row.get("bao_format"))
            ):
                out["assay_field_non_empty"] += 1
            if str(row.get("standard_relation") or "").strip() in allowed_rel:
                out["relation_allowed"] += 1

    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--psi-db", type=Path, required=True)
    ap.add_argument("--activity-context", type=Path, required=True)
    ap.add_argument("--structure-evidence", type=Path, required=True)
    ap.add_argument("--activity-validation", type=Path, required=True)
    ap.add_argument("--structure-validation", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    for p in [
        args.psi_db,
        args.activity_context,
        args.structure_evidence,
        args.activity_validation,
        args.structure_validation,
    ]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    activity_stats = scan_activity(args.activity_context)

    conn = sqlite3.connect(f"file:{args.psi_db.resolve()}?mode=ro", uri=True)
    total_edges = int(conn.execute("SELECT COUNT(*) FROM psi_edges_v1").fetchone()[0])

    activity_load = load_rows_into_temp(conn, args.activity_context, "tmp_activity")
    structure_load = load_rows_into_temp(conn, args.structure_evidence, "tmp_structure")

    mapped_activity = int(
        conn.execute(
            "SELECT COUNT(*) FROM tmp_activity a JOIN psi_edges_v1 e ON e.edge_id = a.edge_id"
        ).fetchone()[0]
    )
    mapped_structure = int(
        conn.execute(
            "SELECT COUNT(*) FROM tmp_structure s JOIN psi_edges_v1 e ON e.edge_id = s.edge_id"
        ).fetchone()[0]
    )
    structure_distinct_edges = int(conn.execute("SELECT COUNT(DISTINCT edge_id) FROM tmp_structure").fetchone()[0])

    conn.close()

    activity_rows = activity_load["rows_total"]
    structure_rows = structure_load["rows_total"]

    edge_join_rate_activity = safe_rate(mapped_activity, activity_rows)
    edge_join_rate_structure = safe_rate(mapped_structure, structure_rows)
    activity_type_coverage = safe_rate(activity_stats["standard_type_non_empty"], activity_stats["rows_total"])
    assay_field_coverage = safe_rate(activity_stats["assay_field_non_empty"], activity_stats["rows_total"])
    relation_symbol_coverage = safe_rate(activity_stats["relation_allowed"], activity_stats["rows_total"])
    structure_edge_coverage = safe_rate(structure_distinct_edges, total_edges)

    activity_validation = json.loads(args.activity_validation.read_text(encoding="utf-8"))
    structure_validation = json.loads(args.structure_validation.read_text(encoding="utf-8"))
    contracts_pass = bool(activity_validation.get("passed")) and bool(structure_validation.get("passed"))

    gates: List[Gate] = [
        Gate(
            gate_id="edge_id_join_rate_activity",
            passed=edge_join_rate_activity >= THRESH_EDGE_JOIN,
            detail={
                "rate": edge_join_rate_activity,
                "threshold": THRESH_EDGE_JOIN,
                "mapped": mapped_activity,
                "rows": activity_rows,
            },
        ),
        Gate(
            gate_id="edge_id_join_rate_structure",
            passed=edge_join_rate_structure >= THRESH_EDGE_JOIN,
            detail={
                "rate": edge_join_rate_structure,
                "threshold": THRESH_EDGE_JOIN,
                "mapped": mapped_structure,
                "rows": structure_rows,
            },
        ),
        Gate(
            gate_id="activity_type_coverage",
            passed=activity_type_coverage >= THRESH_ACTIVITY_TYPE,
            detail={
                "rate": activity_type_coverage,
                "threshold": THRESH_ACTIVITY_TYPE,
                "non_empty": activity_stats["standard_type_non_empty"],
                "rows": activity_stats["rows_total"],
            },
        ),
        Gate(
            gate_id="assay_field_coverage",
            passed=assay_field_coverage >= THRESH_ASSAY,
            detail={
                "rate": assay_field_coverage,
                "threshold": THRESH_ASSAY,
                "non_empty": activity_stats["assay_field_non_empty"],
                "rows": activity_stats["rows_total"],
            },
        ),
        Gate(
            gate_id="structure_edge_coverage",
            passed=structure_edge_coverage >= THRESH_STRUCTURE,
            detail={
                "rate": structure_edge_coverage,
                "threshold": THRESH_STRUCTURE,
                "distinct_edges": structure_distinct_edges,
                "psi_total_edges": total_edges,
            },
        ),
        Gate(
            gate_id="contracts_validation_pass",
            passed=contracts_pass,
            detail={
                "activity_contract_pass": bool(activity_validation.get("passed")),
                "structure_contract_pass": bool(structure_validation.get("passed")),
            },
        ),
    ]

    passed = all(g.passed for g in gates)

    report = {
        "name": "psi_activity_structure_enrichment_v2.qa",
        "created_at": utc_now_iso(),
        "inputs": {
            "psi_db": str(args.psi_db),
            "activity_context": str(args.activity_context),
            "structure_evidence": str(args.structure_evidence),
            "activity_validation": str(args.activity_validation),
            "structure_validation": str(args.structure_validation),
        },
        "metrics": {
            "psi_total_edges": total_edges,
            "activity_rows": activity_rows,
            "structure_rows": structure_rows,
            "edge_id_join_rate_activity": edge_join_rate_activity,
            "edge_id_join_rate_structure": edge_join_rate_structure,
            "activity_type_coverage": activity_type_coverage,
            "assay_field_coverage": assay_field_coverage,
            "relation_symbol_coverage": relation_symbol_coverage,
            "structure_edge_coverage": structure_edge_coverage,
            "contracts_pass": contracts_pass,
        },
        "gates": [
            {"id": g.gate_id, "passed": g.passed, "detail": g.detail}
            for g in gates
        ],
        "passed": passed,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if passed else "FAIL"
    print(f"[{status}] psi_activity_structure_enrichment_v2 QA -> {args.out}")
    print(
        "[METRICS] "
        f"edge_join_activity={edge_join_rate_activity:.4f}, "
        f"edge_join_structure={edge_join_rate_structure:.4f}, "
        f"activity_type_cov={activity_type_coverage:.4f}, "
        f"assay_cov={assay_field_coverage:.4f}, "
        f"structure_cov={structure_edge_coverage:.4f}, "
        f"contracts={contracts_pass}"
    )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
