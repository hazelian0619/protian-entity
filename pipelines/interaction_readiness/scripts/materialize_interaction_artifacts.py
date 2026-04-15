#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def ensure_dirs(paths: Iterable[Path]) -> None:
    for p in paths:
        p.mkdir(parents=True, exist_ok=True)


def sha1_short(text: str, n: int = 20) -> str:
    return hashlib.sha1(text.encode("utf-8")).hexdigest()[:n]


def read_tsv_header(path: Path) -> List[str]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.reader(f, delimiter="\t")
        return next(r)


def write_tsv(path: Path, fieldnames: List[str], rows: Iterable[Dict[str, str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow({k: (row.get(k, "") or "") for k in fieldnames})
            n += 1
    return n


def normalize_bool_text(v: str, default: bool) -> str:
    t = (v or "").strip().lower()
    if t in {"1", "true", "t", "yes", "y"}:
        return "true"
    if t in {"0", "false", "f", "no", "n"}:
        return "false"
    return "true" if default else "false"


def resolve_source_paths(source_root: Path, m3_sqlite: Optional[Path]) -> Dict[str, Path]:
    candidates = {
        "ppi_edges": [
            source_root / "pipelines/edges_ppi/data/output/edges/edges_ppi_v1.tsv",
            source_root / "data/output/edges/edges_ppi_v1.tsv",
        ],
        "ppi_evidence": [
            source_root / "pipelines/edges_ppi/data/output/evidence/ppi_evidence_v1.tsv",
            source_root / "data/output/evidence/ppi_evidence_v1.tsv",
        ],
        "psi_edges": [
            source_root / "data/output/edges/drug_target_edges_v1.tsv",
        ],
        "psi_evidence": [
            source_root / "data/output/evidence/drug_target_evidence_v1.tsv",
        ],
        "m3_sqlite": [
            m3_sqlite,
            source_root / "data/output/molecules/chembl_m3.sqlite",
            source_root / "dist/molecules/chembl_m3.sqlite",
            Path("/Users/pluviophile/graph/12182/out/m3/chembl_m3.sqlite"),
        ],
    }

    resolved: Dict[str, Path] = {}
    missing: List[str] = []
    for key, arr in candidates.items():
        hit = None
        for p in arr:
            if p is None:
                continue
            pp = Path(p)
            if pp.exists():
                hit = pp.resolve()
                break
        if hit is None:
            missing.append(key)
        else:
            resolved[key] = hit

    if missing:
        raise FileNotFoundError(f"missing source artifacts: {missing}")
    return resolved


def transform_ppi_edges(src: Path, dst: Path, limit: Optional[int]) -> Dict[str, Any]:
    cols = [
        "edge_id",
        "src_type",
        "src_id",
        "dst_type",
        "dst_id",
        "predicate",
        "directed",
        "best_score",
        "source",
        "source_version",
        "fetch_date",
    ]

    def rows() -> Iterable[Dict[str, str]]:
        with src.open("r", encoding="utf-8", newline="") as f:
            r = csv.DictReader(f, delimiter="\t")
            n = 0
            for row in r:
                n += 1
                if limit is not None and n > limit:
                    break
                yield {
                    "edge_id": (row.get("edge_id", "") or "").strip(),
                    "src_type": "Protein",
                    "src_id": (row.get("src_id", "") or "").strip(),
                    "dst_type": "Protein",
                    "dst_id": (row.get("dst_id", "") or "").strip(),
                    "predicate": "interacts_with",
                    "directed": normalize_bool_text(row.get("directed", ""), default=False),
                    "best_score": (row.get("best_score", "") or "").strip(),
                    "source": (row.get("source", "") or "STRING").strip(),
                    "source_version": (row.get("source_version", "") or "unknown").strip(),
                    "fetch_date": (row.get("fetch_date", "") or "").strip(),
                }

    row_count = write_tsv(dst, cols, rows())
    return {"path": str(dst), "rows": row_count, "columns": cols}


def transform_ppi_evidence(src: Path, dst: Path, limit: Optional[int]) -> Dict[str, Any]:
    cols = [
        "evidence_id",
        "edge_id",
        "evidence_type",
        "method",
        "score",
        "reference",
        "source",
        "source_version",
        "fetch_date",
    ]

    method_filled = 0
    ref_or_method = 0

    def rows() -> Iterable[Dict[str, str]]:
        nonlocal method_filled, ref_or_method
        with src.open("r", encoding="utf-8", newline="") as f:
            r = csv.DictReader(f, delimiter="\t")
            n = 0
            for row in r:
                n += 1
                if limit is not None and n > limit:
                    break
                ev_type = "protein_protein_interaction"
                method = (row.get("evidence_type", "") or "").strip() or "string_combined_score"
                reference = (row.get("reference", "") or "").strip()
                if method:
                    method_filled += 1
                if method or reference:
                    ref_or_method += 1
                evid = (row.get("evidence_id", "") or "").strip() or f"ev_{sha1_short(str(n))}"
                yield {
                    "evidence_id": evid,
                    "edge_id": (row.get("edge_id", "") or "").strip(),
                    "evidence_type": ev_type,
                    "method": method,
                    "score": (row.get("score", "") or "").strip(),
                    "reference": reference,
                    "source": (row.get("source", "") or "STRING").strip(),
                    "source_version": (row.get("source_version", "") or "unknown").strip(),
                    "fetch_date": (row.get("fetch_date", "") or "").strip(),
                }

    row_count = write_tsv(dst, cols, rows())
    return {
        "path": str(dst),
        "rows": row_count,
        "columns": cols,
        "method_non_empty_rate": (method_filled / row_count) if row_count else 0.0,
        "method_or_reference_rate": (ref_or_method / row_count) if row_count else 0.0,
    }


def transform_psi_edges(src: Path, dst: Path, limit: Optional[int]) -> Dict[str, Any]:
    cols = [
        "edge_id",
        "src_type",
        "src_id",
        "dst_type",
        "dst_id",
        "predicate",
        "directed",
        "best_score",
        "source",
        "source_version",
        "fetch_date",
    ]

    def rows() -> Iterable[Dict[str, str]]:
        with src.open("r", encoding="utf-8", newline="") as f:
            r = csv.DictReader(f, delimiter="\t")
            n = 0
            for row in r:
                n += 1
                if limit is not None and n > limit:
                    break
                yield {
                    "edge_id": (row.get("edge_id", "") or "").strip(),
                    "src_type": "Drug",
                    "src_id": (row.get("src_id", "") or "").strip(),
                    "dst_type": "Protein",
                    "dst_id": (row.get("dst_id", "") or "").strip(),
                    "predicate": "targets",
                    "directed": normalize_bool_text(row.get("directed", ""), default=True),
                    "best_score": (row.get("best_score", "") or "").strip(),
                    "source": (row.get("source", "") or "DrugBank").strip(),
                    "source_version": (row.get("source_version", "") or "unknown").strip(),
                    "fetch_date": (row.get("fetch_date", "") or "").strip(),
                }

    row_count = write_tsv(dst, cols, rows())
    return {"path": str(dst), "rows": row_count, "columns": cols}


def build_pair_to_edge_id(edges_path: Path) -> Dict[Tuple[str, str], str]:
    m: Dict[Tuple[str, str], str] = {}
    with edges_path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            pair = ((row.get("src_id", "") or "").strip(), (row.get("dst_id", "") or "").strip())
            eid = (row.get("edge_id", "") or "").strip()
            if pair[0] and pair[1] and eid and pair not in m:
                m[pair] = eid
    return m


def transform_psi_evidence(src: Path, dst: Path, edge_lookup: Dict[Tuple[str, str], str], limit: Optional[int]) -> Dict[str, Any]:
    cols = [
        "evidence_id",
        "edge_id",
        "evidence_type",
        "method",
        "score",
        "reference",
        "source",
        "source_version",
        "fetch_date",
    ]

    missing_edge = 0
    method_or_ref = 0

    def rows() -> Iterable[Dict[str, str]]:
        nonlocal missing_edge, method_or_ref
        with src.open("r", encoding="utf-8", newline="") as f:
            r = csv.DictReader(f, delimiter="\t")
            n = 0
            for row in r:
                n += 1
                if limit is not None and n > limit:
                    break
                drug = (row.get("drugbank_id", "") or "").strip()
                prot = (row.get("uniprot_id", "") or "").strip()
                edge_id = edge_lookup.get((drug, prot), "")
                if not edge_id:
                    missing_edge += 1
                target_role = (row.get("target_role", "") or "").strip()
                action = (row.get("action", "") or "").strip()
                method = action or target_role or "database_record"
                reference = (row.get("reference", "") or "").strip()
                if method or reference:
                    method_or_ref += 1
                base_eid = (row.get("evidence_id", "") or "").strip()
                if base_eid:
                    evid = base_eid
                else:
                    evid = "ev_" + sha1_short(f"{edge_id}|{method}|{reference}|{n}")
                yield {
                    "evidence_id": evid,
                    "edge_id": edge_id,
                    "evidence_type": "drug_target_interaction",
                    "method": method,
                    "score": "",
                    "reference": reference,
                    "source": (row.get("source", "") or "DrugBank").strip(),
                    "source_version": (row.get("source_version", "") or "unknown").strip(),
                    "fetch_date": (row.get("fetch_date", "") or "").strip(),
                }

    row_count = write_tsv(dst, cols, rows())
    return {
        "path": str(dst),
        "rows": row_count,
        "columns": cols,
        "missing_edge_id_rows": missing_edge,
        "edge_join_rate": ((row_count - missing_edge) / row_count) if row_count else 0.0,
        "method_or_reference_rate": (method_or_ref / row_count) if row_count else 0.0,
    }


def export_m3_psi_edges(sqlite_path: Path, dst: Path, limit: Optional[int]) -> Dict[str, Any]:
    con = sqlite3.connect(f"file:{sqlite_path}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row

    table = "psi_edges_v1"
    cols = [r[1] for r in con.execute(f"PRAGMA table_info({table})").fetchall()]
    if not cols:
        raise RuntimeError(f"table {table} not found in {sqlite_path}")

    sql = f"SELECT {', '.join(cols)} FROM {table}"
    params: Tuple[Any, ...] = ()
    if limit is not None:
        sql += " LIMIT ?"
        params = (limit,)

    dst.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with dst.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for row in con.execute(sql, params):
            w.writerow({c: row[c] for c in cols})
            n += 1

    total_rows = con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    con.close()
    return {
        "path": str(dst),
        "rows_written": n,
        "table_total_rows": int(total_rows),
        "columns": cols,
        "source_sqlite": str(sqlite_path),
    }


def write_contracts(repo_root: Path) -> Dict[str, str]:
    contracts: Dict[str, Dict[str, Any]] = {
        "pipelines/edges_ppi/contracts/edges_ppi_v1.json": {
            "name": "edges_ppi_v1",
            "required_columns": [
                "edge_id",
                "src_type",
                "src_id",
                "dst_type",
                "dst_id",
                "predicate",
                "directed",
                "best_score",
                "source",
                "source_version",
                "fetch_date",
            ],
            "rules": [
                {"id": "edge_id_non_empty", "type": "non_empty_rate", "column": "edge_id", "min_rate": 1.0},
                {"id": "edge_id_unique", "type": "unique", "column": "edge_id"},
                {"id": "src_type_protein", "type": "equals_value_rate", "column": "src_type", "value": "Protein", "min_rate": 1.0},
                {"id": "dst_type_protein", "type": "equals_value_rate", "column": "dst_type", "value": "Protein", "min_rate": 1.0},
                {"id": "predicate_non_empty", "type": "non_empty_rate", "column": "predicate", "min_rate": 1.0},
                {"id": "directed_bool", "type": "regex_rate", "column": "directed", "pattern": "^(true|false)$", "min_rate": 1.0},
                {"id": "source_non_empty", "type": "non_empty_rate", "column": "source", "min_rate": 1.0},
                {"id": "source_version_non_empty", "type": "non_empty_rate", "column": "source_version", "min_rate": 1.0},
                {"id": "fetch_date_iso", "type": "regex_rate", "column": "fetch_date", "pattern": "^\\d{4}-\\d{2}-\\d{2}$", "min_rate": 1.0},
            ],
        },
        "pipelines/edges_ppi/contracts/ppi_evidence_v1.json": {
            "name": "ppi_evidence_v1",
            "required_columns": [
                "evidence_id",
                "edge_id",
                "evidence_type",
                "method",
                "score",
                "reference",
                "source",
                "source_version",
                "fetch_date",
            ],
            "rules": [
                {"id": "evidence_id_non_empty", "type": "non_empty_rate", "column": "evidence_id", "min_rate": 1.0},
                {"id": "evidence_id_unique", "type": "unique", "column": "evidence_id"},
                {"id": "edge_id_non_empty", "type": "non_empty_rate", "column": "edge_id", "min_rate": 1.0},
                {"id": "evidence_type_non_empty", "type": "non_empty_rate", "column": "evidence_type", "min_rate": 1.0},
                {"id": "method_non_empty", "type": "non_empty_rate", "column": "method", "min_rate": 1.0},
                {"id": "score_numeric_or_empty", "type": "regex_rate", "column": "score", "pattern": "^(|-?[0-9]+(\\.[0-9]+)?)$", "min_rate": 1.0},
                {"id": "source_non_empty", "type": "non_empty_rate", "column": "source", "min_rate": 1.0},
                {"id": "source_version_non_empty", "type": "non_empty_rate", "column": "source_version", "min_rate": 1.0},
                {"id": "fetch_date_iso", "type": "regex_rate", "column": "fetch_date", "pattern": "^\\d{4}-\\d{2}-\\d{2}$", "min_rate": 1.0},
            ],
        },
        "pipelines/drugbank/contracts/drug_target_edges_v1.json": {
            "name": "drug_target_edges_v1",
            "required_columns": [
                "edge_id",
                "src_type",
                "src_id",
                "dst_type",
                "dst_id",
                "predicate",
                "directed",
                "best_score",
                "source",
                "source_version",
                "fetch_date",
            ],
            "rules": [
                {"id": "edge_id_non_empty", "type": "non_empty_rate", "column": "edge_id", "min_rate": 1.0},
                {"id": "edge_id_unique", "type": "unique", "column": "edge_id"},
                {"id": "src_type_drug", "type": "equals_value_rate", "column": "src_type", "value": "Drug", "min_rate": 1.0},
                {"id": "dst_type_protein", "type": "equals_value_rate", "column": "dst_type", "value": "Protein", "min_rate": 1.0},
                {"id": "predicate_non_empty", "type": "non_empty_rate", "column": "predicate", "min_rate": 1.0},
                {"id": "directed_bool", "type": "regex_rate", "column": "directed", "pattern": "^(true|false)$", "min_rate": 1.0},
                {"id": "source_non_empty", "type": "non_empty_rate", "column": "source", "min_rate": 1.0},
                {"id": "source_version_non_empty", "type": "non_empty_rate", "column": "source_version", "min_rate": 1.0},
                {"id": "fetch_date_iso", "type": "regex_rate", "column": "fetch_date", "pattern": "^\\d{4}-\\d{2}-\\d{2}$", "min_rate": 1.0},
            ],
        },
        "pipelines/drugbank/contracts/drug_target_evidence_v1.json": {
            "name": "drug_target_evidence_v1",
            "required_columns": [
                "evidence_id",
                "edge_id",
                "evidence_type",
                "method",
                "score",
                "reference",
                "source",
                "source_version",
                "fetch_date",
            ],
            "rules": [
                {"id": "evidence_id_non_empty", "type": "non_empty_rate", "column": "evidence_id", "min_rate": 1.0},
                {"id": "evidence_id_unique", "type": "unique", "column": "evidence_id"},
                {"id": "edge_id_non_empty", "type": "non_empty_rate", "column": "edge_id", "min_rate": 1.0},
                {"id": "evidence_type_non_empty", "type": "non_empty_rate", "column": "evidence_type", "min_rate": 1.0},
                {"id": "method_non_empty", "type": "non_empty_rate", "column": "method", "min_rate": 1.0},
                {"id": "score_numeric_or_empty", "type": "regex_rate", "column": "score", "pattern": "^(|-?[0-9]+(\\.[0-9]+)?)$", "min_rate": 1.0},
                {"id": "source_non_empty", "type": "non_empty_rate", "column": "source", "min_rate": 1.0},
                {"id": "source_version_non_empty", "type": "non_empty_rate", "column": "source_version", "min_rate": 1.0},
                {"id": "fetch_date_iso", "type": "regex_rate", "column": "fetch_date", "pattern": "^\\d{4}-\\d{2}-\\d{2}$", "min_rate": 1.0},
            ],
        },
    }

    written = {}
    for rel_path, content in contracts.items():
        p = repo_root / rel_path
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(content, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        written[rel_path] = "written"
    return written


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Materialize PPI/PSI artifacts into current repo with unified schema")
    p.add_argument("--repo-root", type=Path, default=Path("."))
    p.add_argument("--source-root", type=Path, default=Path("../1218"))
    p.add_argument("--m3-sqlite", type=Path, default=None)
    p.add_argument("--mode", choices=["sample", "full"], default="full")
    p.add_argument("--sample-limit", type=int, default=200)
    p.add_argument("--report", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    source_root = args.source_root.resolve()

    limit = args.sample_limit if args.mode == "sample" else None

    output_dirs = [
        repo_root / "data/output/edges",
        repo_root / "data/output/evidence",
        repo_root / "data/output/molecules",
        repo_root / "pipelines/edges_ppi/contracts",
        repo_root / "pipelines/edges_ppi/reports",
        repo_root / "pipelines/drugbank/contracts",
        repo_root / "pipelines/drugbank/reports",
    ]
    ensure_dirs(output_dirs)

    try:
        sources = resolve_source_paths(source_root=source_root, m3_sqlite=args.m3_sqlite)
    except Exception as e:
        blocked = {
            "pipeline": "interaction_artifact_materialization",
            "generated_at_utc": now_utc(),
            "status": "blocked_missing_inputs",
            "source_root": str(source_root),
            "error": str(e),
        }
        write_json(args.report, blocked)
        print(f"[BLOCKED] {e}")
        return 2

    out_ppi_edges = repo_root / "data/output/edges/edges_ppi_v1.tsv"
    out_ppi_evidence = repo_root / "data/output/evidence/ppi_evidence_v1.tsv"
    out_psi_edges = repo_root / "data/output/edges/drug_target_edges_v1.tsv"
    out_psi_evidence = repo_root / "data/output/evidence/drug_target_evidence_v1.tsv"
    out_m3_psi = repo_root / "data/output/molecules/molecules_m3_psi_edges_v1.tsv"

    ppi_edges_info = transform_ppi_edges(sources["ppi_edges"], out_ppi_edges, limit=limit)
    ppi_ev_info = transform_ppi_evidence(sources["ppi_evidence"], out_ppi_evidence, limit=limit)

    psi_edges_info = transform_psi_edges(sources["psi_edges"], out_psi_edges, limit=limit)
    edge_lookup = build_pair_to_edge_id(out_psi_edges)
    psi_ev_info = transform_psi_evidence(sources["psi_evidence"], out_psi_evidence, edge_lookup=edge_lookup, limit=limit)

    m3_info = export_m3_psi_edges(sources["m3_sqlite"], out_m3_psi, limit=limit)

    contract_info = write_contracts(repo_root)

    report = {
        "pipeline": "interaction_artifact_materialization",
        "generated_at_utc": now_utc(),
        "status": "completed",
        "mode": args.mode,
        "source_root": str(source_root),
        "sources": {k: str(v) for k, v in sources.items()},
        "outputs": {
            "ppi_edges": ppi_edges_info,
            "ppi_evidence": ppi_ev_info,
            "psi_edges": psi_edges_info,
            "psi_evidence": psi_ev_info,
            "m3_psi_edges": m3_info,
        },
        "contracts": contract_info,
        "acceptance": {
            "ppi_method_or_reference_rate": ppi_ev_info.get("method_or_reference_rate", 0.0),
            "psi_evidence_edge_join_rate": psi_ev_info.get("edge_join_rate", 0.0),
            "psi_method_or_reference_rate": psi_ev_info.get("method_or_reference_rate", 0.0),
        },
    }
    write_json(args.report, report)

    print(
        "[OK] materialized interaction artifacts "
        f"mode={args.mode} ppi_edges={ppi_edges_info['rows']} psi_edges={psi_edges_info['rows']} "
        f"m3_rows={m3_info['rows_written']} report={args.report}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
