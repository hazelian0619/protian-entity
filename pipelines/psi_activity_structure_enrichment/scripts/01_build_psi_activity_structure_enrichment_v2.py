#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import re
import sqlite3
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


ASSAY_TYPE_DESC = {
    "B": "Binding assay",
    "F": "Functional assay",
    "A": "ADME assay",
    "T": "Toxicity assay",
    "P": "Physicochemical assay",
}

RELATION_ALLOWED = {"=", "<", ">", "<=", ">=", "~"}

SYSTEM_KEYWORDS = [
    "cell-based",
    "cell based",
    "cell-free",
    "cell free",
    "microsome",
    "plasma",
    "serum",
    "lysate",
    "buffer",
    "recombinant",
    "membrane",
    "whole blood",
    "cytosol",
    "homogenate",
]

PH_RE = re.compile(r"\bp\s*h\s*[:=]?\s*([0-9]+(?:\.[0-9]+)?)\b", re.IGNORECASE)
TEMP_RE = re.compile(r"\b(?:temp(?:erature)?\s*[:=]?\s*|at\s*)([0-9]+(?:\.[0-9]+)?)\s*(?:°?\s*[CF])\b", re.IGNORECASE)
TEMP_C_RE = re.compile(r"\b([0-9]+(?:\.[0-9]+)?)\s*°\s*C\b", re.IGNORECASE)


ACTIVITY_HEADER = [
    "edge_id",
    "activity_id",
    "assay_id",
    "doc_id",
    "compound_inchikey",
    "target_uniprot_accession",
    "target_chembl_id",
    "standard_type",
    "standard_relation",
    "standard_relation_raw",
    "standard_value",
    "standard_units",
    "standard_value_nM",
    "pchembl_value",
    "pchembl_value_eff",
    "assay_type",
    "assay_type_desc",
    "assay_description",
    "assay_context",
    "assay_confidence_score",
    "bao_format",
    "condition_pH",
    "condition_temperature_c",
    "condition_system",
    "condition_context",
    "data_validity_comment",
    "activity_comment",
    "doi",
    "pubmed_id",
    "n_evidence",
    "n_docs",
    "evidence_score_max",
    "source",
    "source_version",
    "fetch_date",
]


STRUCTURE_HEADER = [
    "edge_id",
    "activity_id",
    "compound_inchikey",
    "target_uniprot_accession",
    "target_chembl_id",
    "pdb_id",
    "pdb_experimental_method",
    "pdb_resolution",
    "pdb_release_date",
    "pdb_ligand_count",
    "structure_evidence_type",
    "complex_match_level",
    "standard_type",
    "standard_relation",
    "standard_value_nM",
    "pchembl_value_eff",
    "structure_affinity_score",
    "source",
    "source_version",
    "fetch_date",
]


def utc_now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def utc_today() -> str:
    return dt.datetime.now(dt.timezone.utc).date().isoformat()


def non_empty(v: object) -> bool:
    return str(v or "").strip() not in {"", "NA", "N/A", "None", "null"}


def parse_float(v: object) -> Optional[float]:
    if not non_empty(v):
        return None
    try:
        return float(str(v).strip())
    except Exception:
        return None


def normalize_relation(raw: object) -> str:
    s = str(raw or "").strip()
    if s in RELATION_ALLOWED:
        return s

    mapping = {
        "==": "=",
        "~=": "~",
        "≈": "~",
        "<=": "<=",
        "=<": "<=",
        ">=": ">=",
        "=>": ">=",
    }
    if s in mapping:
        return mapping[s]

    if not s:
        return "~"

    if "<" in s and "=" in s:
        return "<="
    if ">" in s and "=" in s:
        return ">="
    if "<" in s:
        return "<"
    if ">" in s:
        return ">"
    return "~"


def normalize_assay_type(raw: object) -> str:
    s = str(raw or "").strip().upper()
    if s in ASSAY_TYPE_DESC:
        return s
    return "UNK"


def first_sentence(text: str, max_len: int = 180) -> str:
    t = (text or "").strip()
    if not t:
        return ""
    m = re.split(r"[\.;\n]", t, maxsplit=1)
    x = m[0].strip() if m else t
    return x[:max_len].strip()


def extract_conditions(assay_desc: str) -> Tuple[str, str, str, str]:
    txt = (assay_desc or "").strip()
    if not txt:
        return "", "", "", ""

    ph = ""
    m_ph = PH_RE.search(txt)
    if m_ph:
        ph = m_ph.group(1)

    temp = ""
    m_temp = TEMP_RE.search(txt)
    if not m_temp:
        m_temp = TEMP_C_RE.search(txt)
    if m_temp:
        temp = m_temp.group(1)

    lower = txt.lower()
    systems: List[str] = []
    for kw in SYSTEM_KEYWORDS:
        if kw in lower:
            systems.append(kw)
    system = ";".join(sorted(set(systems)))

    parts = []
    if ph:
        parts.append(f"pH={ph}")
    if temp:
        parts.append(f"temperature_c={temp}")
    if system:
        parts.append(f"system={system}")
    ctx = "; ".join(parts)

    return ph, temp, system, ctx


@dataclass
class PDBMeta:
    pdb_id: str
    experimental_method: str
    resolution: str
    release_date: str
    ligand_count: int


def pdb_rank(meta: PDBMeta) -> Tuple[int, float, int, str]:
    # rank: high ligand_count, low resolution, newer date
    res = parse_float(meta.resolution)
    res_rank = res if res is not None else 99.0
    date_rank = int(meta.release_date.replace("-", "")) if re.match(r"^\d{4}-\d{2}-\d{2}$", meta.release_date) else 0
    return (-meta.ligand_count, res_rank, -date_rank, meta.pdb_id)


def load_best_pdb_map(path: Path, min_ligand_count: int) -> Tuple[Dict[str, PDBMeta], Dict[str, int]]:
    best: Dict[str, PDBMeta] = {}
    stats = {
        "rows_total": 0,
        "rows_with_uniprot": 0,
        "rows_ligand_filter_kept": 0,
        "unique_uniprot_with_structure": 0,
    }

    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"pdb_id", "uniprot_id", "experimental_method", "resolution", "release_date", "ligand_count"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] protein pdb table missing columns: {sorted(missing)}")

        for row in reader:
            stats["rows_total"] += 1
            uniprot = str(row.get("uniprot_id") or "").strip()
            if not uniprot:
                continue
            stats["rows_with_uniprot"] += 1

            ligand_count = int(parse_float(row.get("ligand_count")) or 0)
            if ligand_count < min_ligand_count:
                continue

            stats["rows_ligand_filter_kept"] += 1
            meta = PDBMeta(
                pdb_id=str(row.get("pdb_id") or "").strip().upper(),
                experimental_method=str(row.get("experimental_method") or "").strip(),
                resolution=str(row.get("resolution") or "").strip(),
                release_date=str(row.get("release_date") or "").strip(),
                ligand_count=ligand_count,
            )
            if not meta.pdb_id:
                continue

            old = best.get(uniprot)
            if old is None or pdb_rank(meta) < pdb_rank(old):
                best[uniprot] = meta

    stats["unique_uniprot_with_structure"] = len(best)
    return best, stats


def get_source_version(conn: sqlite3.Connection) -> str:
    try:
        row = conn.execute(
            """
            SELECT script_version, rules_version
            FROM meta_run
            ORDER BY started_at_utc DESC
            LIMIT 1
            """
        ).fetchone()
        if row:
            sv = str(row[0] or "").strip()
            rv = str(row[1] or "").strip()
            if sv or rv:
                return "|".join([x for x in [sv, rv] if x])
    except Exception:
        pass
    return "chembl_m3_v1"


def as_str(v: object) -> str:
    if v is None:
        return ""
    return str(v)


def format_float(v: Optional[float], digits: int = 6) -> str:
    if v is None:
        return ""
    x = f"{v:.{digits}f}".rstrip("0").rstrip(".")
    return x if x != "-0" else "0"


def structure_affinity_score(pchembl_eff: Optional[float], pdb_resolution: Optional[float], ligand_count: int) -> float:
    affinity_term = 0.0
    if pchembl_eff is not None:
        affinity_term = max(0.0, min(1.0, (pchembl_eff - 5.0) / 4.0))

    if pdb_resolution is None:
        res_term = 0.3
    else:
        res_term = max(0.0, min(1.0, (4.0 - pdb_resolution) / 3.5))

    ligand_term = max(0.0, min(1.0, ligand_count / 10.0))
    score = 0.5 * affinity_term + 0.35 * res_term + 0.15 * ligand_term
    return round(max(0.0, min(1.0, score)), 6)


def iter_join_rows(conn: sqlite3.Connection, max_rows: int) -> Iterable[sqlite3.Row]:
    sql = """
    SELECT
        e.edge_id,
        e.compound_inchikey,
        e.target_uniprot_accession,
        e.target_chembl_id,
        e.standard_type AS edge_standard_type,
        e.n_evidence,
        e.n_docs,
        e.evidence_score_max,
        e.best_activity_id,

        ev.activity_id,
        ev.assay_id,
        ev.doc_id,
        ev.standard_type,
        ev.standard_relation,
        ev.standard_value,
        ev.standard_units,
        ev.standard_value_nM,
        ev.pchembl_value,
        ev.pchembl_value_eff,
        ev.assay_type,
        ev.assay_confidence_score,
        ev.assay_description,
        ev.bao_format,
        ev.data_validity_comment,
        ev.activity_comment,
        ev.doi,
        ev.pubmed_id,
        ev.source
    FROM psi_edges_v1 e
    JOIN psi_evidence_v1 ev
      ON ev.activity_id = e.best_activity_id
    ORDER BY e.edge_id
    """
    params: Sequence[object] = ()
    if max_rows > 0:
        sql += " LIMIT ?"
        params = (max_rows,)

    cur = conn.execute(sql, params)
    while True:
        rows = cur.fetchmany(10000)
        if not rows:
            break
        for row in rows:
            yield row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chembl-m3-db", type=Path, required=True)
    ap.add_argument("--protein-pdb", type=Path, required=True)
    ap.add_argument("--out-activity", type=Path, required=True)
    ap.add_argument("--out-structure", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--fetch-date", default=utc_today())
    ap.add_argument("--source-version", default="")
    ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("--min-ligand-count", type=int, default=1)
    args = ap.parse_args()

    if not args.chembl_m3_db.exists():
        raise SystemExit(f"[ERROR] missing input: {args.chembl_m3_db}")
    if not args.protein_pdb.exists():
        raise SystemExit(f"[ERROR] missing input: {args.protein_pdb}")

    best_pdb_map, pdb_stats = load_best_pdb_map(args.protein_pdb, min_ligand_count=args.min_ligand_count)

    conn = sqlite3.connect(f"file:{args.chembl_m3_db.resolve()}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row

    total_psi_edges = int(conn.execute("SELECT COUNT(*) FROM psi_edges_v1").fetchone()[0])
    source_version = args.source_version.strip() or get_source_version(conn)

    args.out_activity.parent.mkdir(parents=True, exist_ok=True)
    args.out_structure.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)

    tmp_activity = args.out_activity.with_suffix(args.out_activity.suffix + ".tmp")
    tmp_structure = args.out_structure.with_suffix(args.out_structure.suffix + ".tmp")

    counts = {
        "activity_rows": 0,
        "activity_standard_type_non_empty": 0,
        "activity_relation_allowed": 0,
        "activity_assay_fields_non_empty": 0,
        "activity_condition_any_non_empty": 0,
        "structure_rows": 0,
    }

    with (
        tmp_activity.open("w", encoding="utf-8", newline="") as fa,
        tmp_structure.open("w", encoding="utf-8", newline="") as fs,
    ):
        wa = csv.writer(fa, delimiter="\t", lineterminator="\n")
        ws = csv.writer(fs, delimiter="\t", lineterminator="\n")
        wa.writerow(ACTIVITY_HEADER)
        ws.writerow(STRUCTURE_HEADER)

        for row in iter_join_rows(conn, max_rows=args.max_rows):
            edge_id = as_str(row["edge_id"]).strip()
            standard_type = as_str(row["standard_type"] or row["edge_standard_type"]).strip()
            relation_raw = as_str(row["standard_relation"]).strip()
            relation = normalize_relation(relation_raw)

            assay_type = normalize_assay_type(row["assay_type"])
            assay_type_desc = ASSAY_TYPE_DESC.get(assay_type, "Unknown assay type")
            assay_desc = as_str(row["assay_description"]).strip()
            assay_context = f"{assay_type_desc}; {first_sentence(assay_desc)}".strip("; ")
            pH, temp_c, system, cond_ctx = extract_conditions(assay_desc)

            assay_any = any(
                [
                    assay_type != "UNK",
                    non_empty(assay_desc),
                    non_empty(row["bao_format"]),
                ]
            )

            wa.writerow(
                [
                    edge_id,
                    as_str(row["activity_id"]),
                    as_str(row["assay_id"]),
                    as_str(row["doc_id"]),
                    as_str(row["compound_inchikey"]),
                    as_str(row["target_uniprot_accession"]),
                    as_str(row["target_chembl_id"]),
                    standard_type,
                    relation,
                    relation_raw,
                    as_str(row["standard_value"]),
                    as_str(row["standard_units"]),
                    as_str(row["standard_value_nM"]),
                    as_str(row["pchembl_value"]),
                    as_str(row["pchembl_value_eff"]),
                    assay_type,
                    assay_type_desc,
                    assay_desc,
                    assay_context,
                    as_str(row["assay_confidence_score"]),
                    as_str(row["bao_format"]),
                    pH,
                    temp_c,
                    system,
                    cond_ctx,
                    as_str(row["data_validity_comment"]),
                    as_str(row["activity_comment"]),
                    as_str(row["doi"]),
                    as_str(row["pubmed_id"]),
                    as_str(row["n_evidence"]),
                    as_str(row["n_docs"]),
                    as_str(row["evidence_score_max"]),
                    as_str(row["source"]),
                    source_version,
                    args.fetch_date,
                ]
            )

            counts["activity_rows"] += 1
            if non_empty(standard_type):
                counts["activity_standard_type_non_empty"] += 1
            if relation in RELATION_ALLOWED:
                counts["activity_relation_allowed"] += 1
            if assay_any:
                counts["activity_assay_fields_non_empty"] += 1
            if pH or temp_c or system:
                counts["activity_condition_any_non_empty"] += 1

            uniprot = as_str(row["target_uniprot_accession"]).strip()
            pdb = best_pdb_map.get(uniprot)
            if pdb is not None:
                pchembl_eff = parse_float(row["pchembl_value_eff"])
                std_val_nm = parse_float(row["standard_value_nM"])
                score = structure_affinity_score(
                    pchembl_eff=pchembl_eff,
                    pdb_resolution=parse_float(pdb.resolution),
                    ligand_count=pdb.ligand_count,
                )
                ws.writerow(
                    [
                        edge_id,
                        as_str(row["activity_id"]),
                        as_str(row["compound_inchikey"]),
                        uniprot,
                        as_str(row["target_chembl_id"]),
                        pdb.pdb_id,
                        pdb.experimental_method,
                        pdb.resolution,
                        pdb.release_date,
                        str(pdb.ligand_count),
                        "pdb_target_complex_context",
                        "target_level",
                        standard_type,
                        relation,
                        format_float(std_val_nm),
                        format_float(pchembl_eff),
                        format_float(score),
                        "RCSB_PDB+ChEMBL36",
                        source_version,
                        args.fetch_date,
                    ]
                )
                counts["structure_rows"] += 1

    tmp_activity.replace(args.out_activity)
    tmp_structure.replace(args.out_structure)

    conn.close()

    processed_edges = counts["activity_rows"]
    activity_join_rate = 1.0 if processed_edges else 0.0
    activity_type_cov = (counts["activity_standard_type_non_empty"] / processed_edges) if processed_edges else 0.0
    assay_cov = (counts["activity_assay_fields_non_empty"] / processed_edges) if processed_edges else 0.0
    relation_cov = (counts["activity_relation_allowed"] / processed_edges) if processed_edges else 0.0

    structure_cov_processed = (counts["structure_rows"] / processed_edges) if processed_edges else 0.0
    structure_cov_total = (counts["structure_rows"] / total_psi_edges) if total_psi_edges else 0.0

    report = {
        "name": "psi_activity_structure_enrichment_v2.build",
        "created_at": utc_now_iso(),
        "inputs": {
            "chembl_m3_db": str(args.chembl_m3_db),
            "protein_pdb": str(args.protein_pdb),
            "max_rows": args.max_rows,
            "min_ligand_count": args.min_ligand_count,
        },
        "outputs": {
            "activity_context": str(args.out_activity),
            "structure_evidence": str(args.out_structure),
        },
        "counts": {
            "total_psi_edges_in_db": total_psi_edges,
            **counts,
        },
        "coverage": {
            "edge_id_join_rate_activity": activity_join_rate,
            "activity_type_coverage": activity_type_cov,
            "relation_symbol_coverage": relation_cov,
            "assay_field_coverage": assay_cov,
            "structure_edge_coverage_on_processed": structure_cov_processed,
            "structure_edge_coverage_on_total_psi": structure_cov_total,
        },
        "pdb_selection_stats": pdb_stats,
        "source": {
            "source_version": source_version,
            "fetch_date": args.fetch_date,
        },
    }

    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(
        "[OK] psi_activity_structure_enrichment_v2 build "
        f"activity_rows={processed_edges:,} structure_rows={counts['structure_rows']:,} "
        f"assay_cov={assay_cov:.4f} structure_cov_total={structure_cov_total:.4f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
