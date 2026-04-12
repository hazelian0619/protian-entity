#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List


INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
NUM_RE = re.compile(r"^-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?$")
REQ_EVIDENCE_FIELDS = [
    "source_db",
    "assay_type",
    "standard_type",
    "standard_value",
    "standard_unit",
    "normalized_nM",
    "confidence",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def norm(x: object) -> str:
    return str(x or "").strip()


def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def count_baseline_chembl(db_path: Path) -> Dict[str, int]:
    with sqlite3.connect(db_path) as conn:
        q_rows = """
            SELECT COUNT(*)
            FROM psi_evidence_v1
            WHERE standard_type IN ('IC50','Ki','Kd','EC50')
              AND compound_inchikey IS NOT NULL
              AND target_uniprot_accession IS NOT NULL
              AND TRIM(compound_inchikey) != ''
              AND TRIM(target_uniprot_accession) != ''
              AND standard_value_nM IS NOT NULL
              AND CAST(standard_value_nM AS REAL) > 0
        """
        q_edges = """
            SELECT COUNT(*)
            FROM (
                SELECT compound_inchikey, target_uniprot_accession, standard_type
                FROM psi_evidence_v1
                WHERE standard_type IN ('IC50','Ki','Kd','EC50')
                  AND compound_inchikey IS NOT NULL
                  AND target_uniprot_accession IS NOT NULL
                  AND TRIM(compound_inchikey) != ''
                  AND TRIM(target_uniprot_accession) != ''
                  AND standard_value_nM IS NOT NULL
                  AND CAST(standard_value_nM AS REAL) > 0
                GROUP BY compound_inchikey, target_uniprot_accession, standard_type
            )
        """
        return {
            "chembl_baseline_rows": int(conn.execute(q_rows).fetchone()[0]),
            "chembl_baseline_edges": int(conn.execute(q_edges).fetchone()[0]),
        }


def summarize_evidence(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        fields = set(r.fieldnames)
        missing = [c for c in REQ_EVIDENCE_FIELDS if c not in fields]
        if missing:
            raise SystemExit(f"[ERROR] evidence missing required columns for QA: {missing}")

        rows = 0
        source_rows = Counter()
        required_non_empty = Counter()
        bad_inchikey = 0
        bad_nm = 0
        bad_conf = 0
        chembl_rows = 0
        for row in r:
            rows += 1
            source = norm(row.get("source_db", "")).lower()
            source_rows[source] += 1
            if source == "chembl":
                chembl_rows += 1

            ik = norm(row.get("compound_inchikey", "")).upper()
            if not INCHIKEY_RE.match(ik):
                bad_inchikey += 1

            nm = norm(row.get("normalized_nM", ""))
            if not NUM_RE.match(nm):
                bad_nm += 1
            else:
                try:
                    if float(nm) <= 0:
                        bad_nm += 1
                except Exception:
                    bad_nm += 1

            conf = norm(row.get("confidence", ""))
            if conf not in {"high", "medium", "low"}:
                bad_conf += 1

            for c in REQ_EVIDENCE_FIELDS:
                if norm(row.get(c, "")):
                    required_non_empty[c] += 1

    return {
        "rows": rows,
        "source_rows": dict(source_rows),
        "required_non_empty": {
            c: (required_non_empty[c] / rows if rows else 0.0) for c in REQ_EVIDENCE_FIELDS
        },
        "bad_inchikey": bad_inchikey,
        "bad_normalized_nM": bad_nm,
        "bad_confidence": bad_conf,
        "chembl_rows": chembl_rows,
    }


def summarize_edges(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        rows = 0
        chembl_edges = 0
        bad_best_nm = 0
        for row in r:
            rows += 1
            srcs = {x.strip().lower() for x in norm(row.get("source_dbs", "")).split(";") if x.strip()}
            if "chembl" in srcs:
                chembl_edges += 1
            nm = norm(row.get("best_normalized_nM", ""))
            if not NUM_RE.match(nm):
                bad_best_nm += 1
                continue
            try:
                if float(nm) <= 0:
                    bad_best_nm += 1
            except Exception:
                bad_best_nm += 1
    return {
        "rows": rows,
        "chembl_edges": chembl_edges,
        "bad_best_normalized_nM": bad_best_nm,
    }


def summarize_conflict(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"[ERROR] missing conflict header: {path}")
        required_cols = {
            "compound_inchikey",
            "target_uniprot_accession",
            "standard_type",
            "source_count",
            "fold_change",
            "conflict_flag",
        }
        missing = sorted(required_cols - set(r.fieldnames))
        if missing:
            raise SystemExit(f"[ERROR] conflict audit missing columns: {missing}")
        rows = 0
        flagged = 0
        for row in r:
            rows += 1
            if norm(row.get("conflict_flag")) == "1":
                flagged += 1
    return {"rows": rows, "flagged_rows": flagged}


def gate(gid: str, passed: bool, detail: Dict[str, Any]) -> Dict[str, Any]:
    return {"id": gid, "passed": bool(passed), "detail": detail}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chembl-m3-db", type=Path, required=True)
    ap.add_argument("--evidence", type=Path, required=True)
    ap.add_argument("--edges", type=Path, required=True)
    ap.add_argument("--conflict-audit-tsv", type=Path, required=True)
    ap.add_argument("--build-report", type=Path, required=True)
    ap.add_argument("--evidence-validation", type=Path, required=True)
    ap.add_argument("--edges-validation", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--require-strict", type=int, default=1)
    args = ap.parse_args()

    for p in [
        args.chembl_m3_db,
        args.evidence,
        args.edges,
        args.conflict_audit_tsv,
        args.build_report,
        args.evidence_validation,
        args.edges_validation,
    ]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    strict = bool(args.require_strict)
    baseline = count_baseline_chembl(args.chembl_m3_db)
    ev = summarize_evidence(args.evidence)
    ed = summarize_edges(args.edges)
    cf = summarize_conflict(args.conflict_audit_tsv)
    build = load_json(args.build_report)
    ev_val = load_json(args.evidence_validation)
    ed_val = load_json(args.edges_validation)

    source_rows = ev["source_rows"]
    source_count = sum(1 for _, v in source_rows.items() if v > 0)
    bindingdb_used = source_rows.get("bindingdb", 0) > 0
    pdbbind_used = source_rows.get("pdbbind", 0) > 0

    required_field_ok = all(rate >= 0.99 for rate in ev["required_non_empty"].values())
    chembl_rows_ok = (ev["chembl_rows"] >= baseline["chembl_baseline_rows"]) if strict else True
    chembl_edges_ok = (ed["chembl_edges"] >= baseline["chembl_baseline_edges"]) if strict else True
    contract_ok = bool(ev_val.get("passed")) and bool(ed_val.get("passed"))

    gates: List[Dict[str, Any]] = [
        gate(
            "contracts_pass",
            contract_ok,
            {
                "evidence_contract_passed": bool(ev_val.get("passed")),
                "edges_contract_passed": bool(ed_val.get("passed")),
            },
        ),
        gate(
            "required_fields_non_empty_rate",
            required_field_ok,
            {
                "required_field_rates": ev["required_non_empty"],
                "threshold": 0.99,
            },
        ),
        gate(
            "normalized_nM_and_confidence_valid",
            ev["bad_normalized_nM"] == 0 and ev["bad_confidence"] == 0 and ev["bad_inchikey"] == 0 and ed["bad_best_normalized_nM"] == 0,
            {
                "evidence_bad_normalized_nM": ev["bad_normalized_nM"],
                "evidence_bad_confidence": ev["bad_confidence"],
                "evidence_bad_inchikey": ev["bad_inchikey"],
                "edges_bad_best_normalized_nM": ed["bad_best_normalized_nM"],
            },
        ),
        gate(
            "multi_source_present",
            source_count >= 2 if strict else source_count >= 1,
            {
                "source_rows": source_rows,
                "source_count_nonzero": source_count,
                "bindingdb_used": bindingdb_used,
                "pdbbind_used": pdbbind_used,
                "strict_mode": strict,
            },
        ),
        gate(
            "conflict_audit_exists_and_parseable",
            cf["rows"] >= 0,
            {
                "audit_rows": cf["rows"],
                "flagged_rows": cf["flagged_rows"],
            },
        ),
        gate(
            "chembl_hard_gate_not_degraded",
            chembl_rows_ok and chembl_edges_ok,
            {
                "strict_mode": strict,
                "chembl_rows_v2": ev["chembl_rows"],
                "chembl_rows_baseline": baseline["chembl_baseline_rows"],
                "chembl_edges_v2": ed["chembl_edges"],
                "chembl_edges_baseline": baseline["chembl_baseline_edges"],
                "build_metrics": build.get("metrics", {}),
            },
        ),
    ]

    passed = all(g["passed"] for g in gates)
    out = {
        "name": "molecule_activity_fusion_v2.qa",
        "created_at": utc_now(),
        "strict_mode": strict,
        "inputs": {
            "chembl_m3_db": str(args.chembl_m3_db),
            "evidence": str(args.evidence),
            "edges": str(args.edges),
            "conflict_audit_tsv": str(args.conflict_audit_tsv),
            "build_report": str(args.build_report),
            "evidence_validation": str(args.evidence_validation),
            "edges_validation": str(args.edges_validation),
        },
        "metrics": {
            "baseline": baseline,
            "evidence": ev,
            "edges": ed,
            "conflict_audit": cf,
        },
        "gates": gates,
        "passed": passed,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"[{'PASS' if passed else 'FAIL'}] molecule_activity_fusion_v2 QA")
    for g in gates:
        print(f"  - {'PASS' if g['passed'] else 'FAIL'} {g['id']}")
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
