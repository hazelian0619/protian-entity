#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Set


THRESH_EDGE_MAPPING = 0.98
THRESH_METHOD_NON_EMPTY = 0.95
THRESH_PMID_DOI = 0.80


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _safe_rate(num: int, den: int) -> float:
    if den <= 0:
        return 1.0
    return num / den


def _load_edge_ids(path: Path) -> Set[str]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if not r.fieldnames or "edge_id" not in r.fieldnames:
            raise SystemExit(f"[ERROR] missing edge_id column: {path}")
        return {(row.get("edge_id") or "").strip() for row in r if (row.get("edge_id") or "").strip()}


@dataclass
class Gate:
    gate_id: str
    passed: bool
    detail: Dict[str, Any]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", type=Path, required=True)
    ap.add_argument("--method-context", type=Path, required=True)
    ap.add_argument("--function-context", type=Path, required=True)
    ap.add_argument("--method-validation", type=Path, required=True)
    ap.add_argument("--function-validation", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    for p in [args.edges, args.method_context, args.function_context, args.method_validation, args.function_validation]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    edge_ids = _load_edge_ids(args.edges)

    method_rows = 0
    method_edge_mapped = 0
    method_non_empty = 0
    pmid_or_doi = 0

    with args.method_context.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {"edge_id", "method", "pmid", "doi"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] method context missing columns: {sorted(missing)}")

        for row in r:
            method_rows += 1
            eid = (row.get("edge_id") or "").strip()
            if eid in edge_ids:
                method_edge_mapped += 1
            if (row.get("method") or "").strip():
                method_non_empty += 1
            if (row.get("pmid") or "").strip() or (row.get("doi") or "").strip():
                pmid_or_doi += 1

    function_rows = 0
    function_edge_mapped = 0

    with args.function_context.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {"edge_id"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] function context missing columns: {sorted(missing)}")

        for row in r:
            function_rows += 1
            eid = (row.get("edge_id") or "").strip()
            if eid in edge_ids:
                function_edge_mapped += 1

    method_map_rate = _safe_rate(method_edge_mapped, method_rows)
    function_map_rate = _safe_rate(function_edge_mapped, function_rows)
    method_non_empty_rate = _safe_rate(method_non_empty, method_rows)
    pmid_or_doi_rate = _safe_rate(pmid_or_doi, method_rows)

    method_validation = json.loads(args.method_validation.read_text(encoding="utf-8"))
    function_validation = json.loads(args.function_validation.read_text(encoding="utf-8"))
    contract_pass = bool(method_validation.get("passed")) and bool(function_validation.get("passed"))

    gates: List[Gate] = [
        Gate(
            gate_id="edge_id_mapping_coverage_method_context",
            passed=method_map_rate >= THRESH_EDGE_MAPPING,
            detail={
                "rate": method_map_rate,
                "threshold": THRESH_EDGE_MAPPING,
                "mapped": method_edge_mapped,
                "rows": method_rows,
            },
        ),
        Gate(
            gate_id="edge_id_mapping_coverage_function_context",
            passed=function_map_rate >= THRESH_EDGE_MAPPING,
            detail={
                "rate": function_map_rate,
                "threshold": THRESH_EDGE_MAPPING,
                "mapped": function_edge_mapped,
                "rows": function_rows,
            },
        ),
        Gate(
            gate_id="method_non_empty_rate",
            passed=method_non_empty_rate >= THRESH_METHOD_NON_EMPTY,
            detail={
                "rate": method_non_empty_rate,
                "threshold": THRESH_METHOD_NON_EMPTY,
                "non_empty": method_non_empty,
                "rows": method_rows,
            },
        ),
        Gate(
            gate_id="pmid_or_doi_coverage_rate",
            passed=pmid_or_doi_rate >= THRESH_PMID_DOI,
            detail={
                "rate": pmid_or_doi_rate,
                "threshold": THRESH_PMID_DOI,
                "covered": pmid_or_doi,
                "rows": method_rows,
            },
        ),
        Gate(
            gate_id="contracts_validation_pass",
            passed=contract_pass,
            detail={
                "method_contract_pass": bool(method_validation.get("passed")),
                "function_contract_pass": bool(function_validation.get("passed")),
            },
        ),
    ]

    passed = all(g.passed for g in gates)

    report: Dict[str, Any] = {
        "name": "ppi_semantic_enrichment_v2",
        "created_at": utc_now_iso(),
        "inputs": {
            "edges": str(args.edges),
            "method_context": str(args.method_context),
            "function_context": str(args.function_context),
            "method_validation": str(args.method_validation),
            "function_validation": str(args.function_validation),
        },
        "metrics": {
            "edge_count": len(edge_ids),
            "method_rows": method_rows,
            "function_rows": function_rows,
            "method_edge_mapping_rate": method_map_rate,
            "function_edge_mapping_rate": function_map_rate,
            "method_non_empty_rate": method_non_empty_rate,
            "pmid_or_doi_rate": pmid_or_doi_rate,
            "contract_pass": contract_pass,
        },
        "gates": [
            {
                "id": g.gate_id,
                "passed": g.passed,
                "detail": g.detail,
            }
            for g in gates
        ],
        "passed": passed,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if passed else "FAIL"
    print(f"[{status}] ppi_semantic_enrichment_v2 QA -> {args.out}")
    print(
        "[METRICS] "
        f"edge_map_method={method_map_rate:.4f}, "
        f"edge_map_function={function_map_rate:.4f}, "
        f"method_non_empty={method_non_empty_rate:.4f}, "
        f"pmid_or_doi={pmid_or_doi_rate:.4f}, "
        f"contracts={contract_pass}"
    )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
