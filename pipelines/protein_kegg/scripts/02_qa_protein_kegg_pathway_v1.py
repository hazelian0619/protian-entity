#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Set


KEGG_PATHWAY_RE = re.compile(r"^hsa\d{5}$")
KEGG_GENE_RE = re.compile(r"^hsa:\d+$")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _load_master_uniprots(path: Path) -> Set[str]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if not r.fieldnames or "uniprot_id" not in r.fieldnames:
            raise SystemExit(f"[ERROR] master table missing uniprot_id: {path}")
        return {
            (row.get("uniprot_id") or "").strip()
            for row in r
            if (row.get("uniprot_id") or "").strip()
        }


def _load_reactome_uniprots(path: Path) -> Set[str]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if not r.fieldnames or "uniprot_id" not in r.fieldnames:
            raise SystemExit(f"[ERROR] Reactome table missing uniprot_id: {path}")
        return {
            (row.get("uniprot_id") or "").strip()
            for row in r
            if (row.get("uniprot_id") or "").strip()
        }


@dataclass
class Gate:
    gate_id: str
    passed: bool
    detail: Dict[str, Any]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kegg", type=Path, required=True)
    ap.add_argument("--master", type=Path, required=True)
    ap.add_argument("--reactome", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    if not args.kegg.exists():
        raise SystemExit(f"[ERROR] missing KEGG output: {args.kegg}")
    if not args.master.exists():
        raise SystemExit(f"[ERROR] missing master table: {args.master}")
    if not args.reactome.exists():
        raise SystemExit(f"[ERROR] missing Reactome table: {args.reactome}")

    master_uniprots = _load_master_uniprots(args.master)
    reactome_uniprots = _load_reactome_uniprots(args.reactome)

    row_count = 0
    missing_fk = 0
    bad_pathway_format = 0
    bad_kegg_gene_format = 0
    empty_pathway_name = 0

    kegg_uniprots: Set[str] = set()
    kegg_pathways: Set[str] = set()
    kegg_genes: Set[str] = set()

    with args.kegg.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {
            "uniprot_id",
            "kegg_gene_id",
            "kegg_pathway_id",
            "pathway_name",
            "source_version",
            "fetch_date",
        }
        missing = required - set(r.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] KEGG output schema missing columns: {sorted(missing)}")

        for row in r:
            row_count += 1
            uniprot_id = (row.get("uniprot_id") or "").strip()
            kegg_gene_id = (row.get("kegg_gene_id") or "").strip()
            kegg_pathway_id = (row.get("kegg_pathway_id") or "").strip()
            pathway_name = (row.get("pathway_name") or "").strip()

            if uniprot_id not in master_uniprots:
                missing_fk += 1

            if not KEGG_GENE_RE.match(kegg_gene_id):
                bad_kegg_gene_format += 1
            if not KEGG_PATHWAY_RE.match(kegg_pathway_id):
                bad_pathway_format += 1
            if not pathway_name:
                empty_pathway_name += 1

            if uniprot_id:
                kegg_uniprots.add(uniprot_id)
            if kegg_pathway_id:
                kegg_pathways.add(kegg_pathway_id)
            if kegg_gene_id:
                kegg_genes.add(kegg_gene_id)

    join_rate = 1.0 if row_count == 0 else (row_count - missing_fk) / row_count
    mapped_uniprot_rate = (len(kegg_uniprots) / len(master_uniprots)) if master_uniprots else 0.0

    overlap = kegg_uniprots & reactome_uniprots
    kegg_only = kegg_uniprots - reactome_uniprots
    reactome_only = reactome_uniprots - kegg_uniprots
    union = kegg_uniprots | reactome_uniprots
    jaccard = (len(overlap) / len(union)) if union else 1.0
    neither = master_uniprots - union

    gates: List[Gate] = [
        Gate(
            gate_id="join_to_master_fk",
            passed=(missing_fk == 0),
            detail={
                "missing_rows": missing_fk,
                "join_rate": join_rate,
            },
        ),
        Gate(
            gate_id="kegg_pathway_id_format",
            passed=(bad_pathway_format == 0),
            detail={"bad_rows": bad_pathway_format, "pattern": KEGG_PATHWAY_RE.pattern},
        ),
        Gate(
            gate_id="kegg_gene_id_format",
            passed=(bad_kegg_gene_format == 0),
            detail={"bad_rows": bad_kegg_gene_format, "pattern": KEGG_GENE_RE.pattern},
        ),
        Gate(
            gate_id="pathway_name_non_empty",
            passed=(empty_pathway_name == 0),
            detail={"bad_rows": empty_pathway_name},
        ),
    ]

    report: Dict[str, Any] = {
        "name": "protein_kegg_pathway_v1.qa",
        "created_at": utc_now(),
        "inputs": {
            "kegg": str(args.kegg),
            "master": str(args.master),
            "reactome": str(args.reactome),
        },
        "metrics": {
            "kegg_rows": row_count,
            "kegg_unique_uniprots": len(kegg_uniprots),
            "kegg_unique_kegg_genes": len(kegg_genes),
            "kegg_unique_pathways": len(kegg_pathways),
            "master_unique_uniprots": len(master_uniprots),
            "reactome_unique_uniprots": len(reactome_uniprots),
            "join_rate_to_master": join_rate,
            "mapped_uniprot_rate_vs_master": mapped_uniprot_rate,
        },
        "complementarity_vs_reactome": {
            "overlap_uniprots": len(overlap),
            "kegg_only_uniprots": len(kegg_only),
            "reactome_only_uniprots": len(reactome_only),
            "neither_uniprots": len(neither),
            "union_uniprots": len(union),
            "jaccard_uniprots": jaccard,
        },
        "gates": [
            {
                "gate_id": g.gate_id,
                "passed": g.passed,
                "detail": g.detail,
            }
            for g in gates
        ],
        "passed": all(g.passed for g in gates),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"[OK] QA report -> {args.out}")
    print(f"[OK] QA passed={report['passed']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
