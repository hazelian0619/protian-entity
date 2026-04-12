#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def resolve_reference_path(raw: str) -> Path:
    ref = (raw or "").strip()
    if "#" in ref:
        ref = ref.split("#", 1)[0]
    p = Path(ref)
    if not p.is_absolute():
        p = (Path.cwd() / p).resolve()
    return p


def summarize_existence(rows: Iterable[Dict[str, str]], column: str) -> Dict[str, object]:
    refs = [(r.get(column, "") or "").strip() for r in rows]
    refs = [r for r in refs if r]
    unique_refs = sorted(set(refs))
    missing: List[str] = []
    for ref in unique_refs:
        if not resolve_reference_path(ref).exists():
            missing.append(ref)
    existing = len(unique_refs) - len(missing)
    rate = (existing / len(unique_refs)) if unique_refs else 0.0
    return {
        "checked_column": column,
        "total_rows": len(refs),
        "unique_references": len(unique_refs),
        "existing_references": existing,
        "missing_references": len(missing),
        "existing_rate": rate,
        "missing_examples": missing[:20],
    }


def write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    p = argparse.ArgumentParser(description="QA for structure model file references")
    p.add_argument("--cov-table", type=Path, required=True)
    p.add_argument("--pred-table", type=Path, required=True)
    p.add_argument(
        "--out",
        type=Path,
        default=Path("pipelines/rna_structure_models/reports/rna_structure_models_v1.file_existence.json"),
    )
    args = p.parse_args()

    cov_rows = read_tsv(args.cov_table)
    pred_rows = read_tsv(args.pred_table)

    cov_check = summarize_existence(cov_rows, "cm_file")
    pred_check = summarize_existence(pred_rows, "structure_file")

    report = {
        "pipeline": "rna_structure_models_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": {
            "cov_table": str(args.cov_table),
            "pred_table": str(args.pred_table),
        },
        "covariance": cov_check,
        "predicted": pred_check,
        "acceptance": {
            "covariance_file_refs_resolved": cov_check["existing_rate"] == 1.0 and cov_check["unique_references"] > 0,
            "predicted_file_refs_resolved": pred_check["existing_rate"] == 1.0 and pred_check["unique_references"] > 0,
        },
    }

    write_json(args.out, report)

    passed = report["acceptance"]["covariance_file_refs_resolved"] and report["acceptance"]["predicted_file_refs_resolved"]
    print(
        f"[{'PASS' if passed else 'FAIL'}] cov_ref_rate={cov_check['existing_rate']:.4f} "
        f"pred_ref_rate={pred_check['existing_rate']:.4f}"
    )
    print(f"[OK] report -> {args.out}")
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
