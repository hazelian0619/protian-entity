#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

LEVELS = {"high", "medium", "low"}
FETCH_DATE = date.today().isoformat()

EXPR_OUTPUT_COLUMNS = [
    "rna_id",
    "biosample_tissue",
    "biosample_cell_type",
    "expression_value",
    "expression_unit",
    "assay",
    "sample_id",
    "experiment_id",
    "source_dataset",
    "source_file",
    "evidence_level",
    "source",
    "source_version",
    "fetch_date",
]

RBP_OUTPUT_COLUMNS = [
    "rna_id",
    "rbp_symbol",
    "biosample_cell_type",
    "assay",
    "peak_count",
    "score",
    "evalue",
    "sample_id",
    "experiment_id",
    "source_dataset",
    "source_file",
    "evidence_level",
    "source",
    "source_version",
    "fetch_date",
]


def open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def file_exists_with_alt(path: Path) -> Optional[Path]:
    if path.exists():
        return path
    if path.suffix == ".gz":
        alt = Path(str(path)[:-3])
        if alt.exists():
            return alt
    else:
        alt = Path(str(path) + ".gz")
        if alt.exists():
            return alt
    return None


def detect_delimiter(path: Path) -> str:
    with open_maybe_gzip(path) as f:
        head = f.read(4096)
    tab = head.count("\t")
    comma = head.count(",")
    return "\t" if tab >= comma else ","


def normalize(s: str) -> str:
    return (s or "").strip()


def lower(s: str) -> str:
    return normalize(s).lower()


def to_level(s: str, default: str = "low") -> str:
    x = lower(s)
    return x if x in LEVELS else default


def numeric_or_empty(v: str) -> str:
    t = normalize(v)
    if t == "":
        return ""
    try:
        float(t)
    except Exception:
        return ""
    return t


def first_existing(row: Dict[str, str], aliases: List[str], default: str = "") -> str:
    for a in aliases:
        if a in row and normalize(row[a]) != "":
            return normalize(row[a])
    return default


def read_rna_master_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None or "rna_id" not in set(reader.fieldnames):
            raise SystemExit("[ERROR] rna_master 缺少 rna_id 列")
        for row in reader:
            rid = normalize(row.get("rna_id", ""))
            if rid:
                ids.add(rid)
    return ids


def parse_expression_rows(path: Path, limit: Optional[int]) -> List[Dict[str, str]]:
    delim = detect_delimiter(path)
    out: List[Dict[str, str]] = []
    with open_maybe_gzip(path) as f:
        reader = csv.DictReader(f, delimiter=delim)
        if reader.fieldnames is None:
            raise SystemExit("[ERROR] expression 输入缺少表头")
        for row in reader:
            rid = first_existing(
                row,
                ["rna_id", "RNA_ID", "target_rna_id", "target_id", "transcript_id", "urs_id"],
            )
            if rid == "":
                continue
            expr_value = numeric_or_empty(
                first_existing(
                    row,
                    ["expression_value", "value", "tpm", "fpkm", "rpm", "signal", "quant_value"],
                )
            )
            out.append(
                {
                    "rna_id": rid,
                    "biosample_tissue": first_existing(row, ["biosample_tissue", "tissue", "organ"], "unknown"),
                    "biosample_cell_type": first_existing(row, ["biosample_cell_type", "cell_type", "biosample"], "unknown"),
                    "expression_value": expr_value,
                    "expression_unit": first_existing(row, ["expression_unit", "unit"], "TPM"),
                    "assay": first_existing(row, ["assay", "assay_type"], "RNA-seq"),
                    "sample_id": first_existing(row, ["sample_id", "biosample_id", "sample"], ""),
                    "experiment_id": first_existing(row, ["experiment_id", "accession", "experiment"], ""),
                    "source_dataset": first_existing(row, ["source_dataset", "dataset", "project"], "ENCODE"),
                    "source_file": first_existing(row, ["source_file", "file_accession", "file"], path.name),
                    "evidence_level": to_level(first_existing(row, ["evidence_level", "confidence", "level"], "medium"), "medium"),
                    "source": "ENCODE",
                }
            )
            if limit is not None and len(out) >= limit:
                break
    return out


def parse_rbp_rows(path: Path, limit: Optional[int]) -> List[Dict[str, str]]:
    delim = detect_delimiter(path)
    out: List[Dict[str, str]] = []
    with open_maybe_gzip(path) as f:
        reader = csv.DictReader(f, delimiter=delim)
        if reader.fieldnames is None:
            raise SystemExit("[ERROR] rbp 输入缺少表头")
        for row in reader:
            rid = first_existing(
                row,
                ["rna_id", "RNA_ID", "target_rna_id", "target_id", "transcript_id", "urs_id"],
            )
            if rid == "":
                continue
            out.append(
                {
                    "rna_id": rid,
                    "rbp_symbol": first_existing(row, ["rbp_symbol", "rbp", "protein", "gene_symbol"], "unknown"),
                    "biosample_cell_type": first_existing(row, ["biosample_cell_type", "cell_type", "biosample"], "unknown"),
                    "assay": first_existing(row, ["assay", "assay_type"], "eCLIP"),
                    "peak_count": first_existing(row, ["peak_count", "peaks", "site_count"], ""),
                    "score": numeric_or_empty(first_existing(row, ["score", "signal", "confidence_score"], "")),
                    "evalue": numeric_or_empty(first_existing(row, ["evalue", "qvalue", "pvalue"], "")),
                    "sample_id": first_existing(row, ["sample_id", "biosample_id", "sample"], ""),
                    "experiment_id": first_existing(row, ["experiment_id", "accession", "experiment"], ""),
                    "source_dataset": first_existing(row, ["source_dataset", "dataset", "project"], "ENCODE"),
                    "source_file": first_existing(row, ["source_file", "file_accession", "file"], path.name),
                    "evidence_level": to_level(first_existing(row, ["evidence_level", "confidence", "level"], "medium"), "medium"),
                    "source": "ENCODE",
                }
            )
            if limit is not None and len(out) >= limit:
                break
    return out


def write_tsv(path: Path, rows: List[Dict[str, str]], columns: List[str], source_version: str, fetch_date: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            row = dict(row)
            row["source_version"] = source_version
            row["fetch_date"] = fetch_date
            writer.writerow({k: normalize(str(row.get(k, ""))) for k in columns})


def completeness_rate(rows: List[Dict[str, str]], required_cols: List[str]) -> float:
    if not rows:
        return 0.0
    total_cells = len(rows) * len(required_cols)
    ok = 0
    for row in rows:
        for col in required_cols:
            if normalize(row.get(col, "")) != "":
                ok += 1
    return ok / total_cells if total_cells else 0.0


def level_non_empty_rate(rows: List[Dict[str, str]]) -> float:
    if not rows:
        return 0.0
    non_empty = sum(1 for r in rows if lower(r.get("evidence_level", "")) in LEVELS)
    return non_empty / len(rows)


def join_rate(rows: List[Dict[str, str]], master_ids: set[str]) -> Tuple[float, int, int]:
    if not rows:
        return 0.0, 0, 0
    matched = sum(1 for r in rows if normalize(r.get("rna_id", "")) in master_ids)
    total = len(rows)
    return (matched / total if total else 0.0), matched, total


def write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def make_download_checklist(rna_master: Path, expr_path: Path, rbp_path: Path) -> List[Dict[str, str]]:
    return [
        {
            "file": "RNA 主表",
            "required_path": str(rna_master),
            "suggested_source": "由 pipelines/rna/run.sh 产出",
            "sha256_command": f"sha256sum {rna_master}",
        },
        {
            "file": "ENCODE 表达汇总文件（TSV/CSV，支持 .gz）",
            "required_path": str(expr_path),
            "suggested_source": "ENCODE Portal（RNA-seq quantification）",
            "sha256_command": f"sha256sum {expr_path}",
        },
        {
            "file": "ENCODE RBP 位点文件（TSV/CSV，支持 .gz）",
            "required_path": str(rbp_path),
            "suggested_source": "ENCODE Portal（eCLIP peaks）",
            "sha256_command": f"sha256sum {rbp_path}",
        },
        {
            "file": "ENCODE 元数据清单",
            "required_path": "data/raw/rna/encode/metadata.tsv",
            "suggested_source": "ENCODE search export",
            "sha256_command": "sha256sum data/raw/rna/encode/metadata.tsv",
        },
    ]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build RNA expression + RBP evidence tables")
    p.add_argument("--rna-master", type=Path, default=Path("data/output/rna_master_v1.tsv"))
    p.add_argument("--expression-input", type=Path, default=Path("data/raw/rna/encode/expression_evidence.tsv.gz"))
    p.add_argument("--rbp-input", type=Path, default=Path("data/raw/rna/encode/rbp_sites.tsv.gz"))
    p.add_argument("--expression-output", type=Path, default=Path("data/output/rna_expression_evidence_v1.tsv"))
    p.add_argument("--rbp-output", type=Path, default=Path("data/output/rna_rbp_sites_v1.tsv"))
    p.add_argument("--report", type=Path, default=Path("pipelines/rna_expression_rbp/reports/rna_expression_rbp_v1.metrics.json"))
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--source-version", default="ENCODE:rolling")
    p.add_argument("--fetch-date", default=FETCH_DATE)
    p.add_argument("--check-inputs", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    expr_path = file_exists_with_alt(args.expression_input)
    rbp_path = file_exists_with_alt(args.rbp_input)

    base = {
        "pipeline": "rna_expression_rbp_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": {
            "rna_master": str(args.rna_master),
            "expression_input": str(expr_path) if expr_path else str(args.expression_input),
            "rbp_input": str(rbp_path) if rbp_path else str(args.rbp_input),
        },
    }

    missing = []
    if not args.rna_master.exists():
        missing.append(str(args.rna_master))
    if expr_path is None:
        missing.append(str(args.expression_input))
    if rbp_path is None:
        missing.append(str(args.rbp_input))

    if missing:
        blocked = {
            **base,
            "status": "blocked_missing_inputs",
            "missing_required": missing,
            "download_checklist": make_download_checklist(args.rna_master, args.expression_input, args.rbp_input),
            "message": "缺外部数据，按协作约束已停止执行。",
        }
        write_json(args.report, blocked)
        print(f"[BLOCKED] missing required inputs: {missing}")
        print(f"[BLOCKED] report -> {args.report}")
        return 2

    master_ids = read_rna_master_ids(args.rna_master)

    if args.check_inputs:
        ready = {
            **base,
            "status": "inputs_ready",
            "rna_master_ids": len(master_ids),
            "message": "输入检查通过。",
        }
        write_json(args.report, ready)
        print(f"[OK] inputs ready -> {args.report}")
        return 0

    expr_rows = parse_expression_rows(expr_path, args.limit)
    rbp_rows = parse_rbp_rows(rbp_path, args.limit)

    write_tsv(args.expression_output, expr_rows, EXPR_OUTPUT_COLUMNS, args.source_version, args.fetch_date)
    write_tsv(args.rbp_output, rbp_rows, RBP_OUTPUT_COLUMNS, args.source_version, args.fetch_date)

    expr_join_rate, expr_matched, expr_total = join_rate(expr_rows, master_ids)
    rbp_join_rate, rbp_matched, rbp_total = join_rate(rbp_rows, master_ids)

    expr_metadata_cols = ["biosample_tissue", "biosample_cell_type", "assay", "sample_id", "experiment_id", "source_dataset", "source_file"]
    rbp_metadata_cols = ["rbp_symbol", "biosample_cell_type", "assay", "sample_id", "experiment_id", "source_dataset", "source_file"]

    expr_level_dist = Counter(lower(r.get("evidence_level", "")) for r in expr_rows)
    rbp_level_dist = Counter(lower(r.get("evidence_level", "")) for r in rbp_rows)

    report = {
        **base,
        "status": "completed",
        "output": {
            "expression": str(args.expression_output),
            "rbp": str(args.rbp_output),
        },
        "row_count": {
            "expression_rows": len(expr_rows),
            "rbp_rows": len(rbp_rows),
        },
        "acceptance": {
            "rna_id_join_rate_expression": expr_join_rate,
            "rna_id_join_rate_rbp": rbp_join_rate,
            "metadata_completeness_expression": completeness_rate(expr_rows, expr_metadata_cols),
            "metadata_completeness_rbp": completeness_rate(rbp_rows, rbp_metadata_cols),
            "evidence_level_non_empty_rate_expression": level_non_empty_rate(expr_rows),
            "evidence_level_non_empty_rate_rbp": level_non_empty_rate(rbp_rows),
        },
        "details": {
            "expression": {
                "join": {"matched": expr_matched, "total": expr_total, "rate": expr_join_rate},
                "evidence_level_distribution": dict(expr_level_dist),
            },
            "rbp": {
                "join": {"matched": rbp_matched, "total": rbp_total, "rate": rbp_join_rate},
                "evidence_level_distribution": dict(rbp_level_dist),
            },
        },
    }

    write_json(args.report, report)

    print(
        f"[OK] expression_rows={len(expr_rows)} rbp_rows={len(rbp_rows)} "
        f"expr_join={expr_join_rate:.4f} rbp_join={rbp_join_rate:.4f}"
    )
    print(f"[OK] report -> {args.report}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
