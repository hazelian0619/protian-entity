#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import pandas as pd

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "Biopython is required. Install with: python3 -m pip install biopython"
    ) from exc


INPUT_COLUMNS = ["uniprot_id", "sequence", "sequence_len"]
OUTPUT_COLUMNS = [
    "uniprot_id",
    "sequence_len",
    "mass_recalc",
    "isoelectric_point",
    "gravy",
    "aromaticity",
    "instability_index",
]
NUMERIC_COLUMNS = [
    "sequence_len",
    "mass_recalc",
    "isoelectric_point",
    "gravy",
    "aromaticity",
    "instability_index",
]

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
AA_SUBSTITUTIONS = {
    "U": "C",  # Selenocysteine -> Cysteine (closest canonical fallback)
    "O": "K",  # Pyrrolysine -> Lysine (closest canonical fallback)
}


def sanitize_sequence(sequence: str) -> str:
    return str(sequence).strip().upper()


def prepare_sequence_for_analysis(sequence: str) -> Tuple[str, List[str]]:
    seq = sanitize_sequence(sequence)
    substitutions_used: List[str] = []
    normalized_chars: List[str] = []
    for ch in seq:
        if ch in CANONICAL_AA:
            normalized_chars.append(ch)
            continue
        if ch in AA_SUBSTITUTIONS:
            normalized_chars.append(AA_SUBSTITUTIONS[ch])
            substitutions_used.append(f"{ch}->{AA_SUBSTITUTIONS[ch]}")
            continue
        raise ValueError(f"Unsupported amino acid '{ch}'")
    return "".join(normalized_chars), substitutions_used


def compute_physchem(sequence: str) -> Dict[str, float]:
    seq = sanitize_sequence(sequence)
    seq_for_analysis, _ = prepare_sequence_for_analysis(seq)
    analyzer = ProteinAnalysis(seq_for_analysis)
    return {
        "sequence_len": float(len(seq)),
        "mass_recalc": float(analyzer.molecular_weight()),
        "isoelectric_point": float(analyzer.isoelectric_point()),
        "gravy": float(analyzer.gravy()),
        "aromaticity": float(analyzer.aromaticity()),
        "instability_index": float(analyzer.instability_index()),
    }


def calculate_numeric_parse_rates(df: pd.DataFrame) -> Dict[str, float]:
    rates: Dict[str, float] = {}
    for col in NUMERIC_COLUMNS:
        parsed = pd.to_numeric(df[col], errors="coerce")
        rates[col] = float(parsed.notna().mean())
    return rates


def build_physchem_table(input_df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    substitution_counts: Dict[str, int] = {}
    error_rows: List[Dict[str, str]] = []

    for row in input_df.itertuples(index=False):
        uniprot_id = str(row.uniprot_id)
        seq_raw = sanitize_sequence(row.sequence)
        try:
            seq_for_analysis, substitutions = prepare_sequence_for_analysis(seq_raw)
            analyzer = ProteinAnalysis(seq_for_analysis)
            for item in substitutions:
                substitution_counts[item] = substitution_counts.get(item, 0) + 1
            records.append(
                {
                    "uniprot_id": uniprot_id,
                    "sequence_len": int(len(seq_raw)),
                    "mass_recalc": float(analyzer.molecular_weight()),
                    "isoelectric_point": float(analyzer.isoelectric_point()),
                    "gravy": float(analyzer.gravy()),
                    "aromaticity": float(analyzer.aromaticity()),
                    "instability_index": float(analyzer.instability_index()),
                }
            )
        except Exception as exc:  # pragma: no cover
            error_rows.append({"uniprot_id": uniprot_id, "error": str(exc)})

    out_df = pd.DataFrame(records, columns=OUTPUT_COLUMNS)

    metadata = {
        "substitution_counts": substitution_counts,
        "error_rows": error_rows,
    }
    return out_df, metadata


def build_quality_report(
    input_df: pd.DataFrame,
    output_df: pd.DataFrame,
    build_meta: Dict[str, Any],
    input_path: Path,
    output_path: Path,
) -> Dict[str, Any]:
    row_count_input = int(len(input_df))
    row_count_output = int(len(output_df))

    input_ids = set(input_df["uniprot_id"].astype(str))
    output_ids = set(output_df["uniprot_id"].astype(str))
    covered_ids = input_ids & output_ids
    missing_ids = sorted(input_ids - output_ids)

    coverage_rate = (len(covered_ids) / len(input_ids)) if input_ids else 0.0

    merged = input_df[["uniprot_id", "sequence_len"]].merge(
        output_df[["uniprot_id", "sequence_len"]],
        on="uniprot_id",
        how="left",
        suffixes=("_master", "_recalc"),
    )

    mismatch_mask = merged["sequence_len_master"].astype("Int64") != merged[
        "sequence_len_recalc"
    ].astype("Int64")
    mismatch_df = merged[mismatch_mask]
    mismatch_count = int(len(mismatch_df))
    mismatch_rate = (mismatch_count / row_count_input) if row_count_input else 0.0

    numeric_rates = calculate_numeric_parse_rates(output_df)
    numeric_100 = all(abs(rate - 1.0) < 1e-12 for rate in numeric_rates.values())

    if mismatch_count == 0:
        diff_explanation = (
            "No sequence_len difference found. sequence_len was recalculated as len(sequence) and "
            "fully matches protein_master_v6_clean.tsv."
        )
    else:
        diff_explanation = (
            "Differences indicate sequence normalization issues (e.g., whitespace/hidden chars) or upstream "
            "sequence_len mismatch in master table; inspect samples in this report."
        )

    report = {
        "pipeline": "protein_physchem_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(input_path),
        "output": str(output_path),
        "row_count": {
            "input": row_count_input,
            "output": row_count_output,
        },
        "coverage": {
            "expected_uniprot_ids": len(input_ids),
            "covered_uniprot_ids": len(covered_ids),
            "coverage_rate": coverage_rate,
            "missing_uniprot_ids_top20": missing_ids[:20],
        },
        "numeric_parse_rates": numeric_rates,
        "sequence_len_diff": {
            "compared_rows": row_count_input,
            "mismatch_rows": mismatch_count,
            "mismatch_rate": mismatch_rate,
            "explanation": diff_explanation,
            "samples_top20": mismatch_df.head(20).to_dict(orient="records"),
        },
        "special_residue_substitution": {
            "mapping": AA_SUBSTITUTIONS,
            "counts": build_meta.get("substitution_counts", {}),
            "note": "U/O are substituted before ProtParam calculation to keep 100% parseability.",
        },
        "errors": {
            "error_count": len(build_meta.get("error_rows", [])),
            "samples_top20": build_meta.get("error_rows", [])[:20],
        },
        "acceptance": {
            "uniprot_id_coverage_100": abs(coverage_rate - 1.0) < 1e-12,
            "numeric_parse_100": numeric_100,
            "sequence_len_diff_reported": True,
        },
    }
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build protein physicochemical feature table.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/processed/protein_master_v6_clean.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/output/protein/protein_physchem_v1.tsv"),
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("pipelines/protein_physchem/reports/protein_physchem_v1.metrics.json"),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Use first N rows for dry run/sample test.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if not args.input.exists():
        raise SystemExit(f"Input file not found: {args.input}")

    input_df = pd.read_csv(args.input, sep="\t", dtype={"uniprot_id": "string"})
    missing_input_cols = [c for c in INPUT_COLUMNS if c not in input_df.columns]
    if missing_input_cols:
        raise SystemExit(f"Missing required input columns: {missing_input_cols}")

    input_df = input_df[INPUT_COLUMNS].copy()
    input_df["sequence"] = input_df["sequence"].astype(str)

    if args.limit is not None:
        input_df = input_df.head(args.limit).copy()

    output_df, build_meta = build_physchem_table(input_df)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    output_df.to_csv(args.output, sep="\t", index=False, float_format="%.6f")

    report = build_quality_report(
        input_df=input_df,
        output_df=output_df,
        build_meta=build_meta,
        input_path=args.input,
        output_path=args.output,
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    ok = report["acceptance"]["uniprot_id_coverage_100"] and report["acceptance"]["numeric_parse_100"]
    status = "PASS" if ok else "FAIL"
    print(
        f"[{status}] protein_physchem_v1 rows={report['row_count']['output']} "
        f"coverage={report['coverage']['coverage_rate']:.4f} "
        f"numeric_parse_100={report['acceptance']['numeric_parse_100']}"
    )
    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
