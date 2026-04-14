#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd

try:
    import Bio
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError as exc:  # pragma: no cover
    raise SystemExit("Biopython is required. Install with: python3 -m pip install biopython") from exc


INPUT_COLUMNS = ["uniprot_id", "sequence", "sequence_len", "source", "fetch_date"]
OUTPUT_COLUMNS = [
    "uniprot_id",
    "sequence_len",
    "mass_recalc",
    "isoelectric_point",
    "gravy",
    "aromaticity",
    "instability_index",
    "extinction_coefficient",
    "extinction_coefficient_reduced",
    "aliphatic_index",
    "source",
    "source_version",
    "fetch_date",
]

CORE_COMPARE_COLUMNS = [
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
    "extinction_coefficient",
    "extinction_coefficient_reduced",
    "aliphatic_index",
]

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
AA_SUBSTITUTIONS = {
    "U": "C",
    "O": "K",
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


def compute_aliphatic_index(sequence: str) -> float:
    if not sequence:
        return float("nan")
    length = len(sequence)
    a = sequence.count("A") / length * 100.0
    v = sequence.count("V") / length * 100.0
    i = sequence.count("I") / length * 100.0
    l = sequence.count("L") / length * 100.0
    return float(a + 2.9 * v + 3.9 * (i + l))


def _derive_source_version(source_versions_path: Path) -> str:
    if not source_versions_path.exists():
        return "protein_master_v6_clean.tsv|source_version=unknown"

    try:
        sv = pd.read_csv(source_versions_path, sep="\t", dtype=str)
    except Exception:
        return "protein_master_v6_clean.tsv|source_version=unknown"

    if "dataset" not in sv.columns or "source_version" not in sv.columns:
        return "protein_master_v6_clean.tsv|source_version=unknown"

    matched = sv.loc[sv["dataset"] == "protein_master_v6_clean.tsv", "source_version"]
    if matched.empty:
        return "protein_master_v6_clean.tsv|source_version=unknown"

    value = str(matched.iloc[0]).strip()
    return value or "protein_master_v6_clean.tsv|source_version=unknown"


def build_extended_table(input_df: pd.DataFrame, source_version: str) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    substitution_counts: Dict[str, int] = {}
    error_rows: List[Dict[str, str]] = []

    for row in input_df.itertuples(index=False):
        uniprot_id = str(row.uniprot_id)
        seq_raw = sanitize_sequence(row.sequence)
        source = str(row.source)
        fetch_date = str(row.fetch_date)

        try:
            seq_for_analysis, substitutions = prepare_sequence_for_analysis(seq_raw)
            analyzer = ProteinAnalysis(seq_for_analysis)
            extinction_reduced, extinction_oxidized = analyzer.molar_extinction_coefficient()
            aliphatic_index = compute_aliphatic_index(seq_for_analysis)

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
                    "extinction_coefficient": float(extinction_oxidized),
                    "extinction_coefficient_reduced": float(extinction_reduced),
                    "aliphatic_index": float(aliphatic_index),
                    "source": source,
                    "source_version": source_version,
                    "fetch_date": fetch_date,
                }
            )
        except Exception as exc:  # pragma: no cover
            error_rows.append({"uniprot_id": uniprot_id, "error": str(exc)})

    out_df = pd.DataFrame(records, columns=OUTPUT_COLUMNS)
    meta = {
        "substitution_counts": substitution_counts,
        "error_rows": error_rows,
    }
    return out_df, meta


def calculate_numeric_parse_rates(df: pd.DataFrame) -> Dict[str, float]:
    rates: Dict[str, float] = {}
    for col in NUMERIC_COLUMNS:
        parsed = pd.to_numeric(df[col], errors="coerce")
        rates[col] = float(parsed.notna().mean())
    return rates


def compare_with_core(core_path: Path, extended_df: pd.DataFrame) -> Dict[str, Any]:
    if not core_path.exists():
        return {
            "available": False,
            "core_path": str(core_path),
            "message": "core protein_physchem_v1 table missing; consistency check skipped",
        }

    core_df = pd.read_csv(core_path, sep="\t", dtype=str)
    missing_cols = [c for c in ["uniprot_id", *CORE_COMPARE_COLUMNS] if c not in core_df.columns]
    if missing_cols:
        return {
            "available": False,
            "core_path": str(core_path),
            "message": f"core table missing columns: {missing_cols}",
        }

    core_trimmed = core_df[["uniprot_id", *CORE_COMPARE_COLUMNS]].copy()
    ext_trimmed = extended_df[["uniprot_id", *CORE_COMPARE_COLUMNS]].copy()

    merged = core_trimmed.merge(ext_trimmed, on="uniprot_id", how="inner", suffixes=("_core", "_extended"))

    mismatch_rows = []
    mismatch_count_by_column: Dict[str, int] = {}
    for col in CORE_COMPARE_COLUMNS:
        core_num = pd.to_numeric(merged[f"{col}_core"], errors="coerce")
        ext_num = pd.to_numeric(merged[f"{col}_extended"], errors="coerce")
        mismatch_mask = (core_num - ext_num).abs() > 1e-6
        mismatch_count = int(mismatch_mask.sum())
        mismatch_count_by_column[col] = mismatch_count

    if len(merged) > 0:
        aggregate_mask = pd.Series(False, index=merged.index)
        for col in CORE_COMPARE_COLUMNS:
            core_num = pd.to_numeric(merged[f"{col}_core"], errors="coerce")
            ext_num = pd.to_numeric(merged[f"{col}_extended"], errors="coerce")
            aggregate_mask = aggregate_mask | ((core_num - ext_num).abs() > 1e-6)
        mismatch_rows_df = merged[aggregate_mask]
        mismatch_rows = mismatch_rows_df.head(20).to_dict(orient="records")
        mismatch_row_count = int(len(mismatch_rows_df))
    else:
        mismatch_row_count = 0

    return {
        "available": True,
        "core_path": str(core_path),
        "joined_rows": int(len(merged)),
        "mismatch_rows": mismatch_row_count,
        "mismatch_count_by_column": mismatch_count_by_column,
        "samples_top20": mismatch_rows,
    }


def build_qa_report(
    input_df: pd.DataFrame,
    output_df: pd.DataFrame,
    build_meta: Dict[str, Any],
    input_path: Path,
    output_path: Path,
    core_compare: Dict[str, Any],
) -> Dict[str, Any]:
    row_count_input = int(len(input_df))
    row_count_output = int(len(output_df))

    input_ids = set(input_df["uniprot_id"].astype(str))
    output_ids = set(output_df["uniprot_id"].astype(str))

    coverage_rate = (len(output_ids & input_ids) / len(input_ids)) if input_ids else 0.0

    numeric_rates = calculate_numeric_parse_rates(output_df)
    numeric_100 = all(abs(rate - 1.0) < 1e-12 for rate in numeric_rates.values())

    provenance_rates = {
        "source_non_empty_rate": float(output_df["source"].astype(str).str.strip().ne("").mean()) if len(output_df) else 0.0,
        "source_version_non_empty_rate": float(output_df["source_version"].astype(str).str.strip().ne("").mean()) if len(output_df) else 0.0,
        "fetch_date_non_empty_rate": float(output_df["fetch_date"].astype(str).str.strip().ne("").mean()) if len(output_df) else 0.0,
        "fetch_date_format_rate": float(
            output_df["fetch_date"].astype(str).str.match(r"^\d{4}-\d{2}-\d{2}$", na=False).mean()
        )
        if len(output_df)
        else 0.0,
    }

    compare_ok = True
    if core_compare.get("available"):
        compare_ok = core_compare.get("mismatch_rows", 0) == 0

    qa = {
        "pipeline": "protein_physchem_extended_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(input_path),
        "output": str(output_path),
        "calculation_definition": {
            "extinction_coefficient": "Molar extinction coefficient at 280nm (M^-1 cm^-1), oxidized cystine form.",
            "extinction_coefficient_reduced": "Molar extinction coefficient at 280nm (M^-1 cm^-1), reduced cysteine form.",
            "aliphatic_index": "Ikai aliphatic index (unitless): X(Ala)+2.9*X(Val)+3.9*(X(Ile)+X(Leu)), X is mole percent.",
            "instability_index": "Biopython ProteinAnalysis.instability_index (unitless).",
            "gravy": "Biopython ProteinAnalysis.gravy (Kyte-Doolittle hydropathicity).",
            "isoelectric_point": "Biopython ProteinAnalysis.isoelectric_point (pH).",
        },
        "unit": {
            "extinction_coefficient": "M^-1 cm^-1",
            "extinction_coefficient_reduced": "M^-1 cm^-1",
            "aliphatic_index": "unitless",
            "instability_index": "unitless",
            "gravy": "unitless",
            "isoelectric_point": "pH",
        },
        "row_count": {
            "input": row_count_input,
            "output": row_count_output,
        },
        "coverage": {
            "expected_uniprot_ids": len(input_ids),
            "covered_uniprot_ids": len(output_ids & input_ids),
            "coverage_rate": coverage_rate,
            "missing_uniprot_ids_top20": sorted(input_ids - output_ids)[:20],
        },
        "numeric_parse_rates": numeric_rates,
        "provenance_rates": provenance_rates,
        "special_residue_substitution": {
            "mapping": AA_SUBSTITUTIONS,
            "counts": build_meta.get("substitution_counts", {}),
            "note": "U/O are substituted before ProtParam calculations to keep 100% parseability.",
        },
        "errors": {
            "error_count": len(build_meta.get("error_rows", [])),
            "samples_top20": build_meta.get("error_rows", [])[:20],
        },
        "core_alignment": core_compare,
        "engine": {
            "biopython_version": getattr(Bio, "__version__", "unknown"),
        },
    }

    qa["acceptance"] = {
        "uniprot_id_coverage_100": abs(coverage_rate - 1.0) < 1e-12,
        "numeric_parse_100": numeric_100,
        "provenance_complete": all(abs(v - 1.0) < 1e-12 for v in provenance_rates.values()),
        "error_rows_zero": qa["errors"]["error_count"] == 0,
        "core_overlap_aligned": compare_ok,
    }
    qa["passed"] = all(qa["acceptance"].values())
    return qa


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Protein physicochemical extended feature table.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/processed/protein_master_v6_clean.tsv"),
    )
    parser.add_argument(
        "--source-versions",
        type=Path,
        default=Path("data/output/protein/protein_source_versions_v1.tsv"),
    )
    parser.add_argument(
        "--core-table",
        type=Path,
        default=Path("data/output/protein/protein_physchem_v1.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/output/protein/protein_physchem_extended_v1.tsv"),
    )
    parser.add_argument(
        "--qa",
        type=Path,
        default=Path("pipelines/protein_physchem_extended/reports/protein_physchem_extended_v1.qa.json"),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Use first N rows for smoke tests.",
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

    source_version = _derive_source_version(args.source_versions)
    output_df, build_meta = build_extended_table(input_df, source_version)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    output_df.to_csv(args.output, sep="\t", index=False, float_format="%.6f")

    core_compare = compare_with_core(args.core_table, output_df)
    qa = build_qa_report(
        input_df=input_df,
        output_df=output_df,
        build_meta=build_meta,
        input_path=args.input,
        output_path=args.output,
        core_compare=core_compare,
    )

    args.qa.parent.mkdir(parents=True, exist_ok=True)
    args.qa.write_text(json.dumps(qa, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status = "PASS" if qa["passed"] else "FAIL"
    print(
        f"[{status}] protein_physchem_extended_v1 rows={qa['row_count']['output']} "
        f"coverage={qa['coverage']['coverage_rate']:.4f} core_overlap_aligned={qa['acceptance']['core_overlap_aligned']}"
    )

    return 0 if qa["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
