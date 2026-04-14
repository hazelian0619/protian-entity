#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd


MASTER_REQUIRED_COLUMNS = [
    "uniprot_id",
    "isoforms",
    "source",
    "fetch_date",
    "date_modified",
]

ISOFORM_COLUMNS = [
    "isoform_record_id",
    "canonical_uniprot_id",
    "isoform_uniprot_id",
    "boundary_class",
    "boundary_reason",
    "mapping_scope",
    "relation_type",
    "event_types",
    "named_isoforms_declared",
    "isoform_name_token",
    "isoform_synonyms",
    "isoid_raw",
    "sequence_status_raw",
    "sequence_status_class",
    "is_displayed_sequence",
    "isoform_block_raw",
    "source",
    "source_version",
    "fetch_date",
]

ISOFORM_MAP_COLUMNS = [
    "map_id",
    "canonical_uniprot_id",
    "isoform_uniprot_id",
    "boundary_class",
    "relation_type",
    "mapping_scope",
    "sequence_status_class",
    "source",
    "source_version",
    "fetch_date",
]

ISOID_RE = re.compile(r"IsoId=([^;]+);")
EVENT_RE = re.compile(r"Event=([^;]+);")
DECLARED_RE = re.compile(r"Named isoforms=([0-9]+);")
NAME_BLOCK_RE = re.compile(r"Name=([^;]+);(.*?)(?=(?:Name=[^;]+;)|$)", re.DOTALL)
PREFIX_RE = re.compile(r"^ALTERNATIVE PRODUCTS:\s*", re.IGNORECASE)


@dataclass
class ParseIssue:
    canonical_uniprot_id: str
    reason: str
    raw_text: str


def _clean_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    return str(value).strip()


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

    return _clean_text(matched.iloc[0]) or "protein_master_v6_clean.tsv|source_version=unknown"


def _extract_first_field(block: str, field: str) -> str:
    m = re.search(rf"{re.escape(field)}=([^;]+);", block)
    return _clean_text(m.group(1)) if m else ""


def _split_multi_values(value: str) -> List[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def _classify_sequence_status(sequence_status_raw: str) -> Tuple[str, bool]:
    raw = _clean_text(sequence_status_raw)
    if raw.lower() == "displayed":
        return "displayed", True
    if raw.startswith("VSP_") or "VSP_" in raw:
        return "variant_vsp", False
    if raw == "":
        return "missing", False
    return "other", False


def _mapping_scope(canonical_uniprot_id: str, isoform_uniprot_id: str) -> str:
    return "self_accession" if isoform_uniprot_id.startswith(f"{canonical_uniprot_id}-") else "cross_accession"


def _relation_type(scope: str, is_displayed_sequence: bool) -> str:
    if scope == "self_accession" and is_displayed_sequence:
        return "canonical_displayed"
    if scope == "self_accession" and not is_displayed_sequence:
        return "canonical_isoform"
    if scope == "cross_accession" and is_displayed_sequence:
        return "cross_accession_displayed"
    return "cross_accession_isoform"


def parse_isoform_records(
    canonical_uniprot_id: str,
    isoforms_raw: str,
    source: str,
    source_version: str,
    fetch_date: str,
) -> Tuple[List[Dict[str, object]], List[ParseIssue], Dict[str, object]]:
    payload = PREFIX_RE.sub("", isoforms_raw).strip()
    if not payload:
        return [], [], {
            "declared_named_isoforms": None,
            "event_types": "",
            "isoid_token_count": 0,
            "blocks_with_name": 0,
            "rows_with_isoid": False,
        }

    event_types = sorted({_clean_text(x) for x in EVENT_RE.findall(payload) if _clean_text(x)})
    event_types_str = "|".join(event_types)

    declared_named_isoforms: Optional[int] = None
    declared_match = DECLARED_RE.search(payload)
    if declared_match:
        declared_named_isoforms = int(declared_match.group(1))

    name_blocks = list(NAME_BLOCK_RE.finditer(payload))
    records: List[Dict[str, object]] = []
    issues: List[ParseIssue] = []

    for block_match in name_blocks:
        isoform_name_token = _clean_text(block_match.group(1))
        block_tail = block_match.group(2)
        isoform_block_raw = f"Name={isoform_name_token};{block_tail}".strip()

        synonyms = _extract_first_field(block_tail, "Synonyms")
        isoid_raw = _extract_first_field(block_tail, "IsoId")
        sequence_status_raw = _extract_first_field(block_tail, "Sequence")

        if not isoid_raw:
            issues.append(
                ParseIssue(
                    canonical_uniprot_id=canonical_uniprot_id,
                    reason="missing_isoid_in_name_block",
                    raw_text=isoform_block_raw[:400],
                )
            )
            continue

        sequence_status_class, is_displayed = _classify_sequence_status(sequence_status_raw)
        for isoform_uniprot_id in _split_multi_values(isoid_raw):
            scope = _mapping_scope(canonical_uniprot_id, isoform_uniprot_id)
            relation_type = _relation_type(scope, is_displayed)
            boundary_class = "canonical" if is_displayed else "isoform"
            boundary_reason = "sequence_status_displayed" if is_displayed else "sequence_status_non_displayed"
            isoform_record_id = f"{canonical_uniprot_id}|{isoform_uniprot_id}"

            records.append(
                {
                    "isoform_record_id": isoform_record_id,
                    "canonical_uniprot_id": canonical_uniprot_id,
                    "isoform_uniprot_id": isoform_uniprot_id,
                    "boundary_class": boundary_class,
                    "boundary_reason": boundary_reason,
                    "mapping_scope": scope,
                    "relation_type": relation_type,
                    "event_types": event_types_str,
                    "named_isoforms_declared": declared_named_isoforms,
                    "isoform_name_token": isoform_name_token,
                    "isoform_synonyms": synonyms,
                    "isoid_raw": isoid_raw,
                    "sequence_status_raw": sequence_status_raw,
                    "sequence_status_class": sequence_status_class,
                    "is_displayed_sequence": "TRUE" if is_displayed else "FALSE",
                    "isoform_block_raw": isoform_block_raw,
                    "source": source,
                    "source_version": source_version,
                    "fetch_date": fetch_date,
                }
            )

    meta = {
        "declared_named_isoforms": declared_named_isoforms,
        "event_types": event_types_str,
        "isoid_token_count": len(ISOID_RE.findall(payload)),
        "blocks_with_name": len(name_blocks),
        "rows_with_isoid": len(records) > 0,
    }
    return records, issues, meta


def build_isoform_tables(master_df: pd.DataFrame, source_version: str) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    isoform_records: List[Dict[str, object]] = []
    issues: List[ParseIssue] = []

    rows_with_isoid = 0
    placeholder_rows = 0
    declared_mismatch_rows = 0

    for row in master_df.itertuples(index=False):
        canonical_uniprot_id = _clean_text(row.uniprot_id)
        isoforms_raw = _clean_text(row.isoforms)
        source = _clean_text(row.source)
        fetch_date = _clean_text(row.fetch_date)

        if isoforms_raw.upper().startswith("ALTERNATIVE PRODUCTS") and "IsoId=" not in isoforms_raw:
            placeholder_rows += 1

        records, parse_issues, meta = parse_isoform_records(
            canonical_uniprot_id=canonical_uniprot_id,
            isoforms_raw=isoforms_raw,
            source=source,
            source_version=source_version,
            fetch_date=fetch_date,
        )

        if meta["rows_with_isoid"]:
            rows_with_isoid += 1

        declared = meta["declared_named_isoforms"]
        if declared is not None and meta["blocks_with_name"] and declared != meta["blocks_with_name"]:
            declared_mismatch_rows += 1

        isoform_records.extend(records)
        issues.extend(parse_issues)

    isoform_df = pd.DataFrame(isoform_records, columns=ISOFORM_COLUMNS)
    if not isoform_df.empty:
        isoform_df = isoform_df.sort_values(["canonical_uniprot_id", "isoform_uniprot_id"], kind="stable")

    isoform_map_df = pd.DataFrame(
        [
            {
                "map_id": r["isoform_record_id"],
                "canonical_uniprot_id": r["canonical_uniprot_id"],
                "isoform_uniprot_id": r["isoform_uniprot_id"],
                "boundary_class": r["boundary_class"],
                "relation_type": r["relation_type"],
                "mapping_scope": r["mapping_scope"],
                "sequence_status_class": r["sequence_status_class"],
                "source": r["source"],
                "source_version": r["source_version"],
                "fetch_date": r["fetch_date"],
            }
            for r in isoform_records
        ],
        columns=ISOFORM_MAP_COLUMNS,
    )
    if not isoform_map_df.empty:
        isoform_map_df = isoform_map_df.sort_values(["canonical_uniprot_id", "isoform_uniprot_id"], kind="stable")

    meta = {
        "rows_with_isoid": rows_with_isoid,
        "placeholder_rows": placeholder_rows,
        "declared_mismatch_rows": declared_mismatch_rows,
        "parse_issues": issues,
    }
    return isoform_df, isoform_map_df, meta


def _non_empty_rate(df: pd.DataFrame, column: str) -> float:
    if len(df) == 0:
        return 0.0
    return float(df[column].fillna("").astype(str).str.strip().ne("").mean())


def _duplicates(df: pd.DataFrame, column: str) -> int:
    if df.empty:
        return 0
    return int(df.duplicated(subset=[column]).sum())


def build_qa_reports(
    master_df: pd.DataFrame,
    isoform_df: pd.DataFrame,
    isoform_map_df: pd.DataFrame,
    meta: Dict[str, object],
    input_path: Path,
    isoform_output_path: Path,
    map_output_path: Path,
) -> Tuple[Dict[str, object], Dict[str, object]]:
    parse_issues: List[ParseIssue] = meta.get("parse_issues", [])

    canonical_total = int(len(master_df))
    canonical_with_isoid = int(meta.get("rows_with_isoid", 0))

    canonical_with_output = int(isoform_df["canonical_uniprot_id"].nunique()) if not isoform_df.empty else 0
    coverage_rate = (canonical_with_output / canonical_with_isoid) if canonical_with_isoid else 1.0

    isoform_record_duplicates = _duplicates(isoform_df, "isoform_record_id")
    map_id_duplicates = _duplicates(isoform_map_df, "map_id")

    relation_dist = (
        isoform_map_df["relation_type"].value_counts(dropna=False).to_dict() if not isoform_map_df.empty else {}
    )
    boundary_dist = (
        isoform_df["boundary_class"].value_counts(dropna=False).to_dict() if not isoform_df.empty else {}
    )

    cross_accession_rows = (
        int((isoform_map_df["mapping_scope"] == "cross_accession").sum()) if not isoform_map_df.empty else 0
    )

    isoform_qa = {
        "pipeline": "protein_isoform_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(input_path),
        "output": str(isoform_output_path),
        "scope": {
            "canonical_definition": "Rows in protein_master_v6_clean.tsv keyed by uniprot_id (no '-' suffix).",
            "isoform_definition": "IsoId tokens parsed from master.isoforms field; one row per canonical_uniprot_id x isoform_uniprot_id.",
            "boundary_rule": "sequence_status_raw=Displayed -> boundary_class=canonical; other sequence_status -> boundary_class=isoform.",
            "mapping_rule": "isoform_uniprot_id startswith canonical_uniprot_id + '-' => self_accession, else cross_accession.",
            "missing_strategy": "Rows without IsoId (blank/placeholder) are not emitted; counts are disclosed in QA under missing_isoform_annotation.",
        },
        "row_count": {
            "canonical_input": canonical_total,
            "canonical_with_isoid_input": canonical_with_isoid,
            "isoform_output": int(len(isoform_df)),
            "isoform_map_output": int(len(isoform_map_df)),
        },
        "coverage": {
            "canonical_with_isoid_output": canonical_with_output,
            "coverage_rate": coverage_rate,
        },
        "missing_isoform_annotation": {
            "placeholder_rows": int(meta.get("placeholder_rows", 0)),
            "null_or_empty_rows": int(
                master_df["isoforms"].fillna("").astype(str).str.strip().eq("").sum()
            ),
            "total_without_isoid": int(canonical_total - canonical_with_isoid),
        },
        "distribution": {
            "boundary_class": boundary_dist,
            "relation_type": relation_dist,
            "cross_accession_rows": cross_accession_rows,
            "cross_accession_rate": (cross_accession_rows / len(isoform_map_df)) if len(isoform_map_df) else 0.0,
        },
        "quality": {
            "isoform_record_id_duplicates": isoform_record_duplicates,
            "map_id_duplicates": map_id_duplicates,
            "source_non_empty_rate": _non_empty_rate(isoform_df, "source") if not isoform_df.empty else 0.0,
            "source_version_non_empty_rate": _non_empty_rate(isoform_df, "source_version") if not isoform_df.empty else 0.0,
            "fetch_date_non_empty_rate": _non_empty_rate(isoform_df, "fetch_date") if not isoform_df.empty else 0.0,
            "declared_named_isoforms_mismatch_rows": int(meta.get("declared_mismatch_rows", 0)),
            "parse_issue_count": len(parse_issues),
            "parse_issue_samples_top20": [
                {
                    "canonical_uniprot_id": x.canonical_uniprot_id,
                    "reason": x.reason,
                    "raw_text": x.raw_text,
                }
                for x in parse_issues[:20]
            ],
        },
    }

    isoform_qa["acceptance"] = {
        "coverage_100_on_rows_with_isoid": abs(coverage_rate - 1.0) < 1e-12,
        "isoform_record_id_unique": isoform_record_duplicates == 0,
        "map_id_unique": map_id_duplicates == 0,
        "provenance_complete": all(
            abs(isoform_qa["quality"][k] - 1.0) < 1e-12
            for k in [
                "source_non_empty_rate",
                "source_version_non_empty_rate",
                "fetch_date_non_empty_rate",
            ]
        ),
        "parse_issue_zero": len(parse_issues) == 0,
    }
    isoform_qa["passed"] = all(isoform_qa["acceptance"].values())

    isoform_map_qa = {
        "pipeline": "protein_isoform_map_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(input_path),
        "output": str(map_output_path),
        "row_count": int(len(isoform_map_df)),
        "distribution": {
            "relation_type": relation_dist,
            "mapping_scope": (
                isoform_map_df["mapping_scope"].value_counts(dropna=False).to_dict()
                if not isoform_map_df.empty
                else {}
            ),
            "sequence_status_class": (
                isoform_map_df["sequence_status_class"].value_counts(dropna=False).to_dict()
                if not isoform_map_df.empty
                else {}
            ),
        },
        "quality": {
            "map_id_duplicates": map_id_duplicates,
            "canonical_uniprot_non_empty_rate": _non_empty_rate(isoform_map_df, "canonical_uniprot_id")
            if not isoform_map_df.empty
            else 0.0,
            "isoform_uniprot_non_empty_rate": _non_empty_rate(isoform_map_df, "isoform_uniprot_id")
            if not isoform_map_df.empty
            else 0.0,
            "source_non_empty_rate": _non_empty_rate(isoform_map_df, "source") if not isoform_map_df.empty else 0.0,
            "source_version_non_empty_rate": _non_empty_rate(isoform_map_df, "source_version")
            if not isoform_map_df.empty
            else 0.0,
            "fetch_date_non_empty_rate": _non_empty_rate(isoform_map_df, "fetch_date")
            if not isoform_map_df.empty
            else 0.0,
        },
    }

    isoform_map_qa["acceptance"] = {
        "map_id_unique": map_id_duplicates == 0,
        "canonical_uniprot_complete": abs(isoform_map_qa["quality"]["canonical_uniprot_non_empty_rate"] - 1.0) < 1e-12,
        "isoform_uniprot_complete": abs(isoform_map_qa["quality"]["isoform_uniprot_non_empty_rate"] - 1.0) < 1e-12,
        "provenance_complete": all(
            abs(isoform_map_qa["quality"][k] - 1.0) < 1e-12
            for k in ["source_non_empty_rate", "source_version_non_empty_rate", "fetch_date_non_empty_rate"]
        ),
    }
    isoform_map_qa["passed"] = all(isoform_map_qa["acceptance"].values())

    return isoform_qa, isoform_map_qa


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Protein isoform normalized layer and mapping table.")
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
        "--output-isoform",
        type=Path,
        default=Path("data/output/protein/protein_isoform_v1.tsv"),
    )
    parser.add_argument(
        "--output-map",
        type=Path,
        default=Path("data/output/protein/protein_isoform_map_v1.tsv"),
    )
    parser.add_argument(
        "--qa-isoform",
        type=Path,
        default=Path("pipelines/protein_isoform/reports/protein_isoform_v1.qa.json"),
    )
    parser.add_argument(
        "--qa-map",
        type=Path,
        default=Path("pipelines/protein_isoform/reports/protein_isoform_map_v1.qa.json"),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Use only first N input rows (for smoke tests).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if not args.input.exists():
        raise SystemExit(f"Input file not found: {args.input}")

    master_df = pd.read_csv(args.input, sep="\t", dtype=str)
    missing_cols = [c for c in MASTER_REQUIRED_COLUMNS if c not in master_df.columns]
    if missing_cols:
        raise SystemExit(f"Missing required input columns: {missing_cols}")

    if args.limit is not None:
        master_df = master_df.head(args.limit).copy()

    source_version = _derive_source_version(args.source_versions)
    isoform_df, isoform_map_df, meta = build_isoform_tables(master_df, source_version)

    args.output_isoform.parent.mkdir(parents=True, exist_ok=True)
    args.output_map.parent.mkdir(parents=True, exist_ok=True)
    isoform_df.to_csv(args.output_isoform, sep="\t", index=False)
    isoform_map_df.to_csv(args.output_map, sep="\t", index=False)

    qa_isoform, qa_map = build_qa_reports(
        master_df=master_df,
        isoform_df=isoform_df,
        isoform_map_df=isoform_map_df,
        meta=meta,
        input_path=args.input,
        isoform_output_path=args.output_isoform,
        map_output_path=args.output_map,
    )

    args.qa_isoform.parent.mkdir(parents=True, exist_ok=True)
    args.qa_map.parent.mkdir(parents=True, exist_ok=True)
    args.qa_isoform.write_text(json.dumps(qa_isoform, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    args.qa_map.write_text(json.dumps(qa_map, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    status_isoform = "PASS" if qa_isoform["passed"] else "FAIL"
    status_map = "PASS" if qa_map["passed"] else "FAIL"
    print(
        f"[{status_isoform}] protein_isoform_v1 rows={len(isoform_df)} coverage={qa_isoform['coverage']['coverage_rate']:.4f} "
        f"parse_issues={qa_isoform['quality']['parse_issue_count']}"
    )
    print(
        f"[{status_map}] protein_isoform_map_v1 rows={len(isoform_map_df)} unique={qa_map['acceptance']['map_id_unique']}"
    )

    return 0 if qa_isoform["passed"] and qa_map["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
