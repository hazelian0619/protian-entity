#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
ZINC_ID_RE = re.compile(r"^ZINC\d+$")

BACKFILL_SOURCE_TAG = "molecule_structure_identifiers_v1.canonical_smiles"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def utc_today() -> str:
    return datetime.now(timezone.utc).date().isoformat()


def normalize(x: str) -> str:
    return (x or "").strip()


def split_multi(x: str) -> List[str]:
    text = normalize(x)
    if not text:
        return []
    return [p.strip() for p in text.split(";") if p.strip()]


def append_chain(existing: str, token: str) -> str:
    token = normalize(token)
    values = split_multi(existing)
    if not token:
        return ";".join(values)
    if token not in values:
        values.append(token)
    return ";".join(values)


def valid_inchikey(x: str) -> bool:
    return bool(INCHIKEY_RE.match(normalize(x).upper()))


def valid_zinc_id(x: str) -> bool:
    return bool(ZINC_ID_RE.match(normalize(x).upper()))


def read_tsv(path: Path, max_rows: Optional[int] = None) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        rows: List[Dict[str, str]] = []
        for i, row in enumerate(reader, start=1):
            rows.append(row)
            if max_rows is not None and i >= max_rows:
                break
    return list(reader.fieldnames), rows


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Dict[str, str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    n = 0
    with tmp.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=list(header), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in header})
            n += 1
    tmp.replace(path)
    return n


def non_empty_rate(rows: Sequence[Dict[str, str]], col: str) -> Dict[str, Any]:
    total = len(rows)
    non_empty = sum(1 for r in rows if normalize(r.get(col, "")))
    return {
        "non_empty": non_empty,
        "total": total,
        "rate": (non_empty / total) if total else 0.0,
    }


def key_mapping(row: Dict[str, str]) -> Tuple[str, str, str, str, str]:
    return (
        normalize(row.get("inchikey", "")).upper(),
        normalize(row.get("zinc_id", "")).upper(),
        normalize(row.get("supplier_name", "")),
        normalize(row.get("supplier_code", "")),
        normalize(row.get("catalog_tier", "")),
    )


def load_structure_smiles(path: Path) -> Tuple[Dict[str, str], Dict[str, Any]]:
    cols, rows = read_tsv(path)
    if "inchikey" not in cols or "canonical_smiles" not in cols:
        raise SystemExit(f"[ERROR] structure identifiers missing columns inchikey/canonical_smiles: {path}")

    smiles_by_ik: Dict[str, str] = {}
    duplicate_inchikey = 0
    conflicting_smiles = 0
    bad_inchikey_rows = 0
    structure_fetch_dates: Set[str] = set()

    for r in rows:
        ik = normalize(r.get("inchikey", "")).upper()
        smiles = normalize(r.get("canonical_smiles", ""))
        fd = normalize(r.get("fetch_date", ""))
        if fd:
            structure_fetch_dates.add(fd)

        if not valid_inchikey(ik):
            bad_inchikey_rows += 1
            continue

        if ik in smiles_by_ik:
            duplicate_inchikey += 1
            if smiles and smiles_by_ik[ik] and smiles != smiles_by_ik[ik]:
                conflicting_smiles += 1
            if not smiles_by_ik[ik] and smiles:
                smiles_by_ik[ik] = smiles
            continue

        smiles_by_ik[ik] = smiles

    stats = {
        "rows": len(rows),
        "valid_inchikey_rows": len(smiles_by_ik),
        "bad_inchikey_rows": bad_inchikey_rows,
        "duplicate_inchikey_rows": duplicate_inchikey,
        "conflicting_smiles_for_same_inchikey": conflicting_smiles,
        "structure_fetch_dates": sorted(structure_fetch_dates),
        "canonical_smiles_non_empty": sum(1 for v in smiles_by_ik.values() if normalize(v)),
    }
    return smiles_by_ik, stats


def build_conflict_metrics(rows: Sequence[Dict[str, str]]) -> Dict[str, Any]:
    inchikey_to_zinc: Dict[str, Set[str]] = defaultdict(set)
    inchikey_to_vendor: Dict[str, Set[str]] = defaultdict(set)
    zinc_to_inchikey: Dict[str, Set[str]] = defaultdict(set)

    for r in rows:
        ik = normalize(r.get("inchikey", "")).upper()
        zid = normalize(r.get("zinc_id", "")).upper()
        vendor = normalize(r.get("supplier_name", ""))
        if ik:
            inchikey_to_zinc[ik].add(zid)
            inchikey_to_vendor[ik].add(vendor)
        if zid:
            zinc_to_inchikey[zid].add(ik)

    multi_zinc_inchikey = {ik: len(v) for ik, v in inchikey_to_zinc.items() if len(v) > 1}
    multi_vendor_inchikey = {ik: len(v) for ik, v in inchikey_to_vendor.items() if len(v) > 1}
    multi_inchikey_zinc = {z: len(v) for z, v in zinc_to_inchikey.items() if len(v) > 1}

    samples = {
        "inchikey_with_multi_zinc": [
            {
                "inchikey": ik,
                "zinc_count": len(inchikey_to_zinc[ik]),
                "zinc_ids": sorted(inchikey_to_zinc[ik])[:10],
                "vendors": sorted(inchikey_to_vendor[ik])[:10],
            }
            for ik in sorted(multi_zinc_inchikey.keys())[:100]
        ],
        "zinc_with_multi_inchikey": [
            {
                "zinc_id": zid,
                "inchikey_count": len(zinc_to_inchikey[zid]),
                "inchikey_sample": sorted(zinc_to_inchikey[zid])[:10],
            }
            for zid in sorted(multi_inchikey_zinc.keys())[:100]
        ],
    }

    metrics = {
        "rows": len(rows),
        "inchikey_with_multi_zinc": len(multi_zinc_inchikey),
        "inchikey_with_multi_vendor": len(multi_vendor_inchikey),
        "zinc_with_multi_inchikey": len(multi_inchikey_zinc),
    }

    return {"metrics": metrics, "samples": samples}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--zinc-v1", type=Path, required=True)
    ap.add_argument("--structure-identifiers", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--coverage-report", type=Path, required=True)
    ap.add_argument("--conflict-report", type=Path, required=True)
    ap.add_argument("--max-rows", type=int, default=None)
    args = ap.parse_args()

    for p in [args.zinc_v1, args.structure_identifiers]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    created_at = utc_now()
    fetch_date = utc_today()

    v1_header, v1_rows = read_tsv(args.zinc_v1, max_rows=args.max_rows)
    required_cols = {
        "inchikey",
        "zinc_id",
        "supplier_code",
        "supplier_name",
        "smiles",
        "catalog_tier",
        "mapping_method",
        "source",
        "source_version",
        "fetch_date",
    }
    missing = sorted(list(required_cols - set(v1_header)))
    if missing:
        raise SystemExit(f"[ERROR] zinc v1 missing required columns: {missing}")

    smiles_by_ik, structure_stats = load_structure_smiles(args.structure_identifiers)

    out_rows: List[Dict[str, str]] = []
    invalid_v1_rows = 0
    invalid_zinc_rows = 0

    backfilled_rows = 0
    retained_v1_smiles_rows = 0
    still_missing_smiles_rows = 0
    fallback_candidate_rows = 0
    source_appended_rows = 0
    source_version_appended_rows = 0
    prefilled_but_structure_diff_rows = 0

    for row in v1_rows:
        out = dict(row)
        ik = normalize(out.get("inchikey", "")).upper()
        zid = normalize(out.get("zinc_id", "")).upper()
        out["inchikey"] = ik
        out["zinc_id"] = zid

        if not valid_inchikey(ik):
            invalid_v1_rows += 1
        if not valid_zinc_id(zid):
            invalid_zinc_rows += 1

        v1_smiles = normalize(out.get("smiles", ""))
        fallback_smiles = normalize(smiles_by_ik.get(ik, ""))

        out["smiles_original"] = v1_smiles

        if v1_smiles:
            retained_v1_smiles_rows += 1
            out["smiles_source"] = "zinc_source_file"
            out["smiles_backfill_flag"] = "0"
            out["smiles_backfill_source"] = ""
            if fallback_smiles and fallback_smiles != v1_smiles:
                prefilled_but_structure_diff_rows += 1
        else:
            if fallback_smiles:
                fallback_candidate_rows += 1
                out["smiles"] = fallback_smiles
                out["smiles_source"] = BACKFILL_SOURCE_TAG
                out["smiles_backfill_flag"] = "1"
                out["smiles_backfill_source"] = BACKFILL_SOURCE_TAG
                backfilled_rows += 1

                source_before = normalize(out.get("source", ""))
                source_after = append_chain(source_before, BACKFILL_SOURCE_TAG)
                if source_after != source_before:
                    source_appended_rows += 1
                out["source"] = source_after

                srcv_before = normalize(out.get("source_version", ""))
                srcv_after = append_chain(srcv_before, "molecule_structure_identifiers_v1.tsv")
                if srcv_after != srcv_before:
                    source_version_appended_rows += 1
                out["source_version"] = srcv_after
            else:
                still_missing_smiles_rows += 1
                out["smiles_source"] = ""
                out["smiles_backfill_flag"] = "0"
                out["smiles_backfill_source"] = ""

        out["fetch_date"] = fetch_date
        out_rows.append(out)

    header = list(v1_header)
    for col in ["smiles_original", "smiles_source", "smiles_backfill_flag", "smiles_backfill_source"]:
        if col not in header:
            # place after smiles for readability
            if col == "smiles_original" and "smiles" in header:
                header.insert(header.index("smiles") + 1, col)
            elif col == "smiles_source" and "smiles_original" in header:
                header.insert(header.index("smiles_original") + 1, col)
            elif col == "smiles_backfill_flag" and "smiles_source" in header:
                header.insert(header.index("smiles_source") + 1, col)
            elif col == "smiles_backfill_source" and "smiles_backfill_flag" in header:
                header.insert(header.index("smiles_backfill_flag") + 1, col)
            else:
                header.append(col)

    rows_written = write_tsv(args.out, header=header, rows=out_rows)

    v1_smiles_rate = non_empty_rate(v1_rows, "smiles")
    v2_smiles_rate = non_empty_rate(out_rows, "smiles")
    smiles_delta = v2_smiles_rate["rate"] - v1_smiles_rate["rate"]

    source_rate = non_empty_rate(out_rows, "source")
    source_version_rate = non_empty_rate(out_rows, "source_version")
    fetch_date_rate = non_empty_rate(out_rows, "fetch_date")

    v1_keys = {key_mapping(r) for r in v1_rows}
    v2_keys = {key_mapping(r) for r in out_rows}

    conflict_v1 = build_conflict_metrics(v1_rows)
    conflict_v2 = build_conflict_metrics(out_rows)

    conflict_report = {
        "name": "molecule_zinc_xref_v2.conflict_audit",
        "created_at": created_at,
        "comparison": {
            "v1_metrics": conflict_v1["metrics"],
            "v2_metrics": conflict_v2["metrics"],
            "mapping_key_set_equal": v1_keys == v2_keys,
        },
        "metrics": conflict_v2["metrics"],
        "samples": conflict_v2["samples"],
    }
    args.conflict_report.parent.mkdir(parents=True, exist_ok=True)
    args.conflict_report.write_text(json.dumps(conflict_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    coverage_report = {
        "name": "molecule_zinc_xref_v2.coverage",
        "created_at": created_at,
        "rows": rows_written,
        "smiles_non_empty_rate_v1": v1_smiles_rate,
        "smiles_non_empty_rate_v2": v2_smiles_rate,
        "smiles_non_empty_delta": smiles_delta,
        "smiles_non_empty_delta_pp": smiles_delta * 100.0,
        "backfill": {
            "fallback_candidates": fallback_candidate_rows,
            "rows_backfilled": backfilled_rows,
            "rows_still_missing_smiles": still_missing_smiles_rows,
            "rows_prefilled_preserved": retained_v1_smiles_rows,
            "prefilled_rows_with_structure_smiles_diff": prefilled_but_structure_diff_rows,
        },
        "non_empty_rates": {
            "source": source_rate,
            "source_version": source_version_rate,
            "fetch_date": fetch_date_rate,
        },
    }
    args.coverage_report.parent.mkdir(parents=True, exist_ok=True)
    args.coverage_report.write_text(json.dumps(coverage_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    build_report = {
        "name": "molecule_zinc_xref_v2.build",
        "created_at": created_at,
        "sample_mode": args.max_rows is not None,
        "max_rows": args.max_rows,
        "inputs": {
            "zinc_v1": str(args.zinc_v1),
            "structure_identifiers": str(args.structure_identifiers),
        },
        "output": str(args.out),
        "metrics": {
            "rows_input": len(v1_rows),
            "rows_written": rows_written,
            "invalid_v1_inchikey_rows": invalid_v1_rows,
            "invalid_v1_zinc_id_rows": invalid_zinc_rows,
            "smiles_rate_v1": v1_smiles_rate["rate"],
            "smiles_rate_v2": v2_smiles_rate["rate"],
            "smiles_rate_delta": smiles_delta,
            "smiles_rate_delta_pp": smiles_delta * 100.0,
            "rows_backfilled": backfilled_rows,
            "rows_still_missing_smiles": still_missing_smiles_rows,
            "rows_prefilled_preserved": retained_v1_smiles_rows,
            "source_appended_rows": source_appended_rows,
            "source_version_appended_rows": source_version_appended_rows,
            "mapping_key_set_equal_to_v1": (v1_keys == v2_keys),
            "source_non_empty_rate": source_rate["rate"],
            "source_version_non_empty_rate": source_version_rate["rate"],
            "fetch_date_non_empty_rate": fetch_date_rate["rate"],
        },
        "structure_stats": structure_stats,
        "coverage_report": str(args.coverage_report),
        "conflict_report": str(args.conflict_report),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(build_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(
        "[OK] molecule_zinc_xref_v2 rows={} smiles_rate {:.4f} -> {:.4f} (delta_pp={:.2f})".format(
            rows_written,
            v1_smiles_rate["rate"],
            v2_smiles_rate["rate"],
            smiles_delta * 100.0,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
