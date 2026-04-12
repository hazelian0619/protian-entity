#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

FETCH_DATE = date.today().isoformat()

TRNA_MT_SYMBOL_TO_ANTICODON = {
    "MT-TF": "GAA",
    "MT-TV": "TAC",
    "MT-TL1": "TAA",
    "MT-TL2": "TAG",
    "MT-TI": "GAT",
    "MT-TQ": "TTG",
    "MT-TM": "CAT",
    "MT-TW": "TCA",
    "MT-TA": "TGC",
    "MT-TN": "GTT",
    "MT-TC": "GCA",
    "MT-TY": "GTA",
    "MT-TS1": "GCT",
    "MT-TS2": "TGA",
    "MT-TD": "GTC",
    "MT-TK": "TTT",
    "MT-TG": "TCC",
    "MT-TR": "TCG",
    "MT-TH": "GTG",
    "MT-TE": "TTC",
    "MT-TT": "TGT",
    "MT-TP": "TGG",
}

TRNA_AA_RE = re.compile(r"tRNA-([A-Za-z]+)-([ACGTU]{3})")
TRNA_LABEL_RE = re.compile(r"(tRNA-[A-Za-z]+-([ACGTU]{3})(?:-[0-9]+(?:-[0-9]+)?)?)")
TRNAI_RE = re.compile(r"TRNA[A-Z]-([ACGTU]{3})", re.IGNORECASE)
HGNC_STYLE_RE = re.compile(r"TRN[A-Z]-([ACGTU]{3})", re.IGNORECASE)
MT_SYMBOL_RE = re.compile(r"\b(MT-T[A-Z0-9]{1,3})\b", re.IGNORECASE)


def normalize(v: str) -> str:
    return (v or "").strip()


def open_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def resolve_with_alt(path: Path) -> Optional[Path]:
    if path.exists():
        return path
    if path.suffix == ".gz":
        alt = Path(str(path)[:-3])
        return alt if alt.exists() else None
    alt = Path(str(path) + ".gz")
    return alt if alt.exists() else None


def parse_trna_annotation(external_id: str, symbol_or_aux: str) -> Tuple[str, str, str, str]:
    aa = ""
    anticodon = ""
    label = ""
    for text in [normalize(external_id), normalize(symbol_or_aux)]:
        if text == "":
            continue

        aa_m = TRNA_AA_RE.search(text)
        if aa_m:
            aa = aa or aa_m.group(1)
            anticodon = anticodon or aa_m.group(2).upper()

        label_m = TRNA_LABEL_RE.search(text)
        if label_m:
            label = label or label_m.group(1)
            anticodon = anticodon or label_m.group(2).upper()

        if anticodon == "":
            tri = TRNAI_RE.search(text)
            if tri:
                anticodon = tri.group(1).upper()

        if anticodon == "":
            hg = HGNC_STYLE_RE.search(text)
            if hg:
                anticodon = hg.group(1).upper()

    anticodon_rna = anticodon.replace("T", "U") if anticodon else ""
    return aa, anticodon, anticodon_rna, label


def extract_mt_symbol(*texts: str) -> str:
    for t in texts:
        m = MT_SYMBOL_RE.search(normalize(t).upper())
        if m:
            return m.group(1).upper()
    return ""


def read_tsv(path: Path) -> Tuple[List[Dict[str, str]], List[str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        cols = list(r.fieldnames)
        rows = [dict(row) for row in r]
    return rows, cols


def write_tsv(path: Path, rows: List[Dict[str, str]], columns: List[str], source_version: str, fetch_date: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        w.writeheader()
        for row in rows:
            row = dict(row)
            row["source_version"] = source_version
            row["fetch_date"] = fetch_date
            w.writerow({c: normalize(str(row.get(c, ""))) for c in columns})


def write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def column_non_empty_rate(rows: List[Dict[str, str]], col: str) -> float:
    if not rows:
        return 0.0
    return sum(1 for r in rows if normalize(r.get(col, "")) != "") / len(rows)


def build_candidate_map(id_mapping_path: Path, target_urs: set[str], taxon_id: str) -> Tuple[Dict[str, Dict], Dict[str, int]]:
    cand_by_urs: Dict[str, Dict[str, Dict]] = defaultdict(lambda: defaultdict(lambda: {"score": 0.0, "lines": 0, "dbs": set(), "methods": set(), "aa": Counter()}))
    stats = {
        "idmap_lines_scanned": 0,
        "taxon_pass": 0,
        "target_rows": 0,
    }

    with open_maybe_gz(id_mapping_path) as f:
        for line in f:
            stats["idmap_lines_scanned"] += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue

            urs_raw = normalize(parts[0])
            db = normalize(parts[1]).upper()
            external_id = normalize(parts[2])
            tx = normalize(parts[3])
            rna_type = normalize(parts[4])
            aux = normalize(parts[5]) if len(parts) > 5 else ""

            if tx != taxon_id or rna_type not in {"tRNA", "Mt_tRNA"}:
                continue
            stats["taxon_pass"] += 1

            urs = urs_raw if urs_raw.endswith(f"_{taxon_id}") else f"{urs_raw}_{taxon_id}"
            if urs not in target_urs:
                continue
            stats["target_rows"] += 1

            aa, anti, _, _ = parse_trna_annotation(external_id, aux)
            if anti:
                method = "current_parser"
                score = 3.0 if db == "GTRNADB" else 2.0
                rec = cand_by_urs[urs][anti]
                rec["score"] += score
                rec["lines"] += 1
                rec["dbs"].add(db)
                rec["methods"].add(method)
                if aa:
                    rec["aa"][aa] += 1

            mt_symbol = extract_mt_symbol(external_id, aux)
            if mt_symbol and mt_symbol in TRNA_MT_SYMBOL_TO_ANTICODON:
                anti = TRNA_MT_SYMBOL_TO_ANTICODON[mt_symbol]
                rec = cand_by_urs[urs][anti]
                rec["score"] += 0.7
                rec["lines"] += 1
                rec["dbs"].add(db)
                rec["methods"].add("mt_symbol_map")

    return cand_by_urs, stats


def choose_balanced(candidates: Dict[str, Dict]) -> Tuple[str, str, str, str, List[Dict]]:
    if not candidates:
        return "", "none", "", "no_candidate", []

    ranked = []
    for anti, rec in candidates.items():
        methods = sorted(rec["methods"])
        dbs = sorted(rec["dbs"])
        score = float(rec["score"])
        lines = int(rec["lines"])
        aa = ""
        if rec["aa"]:
            aa = rec["aa"].most_common(1)[0][0]
        ranked.append({
            "anti": anti,
            "score": score,
            "lines": lines,
            "db_count": len(dbs),
            "dbs": dbs,
            "methods": methods,
            "aa": aa,
        })

    ranked.sort(key=lambda x: (-x["score"], -x["lines"], -x["db_count"], x["anti"]))

    # Balanced low-confidence gate for mt-only candidates
    filtered = []
    for item in ranked:
        mt_only = set(item["methods"]) == {"mt_symbol_map"}
        if mt_only and item["lines"] < 2 and item["db_count"] < 2:
            continue
        filtered.append(item)

    if not filtered:
        return "", "none", "", "low_conf_filtered", ranked

    if len(filtered) == 1:
        top = filtered[0]
        source = "gtrnadb" if "current_parser" in top["methods"] and "GTRNADB" in top["dbs"] else "symbol_inference"
        confidence = "high" if top["score"] >= 2.0 else "medium_low"
        return top["anti"], source, confidence, "single_or_high_conf", ranked

    top = filtered[0]
    second = filtered[1]
    if (top["score"] - second["score"]) >= 1.0:
        source = "gtrnadb" if "GTRNADB" in top["dbs"] and "current_parser" in top["methods"] else "symbol_inference"
        confidence = "high" if top["score"] >= 3.0 else "medium"
        return top["anti"], source, confidence, "weighted_conflict_resolve", ranked

    return "", "none", "", "conflict_unresolved", ranked


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build balanced optimized tRNA anticodon table v2")
    p.add_argument("--trna-v1", type=Path, default=Path("data/output/rna_trna_features_v1.tsv"))
    p.add_argument("--id-mapping", type=Path, default=Path("data/raw/rna/rnacentral/id_mapping.tsv.gz"))
    p.add_argument("--taxon-id", default="9606")
    p.add_argument("--output", type=Path, default=Path("data/output/rna_trna_features_v2.tsv"))
    p.add_argument("--conflict-audit", type=Path, default=Path("pipelines/rna_type_features/reports/rna_trna_anticodon_conflicts_v2.tsv"))
    p.add_argument("--report", type=Path, default=Path("pipelines/rna_type_features/reports/rna_trna_features_v2.metrics.json"))
    p.add_argument("--source-version", default="RNAcentral:id_mapping.tsv;logic:balanced_v1")
    p.add_argument("--fetch-date", default=FETCH_DATE)
    p.add_argument("--check-inputs", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    idmap_path = resolve_with_alt(args.id_mapping)

    missing = []
    if not args.trna_v1.exists():
        missing.append(str(args.trna_v1))
    if idmap_path is None:
        missing.append(str(args.id_mapping))

    if missing:
        write_json(
            args.report,
            {
                "pipeline": "rna_trna_features_v2",
                "status": "blocked_missing_inputs",
                "missing_required": missing,
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            },
        )
        print(f"[BLOCKED] missing required inputs: {missing}")
        return 2

    if args.check_inputs:
        write_json(
            args.report,
            {
                "pipeline": "rna_trna_features_v2",
                "status": "inputs_ready",
                "inputs": {
                    "trna_v1": str(args.trna_v1),
                    "id_mapping": str(idmap_path),
                },
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            },
        )
        print(f"[OK] inputs ready -> {args.report}")
        return 0

    rows, columns = read_tsv(args.trna_v1)
    if "anticodon_confidence" not in columns:
        columns = columns + ["anticodon_confidence"]

    baseline_rate = column_non_empty_rate(rows, "anticodon")
    target_urs = {normalize(r.get("urs_id", "")) for r in rows if normalize(r.get("anticodon", "")) == ""}

    cand_by_urs, scan_stats = build_candidate_map(idmap_path, target_urs, args.taxon_id)

    status_counter = Counter()
    source_counter = Counter(r.get("anticodon_source", "none") for r in rows)
    improved = 0
    unchanged_non_empty = sum(1 for r in rows if normalize(r.get("anticodon", "")) != "")

    conflict_rows: List[Dict[str, str]] = []

    for row in rows:
        if normalize(row.get("anticodon", "")) != "":
            row["anticodon_confidence"] = row.get("anticodon_confidence", "high") or "high"
            continue

        urs = normalize(row.get("urs_id", ""))
        anti, anti_source, conf, status, ranked = choose_balanced(cand_by_urs.get(urs, {}))
        status_counter[status] += 1

        if anti:
            row["anticodon"] = anti
            row["anticodon_rna"] = anti.replace("T", "U")
            row["anticodon_source"] = anti_source
            row["anticodon_confidence"] = conf
            improved += 1
            source_counter[anti_source] += 1
        else:
            row["anticodon_confidence"] = ""
            if status == "conflict_unresolved":
                # keep top 5 for audit
                top = ranked[:5]
                conflict_rows.append(
                    {
                        "urs_id": urs,
                        "candidate_count": str(len(ranked)),
                        "top_candidates_json": json.dumps(top, ensure_ascii=False),
                    }
                )

    v2_rate = column_non_empty_rate(rows, "anticodon")

    # conflict audit
    args.conflict_audit.parent.mkdir(parents=True, exist_ok=True)
    with args.conflict_audit.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["urs_id", "candidate_count", "top_candidates_json"], delimiter="\t")
        w.writeheader()
        w.writerows(conflict_rows)

    write_tsv(args.output, rows, columns, args.source_version, args.fetch_date)

    report = {
        "pipeline": "rna_trna_features_v2",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "trna_v1": str(args.trna_v1),
            "id_mapping": str(idmap_path),
        },
        "outputs": {
            "trna_v2": str(args.output),
            "conflict_audit": str(args.conflict_audit),
        },
        "row_count": {
            "rows_total": len(rows),
            "rows_non_empty_baseline": unchanged_non_empty,
            "rows_newly_filled": improved,
            "rows_conflict_unresolved": len(conflict_rows),
        },
        "rates": {
            "anticodon_non_empty_rate_baseline": baseline_rate,
            "anticodon_non_empty_rate_v2": v2_rate,
            "delta": v2_rate - baseline_rate,
        },
        "distribution": {
            "status": dict(status_counter),
            "anticodon_source": dict(source_counter),
        },
        "scan_stats": scan_stats,
        "acceptance": {
            "non_empty_rate_not_degrade": v2_rate >= baseline_rate,
            "target_0_115_met": v2_rate >= 0.115,
            "conflict_audit_generated": True,
        },
    }
    write_json(args.report, report)

    print(
        f"[OK] rows={len(rows)} baseline={baseline_rate:.6f} v2={v2_rate:.6f} "
        f"newly_filled={improved} unresolved_conflicts={len(conflict_rows)}"
    )
    print(f"[OK] output -> {args.output}")
    print(f"[OK] report -> {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
