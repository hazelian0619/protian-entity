#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def open_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


@dataclass
class MasterRecord:
    rna_id: str
    rna_type: str
    taxon_id: str
    enst_base_from_rna_id: str
    ensembl_transcript_id_versioned: str
    ensembl_transcript_id_base: str
    enst_base_from_rna_name: str


def resolve_master_path(explicit: Optional[Path]) -> Path:
    if explicit is not None:
        if explicit.exists():
            return explicit
        raise SystemExit(f"[ERROR] --master not found: {explicit}")

    candidates = [
        Path("data/output/rna_master_v1.tsv"),
        Path("data/output/rna_master_mrna_v1.fixed.tsv"),
        Path("data/output/rna_master_mrna_v1.tsv"),
        Path("dist/rna-l1-v1/rna_master_mrna_v1.tsv.gz"),
    ]
    for c in candidates:
        if c.exists():
            return c

    raise SystemExit(
        "[ERROR] missing mRNA master table. Checked: "
        + ", ".join(str(c) for c in candidates)
    )


def read_master_records(
    master_path: Path,
    taxon_id: str,
    max_master_rows: Optional[int],
) -> Tuple[List[MasterRecord], Dict[str, int]]:
    stats = {
        "master_rows_scanned": 0,
        "master_mrna_rows": 0,
        "master_mrna_target_taxon_rows": 0,
        "master_mrna_enst_rows": 0,
    }

    records: List[MasterRecord] = []

    with open_maybe_gz(master_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header in: {master_path}")

        required = {"rna_id", "rna_type", "taxon_id"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise SystemExit(f"[ERROR] missing required columns in master: {sorted(missing)}")

        for row in reader:
            stats["master_rows_scanned"] += 1

            rna_type = (row.get("rna_type") or "").strip().lower()
            if rna_type != "mrna":
                continue
            stats["master_mrna_rows"] += 1

            tx = (row.get("taxon_id") or "").strip()
            if tx != taxon_id:
                continue
            stats["master_mrna_target_taxon_rows"] += 1

            rna_id = (row.get("rna_id") or "").strip()
            if rna_id.startswith("ENST") and rna_id.endswith(f"_{taxon_id}"):
                enst_base_from_rna_id = rna_id[: -(len(taxon_id) + 1)]
                stats["master_mrna_enst_rows"] += 1
            else:
                enst_base_from_rna_id = ""

            et_ver = (row.get("ensembl_transcript_id") or "").strip()
            et_base = et_ver.split(".")[0] if et_ver.startswith("ENST") else ""

            rna_name = (row.get("rna_name") or "").strip()
            rna_name_base = rna_name.split(".")[0] if rna_name.startswith("ENST") else ""

            records.append(
                MasterRecord(
                    rna_id=rna_id,
                    rna_type="mrna",
                    taxon_id=tx,
                    enst_base_from_rna_id=enst_base_from_rna_id,
                    ensembl_transcript_id_versioned=et_ver if et_ver.startswith("ENST") else "",
                    ensembl_transcript_id_base=et_base,
                    enst_base_from_rna_name=rna_name_base,
                )
            )

            if max_master_rows is not None and len(records) >= max_master_rows:
                break

    return records, stats


def read_id_mapping(
    id_mapping_path: Path,
    taxon_id: str,
    allowed_dbs: Tuple[str, ...],
    max_idmap_lines: Optional[int],
) -> Tuple[Dict[str, Counter], Dict[str, Counter], Dict[str, int]]:
    by_base: Dict[str, Counter] = defaultdict(Counter)
    by_versioned: Dict[str, Counter] = defaultdict(Counter)

    stats = {
        "idmap_lines_scanned": 0,
        "idmap_pass_taxid": 0,
        "idmap_pass_db": 0,
        "idmap_pass_enst": 0,
    }

    with open_maybe_gz(id_mapping_path) as f:
        for line in f:
            stats["idmap_lines_scanned"] += 1
            if max_idmap_lines is not None and stats["idmap_lines_scanned"] > max_idmap_lines:
                break

            if line.strip() == "":
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue

            urs_raw = parts[0].strip()
            db = parts[1].strip().upper()
            ext_id = parts[2].strip()
            tx = parts[3].strip()

            if tx != taxon_id:
                continue
            stats["idmap_pass_taxid"] += 1

            if db not in allowed_dbs:
                continue
            stats["idmap_pass_db"] += 1

            if not ext_id.startswith("ENST"):
                continue
            stats["idmap_pass_enst"] += 1

            urs_norm = urs_raw if urs_raw.endswith(f"_{taxon_id}") else f"{urs_raw}_{taxon_id}"
            ext_base = ext_id.split(".")[0]

            by_base[ext_base][urs_norm] += 1
            by_versioned[ext_id][urs_norm] += 1

    return by_base, by_versioned, stats


def choose_best_urs(counter: Counter) -> Tuple[str, List[Tuple[str, int]]]:
    ranked = sorted(counter.items(), key=lambda kv: (-kv[1], kv[0]))
    best_urs = ranked[0][0]
    return best_urs, ranked


def manual_download_plan(min_coverage: float, observed_coverage: float) -> Dict:
    return {
        "reason": "RNAcentral id_mapping.tsv alone cannot reach ENST↔URS coverage gate for current mRNA universe.",
        "current_coverage": observed_coverage,
        "min_coverage_required": min_coverage,
        "download_checklist": [
            {
                "name": "Ensembl transcript xref export (GRCh38 / Ensembl release aligned to mRNA master)",
                "version": "match mRNA master release (e.g., Ensembl 113)",
                "suggested_filename": "ensembl_transcript_xref_grch38_r113.tsv.gz",
                "required_columns": [
                    "ensembl_transcript_id",
                    "refseq_mrna",
                    "insdc_mrna",
                    "ccds_id",
                ],
                "purpose": "Bridge ENST to RefSeq/INSDC accessions so RNAcentral REFSEQ/ENA mappings can be joined.",
            },
            {
                "name": "NCBI gene2refseq",
                "version": "latest compatible snapshot",
                "suggested_filename": "gene2refseq.gz",
                "required_columns": [
                    "tax_id",
                    "RNA_nucleotide_accession.version",
                    "GeneID",
                ],
                "purpose": "Fallback ENST→GeneID→RefSeq bridge for transcripts missing direct xref.",
            },
            {
                "name": "RNAcentral release notes / checksum file",
                "version": "same as local id_mapping.tsv",
                "suggested_filename": "CHECKSUMS",
                "required_columns": [],
                "purpose": "Verify integrity/version alignment of id_mapping and auxiliary files.",
            },
        ],
        "place_under": "data/raw/rna/aux_xref/",
        "checksum": {
            "method": "sha256",
            "commands": [
                "mkdir -p data/raw/rna/aux_xref",
                "cd data/raw/rna/aux_xref",
                "shasum -a 256 ensembl_transcript_xref_grch38_r113.tsv.gz gene2refseq.gz > SHA256SUMS.txt",
                "shasum -a 256 -c SHA256SUMS.txt",
            ],
        },
    }


def write_tsv(path: Path, header: List[str], rows: Iterable[List[str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    n = 0
    with tmp.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for row in rows:
            w.writerow(row)
            n += 1
    tmp.replace(path)
    return n


def write_json(path: Path, obj: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", type=Path, default=None)
    ap.add_argument("--id-mapping", type=Path, default=Path("data/raw/rna/rnacentral/id_mapping.tsv"))
    ap.add_argument("--taxon-id", default="9606")
    ap.add_argument("--out-xref", type=Path, required=True)
    ap.add_argument("--coverage-report", type=Path, required=True)
    ap.add_argument("--conflicts-report", type=Path, required=True)
    ap.add_argument("--manual-download-report", type=Path, required=True)
    ap.add_argument("--fetch-date", default=date.today().isoformat())
    ap.add_argument("--source-version", default="RNAcentral:25")
    ap.add_argument("--source", default="RNAcentral:id_mapping.tsv")
    ap.add_argument("--min-coverage", type=float, default=0.70)
    ap.add_argument("--target-coverage", type=float, default=0.90)
    ap.add_argument("--max-master-rows", type=int, default=None)
    ap.add_argument("--max-idmap-lines", type=int, default=None)
    args = ap.parse_args()

    master_path = resolve_master_path(args.master)
    if not args.id_mapping.exists():
        raise SystemExit(f"[ERROR] missing id_mapping: {args.id_mapping}")

    allowed_dbs = ("ENSEMBL", "ENSEMBL_GENCODE", "GENCODE")

    records, master_stats = read_master_records(
        master_path=master_path,
        taxon_id=args.taxon_id,
        max_master_rows=args.max_master_rows,
    )

    by_base, by_versioned, idmap_stats = read_id_mapping(
        id_mapping_path=args.id_mapping,
        taxon_id=args.taxon_id,
        allowed_dbs=allowed_dbs,
        max_idmap_lines=args.max_idmap_lines,
    )

    strategy_counts: Counter = Counter()
    unresolved_examples: List[str] = []
    one_to_many_conflicts: List[Dict] = []

    xref_rows: List[List[str]] = []

    for rec in records:
        if rec.enst_base_from_rna_id == "":
            if len(unresolved_examples) < 20:
                unresolved_examples.append(rec.rna_id)
            continue

        strategy_chain = [
            ("exact_rna_id_base", rec.enst_base_from_rna_id, by_base),
            ("fallback_ensembl_transcript_version", rec.ensembl_transcript_id_versioned, by_versioned),
            ("fallback_ensembl_transcript_base", rec.ensembl_transcript_id_base, by_base),
            ("fallback_rna_name_base", rec.enst_base_from_rna_name, by_base),
        ]

        matched = False
        for strategy, key, mapping in strategy_chain:
            if key == "":
                continue
            urs_counter = mapping.get(key)
            if not urs_counter:
                continue

            chosen_urs, ranked = choose_best_urs(urs_counter)

            if len(ranked) > 1:
                one_to_many_conflicts.append(
                    {
                        "rna_id": rec.rna_id,
                        "strategy": strategy,
                        "candidate_count": len(ranked),
                        "candidates": [{"urs": urs, "votes": votes} for urs, votes in ranked[:20]],
                    }
                )

            strategy_counts[strategy] += 1
            xref_rows.append(
                [
                    rec.rna_id,
                    rec.rna_type,
                    rec.taxon_id,
                    "RNAcentral",
                    chosen_urs,
                    "urs_id",
                    strategy,
                    args.source,
                    args.fetch_date,
                    args.source_version,
                ]
            )
            matched = True
            break

        if not matched and len(unresolved_examples) < 20:
            unresolved_examples.append(rec.rna_id)

    xref_rows.sort(key=lambda r: r[0])

    header = [
        "rna_id",
        "rna_type",
        "taxon_id",
        "xref_db",
        "xref_id",
        "xref_level",
        "match_strategy",
        "source",
        "fetch_date",
        "source_version",
    ]
    row_count = write_tsv(args.out_xref, header=header, rows=xref_rows)

    reverse = defaultdict(list)
    for row in xref_rows:
        reverse[row[4]].append(row[0])

    many_to_one = [
        {"xref_id": urs, "rna_ids": sorted(ids), "count": len(ids)}
        for urs, ids in reverse.items()
        if len(ids) > 1
    ]
    many_to_one.sort(key=lambda x: (-x["count"], x["xref_id"]))

    conflicts_report = {
        "name": "mrna_enst_urs_conflicts_v2",
        "generated_at": utc_now_iso(),
        "one_to_many_conflicts": {
            "count": len(one_to_many_conflicts),
            "items": one_to_many_conflicts,
        },
        "many_to_one_conflicts": {
            "count": len(many_to_one),
            "items": many_to_one,
        },
    }
    write_json(args.conflicts_report, conflicts_report)

    total_enst = master_stats["master_mrna_enst_rows"]
    mapped = row_count
    coverage = (mapped / total_enst) if total_enst else 0.0

    met_minimum = coverage >= args.min_coverage
    met_target = coverage >= args.target_coverage

    manual_plan = None
    if not met_minimum:
        manual_plan = manual_download_plan(min_coverage=args.min_coverage, observed_coverage=coverage)
        write_json(args.manual_download_report, manual_plan)
    elif args.manual_download_report.exists():
        args.manual_download_report.unlink()

    coverage_report = {
        "name": "mrna_enst_urs_coverage_v2",
        "generated_at": utc_now_iso(),
        "inputs": {
            "master": str(master_path),
            "id_mapping": str(args.id_mapping),
            "taxon_id": args.taxon_id,
            "allowed_dbs": list(allowed_dbs),
            "sample_mode": bool(args.max_master_rows is not None or args.max_idmap_lines is not None),
            "max_master_rows": args.max_master_rows,
            "max_idmap_lines": args.max_idmap_lines,
        },
        "counts": {
            **master_stats,
            **idmap_stats,
            "xref_rows": row_count,
            "mapped_mrna_enst_rows": mapped,
            "unmapped_mrna_enst_rows": max(total_enst - mapped, 0),
        },
        "coverage": {
            "mapped": mapped,
            "total": total_enst,
            "rate": coverage,
            "percent": round(coverage * 100.0, 4),
            "strategy_counts": dict(strategy_counts),
            "unresolved_examples": unresolved_examples,
        },
        "thresholds": {
            "minimum": args.min_coverage,
            "target": args.target_coverage,
            "met_minimum": met_minimum,
            "met_target": met_target,
        },
        "gates": {
            "status": "PASS" if met_minimum else "LOW_COVERAGE_KPI_NOT_MET",
            "reason": (
                "coverage meets minimum gate"
                if met_minimum
                else "coverage below KPI minimum; auxiliary xref files recommended"
            ),
            "recommended_action": (
                "none"
                if met_minimum
                else "download auxiliary xref files and improve mapping in follow-up run"
            ),
        },
        "conflicts_report": str(args.conflicts_report),
        "manual_download_report": str(args.manual_download_report) if manual_plan else None,
    }
    write_json(args.coverage_report, coverage_report)

    print(f"[OK] wrote xref: {args.out_xref} rows={row_count}")
    print(f"[OK] wrote coverage report: {args.coverage_report}")
    print(f"[OK] wrote conflicts report: {args.conflicts_report}")
    print(f"[INFO] coverage={coverage * 100:.4f}% (min={args.min_coverage * 100:.2f}%)")
    if manual_plan is not None:
        print(f"[WARN] coverage below gate; manual download checklist: {args.manual_download_report}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
