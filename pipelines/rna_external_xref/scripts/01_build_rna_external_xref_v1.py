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
from typing import Dict, Iterable, List, Optional, Set, Tuple


TARGET_TAXON = "9606"
DIRECT_ENST_DBS = {"ENSEMBL", "ENSEMBL_GENCODE", "GENCODE"}
TARGET_XREF_DBS = {"ENSEMBL", "ENSEMBL_GENCODE", "GENCODE", "REFSEQ", "MIRBASE", "RNACENTRAL"}


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def open_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def write_json(path: Path, obj: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


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


def normalize_urs(urs_raw: str, taxon_id: str) -> str:
    urs = (urs_raw or "").strip()
    if urs == "":
        return ""
    if urs.endswith(f"_{taxon_id}"):
        return urs
    return f"{urs}_{taxon_id}"


def normalize_enst(value: str) -> Tuple[str, str]:
    s = (value or "").strip()
    if not s.startswith("ENST"):
        return "", ""
    base = s.split(".")[0]
    return s, base


def normalize_refseq(value: str) -> Tuple[str, str]:
    s = (value or "").strip()
    if s == "":
        return "", ""
    token = s.split(":", 1)[0].strip()
    base = token.split(".")[0]
    return token, base


def classify_xref_type(xref_db: str, xref_id: str) -> str:
    db = (xref_db or "").upper()
    xid = (xref_id or "").strip()
    if db == "RNACENTRAL" or xid.startswith("URS"):
        return "urs_id"
    if db in {"ENSEMBL", "ENSEMBL_GENCODE", "GENCODE"} or xid.startswith("ENST"):
        return "transcript_id"
    if db == "REFSEQ" or xid.startswith(("NM_", "NR_", "XM_", "XR_", "NP_", "YP_")):
        return "refseq_accession"
    if db == "MIRBASE" or xid.startswith(("MIMAT", "MI")):
        return "mirbase_id"
    return "accession"


@dataclass
class MasterRecord:
    rna_id: str
    rna_type: str
    rna_name: str
    ensembl_transcript_id: str
    mirbase_id: str


class AuxParseError(RuntimeError):
    pass


def _header_index(header: List[str]) -> Dict[str, int]:
    return {h.strip().lower(): i for i, h in enumerate(header)}


def _pick_col(idx: Dict[str, int], candidates: List[str]) -> Optional[int]:
    for c in candidates:
        if c in idx:
            return idx[c]
    return None


def load_aux_enst_to_refseq(aux_dir: Path, max_rows_per_file: Optional[int] = None) -> Tuple[Dict[str, Set[str]], Dict]:
    if not aux_dir.exists() or not aux_dir.is_dir():
        raise AuxParseError(f"aux_xref directory not found: {aux_dir}")

    files = sorted(
        [
            p
            for p in aux_dir.iterdir()
            if p.is_file() and any(str(p).endswith(ext) for ext in (".tsv", ".tsv.gz", ".txt", ".txt.gz"))
        ]
    )
    if not files:
        raise AuxParseError(f"no aux_xref mapping files found in: {aux_dir}")

    enst_to_refseq: Dict[str, Set[str]] = defaultdict(set)
    file_reports: List[Dict] = []

    enst_cols = [
        "ensembl_transcript_id",
        "enst",
        "enst_id",
        "transcript_id",
        "transcript stable id",
        "transcript stable id version",
    ]
    refseq_cols = [
        "refseq_mrna",
        "refseq_accession",
        "refseq mrna id",
        "refseq mrna predicted id",
        "refseq ncrna id",
        "refseq ncrna predicted id",
        "rna_nucleotide_accession.version",
        "rna_nucleotide_accession",
        "insdc_mrna",
        "ena_mrna",
        "mrna_accession",
    ]

    for fp in files:
        parsed_rows = 0
        mapped_rows = 0
        skipped_rows = 0
        used = False

        with open_maybe_gz(fp) as f:
            first = f.readline()
            if not first:
                file_reports.append(
                    {
                        "file": str(fp),
                        "status": "empty",
                        "parsed_rows": 0,
                        "mapped_rows": 0,
                    }
                )
                continue

            header = first.rstrip("\n\r").split("\t")
            idx = _header_index(header)
            enst_i = _pick_col(idx, enst_cols)
            ref_i = _pick_col(idx, refseq_cols)

            if enst_i is None or ref_i is None:
                file_reports.append(
                    {
                        "file": str(fp),
                        "status": "unsupported_header",
                        "parsed_rows": 0,
                        "mapped_rows": 0,
                        "required_columns": {
                            "ensembl": enst_cols,
                            "refseq": refseq_cols,
                        },
                        "header": header,
                    }
                )
                continue

            used = True
            for line in f:
                if line.strip() == "":
                    continue
                parsed_rows += 1
                parts = line.rstrip("\n\r").split("\t")
                if len(parts) <= max(enst_i, ref_i):
                    skipped_rows += 1
                    continue

                enst_ver, enst_base = normalize_enst(parts[enst_i])
                ref_ver, ref_base = normalize_refseq(parts[ref_i])
                if enst_base == "" or (ref_ver == "" and ref_base == ""):
                    skipped_rows += 1
                    continue

                # Keep both versioned + base token to maximize join hit rate.
                if ref_ver:
                    enst_to_refseq[enst_base].add(ref_ver)
                if ref_base:
                    enst_to_refseq[enst_base].add(ref_base)
                if enst_ver:
                    enst_to_refseq[enst_ver].add(ref_ver or ref_base)

                mapped_rows += 1
                if max_rows_per_file is not None and parsed_rows >= max_rows_per_file:
                    break

        file_reports.append(
            {
                "file": str(fp),
                "status": "used" if used else "unused",
                "parsed_rows": parsed_rows,
                "mapped_rows": mapped_rows,
                "skipped_rows": skipped_rows,
            }
        )

    used_files = [r for r in file_reports if r.get("status") == "used"]
    if not used_files:
        raise AuxParseError(
            "no usable aux_xref file: expected at least one file containing ENST + RefSeq-like columns"
        )

    summary = {
        "aux_dir": str(aux_dir),
        "files_total": len(files),
        "files_used": len(used_files),
        "enst_keys": len(enst_to_refseq),
        "enst_refseq_pairs": int(sum(len(v) for v in enst_to_refseq.values())),
        "files": file_reports,
    }
    return enst_to_refseq, summary


def read_master(master_path: Path, max_rows: Optional[int] = None) -> Tuple[List[MasterRecord], Dict]:
    if not master_path.exists():
        raise SystemExit(f"[ERROR] missing master: {master_path}")

    out: List[MasterRecord] = []
    stats = Counter()
    required = {"rna_id", "rna_type", "rna_name", "ensembl_transcript_id", "mirbase_id"}

    with master_path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        missing = required - set(r.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] master missing required columns: {sorted(missing)}")

        for row in r:
            stats["master_rows_scanned"] += 1
            rec = MasterRecord(
                rna_id=(row.get("rna_id") or "").strip(),
                rna_type=(row.get("rna_type") or "").strip().lower(),
                rna_name=(row.get("rna_name") or "").strip(),
                ensembl_transcript_id=(row.get("ensembl_transcript_id") or "").strip(),
                mirbase_id=(row.get("mirbase_id") or "").strip(),
            )
            if rec.rna_id == "":
                stats["master_empty_rna_id"] += 1
                continue
            out.append(rec)
            stats[f"master_type_{rec.rna_type or 'unknown'}"] += 1
            if max_rows is not None and len(out) >= max_rows:
                break

    stats["master_rows_used"] = len(out)
    return out, dict(stats)


def choose_best(counter: Counter) -> Tuple[str, List[Tuple[str, int]]]:
    ranked = sorted(counter.items(), key=lambda kv: (-kv[1], kv[0]))
    return ranked[0][0], ranked


def build_mapping_index(master_records: List[MasterRecord]) -> Tuple[Dict[str, Set[str]], Set[str], Set[str], Dict[str, str]]:
    enst_key_to_rna_ids: Dict[str, Set[str]] = defaultdict(set)
    enst_keys_needed: Set[str] = set()
    urs_seed_to_rna_id: Dict[str, str] = {}
    master_rna_to_type: Dict[str, str] = {}

    for rec in master_records:
        master_rna_to_type[rec.rna_id] = rec.rna_type

        if rec.rna_id.startswith("URS") and rec.rna_id.endswith(f"_{TARGET_TAXON}"):
            urs_seed_to_rna_id[rec.rna_id] = rec.rna_id

        if rec.rna_id.startswith("ENST") and rec.rna_id.endswith(f"_{TARGET_TAXON}"):
            base = rec.rna_id[: -(len(TARGET_TAXON) + 1)]
            enst_key_to_rna_ids[base].add(rec.rna_id)
            enst_keys_needed.add(base)

        enst_ver, enst_base = normalize_enst(rec.ensembl_transcript_id)
        if enst_ver:
            enst_key_to_rna_ids[enst_ver].add(rec.rna_id)
            enst_keys_needed.add(enst_ver)
        if enst_base:
            enst_key_to_rna_ids[enst_base].add(rec.rna_id)
            enst_keys_needed.add(enst_base)

        name_ver, name_base = normalize_enst(rec.rna_name)
        if name_ver:
            enst_key_to_rna_ids[name_ver].add(rec.rna_id)
            enst_keys_needed.add(name_ver)
        if name_base:
            enst_key_to_rna_ids[name_base].add(rec.rna_id)
            enst_keys_needed.add(name_base)

    return enst_key_to_rna_ids, enst_keys_needed, set(urs_seed_to_rna_id.keys()), master_rna_to_type


def build_rna_to_urs(
    master_records: List[MasterRecord],
    enst_to_urs: Dict[str, Counter],
    refseq_to_urs: Dict[str, Counter],
    aux_enst_to_refseq: Dict[str, Set[str]],
) -> Tuple[Dict[str, str], Dict]:
    rna_to_urs: Dict[str, str] = {}
    strategy_counts = Counter()
    unresolved_examples: List[str] = []
    direct_conflicts: List[Dict] = []

    # Record ENST keys with multi-URS ambiguity (report only)
    for key, counter in enst_to_urs.items():
        if len(counter) > 1 and key.startswith("ENST"):
            _, ranked = choose_best(counter)
            direct_conflicts.append(
                {
                    "enst_key": key,
                    "candidate_count": len(ranked),
                    "candidates": [{"urs": u, "votes": c} for u, c in ranked[:10]],
                }
            )

    for rec in master_records:
        rid = rec.rna_id
        if rid.startswith("URS") and rid.endswith(f"_{TARGET_TAXON}"):
            rna_to_urs[rid] = rid
            strategy_counts["legacy_urs"] += 1
            continue

        if not rid.startswith("ENST"):
            strategy_counts["unsupported_rna_id"] += 1
            if len(unresolved_examples) < 30:
                unresolved_examples.append(rid)
            continue

        candidates: List[Tuple[str, str]] = []

        enst_from_rid = rid[: -(len(TARGET_TAXON) + 1)] if rid.endswith(f"_{TARGET_TAXON}") else ""
        enst_ver, enst_base = normalize_enst(rec.ensembl_transcript_id)
        name_ver, name_base = normalize_enst(rec.rna_name)

        if enst_from_rid:
            candidates.append(("rna_id_base", enst_from_rid))
        if enst_ver:
            candidates.append(("ensembl_transcript_version", enst_ver))
        if enst_base:
            candidates.append(("ensembl_transcript_base", enst_base))
        if name_ver:
            candidates.append(("rna_name_version", name_ver))
        if name_base:
            candidates.append(("rna_name_base", name_base))

        chosen = ""
        used_strategy = ""
        for strategy, key in candidates:
            counter = enst_to_urs.get(key)
            if not counter:
                continue
            chosen, _ = choose_best(counter)
            used_strategy = f"direct_{strategy}"
            break

        if chosen == "":
            # aux fallback: ENST -> RefSeq -> URS
            aux_keys = []
            if enst_from_rid:
                aux_keys.append(enst_from_rid)
            if enst_ver:
                aux_keys.append(enst_ver)
            if enst_base:
                aux_keys.append(enst_base)
            if name_ver:
                aux_keys.append(name_ver)
            if name_base:
                aux_keys.append(name_base)

            bridge_counter: Counter = Counter()
            for k in aux_keys:
                for ref in aux_enst_to_refseq.get(k, set()):
                    ref_ver, ref_base = normalize_refseq(ref)
                    if ref_ver:
                        bridge_counter.update(refseq_to_urs.get(ref_ver, Counter()))
                    if ref_base and ref_base != ref_ver:
                        bridge_counter.update(refseq_to_urs.get(ref_base, Counter()))

            if bridge_counter:
                chosen, _ = choose_best(bridge_counter)
                used_strategy = "aux_refseq_bridge"

        if chosen:
            rna_to_urs[rid] = chosen
            strategy_counts[used_strategy or "unknown"] += 1
        else:
            strategy_counts["unmapped"] += 1
            if len(unresolved_examples) < 30:
                unresolved_examples.append(rid)

    report = {
        "rna_to_urs_assigned": len(rna_to_urs),
        "strategy_counts": dict(strategy_counts),
        "direct_enst_ambiguous_keys": {
            "count": len(direct_conflicts),
            "items": direct_conflicts[:200],
        },
        "unresolved_examples": unresolved_examples,
    }
    return rna_to_urs, report


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", type=Path, required=True)
    ap.add_argument("--id-mapping", type=Path, required=True)
    ap.add_argument("--aux-dir", type=Path, required=True)
    ap.add_argument("--out-table", type=Path, required=True)
    ap.add_argument("--build-report", type=Path, required=True)
    ap.add_argument("--coverage-report", type=Path, required=True)
    ap.add_argument("--conflicts-report", type=Path, required=True)
    ap.add_argument("--fetch-date", default=date.today().isoformat())
    ap.add_argument("--source", default="RNAcentral:id_mapping.tsv")
    ap.add_argument("--source-version", default="RNAcentral:25")
    ap.add_argument("--min-backlink-rate", type=float, default=0.99)
    ap.add_argument("--max-master-rows", type=int, default=None)
    ap.add_argument("--max-idmap-lines", type=int, default=None)
    ap.add_argument("--max-aux-rows-per-file", type=int, default=None)
    args = ap.parse_args()

    if not args.id_mapping.exists():
        raise SystemExit(f"[ERROR] missing id_mapping: {args.id_mapping}")

    master_records, master_stats = read_master(args.master, max_rows=args.max_master_rows)
    if not master_records:
        raise SystemExit("[ERROR] empty master input")

    try:
        aux_enst_to_refseq, aux_summary = load_aux_enst_to_refseq(
            args.aux_dir,
            max_rows_per_file=args.max_aux_rows_per_file,
        )
    except AuxParseError as e:
        raise SystemExit(f"[ERROR] {e}")

    enst_key_to_rna_ids, enst_keys_needed, seed_urs_set, master_rna_to_type = build_mapping_index(master_records)

    enst_to_urs: Dict[str, Counter] = defaultdict(Counter)
    refseq_to_urs: Dict[str, Counter] = defaultdict(Counter)
    idmap_stats = Counter()

    # Pass 1: build ENST->URS and REFSEQ->URS mapping indexes.
    with open_maybe_gz(args.id_mapping) as f:
        for line in f:
            idmap_stats["idmap_lines_scanned"] += 1
            if args.max_idmap_lines is not None and idmap_stats["idmap_lines_scanned"] > args.max_idmap_lines:
                break
            if line.strip() == "":
                continue

            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 4:
                idmap_stats["idmap_lines_bad"] += 1
                continue

            urs_raw, db_raw, ext_raw, taxon = parts[0].strip(), parts[1].strip().upper(), parts[2].strip(), parts[3].strip()
            if taxon != TARGET_TAXON:
                continue
            idmap_stats["idmap_lines_taxon_pass"] += 1

            urs = normalize_urs(urs_raw, TARGET_TAXON)
            if urs == "" or ext_raw == "":
                continue

            if db_raw in DIRECT_ENST_DBS and ext_raw.startswith("ENST"):
                idmap_stats["idmap_direct_enst_rows"] += 1
                ext_ver, ext_base = normalize_enst(ext_raw)
                if ext_ver and ext_ver in enst_keys_needed:
                    enst_to_urs[ext_ver][urs] += 1
                if ext_base and ext_base in enst_keys_needed:
                    enst_to_urs[ext_base][urs] += 1

            if db_raw == "REFSEQ":
                idmap_stats["idmap_refseq_rows"] += 1
                ref_ver, ref_base = normalize_refseq(ext_raw)
                if ref_ver:
                    refseq_to_urs[ref_ver][urs] += 1
                if ref_base:
                    refseq_to_urs[ref_base][urs] += 1

    rna_to_urs, mapping_report = build_rna_to_urs(
        master_records=master_records,
        enst_to_urs=enst_to_urs,
        refseq_to_urs=refseq_to_urs,
        aux_enst_to_refseq=aux_enst_to_refseq,
    )

    urs_to_rna_ids: Dict[str, Set[str]] = defaultdict(set)
    for rna_id, urs in rna_to_urs.items():
        urs_to_rna_ids[urs].add(rna_id)

    target_urs = set(urs_to_rna_ids.keys()) | seed_urs_set

    # Build xref rows
    row_map: Dict[Tuple[str, str, str, str], Tuple[str, str, str]] = {}
    output_db_counter = Counter()
    output_source_counter = Counter()

    def emit_row(rna_id: str, xref_db: str, xref_id: str, source: str) -> None:
        xdb = (xref_db or "").strip().upper()
        xid = (xref_id or "").strip()
        if rna_id == "" or xdb == "" or xid == "":
            return
        xref_type = classify_xref_type(xdb, xid)
        row_key = (rna_id, xdb, xid, xref_type)
        prev = row_map.get(row_key)
        if prev is None:
            row_map[row_key] = (source, args.source_version, args.fetch_date)
            return

        # Prefer RNAcentral id_mapping source over fallback master seed.
        prev_source, _prev_version, _prev_date = prev
        if prev_source == "rna_master_v1.tsv" and source != "rna_master_v1.tsv":
            row_map[row_key] = (source, args.source_version, args.fetch_date)

    # Seed rows from master to guarantee backlink/master coverage visibility.
    for rec in master_records:
        rid = rec.rna_id
        if rid.startswith("ENST") and rid.endswith(f"_{TARGET_TAXON}"):
            enst_from_rid = rid[: -(len(TARGET_TAXON) + 1)]
            emit_row(rid, "ENSEMBL", enst_from_rid, "rna_master_v1.tsv")
        if rec.mirbase_id:
            emit_row(rid, "MIRBASE", rec.mirbase_id, "rna_master_v1.tsv")
        if rid.startswith("URS") and rid.endswith(f"_{TARGET_TAXON}"):
            emit_row(rid, "RNACENTRAL", rid, "rna_master_v1.tsv")

    # Pass 2: fetch target db xrefs for mapped URS universe.
    with open_maybe_gz(args.id_mapping) as f:
        for line in f:
            idmap_stats["idmap_lines_scanned_pass2"] += 1
            if args.max_idmap_lines is not None and idmap_stats["idmap_lines_scanned_pass2"] > args.max_idmap_lines:
                break
            if line.strip() == "":
                continue

            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 4:
                continue

            urs_raw, db_raw, ext_raw, taxon = parts[0].strip(), parts[1].strip().upper(), parts[2].strip(), parts[3].strip()
            if taxon != TARGET_TAXON:
                continue
            if db_raw not in TARGET_XREF_DBS:
                continue
            if ext_raw.strip() == "":
                continue

            urs = normalize_urs(urs_raw, TARGET_TAXON)
            if urs not in target_urs:
                continue

            # Ensure RNACENTRAL xref_id is normalized URS.
            xref_id = urs if db_raw == "RNACENTRAL" else ext_raw.strip()

            for rid in urs_to_rna_ids.get(urs, set()):
                emit_row(rid, db_raw, xref_id, args.source)

            # For native URS master rows (miRNA), always include row even if not in rna_to_urs map
            if urs in seed_urs_set:
                emit_row(urs, db_raw, xref_id, args.source)

    # Final table
    rows_sorted = sorted(
        [
            (rid, db, xid, xref_type, src, sv, fd)
            for (rid, db, xid, xref_type), (src, sv, fd) in row_map.items()
        ],
        key=lambda x: (x[0], x[1], x[2], x[3]),
    )
    output_db_counter = Counter(row[1] for row in rows_sorted)
    output_source_counter = Counter(row[4] for row in rows_sorted)
    header = ["rna_id", "xref_db", "xref_id", "xref_type", "source", "source_version", "fetch_date"]
    out_count = write_tsv(args.out_table, header=header, rows=[list(r) for r in rows_sorted])

    # QA metrics shared by coverage/conflicts
    master_ids = {r.rna_id for r in master_records}
    out_ids = {r[0] for r in rows_sorted}
    backlink_hits = sum(1 for rid in out_ids if rid in master_ids)
    backlink_rate = (backlink_hits / len(out_ids)) if out_ids else 0.0
    master_covered = len(master_ids & out_ids)
    master_coverage_rate = (master_covered / len(master_ids)) if master_ids else 0.0

    per_db_rna: Dict[str, Set[str]] = defaultdict(set)
    per_db_xid: Dict[str, Set[str]] = defaultdict(set)
    per_db_rna_to_xids: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))
    per_db_xid_to_rnas: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))

    for rid, db, xid, *_rest in rows_sorted:
        per_db_rna[db].add(rid)
        per_db_xid[db].add(xid)
        per_db_rna_to_xids[db][rid].add(xid)
        per_db_xid_to_rnas[db][xid].add(rid)

    db_coverage = {}
    db_conflicts = {}

    for db in sorted(per_db_rna.keys()):
        covered = len(per_db_rna[db])
        unique_xids = len(per_db_xid[db])
        rna_with_multi = sum(1 for _rid, s in per_db_rna_to_xids[db].items() if len(s) > 1)
        xid_with_multi = sum(1 for _xid, s in per_db_xid_to_rnas[db].items() if len(s) > 1)

        db_coverage[db] = {
            "rows": output_db_counter.get(db, 0),
            "covered_rna_ids": covered,
            "master_total_rna_ids": len(master_ids),
            "coverage_rate": (covered / len(master_ids)) if master_ids else 0.0,
            "unique_xref_ids": unique_xids,
        }

        db_conflicts[db] = {
            "covered_rna_ids": covered,
            "unique_xref_ids": unique_xids,
            "rna_to_many_xref_count": rna_with_multi,
            "rna_to_many_conflict_rate": (rna_with_multi / covered) if covered else 0.0,
            "xref_to_many_rna_count": xid_with_multi,
            "xref_to_many_conflict_rate": (xid_with_multi / unique_xids) if unique_xids else 0.0,
            "samples": {
                "rna_to_many": [
                    {"rna_id": rid, "xref_ids": sorted(list(xids))[:10], "count": len(xids)}
                    for rid, xids in sorted(
                        per_db_rna_to_xids[db].items(), key=lambda kv: (-len(kv[1]), kv[0])
                    )
                    if len(xids) > 1
                ][:50],
                "xref_to_many": [
                    {"xref_id": xid, "rna_ids": sorted(list(rids))[:10], "count": len(rids)}
                    for xid, rids in sorted(
                        per_db_xid_to_rnas[db].items(), key=lambda kv: (-len(kv[1]), kv[0])
                    )
                    if len(rids) > 1
                ][:50],
            },
        }

    coverage_report = {
        "name": "rna_external_xref_v1.coverage",
        "generated_at": utc_now_iso(),
        "inputs": {
            "master": str(args.master),
            "id_mapping": str(args.id_mapping),
            "aux_dir": str(args.aux_dir),
            "sample_mode": bool(args.max_master_rows is not None or args.max_idmap_lines is not None),
            "max_master_rows": args.max_master_rows,
            "max_idmap_lines": args.max_idmap_lines,
        },
        "summary": {
            "rows_output": out_count,
            "unique_output_rna_ids": len(out_ids),
            "master_rna_ids": len(master_ids),
            "master_rna_ids_covered": master_covered,
            "master_coverage_rate": master_coverage_rate,
            "backlink_hits": backlink_hits,
            "backlink_rate": backlink_rate,
            "min_backlink_rate_required": args.min_backlink_rate,
            "backlink_passed": backlink_rate >= args.min_backlink_rate,
        },
        "per_xref_db": db_coverage,
    }

    conflicts_report = {
        "name": "rna_external_xref_v1.conflicts",
        "generated_at": utc_now_iso(),
        "per_xref_db": db_conflicts,
    }

    build_report = {
        "name": "rna_external_xref_v1.build",
        "generated_at": utc_now_iso(),
        "inputs": {
            "master": str(args.master),
            "id_mapping": str(args.id_mapping),
            "aux_dir": str(args.aux_dir),
            "source": args.source,
            "source_version": args.source_version,
            "fetch_date": args.fetch_date,
            "sample_mode": bool(args.max_master_rows is not None or args.max_idmap_lines is not None),
            "max_master_rows": args.max_master_rows,
            "max_idmap_lines": args.max_idmap_lines,
        },
        "master_stats": master_stats,
        "aux_summary": aux_summary,
        "id_mapping_stats": dict(idmap_stats),
        "mapping": mapping_report,
        "output": {
            "rows": out_count,
            "db_counts": dict(output_db_counter),
            "source_counts": dict(output_source_counter),
        },
        "gates": {
            "backlink_rate_gate": {
                "required": args.min_backlink_rate,
                "actual": backlink_rate,
                "passed": backlink_rate >= args.min_backlink_rate,
            }
        },
    }

    write_json(args.build_report, build_report)
    write_json(args.coverage_report, coverage_report)
    write_json(args.conflicts_report, conflicts_report)

    print(f"[OK] wrote table: {args.out_table} rows={out_count}")
    print(f"[OK] wrote build report: {args.build_report}")
    print(f"[OK] wrote coverage report: {args.coverage_report}")
    print(f"[OK] wrote conflicts report: {args.conflicts_report}")
    print(
        "[INFO] backlink_rate="
        f"{backlink_rate:.6f} (required>={args.min_backlink_rate:.2f}), "
        f"master_coverage={master_coverage_rate:.6f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
