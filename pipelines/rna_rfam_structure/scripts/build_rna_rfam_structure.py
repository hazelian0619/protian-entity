#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import defaultdict
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

RFAM_RE = re.compile(r"^RF\d{5}$")
FETCH_DATE = date.today().isoformat()

OUTPUT_COLUMNS = [
    "rna_id",
    "rfam_id",
    "rfam_desc",
    "secondary_structure",
    "structure_source",
    "score",
    "evalue",
    "source_version",
    "fetch_date",
]

DOTBRACKET_ALLOWED = set(".()[]{}<>,-_:~")


def open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def resolve_id_mapping(path: Path) -> Optional[Path]:
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


def detect_xref_candidates() -> List[Path]:
    return [
        Path("data/output/rna_xref_mrna_enst_urs_v2.tsv"),
        Path("data/output/rna_xref_mrna_enst_urs_v1.tsv"),
    ]


def read_rna_master(path: Path, limit: Optional[int] = None) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"rna_id", "rna_type"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise SystemExit(f"[ERROR] rna_master 缺少必需列: {sorted(required)}")
        for i, row in enumerate(reader):
            rows.append(
                {
                    "rna_id": str(row.get("rna_id", "")).strip(),
                    "rna_type": str(row.get("rna_type", "")).strip().lower(),
                    "rna_name": str(row.get("rna_name", "")).strip(),
                    "symbol": str(row.get("symbol", "")).strip(),
                }
            )
            if limit is not None and len(rows) >= limit:
                break
    return rows


def load_enst_to_urs(xref_path: Optional[Path]) -> Dict[str, str]:
    if xref_path is None or not xref_path.exists():
        return {}
    mapping: Dict[str, str] = {}
    with xref_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            enst = str(row.get("rna_id", "")).strip()
            urs = str(row.get("xref_id", "")).strip()
            if enst.startswith("ENST") and urs.startswith("URS"):
                mapping[enst] = urs
    return mapping


def canonical_urs_from_rna_id(rna_id: str, enst_to_urs: Dict[str, str]) -> Optional[str]:
    if rna_id.startswith("URS") and rna_id.endswith("_9606"):
        return rna_id
    if rna_id.startswith("ENST"):
        return enst_to_urs.get(rna_id)
    return None


def load_rfam_hits_by_urs(rfam_annotations_path: Path, wanted_urs: set[str]) -> Dict[str, List[Dict[str, str]]]:
    # RNAcentral rfam_annotations format:
    # URS  RFAM_ID  SCORE  EVALUE  SEQ_START  SEQ_STOP  MODEL_START  MODEL_STOP  DESCRIPTION
    if not wanted_urs:
        return {}
    wanted_urs_bare = {u.split("_", 1)[0] for u in wanted_urs} if wanted_urs else set()
    best: Dict[Tuple[str, str], Dict[str, str]] = {}
    with open_maybe_gzip(rfam_annotations_path) as f:
        for line in f:
            if not line.strip():
                continue
            # fast path: only parse full row for wanted URS
            first_tab = line.find("\t")
            if first_tab <= 0:
                continue
            urs_raw = line[:first_tab]
            if wanted_urs_bare and urs_raw not in wanted_urs_bare:
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            urs_raw, rfam_id = parts[0], parts[1]
            if not RFAM_RE.match(rfam_id):
                continue
            urs = urs_raw if urs_raw.endswith("_9606") else f"{urs_raw}_9606"
            if wanted_urs and urs not in wanted_urs:
                continue

            score = parts[2].strip()
            evalue = parts[3].strip()
            desc = parts[8].strip()
            key = (urs, rfam_id)
            old = best.get(key)
            if old is None:
                best[key] = {"rfam_id": rfam_id, "score": score, "evalue": evalue, "rfam_desc": desc}
                continue
            # keep higher score when duplicated
            try:
                old_s = float(old.get("score") or 0.0)
            except Exception:
                old_s = 0.0
            try:
                new_s = float(score or 0.0)
            except Exception:
                new_s = 0.0
            if new_s > old_s:
                best[key] = {"rfam_id": rfam_id, "score": score, "evalue": evalue, "rfam_desc": desc}

    out: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for (urs, _), hit in best.items():
        out[urs].append(hit)
    return out


def load_rfam_family_desc(path: Path) -> Dict[str, str]:
    desc: Dict[str, str] = {}
    with open_maybe_gzip(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            rfam_id = None
            rfam_idx = -1
            for i, col in enumerate(cols):
                c = col.strip()
                if RFAM_RE.match(c):
                    rfam_id = c
                    rfam_idx = i
                    break
            if rfam_id is None:
                continue

            text_candidates = [c.strip() for c in cols[rfam_idx + 1 :] if c.strip()]
            if not text_candidates:
                continue
            desc[rfam_id] = text_candidates[0]
    return desc


def load_rfam_family_name_alias(path: Path) -> Dict[str, List[str]]:
    alias: Dict[str, List[str]] = defaultdict(list)
    with open_maybe_gzip(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            rfam_id = cols[0].strip()
            family_name = cols[1].strip().lower()
            if not RFAM_RE.match(rfam_id) or not family_name:
                continue
            if rfam_id not in alias[family_name]:
                alias[family_name].append(rfam_id)
    return alias


def _rm_tail_letter_after_digit(s: str) -> str:
    return re.sub(r"(?<=\d)[a-z]$", "", s)


def mirna_family_candidates(rna_name: str, symbol: str) -> List[str]:
    cands: List[str] = []

    def add(x: str) -> None:
        t = (x or "").strip().lower()
        if t and t not in cands:
            cands.append(t)

    n = (rna_name or "").strip().lower()
    s = (symbol or "").strip().lower()

    if n:
        add(n)
        n1 = re.sub(r"^[a-z]{3}-", "", n)  # hsa-miR-XXX -> miR-XXX
        add(n1)
        n2 = re.sub(r"-(3p|5p)$", "", n1)  # miR-XXX-3p -> miR-XXX
        add(n2)
        n3 = re.sub(r"-\d+$", "", n2)  # let-7a-2 -> let-7a
        add(n3)
        add(_rm_tail_letter_after_digit(n3))  # mir-320d -> mir-320

    if s:
        add(s)
        m_let = re.match(r"^mirlet(\d+)[a-z]?$", s)
        if m_let:
            add(f"let-{m_let.group(1)}")
        m_mir = re.match(r"^mir(\d+[a-z]?)$", s)
        if m_mir:
            raw = m_mir.group(1)
            add(f"mir-{raw}")
            add(f"mir-{_rm_tail_letter_after_digit(raw)}")

    # extra fallbacks:
    expanded: List[str] = []
    for x in cands:
        expanded.append(x)
        if x.startswith("mir-"):
            trimmed = re.sub(r"-(\d+)$", "", x)
            if trimmed and trimmed not in expanded:
                expanded.append(trimmed)

    # de-duplicate preserving order
    out: List[str] = []
    for x in expanded:
        if x and x not in out:
            out.append(x)
    return out


def heuristic_rfam_hits_for_mirna(
    row: Dict[str, str],
    family_alias: Dict[str, List[str]],
) -> List[Dict[str, str]]:
    hits: List[Dict[str, str]] = []
    if (row.get("rna_type") or "").lower() != "mirna":
        return hits

    candidates = mirna_family_candidates(row.get("rna_name", ""), row.get("symbol", ""))
    for key in candidates:
        rfam_ids = family_alias.get(key, [])
        if not rfam_ids:
            continue
        for rfam_id in rfam_ids:
            hits.append(
                {
                    "rfam_id": rfam_id,
                    "score": "",
                    "evalue": "",
                    "rfam_desc": "",
                    "_hit_source": "rfam_family_name_heuristic",
                }
            )
        break
    return hits


def load_rfam_seed_consensus(path: Path) -> Dict[str, str]:
    consensus: Dict[str, str] = {}
    current_acc: Optional[str] = None
    ss_buf: List[str] = []

    with open_maybe_gzip(path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith("#=GF AC"):
                parts = line.split()
                if len(parts) >= 3 and RFAM_RE.match(parts[2]):
                    current_acc = parts[2]
                    ss_buf = []
            elif line.startswith("#=GC SS_cons") and current_acc:
                parts = line.split(maxsplit=2)
                if len(parts) == 3:
                    ss_buf.append(parts[2].strip())
            elif line.strip() == "//":
                if current_acc and ss_buf:
                    consensus[current_acc] = "".join(ss_buf)
                current_acc = None
                ss_buf = []

    return consensus


def is_valid_dotbracket(s: str) -> bool:
    if s is None:
        return False
    seq = s.strip()
    if seq == "":
        return False
    return set(seq) <= DOTBRACKET_ALLOWED


def make_download_checklist(
    rna_master_path: Path,
    id_mapping_path: Path,
    rfam_annotations_path: Path,
    rfam_family_path: Path,
    rfam_seed_path: Path,
) -> List[Dict[str, str]]:
    return [
        {
            "file": "RNA 主表",
            "required_path": str(rna_master_path),
            "suggested_source": "由 pipelines/rna/run.sh 产出",
            "sha256_command": f"sha256sum {rna_master_path}",
        },
        {
            "file": "RNAcentral id_mapping",
            "required_path": str(id_mapping_path),
            "suggested_source": "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/",
            "sha256_command": f"sha256sum {id_mapping_path}",
        },
        {
            "file": "RNAcentral Rfam annotations",
            "required_path": str(rfam_annotations_path),
            "suggested_source": "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/rfam/rfam_annotations.tsv.gz",
            "sha256_command": f"sha256sum {rfam_annotations_path}",
        },
        {
            "file": "Rfam family metadata",
            "required_path": str(rfam_family_path),
            "suggested_source": "https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/",
            "sha256_command": f"sha256sum {rfam_family_path}",
        },
        {
            "file": "Rfam seed alignment",
            "required_path": str(rfam_seed_path),
            "suggested_source": "https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/",
            "sha256_command": f"sha256sum {rfam_seed_path}",
        },
        {
            "file": "mRNA ENST↔URS 映射（若主表含 ENST）",
            "required_path": "data/output/rna_xref_mrna_enst_urs_v2.tsv（或 v1）",
            "suggested_source": "助手A产物（rna_xref_enst_urs）",
            "sha256_command": "sha256sum data/output/rna_xref_mrna_enst_urs_v2.tsv",
        },
    ]


def write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build RNA Rfam structure mapping table.")
    parser.add_argument("--rna-master", type=Path, default=Path("data/output/rna_master_v1.tsv"))
    parser.add_argument("--id-mapping", type=Path, default=Path("data/raw/rna/rnacentral/id_mapping.tsv.gz"))
    parser.add_argument(
        "--rfam-annotations",
        type=Path,
        default=Path("data/raw/rna/rnacentral/rfam_annotations.tsv.gz"),
    )
    parser.add_argument("--rfam-family", type=Path, default=Path("data/raw/rna/rfam/family.txt.gz"))
    parser.add_argument("--rfam-seed", type=Path, default=Path("data/raw/rna/rfam/Rfam.seed.gz"))
    parser.add_argument("--xref", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=Path("data/output/rna_structure_rfam_v1.tsv"))
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("pipelines/rna_rfam_structure/reports/rna_structure_rfam_v1.metrics.json"),
    )
    parser.add_argument("--limit", type=int, default=None, help="最小样本模式，仅处理前 N 行 rna_master")
    parser.add_argument("--source-version", default="RNAcentral:25;Rfam:15.0")
    parser.add_argument("--fetch-date", default=FETCH_DATE)
    parser.add_argument("--check-inputs", action="store_true", help="仅检查输入可用性并输出阻塞报告")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    id_mapping_path = resolve_id_mapping(args.id_mapping)
    xref_path = args.xref
    if xref_path is None:
        for p in detect_xref_candidates():
            if p.exists():
                xref_path = p
                break

    missing_required: List[str] = []
    if not args.rna_master.exists():
        missing_required.append(str(args.rna_master))
    if id_mapping_path is None:
        missing_required.append(str(args.id_mapping))
    if not args.rfam_annotations.exists():
        missing_required.append(str(args.rfam_annotations))
    if not args.rfam_family.exists():
        missing_required.append(str(args.rfam_family))
    if not args.rfam_seed.exists():
        missing_required.append(str(args.rfam_seed))

    base_report = {
        "pipeline": "rna_rfam_structure_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": {
            "rna_master": str(args.rna_master),
            "id_mapping": str(id_mapping_path) if id_mapping_path else str(args.id_mapping),
            "rfam_annotations": str(args.rfam_annotations),
            "rfam_family": str(args.rfam_family),
            "rfam_seed": str(args.rfam_seed),
            "xref": str(xref_path) if xref_path else None,
        },
    }

    if missing_required:
        blocked = {
            **base_report,
            "status": "blocked_missing_inputs",
            "missing_required": missing_required,
            "download_checklist": make_download_checklist(
                args.rna_master, args.id_mapping, args.rfam_annotations, args.rfam_family, args.rfam_seed
            ),
            "message": "缺外部数据，按协作约束已停止执行。",
        }
        write_json(args.report, blocked)
        print(f"[BLOCKED] missing required inputs: {missing_required}")
        print(f"[BLOCKED] report -> {args.report}")
        return 2

    rna_rows = read_rna_master(args.rna_master, limit=args.limit)
    if len(rna_rows) == 0:
        blocked = {
            **base_report,
            "status": "blocked_empty_rna_master",
            "message": "rna_master 行数为 0，无法继续。",
        }
        write_json(args.report, blocked)
        print("[BLOCKED] rna_master is empty")
        return 2

    enst_count = sum(1 for r in rna_rows if r["rna_id"].startswith("ENST"))
    enst_to_urs = load_enst_to_urs(xref_path)

    if enst_count > 0 and len(enst_to_urs) == 0:
        blocked = {
            **base_report,
            "status": "blocked_missing_enst_urs_xref",
            "message": "rna_master 含 ENST 主键，但未找到 ENST↔URS 映射文件。",
            "enst_rows": enst_count,
            "download_checklist": make_download_checklist(
                args.rna_master, args.id_mapping, args.rfam_annotations, args.rfam_family, args.rfam_seed
            ),
        }
        write_json(args.report, blocked)
        print("[BLOCKED] ENST rows found but xref file missing")
        print(f"[BLOCKED] report -> {args.report}")
        return 2

    if args.check_inputs:
        ok_report = {
            **base_report,
            "status": "inputs_ready",
            "rna_master_rows_checked": len(rna_rows),
            "enst_rows": enst_count,
            "xref_rows": len(enst_to_urs),
            "message": "输入检查通过。",
        }
        write_json(args.report, ok_report)
        print(f"[OK] inputs ready -> {args.report}")
        return 0

    canonical_by_rna: Dict[str, Optional[str]] = {}
    wanted_urs: set[str] = set()
    no_canonical = 0
    for row in rna_rows:
        rid = row["rna_id"]
        urs = canonical_urs_from_rna_id(rid, enst_to_urs)
        canonical_by_rna[rid] = urs
        if urs:
            wanted_urs.add(urs)
        else:
            no_canonical += 1

    urs_to_hits = load_rfam_hits_by_urs(args.rfam_annotations, wanted_urs)
    family_desc = load_rfam_family_desc(args.rfam_family)
    family_alias = load_rfam_family_name_alias(args.rfam_family)
    seed_ss = load_rfam_seed_consensus(args.rfam_seed)

    out_rows: List[Dict[str, str]] = []
    mapped_rna_ids: set[str] = set()
    no_rfam_hits = 0
    mapped_via_rnacentral = 0
    mapped_via_heuristic = 0

    for row in rna_rows:
        rid = row["rna_id"]
        urs = canonical_by_rna.get(rid)
        hits: List[Dict[str, str]] = []
        if urs:
            hits = list(urs_to_hits.get(urs, []))
            for h in hits:
                h["_hit_source"] = "rnacentral_rfam_xref"

        if not hits:
            h_hits = heuristic_rfam_hits_for_mirna(row, family_alias)
            if h_hits:
                hits = h_hits
            else:
                if urs:
                    no_rfam_hits += 1
                continue

        mapped_rna_ids.add(rid)
        for hit in hits:
            rfam_id = hit["rfam_id"]
            ss = seed_ss.get(rfam_id, "")
            hit_source = hit.get("_hit_source", "rnacentral_rfam_xref")
            if hit_source == "rfam_family_name_heuristic":
                mapped_via_heuristic += 1
            else:
                mapped_via_rnacentral += 1

            structure_source = f"{hit_source}+rfam_seed_consensus" if ss else hit_source
            out_rows.append(
                {
                    "rna_id": rid,
                    "rfam_id": rfam_id,
                    "rfam_desc": family_desc.get(rfam_id, "") or hit.get("rfam_desc", ""),
                    "secondary_structure": ss,
                    "structure_source": structure_source,
                    "score": hit.get("score", ""),
                    "evalue": hit.get("evalue", ""),
                    "source_version": args.source_version,
                    "fetch_date": args.fetch_date,
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    total = len(rna_rows)
    mapped = len(mapped_rna_ids)
    overall_cov = (mapped / total) if total else 0.0

    type_total = defaultdict(int)
    type_mapped = defaultdict(int)
    for row in rna_rows:
        t = row["rna_type"] or "unknown"
        type_total[t] += 1
        if row["rna_id"] in mapped_rna_ids:
            type_mapped[t] += 1

    rfam_valid = sum(1 for r in out_rows if RFAM_RE.match(r["rfam_id"]))
    rfam_valid_rate = (rfam_valid / len(out_rows)) if out_rows else 0.0

    ss_non_empty = [r["secondary_structure"] for r in out_rows if r["secondary_structure"].strip()]
    ss_non_empty_rate = (len(ss_non_empty) / len(out_rows)) if out_rows else 0.0
    ss_valid = sum(1 for s in ss_non_empty if is_valid_dotbracket(s))
    ss_valid_rate = (ss_valid / len(ss_non_empty)) if ss_non_empty else 1.0

    metrics = {
        **base_report,
        "status": "completed",
        "output": str(args.output),
        "row_count": {
            "rna_master": total,
            "output_rows": len(out_rows),
            "mapped_unique_rna": mapped,
        },
        "coverage": {
            "overall_mapped_rate": overall_cov,
            "by_rna_type": {
                t: {
                    "total": type_total[t],
                    "mapped": type_mapped[t],
                    "mapped_rate": (type_mapped[t] / type_total[t]) if type_total[t] else 0.0,
                }
                for t in sorted(type_total.keys())
            },
        },
        "quality": {
            "rfam_id_valid_rows": rfam_valid,
            "rfam_id_valid_rate": rfam_valid_rate,
            "secondary_structure_non_empty_rate": ss_non_empty_rate,
            "secondary_structure_valid_rate_non_empty": ss_valid_rate,
        },
        "unmatched_analysis": {
            "no_canonical_urs": no_canonical,
            "no_rfam_hits_after_urs_mapping": no_rfam_hits,
        },
        "mapping_sources": {
            "rnacentral_rfam_xref_hits": mapped_via_rnacentral,
            "rfam_family_name_heuristic_hits": mapped_via_heuristic,
        },
        "acceptance": {
            "rna_id_joinable": True,
            "rfam_id_format_pass": rfam_valid_rate == 1.0 if out_rows else False,
            "coverage_by_type_reported": True,
            "dot_bracket_validity_pass": ss_valid_rate == 1.0,
        },
    }

    write_json(args.report, metrics)

    passed = (
        metrics["acceptance"]["rfam_id_format_pass"]
        and metrics["acceptance"]["dot_bracket_validity_pass"]
        and len(out_rows) > 0
    )
    print(
        f"[{'PASS' if passed else 'WARN'}] rows={len(out_rows)} mapped={mapped}/{total} "
        f"coverage={overall_cov:.4f} rfam_valid={rfam_valid_rate:.4f} ss_valid={ss_valid_rate:.4f}"
    )
    print(f"[OK] output -> {args.output}")
    print(f"[OK] report -> {args.report}")

    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
