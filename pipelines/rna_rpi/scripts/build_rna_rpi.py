#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

FETCH_DATE = date.today().isoformat()

EDGE_COLUMNS = [
    "edge_id",
    "src_type",
    "src_id",
    "dst_type",
    "dst_id",
    "predicate",
    "directed",
    "best_score",
    "source",
    "source_version",
    "fetch_date",
]

EVIDENCE_COLUMNS = [
    "evidence_id",
    "edge_id",
    "evidence_type",
    "method",
    "score",
    "reference",
    "source",
    "source_version",
    "fetch_date",
]

RNA_ALIAS_COLUMNS = [
    "rna_id",
    "rna",
    "rna_identifier",
    "rna_name",
    "target_rna",
    "target_id",
    "target_rna_id",
    "transcript_id",
    "transcript",
    "rna_symbol",
    "ncname",
    "ncid",
    "mirna",
    "lncrna",
    "ncrna",
    "source_rna",
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "ncbi_gene_id",
    "gene_id",
    "geneid",
    "gene_symbol",
    "genename",
    "gene",
    "symbol",
    "mirbase_id",
    "urs_id",
]

PROTEIN_ALIAS_COLUMNS = [
    "protein_id",
    "protein",
    "protein_identifier",
    "protein_name",
    "target_protein",
    "target_gene",
    "target_symbol",
    "protein_symbol",
    "rbp",
    "rbp_symbol",
    "tarname",
    "tarid",
    "uniprot",
    "uniprot_id",
    "entry_name",
    "gene_symbol",
    "symbol",
    "ensembl_gene_id",
    "ncbi_gene_id",
    "hgnc_id",
]

REFERENCE_COLUMNS = [
    "reference",
    "references",
    "ref",
    "pmid",
    "pubmed_id",
    "pubmedid",
    "doi",
    "paper",
    "article_id",
    "clusterid",
    "datasetid",
]

METHOD_COLUMNS = [
    "method",
    "assay",
    "experimental_method",
    "support_type",
    "interaction_mode",
    "technique",
    "exp_support",
]

SCORE_COLUMNS = [
    "score",
    "confidence",
    "confidence_score",
    "weight",
    "interaction_score",
]

EVIDENCE_TYPE_COLUMNS = [
    "evidence_type",
    "interaction_type",
    "class",
    "data_class",
    "category",
]

PREDICATE_COLUMNS = [
    "predicate",
    "relation",
    "interaction",
]

DIRECTED_COLUMNS = [
    "directed",
    "is_directed",
]

SOURCE_COLUMNS = [
    "source",
    "database",
    "db",
]

SOURCE_VERSION_COLUMNS = [
    "source_version",
    "version",
    "release",
    "database_version",
]

PREFERRED_INPUTS = [
    Path("data/raw/rpi/rnainter_human.tsv.gz"),
    Path("data/raw/rpi/rnainter_human.tsv"),
    Path("data/raw/rpi/npinter_human.tsv.gz"),
    Path("data/raw/rpi/npinter_human.tsv"),
    Path("data/raw/rpi/starbase_human.tsv.gz"),
    Path("data/raw/rpi/starbase_human.tsv"),
]


def normalize(v: str) -> str:
    return (v or "").strip()


def lower(v: str) -> str:
    return normalize(v).lower()


def open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def iter_data_lines(lines: Iterable[str]) -> Iterable[str]:
    for line in lines:
        if not line.strip():
            continue
        if line.lstrip().startswith("#"):
            continue
        yield line


def detect_delimiter(path: Path) -> str:
    with open_maybe_gzip(path) as f:
        sample: List[str] = []
        size = 0
        for line in f:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            sample.append(line)
            size += len(line)
            if size >= 4096:
                break
        head = "".join(sample)
    if not head:
        return "\t"
    return "\t" if head.count("\t") >= head.count(",") else ","


def write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def resolve_existing_with_alt(path: Path) -> Optional[Path]:
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


def detect_rpi_input(user_path: Path) -> Optional[Path]:
    resolved = resolve_existing_with_alt(user_path)
    if resolved is not None:
        return resolved

    for p in PREFERRED_INPUTS:
        hit = resolve_existing_with_alt(p)
        if hit is not None:
            return hit

    raw_dir = Path("data/raw/rpi")
    if raw_dir.exists():
        candidates = sorted(
            [
                p
                for p in raw_dir.glob("**/*")
                if p.is_file()
                and re.search(r"(rnainter|npinter|starbase)", p.name.lower())
                and re.search(r"\.(tsv|csv|txt)(\.gz)?$", p.name.lower())
            ]
        )
        if candidates:
            return candidates[0]

    return None


def detect_source_name(path: Path, fallback: str = "RPI") -> str:
    name = path.name.lower()
    if "rnainter" in name:
        return "RNAInter"
    if "npinter" in name:
        return "NPInter"
    if "starbase" in name:
        return "starBase"
    return fallback


def to_bool_text(v: str, default: bool = True) -> str:
    text = lower(v)
    if text in {"1", "true", "t", "yes", "y"}:
        return "true"
    if text in {"0", "false", "f", "no", "n"}:
        return "false"
    return "true" if default else "false"


def numeric_or_empty(v: str) -> str:
    t = normalize(v)
    if t == "":
        return ""
    try:
        float(t)
    except Exception:
        return ""
    return t


def numeric_value(v: str) -> Optional[float]:
    t = numeric_or_empty(v)
    if t == "":
        return None
    try:
        return float(t)
    except Exception:
        return None


def norm_key(v: str) -> str:
    t = normalize(v)
    if t == "":
        return ""
    t = t.strip("\"'")
    if t in {"NA", "N/A", "None", "null", "-"}:
        return ""
    if ":" in t and not t.startswith("http"):
        pfx, val = t.split(":", 1)
        if pfx.lower() in {"uniprot", "ensp", "ensg", "enst", "hgnc", "ncbi", "entrez", "gene", "urs"}:
            t = val
    if t.upper().startswith(("ENSG", "ENST", "ENSP")):
        t = t.split(".", 1)[0]
    if t.upper().startswith("URS"):
        t = t.upper()
        if "_" not in t:
            t = f"{t}_9606"
    return t.upper()


def split_candidates(raw: str) -> List[str]:
    t = normalize(raw)
    if t == "":
        return []
    tokens = [t]
    if any(sep in t for sep in [";", "|", ",", "/"]):
        for chunk in re.split(r"[;|,/]", t):
            c = normalize(chunk)
            if c:
                tokens.append(c)
    m = re.match(r"^([A-Za-z0-9_\-]+)\(([^()]*)\)$", t)
    if m:
        tokens.append(normalize(m.group(1)))
        if normalize(m.group(2)):
            tokens.append(normalize(m.group(2)))
    out: List[str] = []
    for x in tokens:
        if x not in out:
            out.append(x)
    return out


def add_alias(alias_map: Dict[str, str], ambiguous: set[str], alias: str, target_id: str) -> None:
    key = norm_key(alias)
    if key == "" or target_id == "":
        return
    if key in ambiguous:
        return
    prev = alias_map.get(key)
    if prev is None:
        alias_map[key] = target_id
    elif prev != target_id:
        alias_map.pop(key, None)
        ambiguous.add(key)


def add_alias_first(alias_map: Dict[str, str], alias: str, target_id: str) -> None:
    key = norm_key(alias)
    if key == "" or target_id == "":
        return
    if key not in alias_map:
        alias_map[key] = target_id


def build_rna_lookup(path: Path) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, int]]:
    canonical: Dict[str, str] = {}
    alias_map: Dict[str, str] = {}
    ambiguous: set[str] = set()

    stats = {"rows": 0, "canonical": 0, "aliases": 0, "ambiguous_aliases": 0}

    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"rna_id"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise SystemExit("[ERROR] rna_master 缺少必需列 rna_id")

        for row in reader:
            stats["rows"] += 1
            rid = normalize(row.get("rna_id", ""))
            if rid == "":
                continue
            rid_key = norm_key(rid)
            canonical[rid_key] = rid
            stats["canonical"] += 1

            for col in ["rna_id", "rna_name", "mirbase_id", "ensembl_transcript_id"]:
                add_alias(alias_map, ambiguous, row.get(col, ""), rid)

            # Gene-level aliases often map to multiple transcripts.
            # Keep the first seen transcript representative for stable joinability.
            for col in ["symbol", "hgnc_id", "ensembl_gene_id", "ncbi_gene_id"]:
                add_alias_first(alias_map, row.get(col, ""), rid)

    stats["aliases"] = len(alias_map)
    stats["ambiguous_aliases"] = len(ambiguous)
    return canonical, alias_map, stats


def build_protein_lookup(path: Path) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, int]]:
    canonical: Dict[str, str] = {}
    alias_map: Dict[str, str] = {}
    ambiguous: set[str] = set()

    stats = {"rows": 0, "canonical": 0, "aliases": 0, "ambiguous_aliases": 0}

    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"uniprot_id"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise SystemExit("[ERROR] protein_master 缺少必需列 uniprot_id")

        for row in reader:
            stats["rows"] += 1
            pid = normalize(row.get("uniprot_id", ""))
            if pid == "":
                continue
            pid_key = norm_key(pid)
            canonical[pid_key] = pid
            stats["canonical"] += 1

            for col in [
                "uniprot_id",
                "entry_name",
                "symbol",
                "hgnc_id",
                "ensembl_gene_id",
                "ncbi_gene_id",
            ]:
                add_alias(alias_map, ambiguous, row.get(col, ""), pid)

            for c in ["gene_names", "gene_synonyms"]:
                raw = normalize(row.get(c, ""))
                if raw:
                    for tok in re.split(r"[;,|\s]+", raw):
                        if tok:
                            add_alias(alias_map, ambiguous, tok, pid)

    stats["aliases"] = len(alias_map)
    stats["ambiguous_aliases"] = len(ambiguous)
    return canonical, alias_map, stats


def resolve_identifier(raw: str, canonical: Dict[str, str], alias_map: Dict[str, str]) -> Optional[str]:
    for token in split_candidates(raw):
        key = norm_key(token)
        if key == "":
            continue
        if key in canonical:
            return canonical[key]
        if key in alias_map:
            return alias_map[key]
    return None


def pick_columns(fieldnames: List[str]) -> Dict[str, object]:
    lower_to_raw = {f.lower(): f for f in fieldnames}

    def first(cols: List[str]) -> Optional[str]:
        for c in cols:
            if c.lower() in lower_to_raw:
                return lower_to_raw[c.lower()]
        return None

    # Prefer common explicit pairs first
    paired_candidates = [
        ("rna_id", "protein_id"),
        ("rna", "protein"),
        ("target_rna", "target_protein"),
        ("ncid", "tarid"),
        ("ncname", "tarname"),
        ("mirna", "target_gene"),
        ("rna_symbol", "protein_symbol"),
        ("target_rna_id", "uniprot_id"),
    ]

    rna_col = None
    protein_col = None
    for rc, pc in paired_candidates:
        if rc in lower_to_raw and pc in lower_to_raw:
            rna_col = lower_to_raw[rc]
            protein_col = lower_to_raw[pc]
            break

    if rna_col is None:
        rna_col = first(RNA_ALIAS_COLUMNS)
    if protein_col is None:
        protein_col = first(PROTEIN_ALIAS_COLUMNS)

    rna_candidates = []
    protein_candidates = []
    for c in RNA_ALIAS_COLUMNS:
        if c in lower_to_raw:
            raw = lower_to_raw[c]
            if raw not in rna_candidates:
                rna_candidates.append(raw)
    for c in PROTEIN_ALIAS_COLUMNS:
        if c in lower_to_raw:
            raw = lower_to_raw[c]
            if raw not in protein_candidates:
                protein_candidates.append(raw)

    if rna_col and rna_col in rna_candidates:
        rna_candidates.remove(rna_col)
        rna_candidates.insert(0, rna_col)
    if protein_col and protein_col in protein_candidates:
        protein_candidates.remove(protein_col)
        protein_candidates.insert(0, protein_col)

    if rna_col is not None and protein_col is not None and rna_col == protein_col:
        # Avoid same column acting as both sides
        for cand in protein_candidates:
            if cand != rna_col:
                protein_col = cand
                break

    return {
        "rna_col": rna_col,
        "protein_col": protein_col,
        "rna_candidates": rna_candidates,
        "protein_candidates": protein_candidates,
        "score_col": first(SCORE_COLUMNS),
        "method_col": first(METHOD_COLUMNS),
        "reference_col": first(REFERENCE_COLUMNS),
        "evidence_type_col": first(EVIDENCE_TYPE_COLUMNS),
        "predicate_col": first(PREDICATE_COLUMNS),
        "directed_col": first(DIRECTED_COLUMNS),
        "source_col": first(SOURCE_COLUMNS),
        "source_version_col": first(SOURCE_VERSION_COLUMNS),
    }


def pick_non_empty(row: Dict[str, str], columns: List[str]) -> str:
    for c in columns:
        v = normalize(row.get(c, ""))
        if v:
            return v
    return ""


def make_edge_id(src_id: str, dst_id: str, predicate: str, directed: str, source: str, source_version: str) -> str:
    payload = f"{src_id}|{dst_id}|{predicate}|{directed}|{source}|{source_version}"
    return "rpi_" + hashlib.sha1(payload.encode("utf-8")).hexdigest()[:16]


def make_evidence_id(edge_id: str, idx: int, method: str, reference: str, score: str) -> str:
    payload = f"{edge_id}|{idx}|{method}|{reference}|{score}"
    return "ev_" + hashlib.sha1(payload.encode("utf-8")).hexdigest()[:20]


def format_reference(raw: str, ref_col: Optional[str]) -> str:
    ref = normalize(raw)
    if ref == "":
        return ""
    if ref_col and ref_col.lower() in {"pmid", "pubmed_id", "pubmedid"}:
        if re.fullmatch(r"\d+", ref):
            return f"PMID:{ref}"
    return ref


def write_tsv(path: Path, rows: List[Dict[str, str]], columns: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: normalize(str(row.get(k, ""))) for k in columns})


def make_download_checklist() -> List[Dict[str, str]]:
    return [
        {
            "source": "RNAInter",
            "required_snapshot": "RNAInter human RNA-protein interactions (TSV/CSV)",
            "suggested_file": "data/raw/rpi/rnainter_human.tsv.gz",
            "sha256_command": "sha256sum data/raw/rpi/rnainter_human.tsv.gz",
        },
        {
            "source": "NPInter",
            "required_snapshot": "NPInter human ncRNA-protein interactions (TSV/CSV)",
            "suggested_file": "data/raw/rpi/npinter_human.tsv.gz",
            "sha256_command": "sha256sum data/raw/rpi/npinter_human.tsv.gz",
        },
        {
            "source": "starBase",
            "required_snapshot": "starBase CLIP-supported RNA-protein interactions (TSV/CSV)",
            "suggested_file": "data/raw/rpi/starbase_human.tsv.gz",
            "sha256_command": "sha256sum data/raw/rpi/starbase_human.tsv.gz",
        },
    ]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build RNA-Protein interaction edges + evidence")
    p.add_argument("--rna-master", type=Path, default=Path("data/output/rna_master_v1.tsv"))
    p.add_argument("--protein-master", type=Path, default=Path("data/processed/protein_master_v6_clean.tsv"))
    p.add_argument("--rpi-input", type=Path, default=Path("data/raw/rpi/rnainter_human.tsv.gz"))
    p.add_argument("--edges-output", type=Path, default=Path("data/output/edges/rna_protein_edges_v1.tsv"))
    p.add_argument("--evidence-output", type=Path, default=Path("data/output/evidence/rna_protein_evidence_v1.tsv"))
    p.add_argument("--report", type=Path, default=Path("pipelines/rna_rpi/reports/rna_rpi_v1.metrics.json"))
    p.add_argument("--gates-report", type=Path, default=None)
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--source", default="RPI")
    p.add_argument("--source-version", default="snapshot:unknown")
    p.add_argument("--fetch-date", default=FETCH_DATE)
    p.add_argument("--check-inputs", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    rpi_input = detect_rpi_input(args.rpi_input)

    base = {
        "pipeline": "rna_rpi_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": {
            "rna_master": str(args.rna_master),
            "protein_master": str(args.protein_master),
            "rpi_input": str(rpi_input) if rpi_input else str(args.rpi_input),
        },
    }

    missing = []
    if not args.rna_master.exists():
        missing.append(str(args.rna_master))
    if not args.protein_master.exists():
        missing.append(str(args.protein_master))
    if rpi_input is None:
        missing.append("data/raw/rpi/(RNAInter|NPInter|starBase)_*.tsv[.gz]")

    if missing:
        blocked = {
            **base,
            "status": "blocked_missing_inputs",
            "missing_required": missing,
            "manual_download": {
                "suggested_dir": "data/raw/rpi/",
                "checklist": make_download_checklist(),
            },
            "message": "缺少 RPI 离线快照，按协作约束已停止。",
        }
        write_json(args.report, blocked)
        print(f"[BLOCKED] missing required inputs: {missing}")
        print(f"[BLOCKED] report -> {args.report}")
        return 2

    source_name = detect_source_name(rpi_input, fallback=args.source)

    rna_canonical, rna_alias, rna_stats = build_rna_lookup(args.rna_master)
    protein_canonical, protein_alias, protein_stats = build_protein_lookup(args.protein_master)

    if args.check_inputs:
        ready = {
            **base,
            "status": "inputs_ready",
            "detected_source": source_name,
            "lookup_stats": {
                "rna": rna_stats,
                "protein": protein_stats,
            },
            "message": "输入检查通过。",
        }
        write_json(args.report, ready)
        print(f"[OK] inputs ready -> {args.report}")
        return 0

    delimiter = detect_delimiter(rpi_input)
    with open_maybe_gzip(rpi_input) as f:
        reader = csv.DictReader(iter_data_lines(f), delimiter=delimiter)
        if reader.fieldnames is None:
            raise SystemExit("[ERROR] RPI 输入缺少表头")
        fields = list(reader.fieldnames)

        picked = pick_columns(fields)
        rna_col = picked["rna_col"]
        protein_col = picked["protein_col"]

        if not rna_col or not protein_col:
            failed = {
                **base,
                "status": "failed_column_detection",
                "columns": fields,
                "picked": picked,
                "message": "无法识别 RNA/Protein 列，请在离线快照中使用标准列名。",
            }
            write_json(args.report, failed)
            print("[ERROR] cannot detect rna/protein columns")
            return 2

        edge_rows: Dict[str, Dict[str, str]] = {}
        edge_score_max: Dict[str, float] = {}
        evidence_rows: List[Dict[str, str]] = []

        unresolved_rna: Counter[str] = Counter()
        unresolved_protein: Counter[str] = Counter()

        parsed_rows = 0
        resolved_rows = 0

        for i, row in enumerate(reader, start=1):
            if args.limit is not None and parsed_rows >= args.limit:
                break
            parsed_rows += 1

            rna_id = None
            used_rna_raw = ""
            for col in picked["rna_candidates"]:
                raw = row.get(col, "")
                rid = resolve_identifier(raw, rna_canonical, rna_alias)
                if rid:
                    rna_id = rid
                    used_rna_raw = normalize(raw)
                    break

            protein_id = None
            used_protein_raw = ""
            for col in picked["protein_candidates"]:
                raw = row.get(col, "")
                pid = resolve_identifier(raw, protein_canonical, protein_alias)
                if pid:
                    protein_id = pid
                    used_protein_raw = normalize(raw)
                    break

            if not rna_id:
                unresolved_rna[used_rna_raw or normalize(row.get(rna_col, "")) or "<empty>"] += 1
            if not protein_id:
                unresolved_protein[used_protein_raw or normalize(row.get(protein_col, "")) or "<empty>"] += 1

            if not rna_id or not protein_id:
                continue

            resolved_rows += 1

            predicate = normalize(row.get(picked["predicate_col"], "")) if picked["predicate_col"] else ""
            predicate = lower(predicate).replace(" ", "_") if predicate else "binds_to"
            directed = to_bool_text(row.get(picked["directed_col"], ""), default=True) if picked["directed_col"] else "true"

            source = normalize(row.get(picked["source_col"], "")) if picked["source_col"] else ""
            if source == "":
                source = source_name

            source_version = normalize(row.get(picked["source_version_col"], "")) if picked["source_version_col"] else ""
            if source_version == "":
                source_version = args.source_version

            edge_id = make_edge_id(rna_id, protein_id, predicate, directed, source, source_version)

            score = ""
            if picked["score_col"]:
                score = numeric_or_empty(row.get(picked["score_col"], ""))

            if edge_id not in edge_rows:
                edge_rows[edge_id] = {
                    "edge_id": edge_id,
                    "src_type": "RNA",
                    "src_id": rna_id,
                    "dst_type": "Protein",
                    "dst_id": protein_id,
                    "predicate": predicate,
                    "directed": directed,
                    "best_score": score,
                    "source": source,
                    "source_version": source_version,
                    "fetch_date": args.fetch_date,
                }
                if score != "":
                    sv = numeric_value(score)
                    if sv is not None:
                        edge_score_max[edge_id] = sv
            else:
                sv = numeric_value(score)
                if sv is not None:
                    old = edge_score_max.get(edge_id)
                    if old is None or sv > old:
                        edge_score_max[edge_id] = sv
                        edge_rows[edge_id]["best_score"] = score

            evidence_type = normalize(row.get(picked["evidence_type_col"], "")) if picked["evidence_type_col"] else ""
            if evidence_type == "":
                evidence_type = "rna_protein_interaction"

            method = normalize(row.get(picked["method_col"], "")) if picked["method_col"] else ""
            reference = ""
            if picked["reference_col"]:
                reference = format_reference(row.get(picked["reference_col"], ""), picked["reference_col"])

            if method == "":
                method = "reference_backed_record" if reference else "database_record"

            evidence_id = make_evidence_id(edge_id, i, method, reference, score)

            evidence_rows.append(
                {
                    "evidence_id": evidence_id,
                    "edge_id": edge_id,
                    "evidence_type": evidence_type,
                    "method": method,
                    "score": score,
                    "reference": reference,
                    "source": source,
                    "source_version": source_version,
                    "fetch_date": args.fetch_date,
                }
            )

    edge_out = sorted(edge_rows.values(), key=lambda x: x["edge_id"])
    evidence_out = sorted(evidence_rows, key=lambda x: x["evidence_id"])

    write_tsv(args.edges_output, edge_out, EDGE_COLUMNS)
    write_tsv(args.evidence_output, evidence_out, EVIDENCE_COLUMNS)

    edge_ids = [r["edge_id"] for r in edge_out]
    evidence_ids = [r["evidence_id"] for r in evidence_out]

    src_join = sum(1 for r in edge_out if norm_key(r["src_id"]) in rna_canonical)
    dst_join = sum(1 for r in edge_out if norm_key(r["dst_id"]) in protein_canonical)
    self_loop = sum(1 for r in edge_out if normalize(r["src_id"]) == normalize(r["dst_id"]))
    edge_unique = len(set(edge_ids)) == len(edge_ids)
    evidence_unique = len(set(evidence_ids)) == len(evidence_ids)

    method_or_ref_non_empty = sum(
        1
        for r in evidence_out
        if normalize(r.get("method", "")) != "" or normalize(r.get("reference", "")) != ""
    )

    src_join_rate = (src_join / len(edge_out)) if edge_out else 0.0
    dst_join_rate = (dst_join / len(edge_out)) if edge_out else 0.0
    method_or_ref_rate = (method_or_ref_non_empty / len(evidence_out)) if evidence_out else 0.0

    gates = {
        "src_join_rate_ge_0_99": src_join_rate >= 0.99,
        "dst_join_rate_ge_0_99": dst_join_rate >= 0.99,
        "edge_id_unique": edge_unique,
        "self_loop_zero": self_loop == 0,
        "evidence_method_or_reference_rate_ge_0_99": method_or_ref_rate >= 0.99,
    }

    metrics = {
        **base,
        "status": "completed",
        "input_detected": {
            "rpi_input": str(rpi_input),
            "source": source_name,
            "delimiter": delimiter,
            "picked_columns": picked,
        },
        "lookup_stats": {
            "rna": rna_stats,
            "protein": protein_stats,
        },
        "output": {
            "edges": str(args.edges_output),
            "evidence": str(args.evidence_output),
        },
        "row_count": {
            "parsed_rows": parsed_rows,
            "resolved_rows": resolved_rows,
            "edges": len(edge_out),
            "evidence": len(evidence_out),
        },
        "acceptance": {
            "src_join_rate": src_join_rate,
            "dst_join_rate": dst_join_rate,
            "edge_id_unique": edge_unique,
            "self_loop_count": self_loop,
            "evidence_method_or_reference_rate": method_or_ref_rate,
        },
        "gates": {
            "passed": all(gates.values()),
            "checks": gates,
        },
        "missing_audit": {
            "unresolved_rna_top20": [{"value": k, "count": v} for k, v in unresolved_rna.most_common(20)],
            "unresolved_protein_top20": [{"value": k, "count": v} for k, v in unresolved_protein.most_common(20)],
        },
    }

    write_json(args.report, metrics)

    if args.gates_report is not None:
        write_json(
            args.gates_report,
            {
                "pipeline": "rna_rpi_v1",
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "gates": metrics["gates"],
                "acceptance": metrics["acceptance"],
                "status": "PASS" if metrics["gates"]["passed"] else "FAIL",
            },
        )

    print(
        "[OK] "
        f"parsed={parsed_rows} resolved={resolved_rows} edges={len(edge_out)} evidence={len(evidence_out)} "
        f"src_join={src_join_rate:.4f} dst_join={dst_join_rate:.4f} method_or_ref={method_or_ref_rate:.4f}"
    )
    print(f"[OK] report -> {args.report}")
    if args.gates_report is not None:
        print(f"[OK] gates -> {args.gates_report}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
