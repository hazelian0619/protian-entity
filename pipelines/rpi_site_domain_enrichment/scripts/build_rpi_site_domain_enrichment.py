#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

FETCH_DATE = date.today().isoformat()

SITE_COLUMNS = [
    "context_id",
    "edge_id",
    "reference",
    "method_type",
    "method_subtype",
    "site_type",
    "chromosome",
    "start",
    "end",
    "strand",
    "transcript_id",
    "gene_id",
    "cell_context",
    "support_count",
    "source",
    "source_version",
    "fetch_date",
]

DOMAIN_COLUMNS = [
    "context_id",
    "edge_id",
    "protein_id",
    "domain_class",
    "interpro_id",
    "pfam_id",
    "domain_name",
    "domain_start",
    "domain_end",
    "source",
    "source_version",
    "fetch_date",
]

FUNCTION_COLUMNS = [
    "context_id",
    "edge_id",
    "protein_id",
    "function_relation",
    "inference_basis",
    "evidence_snippet",
    "source",
    "source_version",
    "fetch_date",
]

METHOD_COLUMNS = [
    "method",
    "experimental_method",
    "assay",
    "technique",
    "interaction_mode",
    "support_type",
    "exp_support",
]

REFERENCE_COLUMNS = [
    "reference",
    "references",
    "ref",
    "clusterid",
    "cluster_id",
    "pmid",
    "pubmed_id",
    "doi",
]

CHROM_COLUMNS = ["chromosome", "chr", "chrom", "seqname"]
START_COLUMNS = ["narrowstart", "start", "start_position", "genome_start", "tx_start", "broadstart"]
END_COLUMNS = ["narrowend", "end", "end_position", "genome_end", "tx_end", "broadend"]
STRAND_COLUMNS = ["strand", "orientation"]
TRANSCRIPT_COLUMNS = ["transcript_id", "target_rna_id", "rna_id", "enst", "ensembl_transcript_id"]
GENE_COLUMNS = ["geneid", "gene_id", "target_id", "ensembl_gene_id"]
CELL_COLUMNS = ["cellline/tissue", "cell_type", "biosample", "tissue"]
SUPPORT_COLUMNS = ["clipexpnum", "totalclipexpnum", "score", "confidence"]

RRM_PFAM = {"PF00076", "PF13893", "PF08777", "PF04059"}
RRM_IPR = {"IPR000504", "IPR035979", "IPR036240"}
KH_PFAM = {"PF00013", "PF07650", "PF02272", "PF13014"}
KH_IPR = {"IPR004088", "IPR003100", "IPR036975"}
CCCH_PFAM = {"PF00642", "PF14608", "PF16817", "PF16589", "PF14442"}
CCCH_IPR = {"IPR000571", "IPR013597", "IPR028919"}

FUNCTION_RELATION_RULES: List[Tuple[str, List[str]]] = [
    (
        "splicing_regulation",
        [
            r"splic",
            r"spliceosome",
            r"exon\s+junction",
            r"\bejc\b",
            r"mRNA\s+processing",
        ],
    ),
    (
        "translation_regulation",
        [
            r"translat",
            r"ribosom",
            r"ires",
            r"cap[-\s]?dependent",
            r"translation\s+initiation",
        ],
    ),
    (
        "transcription_regulation",
        [
            r"transcription",
            r"transcriptional",
            r"gene\s+expression",
            r"rna\s+polymerase",
        ],
    ),
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


def write_tsv(path: Path, rows: List[Dict[str, str]], columns: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow({k: normalize(str(row.get(k, ""))) for k in columns})


def write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    import hashlib

    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def make_id(prefix: str, payload: str) -> str:
    return f"{prefix}_" + hashlib.sha1(payload.encode("utf-8")).hexdigest()[:20]


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


def detect_source_from_name(path: Path) -> str:
    name = path.name.lower()
    if "starbase" in name or "encori" in name:
        return "starBase"
    if "rnainter" in name:
        return "RNAInter"
    if "npinter" in name:
        return "NPInter"
    return "RPI"


def discover_rpi_snapshots(raw_rpi_dir: Path) -> List[Path]:
    if not raw_rpi_dir.exists():
        return []
    out: List[Path] = []
    for p in sorted(raw_rpi_dir.glob("**/*")):
        if not p.is_file():
            continue
        low = p.name.lower()
        if not re.search(r"(starbase|encori|rnainter|npinter)", low):
            continue
        if not re.search(r"\.(tsv|csv|txt)(\.gz)?$", low):
            continue
        out.append(p)
    return out


def numeric_int_or_empty(v: str) -> str:
    t = normalize(v)
    if t == "":
        return ""
    try:
        n = int(float(t))
    except Exception:
        return ""
    return str(n)


def numeric_or_empty(v: str) -> str:
    t = normalize(v)
    if t == "":
        return ""
    try:
        float(t)
    except Exception:
        return ""
    return t


def infer_method_type(method: str, source: str) -> Tuple[str, str]:
    text = f"{method} {source}".lower()
    if any(x in text for x in ["eclip", "par-clip", "parclip", "iclip", "hits-clip", "clip"]):
        if "eclip" in text:
            return "CLIP", "eCLIP"
        if "par-clip" in text or "parclip" in text:
            return "CLIP", "PAR-CLIP"
        if "iclip" in text:
            return "CLIP", "iCLIP"
        if "hits-clip" in text:
            return "CLIP", "HITS-CLIP"
        return "CLIP", "CLIP-seq"
    if "rip" in text:
        return "RIP", "RIP-seq"
    if "chirp" in text:
        return "CHIRP", "CHIRP"
    if re.search(r"\brap\b", text):
        return "RAP", "RAP"
    if "paris" in text:
        return "PARIS", "PARIS"
    if any(x in text for x in ["starbase", "encori"]):
        return "CLIP", "CLIP-seq"
    return "OTHER", "unspecified"


def first_existing(row_lc: Dict[str, str], names: List[str]) -> str:
    for n in names:
        if n in row_lc and normalize(row_lc[n]) != "":
            return normalize(row_lc[n])
    return ""


def parse_site_records_from_snapshot(path: Path, limit: Optional[int] = None) -> Dict[str, List[Dict[str, str]]]:
    source_name = detect_source_from_name(path)
    delim = detect_delimiter(path)
    ref_map: Dict[str, List[Dict[str, str]]] = defaultdict(list)

    with open_maybe_gzip(path) as f:
        reader = csv.DictReader(iter_data_lines(f), delimiter=delim)
        if reader.fieldnames is None:
            return ref_map

        for i, row in enumerate(reader, start=1):
            if limit is not None and i > limit:
                break
            row_lc = {k.lower(): v for k, v in row.items() if k is not None}

            ref = first_existing(row_lc, REFERENCE_COLUMNS)
            if ref == "":
                continue

            raw_method = first_existing(row_lc, METHOD_COLUMNS)
            method_type, method_subtype = infer_method_type(raw_method, source_name)

            chromosome = first_existing(row_lc, CHROM_COLUMNS)
            start = numeric_int_or_empty(first_existing(row_lc, START_COLUMNS))
            end = numeric_int_or_empty(first_existing(row_lc, END_COLUMNS))
            strand = first_existing(row_lc, STRAND_COLUMNS)
            transcript_id = first_existing(row_lc, TRANSCRIPT_COLUMNS)
            gene_id = first_existing(row_lc, GENE_COLUMNS)
            cell_context = first_existing(row_lc, CELL_COLUMNS)
            support_count = numeric_or_empty(first_existing(row_lc, SUPPORT_COLUMNS))

            if chromosome and start and end:
                site_type = "genomic"
            elif transcript_id:
                site_type = "transcript"
            else:
                site_type = "unknown"

            ref_map[ref].append(
                {
                    "reference": ref,
                    "method_type": method_type,
                    "method_subtype": method_subtype,
                    "site_type": site_type,
                    "chromosome": chromosome,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "transcript_id": transcript_id,
                    "gene_id": gene_id,
                    "cell_context": cell_context,
                    "support_count": support_count,
                    "source": source_name,
                }
            )

    return ref_map


def merge_site_maps(site_maps: List[Dict[str, List[Dict[str, str]]]]) -> Dict[str, List[Dict[str, str]]]:
    merged: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for mp in site_maps:
        for ref, rows in mp.items():
            merged[ref].extend(rows)
    return merged


def load_edges(path: Path) -> Tuple[Dict[str, Dict[str, str]], List[str]]:
    edges: Dict[str, Dict[str, str]] = {}
    order: List[str] = []
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None or "edge_id" not in set(r.fieldnames):
            raise SystemExit("[ERROR] edges 输入缺少 edge_id")
        for row in r:
            eid = normalize(row.get("edge_id", ""))
            if not eid:
                continue
            edges[eid] = row
            order.append(eid)
    return edges, order


def load_evidence(path: Path, limit: Optional[int]) -> List[Dict[str, str]]:
    out: List[Dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None or "edge_id" not in set(r.fieldnames):
            raise SystemExit("[ERROR] evidence 输入缺少 edge_id")
        for i, row in enumerate(r, start=1):
            if limit is not None and i > limit:
                break
            out.append(row)
    return out


def load_domains(path: Path) -> Dict[str, List[Dict[str, str]]]:
    out: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None or "uniprot_id" not in set(r.fieldnames):
            raise SystemExit("[ERROR] protein domains 输入缺少 uniprot_id")
        for row in r:
            uid = normalize(row.get("uniprot_id", ""))
            if uid == "":
                continue
            out[uid].append(row)
    return out


def load_protein_annotations(path: Path) -> Dict[str, Dict[str, str]]:
    out: Dict[str, Dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None or "uniprot_id" not in set(r.fieldnames):
            raise SystemExit("[ERROR] protein master 输入缺少 uniprot_id")
        for row in r:
            uid = normalize(row.get("uniprot_id", ""))
            if not uid:
                continue
            text_fields = [
                normalize(row.get("protein_name", "")),
                normalize(row.get("function", "")),
                normalize(row.get("go_biological_process", "")),
                normalize(row.get("go_molecular_function", "")),
                normalize(row.get("keywords", "")),
                normalize(row.get("domains", "")),
            ]
            out[uid] = {
                "source": normalize(row.get("source", "")) or "UniProtKB",
                "source_version": "protein_master_v6_clean",
                "text": " ; ".join([x for x in text_fields if x]),
            }
    return out


def classify_domain(domain_name: str, pfam_id: str, interpro_id: str) -> str:
    p = normalize(pfam_id).upper()
    i = normalize(interpro_id).upper()
    text = f"{domain_name} {p} {i}".lower()

    if p in RRM_PFAM or i in RRM_IPR or "rna recognition motif" in text or re.search(r"\brrm\b", text):
        return "RRM"
    if p in KH_PFAM or i in KH_IPR or "k-homology" in text or re.search(r"\bkh\b", text):
        return "KH"
    if p in CCCH_PFAM or i in CCCH_IPR or "ccch" in text or "c3h" in text:
        return "CCCH"
    return "PFAM"


def infer_function_relations(annotation_text: str) -> List[Tuple[str, str, str]]:
    text = lower(annotation_text)
    if text == "":
        return [("post_transcriptional_regulation", "fallback:no_annotation", "")]

    out: List[Tuple[str, str, str]] = []
    seen: set[str] = set()
    for relation, pats in FUNCTION_RELATION_RULES:
        for p in pats:
            m = re.search(p, text)
            if m:
                if relation in seen:
                    break
                seen.add(relation)
                snippet = annotation_text[max(0, m.start() - 20) : min(len(annotation_text), m.end() + 40)]
                out.append((relation, f"keyword:{p}", normalize(snippet)))
                break

    if not out:
        out.append(("post_transcriptional_regulation", "fallback:rbp_default", ""))

    return out


def rate(n: int, d: int) -> float:
    return (n / d) if d else 0.0


def make_download_checklist() -> List[Dict[str, str]]:
    return [
        {
            "source": "RNAInter",
            "suggested_path": "data/raw/rpi/rnainter_human.tsv.gz",
            "download_hint": "RNAInter human RNA-protein snapshot",
            "sha256_command": "sha256sum data/raw/rpi/rnainter_human.tsv.gz",
        },
        {
            "source": "NPInter",
            "suggested_path": "data/raw/rpi/npinter_human.tsv.gz",
            "download_hint": "NPInter human ncRNA-protein snapshot",
            "sha256_command": "sha256sum data/raw/rpi/npinter_human.tsv.gz",
        },
        {
            "source": "starBase/ENCORI",
            "suggested_path": "data/raw/rpi/starbase_human.tsv",
            "download_hint": "ENCORI API RBPTarget table",
            "sha256_command": "sha256sum data/raw/rpi/starbase_human.tsv",
        },
        {
            "source": "Pfam/InterPro mapping",
            "suggested_path": "data/output/protein/protein_domains_interpro_v1.tsv",
            "download_hint": "from pipelines/protein_domains",
            "sha256_command": "sha256sum data/output/protein/protein_domains_interpro_v1.tsv",
        },
    ]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build RPI site/domain/function enrichment context tables")
    p.add_argument("--edges-input", type=Path, default=Path("data/output/edges/rna_protein_edges_v1.tsv"))
    p.add_argument("--evidence-input", type=Path, default=Path("data/output/evidence/rna_protein_evidence_v1.tsv"))
    p.add_argument("--protein-domains", type=Path, default=Path("data/output/protein/protein_domains_interpro_v1.tsv"))
    p.add_argument("--protein-master", type=Path, default=Path("data/processed/protein_master_v6_clean.tsv"))
    p.add_argument("--raw-rpi-dir", type=Path, default=Path("data/raw/rpi"))
    p.add_argument("--site-output", type=Path, default=Path("data/output/evidence/rpi_site_context_v2.tsv"))
    p.add_argument("--domain-output", type=Path, default=Path("data/output/evidence/rpi_domain_context_v2.tsv"))
    p.add_argument("--function-output", type=Path, default=Path("data/output/evidence/rpi_function_context_v2.tsv"))
    p.add_argument(
        "--report",
        type=Path,
        default=Path("pipelines/rpi_site_domain_enrichment/reports/rpi_site_domain_enrichment_v2.metrics.json"),
    )
    p.add_argument("--gates-report", type=Path, default=None)
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--fetch-date", default=FETCH_DATE)
    p.add_argument("--check-inputs", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    base = {
        "pipeline": "rpi_site_domain_enrichment_v2",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": {
            "edges": str(args.edges_input),
            "evidence": str(args.evidence_input),
            "protein_domains": str(args.protein_domains),
            "protein_master": str(args.protein_master),
            "raw_rpi_dir": str(args.raw_rpi_dir),
        },
        "sample_mode": args.limit is not None,
    }

    missing: List[str] = []
    for p in [args.edges_input, args.evidence_input, args.protein_domains, args.protein_master]:
        if not p.exists():
            missing.append(str(p))

    snapshots = discover_rpi_snapshots(args.raw_rpi_dir)
    snapshot_info = []
    for sp in snapshots:
        snapshot_info.append(
            {
                "path": str(sp),
                "source": detect_source_from_name(sp),
                "sha256": sha256(sp),
                "size_bytes": sp.stat().st_size,
            }
        )

    if args.check_inputs:
        if missing:
            blocked = {
                **base,
                "status": "blocked_missing_inputs",
                "missing_required": missing,
                "snapshot_detected": snapshot_info,
                "manual_download": {
                    "checklist": make_download_checklist(),
                },
                "message": "缺少必需输入，已按协作约束中断。",
            }
            write_json(args.report, blocked)
            print(f"[BLOCKED] missing inputs: {missing}")
            print(f"[BLOCKED] report -> {args.report}")
            return 2

        ready = {
            **base,
            "status": "inputs_ready",
            "snapshot_detected": snapshot_info,
            "message": "输入检查通过。",
        }
        write_json(args.report, ready)
        print(f"[OK] inputs ready -> {args.report}")
        return 0

    if missing:
        blocked = {
            **base,
            "status": "blocked_missing_inputs",
            "missing_required": missing,
            "snapshot_detected": snapshot_info,
            "manual_download": {
                "checklist": make_download_checklist(),
            },
            "message": "缺少必需输入，已按协作约束中断。",
        }
        write_json(args.report, blocked)
        print(f"[BLOCKED] missing inputs: {missing}")
        return 2

    edges, edge_order = load_edges(args.edges_input)
    evidence_rows = load_evidence(args.evidence_input, args.limit)
    domains_by_protein = load_domains(args.protein_domains)
    protein_annotations = load_protein_annotations(args.protein_master)

    if args.limit is None:
        target_edge_ids = edge_order
    else:
        seen_e: set[str] = set()
        target_edge_ids = []
        for row in evidence_rows:
            eid = normalize(row.get("edge_id", ""))
            if eid and eid in edges and eid not in seen_e:
                seen_e.add(eid)
                target_edge_ids.append(eid)

    target_edge_set = set(target_edge_ids)

    site_maps = [parse_site_records_from_snapshot(p) for p in snapshots]
    site_ref_map = merge_site_maps(site_maps)

    site_rows: List[Dict[str, str]] = []
    site_edge_with_coord: set[str] = set()
    unmatched_reference = Counter()
    matched_reference = 0

    for idx, ev in enumerate(evidence_rows, start=1):
        edge_id = normalize(ev.get("edge_id", ""))
        if edge_id == "" or edge_id not in target_edge_set:
            continue

        reference = normalize(ev.get("reference", ""))
        source = normalize(ev.get("source", ""))
        source_version = normalize(ev.get("source_version", "")) or "snapshot:unknown"

        fallback_method_type, fallback_subtype = infer_method_type(ev.get("method", ""), source)
        hit = site_ref_map.get(reference, [])

        if hit:
            matched_reference += 1
            candidate = hit[0]
        else:
            if reference:
                unmatched_reference[reference] += 1
            candidate = {
                "method_type": fallback_method_type,
                "method_subtype": fallback_subtype,
                "site_type": "unknown",
                "chromosome": "",
                "start": "",
                "end": "",
                "strand": "",
                "transcript_id": "",
                "gene_id": "",
                "cell_context": "",
                "support_count": "",
                "source": source,
            }

        method_type = normalize(candidate.get("method_type", "")) or fallback_method_type
        method_subtype = normalize(candidate.get("method_subtype", "")) or fallback_subtype

        row = {
            "context_id": make_id("rpi_site", f"{edge_id}|{reference}|{idx}"),
            "edge_id": edge_id,
            "reference": reference,
            "method_type": method_type,
            "method_subtype": method_subtype,
            "site_type": normalize(candidate.get("site_type", "")) or "unknown",
            "chromosome": normalize(candidate.get("chromosome", "")),
            "start": numeric_int_or_empty(candidate.get("start", "")),
            "end": numeric_int_or_empty(candidate.get("end", "")),
            "strand": normalize(candidate.get("strand", "")),
            "transcript_id": normalize(candidate.get("transcript_id", "")),
            "gene_id": normalize(candidate.get("gene_id", "")),
            "cell_context": normalize(candidate.get("cell_context", "")),
            "support_count": numeric_or_empty(candidate.get("support_count", "")),
            "source": source or normalize(candidate.get("source", "")) or "RPI",
            "source_version": source_version,
            "fetch_date": args.fetch_date,
        }
        site_rows.append(row)

        if row["site_type"] in {"genomic", "transcript"}:
            site_edge_with_coord.add(edge_id)

    # ensure full-edge coverage in site table (one fallback row per uncovered edge)
    missing_site_edges = target_edge_set - {r["edge_id"] for r in site_rows}
    for eid in sorted(missing_site_edges):
        edge = edges[eid]
        mtype, msub = infer_method_type("", edge.get("source", ""))
        site_rows.append(
            {
                "context_id": make_id("rpi_site", f"{eid}|fallback"),
                "edge_id": eid,
                "reference": "",
                "method_type": mtype,
                "method_subtype": msub,
                "site_type": "unknown",
                "chromosome": "",
                "start": "",
                "end": "",
                "strand": "",
                "transcript_id": "",
                "gene_id": "",
                "cell_context": "",
                "support_count": "",
                "source": normalize(edge.get("source", "")) or "RPI",
                "source_version": normalize(edge.get("source_version", "")) or "snapshot:unknown",
                "fetch_date": args.fetch_date,
            }
        )

    domain_rows: List[Dict[str, str]] = []
    domain_edge_covered: set[str] = set()
    domain_class_counter = Counter()

    for eid in target_edge_ids:
        edge = edges[eid]
        pid = normalize(edge.get("dst_id", ""))
        domain_hits = domains_by_protein.get(pid, [])
        if not domain_hits:
            continue

        domain_edge_covered.add(eid)
        seen_domain_key: set[Tuple[str, str, str, str, str]] = set()
        for d in domain_hits:
            interpro_id = normalize(d.get("interpro_id", ""))
            pfam_id = normalize(d.get("pfam_id", ""))
            domain_name = normalize(d.get("entry_name", ""))
            dstart = numeric_int_or_empty(d.get("start", ""))
            dend = numeric_int_or_empty(d.get("end", ""))
            dclass = classify_domain(domain_name, pfam_id, interpro_id)
            key = (dclass, interpro_id, pfam_id, dstart, dend)
            if key in seen_domain_key:
                continue
            seen_domain_key.add(key)
            domain_class_counter[dclass] += 1
            domain_rows.append(
                {
                    "context_id": make_id("rpi_domain", f"{eid}|{pid}|{interpro_id}|{pfam_id}|{dstart}|{dend}"),
                    "edge_id": eid,
                    "protein_id": pid,
                    "domain_class": dclass,
                    "interpro_id": interpro_id,
                    "pfam_id": pfam_id,
                    "domain_name": domain_name,
                    "domain_start": dstart,
                    "domain_end": dend,
                    "source": "InterPro/Pfam",
                    "source_version": normalize(d.get("source_version", "")) or "unknown",
                    "fetch_date": args.fetch_date,
                }
            )

    function_rows: List[Dict[str, str]] = []
    function_relation_counter = Counter()

    for eid in target_edge_ids:
        edge = edges[eid]
        pid = normalize(edge.get("dst_id", ""))
        ann = protein_annotations.get(pid, {"source": "UniProtKB", "source_version": "protein_master_v6_clean", "text": ""})
        rels = infer_function_relations(ann.get("text", ""))
        for rel, basis, snippet in rels:
            function_relation_counter[rel] += 1
            function_rows.append(
                {
                    "context_id": make_id("rpi_function", f"{eid}|{pid}|{rel}|{basis}"),
                    "edge_id": eid,
                    "protein_id": pid,
                    "function_relation": rel,
                    "inference_basis": basis,
                    "evidence_snippet": snippet,
                    "source": ann.get("source", "") or "UniProtKB",
                    "source_version": ann.get("source_version", "") or "protein_master_v6_clean",
                    "fetch_date": args.fetch_date,
                }
            )

    # sort deterministic
    site_rows = sorted(site_rows, key=lambda r: r["context_id"])
    domain_rows = sorted(domain_rows, key=lambda r: r["context_id"])
    function_rows = sorted(function_rows, key=lambda r: r["context_id"])

    write_tsv(args.site_output, site_rows, SITE_COLUMNS)
    write_tsv(args.domain_output, domain_rows, DOMAIN_COLUMNS)
    write_tsv(args.function_output, function_rows, FUNCTION_COLUMNS)

    total_edges = len(target_edge_ids)
    site_edges_present = {r["edge_id"] for r in site_rows}
    domain_edges_present = {r["edge_id"] for r in domain_rows}
    union_edges = site_edge_with_coord | domain_edges_present

    site_method_non_empty = sum(1 for r in site_rows if normalize(r.get("method_type", "")) != "")
    site_method_rate = rate(site_method_non_empty, len(site_rows))

    site_edge_coverage = rate(len(site_edge_with_coord), total_edges)
    domain_edge_coverage = rate(len(domain_edges_present), total_edges)
    site_or_domain_coverage = rate(len(union_edges), total_edges)

    def join_rate(rows: List[Dict[str, str]]) -> float:
        if not rows:
            return 0.0
        matched = sum(1 for r in rows if normalize(r.get("edge_id", "")) in target_edge_set)
        return matched / len(rows)

    site_join_rate = join_rate(site_rows)
    domain_join_rate = join_rate(domain_rows)
    function_join_rate = join_rate(function_rows)
    edge_join_rate = min(site_join_rate, domain_join_rate, function_join_rate) if function_rows else 0.0

    gates = {
        "edge_id_join_rate_ge_0_99": edge_join_rate >= 0.99,
        "method_field_coverage_ge_0_90": site_method_rate >= 0.90,
        "site_or_domain_coverage_ge_0_60": site_or_domain_coverage >= 0.60,
    }

    metrics = {
        **base,
        "status": "completed",
        "row_count": {
            "target_edges": total_edges,
            "evidence_rows_scanned": len(evidence_rows),
            "site_rows": len(site_rows),
            "domain_rows": len(domain_rows),
            "function_rows": len(function_rows),
        },
        "output": {
            "site_context": str(args.site_output),
            "domain_context": str(args.domain_output),
            "function_context": str(args.function_output),
        },
        "coverage": {
            "edge_id_join_rate": edge_join_rate,
            "site_join_rate": site_join_rate,
            "domain_join_rate": domain_join_rate,
            "function_join_rate": function_join_rate,
            "method_field_coverage": site_method_rate,
            "site_edge_coverage": site_edge_coverage,
            "domain_edge_coverage": domain_edge_coverage,
            "site_or_domain_coverage": site_or_domain_coverage,
        },
        "distribution": {
            "site_method_type": dict(sorted(Counter(r["method_type"] for r in site_rows).items())),
            "domain_class": dict(sorted(domain_class_counter.items())),
            "function_relation": dict(sorted(function_relation_counter.items())),
        },
        "lineage": {
            "raw_snapshots": snapshot_info,
            "matched_reference_count": matched_reference,
            "unmatched_reference_top20": [
                {"reference": k, "count": v} for k, v in unmatched_reference.most_common(20)
            ],
        },
        "gates": {
            "passed": all(gates.values()),
            "checks": gates,
        },
    }

    write_json(args.report, metrics)

    if args.gates_report is not None:
        write_json(
            args.gates_report,
            {
                "pipeline": "rpi_site_domain_enrichment_v2",
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "status": "PASS" if metrics["gates"]["passed"] else "FAIL",
                "coverage": metrics["coverage"],
                "checks": metrics["gates"]["checks"],
            },
        )

    print(
        "[OK] "
        f"edges={total_edges} site_rows={len(site_rows)} domain_rows={len(domain_rows)} function_rows={len(function_rows)} "
        f"join={edge_join_rate:.4f} method={site_method_rate:.4f} site_or_domain={site_or_domain_coverage:.4f}"
    )
    print(f"[OK] report -> {args.report}")
    if args.gates_report is not None:
        print(f"[OK] gates -> {args.gates_report}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
