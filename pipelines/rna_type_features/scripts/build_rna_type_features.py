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

LNC_COLUMNS = [
    "rna_id",
    "urs_id",
    "taxon_id",
    "ensembl_transcript_id",
    "ensembl_gene_id",
    "gene_symbol",
    "rna_name",
    "chrom",
    "locus_start",
    "locus_end",
    "strand",
    "biotype",
    "mapping_status",
    "source",
    "source_version",
    "fetch_date",
]

TRNA_COLUMNS = [
    "rna_id",
    "urs_id",
    "taxon_id",
    "trna_label",
    "amino_acid",
    "anticodon",
    "anticodon_rna",
    "anticodon_source",
    "seq_accession",
    "locus_start",
    "locus_end",
    "locus_type",
    "evidence_dbs",
    "source",
    "source_version",
    "fetch_date",
]

RRNA_COLUMNS = [
    "rrna_locus_id",
    "rna_id",
    "urs_id",
    "taxon_id",
    "rrna_type",
    "seq_accession",
    "locus_start",
    "locus_end",
    "strand",
    "source_db",
    "external_id",
    "gene_symbol",
    "source",
    "source_version",
    "fetch_date",
]


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


def parse_gtf_attributes(attr_text: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for key, value in re.findall(r'([A-Za-z0-9_]+)\s+"([^"]*)";', attr_text):
        out[key] = value
    return out


_ENA_LOCUS_RE = re.compile(r"^([^:]+):(\d+)\.\.(\d+):([^:]+)$")
_ENA_COMP_LOCUS_RE = re.compile(r"^([^:]+):complement\((\d+)\.\.(\d+)\):([^:]+)$")
_GTRNADB_LOCUS_RE = re.compile(r"^(?:GTRNADB:)?[^:]+:([^:]+):(\d+)-(\d+)$")


def parse_external_locus(external_id: str) -> Optional[Dict[str, str]]:
    text = normalize(external_id)
    if text == "":
        return None

    m = _ENA_LOCUS_RE.match(text)
    if m:
        return {
            "seq_accession": m.group(1),
            "start": m.group(2),
            "end": m.group(3),
            "strand": ".",
            "feature": m.group(4),
        }

    m = _ENA_COMP_LOCUS_RE.match(text)
    if m:
        return {
            "seq_accession": m.group(1),
            "start": m.group(2),
            "end": m.group(3),
            "strand": "-",
            "feature": m.group(4),
        }

    m = _GTRNADB_LOCUS_RE.match(text)
    if m:
        return {
            "seq_accession": m.group(1),
            "start": m.group(2),
            "end": m.group(3),
            "strand": ".",
            "feature": "tRNA",
        }

    return None


_TRNA_ANTICODON_RE = re.compile(r"(tRNA-[A-Za-z]+-([ACGTU]{3})(?:-[0-9]+(?:-[0-9]+)?)?)")
_TRNA_AA_RE = re.compile(r"tRNA-([A-Za-z]+)-([ACGTU]{3})")
_TRNA_TRNAI_RE = re.compile(r"TRNA[A-Z]-([ACGTU]{3})", re.IGNORECASE)


def parse_trna_annotation(external_id: str, symbol_or_aux: str) -> Tuple[str, str, str, str]:
    aa = ""
    anticodon = ""
    label = ""
    for text in [normalize(external_id), normalize(symbol_or_aux)]:
        if text == "":
            continue

        aa_m = _TRNA_AA_RE.search(text)
        if aa_m:
            aa = aa or aa_m.group(1)
            anticodon = anticodon or aa_m.group(2).upper()

        label_m = _TRNA_ANTICODON_RE.search(text)
        if label_m:
            label = label or label_m.group(1)
            anticodon = anticodon or label_m.group(2).upper()

        if anticodon == "":
            tri = _TRNA_TRNAI_RE.search(text)
            if tri:
                anticodon = tri.group(1).upper()

    anticodon_rna = anticodon.replace("T", "U") if anticodon else ""
    return aa, anticodon, anticodon_rna, label


def choose_best_urs(counter_like: Dict[str, int]) -> str:
    if not counter_like:
        return ""
    return sorted(counter_like.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]


def infer_rrna_class(symbol: str, external_id: str, feature: str) -> str:
    text = " ".join([normalize(symbol), normalize(external_id), normalize(feature)]).upper()
    if "5.8" in text or "5-8" in text:
        return "5.8S_rRNA"
    if "18S" in text:
        return "18S_rRNA"
    if "28S" in text:
        return "28S_rRNA"
    if "16S" in text:
        return "16S_rRNA"
    if "5S" in text:
        return "5S_rRNA"
    return "rRNA"


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
    ok = sum(1 for r in rows if normalize(r.get(col, "")) != "")
    return ok / len(rows)


def read_master_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None or "rna_id" not in set(r.fieldnames):
            raise SystemExit("[ERROR] rna_master_v1.tsv missing rna_id column")
        for row in r:
            rid = normalize(row.get("rna_id", ""))
            if rid:
                ids.add(rid)
    return ids


def scan_id_mapping(
    id_mapping_path: Path,
    taxon_id: str,
    max_lines: Optional[int] = None,
) -> Tuple[Dict[str, Counter], Dict[str, Counter], Dict[str, Dict], List[Dict[str, str]], Dict]:
    enst_version_to_urs: Dict[str, Counter] = defaultdict(Counter)
    enst_base_to_urs: Dict[str, Counter] = defaultdict(Counter)
    trna_by_urs: Dict[str, Dict] = {}
    rrna_loci: List[Dict[str, str]] = []
    rrna_seen = set()

    stats = {
        "idmap_lines_scanned": 0,
        "taxon_pass": 0,
        "lnc_enst_links": 0,
        "trna_rows": 0,
        "rrna_rows": 0,
        "rrna_with_locus": 0,
    }

    with open_maybe_gz(id_mapping_path) as f:
        for line in f:
            stats["idmap_lines_scanned"] += 1
            if max_lines is not None and stats["idmap_lines_scanned"] > max_lines:
                break
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue

            urs_raw = normalize(parts[0])
            db = normalize(parts[1]).upper()
            external_id = normalize(parts[2])
            tx = normalize(parts[3])
            rna_type = normalize(parts[4])
            aux = normalize(parts[5]) if len(parts) > 5 else ""

            if tx != taxon_id:
                continue
            stats["taxon_pass"] += 1

            urs = urs_raw if urs_raw.endswith(f"_{taxon_id}") else f"{urs_raw}_{taxon_id}"

            if rna_type.lower() == "lncrna" and db in {"ENSEMBL", "ENSEMBL_GENCODE"} and external_id.startswith("ENST"):
                enst_version_to_urs[external_id][urs] += 1
                enst_base_to_urs[external_id.split(".")[0]][urs] += 1
                stats["lnc_enst_links"] += 1

            if rna_type in {"tRNA", "Mt_tRNA"}:
                stats["trna_rows"] += 1
                rec = trna_by_urs.get(urs)
                if rec is None:
                    rec = {
                        "rna_id": urs,
                        "urs_id": urs,
                        "taxon_id": taxon_id,
                        "trna_label": "",
                        "amino_acid": "",
                        "anticodon": "",
                        "anticodon_rna": "",
                        "anticodon_source": "none",
                        "seq_accession": "",
                        "locus_start": "",
                        "locus_end": "",
                        "locus_type": rna_type if rna_type == "Mt_tRNA" else "tRNA",
                        "evidence_dbs": set(),
                        "source": "RNAcentral",
                        "_locus_priority": -1,
                    }
                    trna_by_urs[urs] = rec

                rec["evidence_dbs"].add(db)
                if rna_type == "Mt_tRNA":
                    rec["locus_type"] = "Mt_tRNA"

                aa, anti, anti_rna, label = parse_trna_annotation(external_id, aux)
                if label and not rec["trna_label"]:
                    rec["trna_label"] = label
                if aa and not rec["amino_acid"]:
                    rec["amino_acid"] = aa
                if anti and not rec["anticodon"]:
                    rec["anticodon"] = anti
                    rec["anticodon_rna"] = anti_rna
                    rec["anticodon_source"] = "gtrnadb" if db == "GTRNADB" else "symbol_inference"

                if not rec["trna_label"] and aux and not aux.startswith("ENSG"):
                    rec["trna_label"] = aux

                locus = parse_external_locus(external_id)
                if locus:
                    priority = 2 if db == "GTRNADB" else 1
                    if priority > rec["_locus_priority"]:
                        rec["seq_accession"] = locus["seq_accession"]
                        rec["locus_start"] = locus["start"]
                        rec["locus_end"] = locus["end"]
                        rec["_locus_priority"] = priority

            if rna_type == "rRNA":
                stats["rrna_rows"] += 1
                locus = parse_external_locus(external_id)
                if not locus:
                    continue
                stats["rrna_with_locus"] += 1

                key = (urs, db, locus["seq_accession"], locus["start"], locus["end"], locus["strand"])
                if key in rrna_seen:
                    continue
                rrna_seen.add(key)

                raw_id = "|".join([*key, external_id])
                rrna_locus_id = "RRL" + hashlib.sha1(raw_id.encode("utf-8")).hexdigest()[:16].upper()
                rrna_loci.append(
                    {
                        "rrna_locus_id": rrna_locus_id,
                        "rna_id": urs,
                        "urs_id": urs,
                        "taxon_id": taxon_id,
                        "rrna_type": infer_rrna_class(aux, external_id, locus["feature"]),
                        "seq_accession": locus["seq_accession"],
                        "locus_start": locus["start"],
                        "locus_end": locus["end"],
                        "strand": locus["strand"],
                        "source_db": db,
                        "external_id": external_id,
                        "gene_symbol": aux,
                        "source": "RNAcentral",
                    }
                )

    return enst_version_to_urs, enst_base_to_urs, trna_by_urs, rrna_loci, stats


def build_lnc_rows(
    gtf_path: Path,
    taxon_id: str,
    enst_version_to_urs: Dict[str, Counter],
    enst_base_to_urs: Dict[str, Counter],
    max_lines: Optional[int] = None,
) -> Tuple[List[Dict[str, str]], Dict]:
    rows: List[Dict[str, str]] = []
    seen_enst = set()
    stats = {
        "gtf_lines_scanned": 0,
        "gtf_transcript_rows": 0,
        "gtf_lnc_transcripts": 0,
        "mapped_enst_version": 0,
        "mapped_enst_base": 0,
        "fallback_enst": 0,
    }

    with open_maybe_gz(gtf_path) as f:
        for line in f:
            stats["gtf_lines_scanned"] += 1
            if max_lines is not None and stats["gtf_lines_scanned"] > max_lines:
                break
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "transcript":
                continue
            stats["gtf_transcript_rows"] += 1

            attrs = parse_gtf_attributes(parts[8])
            biotype = normalize(attrs.get("transcript_biotype") or attrs.get("gene_biotype"))
            if biotype.lower() != "lncrna":
                continue
            stats["gtf_lnc_transcripts"] += 1

            enst_version = normalize(attrs.get("transcript_id"))
            if not enst_version.startswith("ENST"):
                continue
            enst_base = enst_version.split(".")[0]
            if enst_base in seen_enst:
                continue
            seen_enst.add(enst_base)

            urs = ""
            mapping_status = "fallback_enst"
            c1 = enst_version_to_urs.get(enst_version)
            c2 = enst_base_to_urs.get(enst_base)
            if c1:
                urs = choose_best_urs(c1)
                mapping_status = "mapped_enst_version"
                stats["mapped_enst_version"] += 1
            elif c2:
                urs = choose_best_urs(c2)
                mapping_status = "mapped_enst_base"
                stats["mapped_enst_base"] += 1
            else:
                stats["fallback_enst"] += 1

            rna_id = urs if urs else f"{enst_base}_{taxon_id}"
            rows.append(
                {
                    "rna_id": rna_id,
                    "urs_id": urs,
                    "taxon_id": taxon_id,
                    "ensembl_transcript_id": enst_version,
                    "ensembl_gene_id": normalize(attrs.get("gene_id")),
                    "gene_symbol": normalize(attrs.get("gene_name")),
                    "rna_name": normalize(attrs.get("transcript_name")) or normalize(attrs.get("gene_name")) or enst_base,
                    "chrom": normalize(parts[0]),
                    "locus_start": normalize(parts[3]),
                    "locus_end": normalize(parts[4]),
                    "strand": normalize(parts[6]),
                    "biotype": "lncRNA",
                    "mapping_status": mapping_status,
                    "source": "Ensembl;RNAcentral",
                }
            )

    rows.sort(key=lambda r: (r["chrom"], int(r["locus_start"]), r["ensembl_transcript_id"]))
    return rows, stats


def finalize_trna_rows(trna_by_urs: Dict[str, Dict]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for urs, rec in trna_by_urs.items():
        dbs = sorted(rec.get("evidence_dbs", set()))
        rows.append(
            {
                "rna_id": rec["rna_id"],
                "urs_id": rec["urs_id"],
                "taxon_id": rec["taxon_id"],
                "trna_label": rec["trna_label"] or urs,
                "amino_acid": rec["amino_acid"],
                "anticodon": rec["anticodon"],
                "anticodon_rna": rec["anticodon_rna"],
                "anticodon_source": rec["anticodon_source"],
                "seq_accession": rec["seq_accession"],
                "locus_start": rec["locus_start"],
                "locus_end": rec["locus_end"],
                "locus_type": rec["locus_type"],
                "evidence_dbs": ";".join(dbs),
                "source": rec["source"],
            }
        )
    rows.sort(key=lambda r: r["rna_id"])
    return rows


def build_metrics(
    lnc_rows: List[Dict[str, str]],
    trna_rows: List[Dict[str, str]],
    rrna_rows: List[Dict[str, str]],
    master_ids: set[str],
    idmap_stats: Dict,
    lnc_stats: Dict,
    inputs: Dict[str, str],
    outputs: Dict[str, str],
) -> Dict:
    lnc_map_dist = Counter(r["mapping_status"] for r in lnc_rows)
    trna_anti_source = Counter(r["anticodon_source"] for r in trna_rows)
    trna_locus_type = Counter(r["locus_type"] for r in trna_rows)
    rrna_type_dist = Counter(r["rrna_type"] for r in rrna_rows)
    rrna_db_dist = Counter(r["source_db"] for r in rrna_rows)

    def overlap_rate(rows: List[Dict[str, str]]) -> Dict[str, float]:
        if not rows:
            return {"matched": 0, "total": 0, "rate": 0.0}
        matched = sum(1 for r in rows if normalize(r.get("rna_id", "")) in master_ids)
        total = len(rows)
        return {"matched": matched, "total": total, "rate": matched / total}

    return {
        "pipeline": "rna_type_features_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": inputs,
        "outputs": outputs,
        "row_count": {
            "rna_lnc_entries_v1": len(lnc_rows),
            "rna_trna_features_v1": len(trna_rows),
            "rna_rrna_loci_v1": len(rrna_rows),
        },
        "acceptance_metrics": {
            "lnc_traceable_rate_ensembl_transcript_id": column_non_empty_rate(lnc_rows, "ensembl_transcript_id"),
            "trna_traceable_rate_rna_id": column_non_empty_rate(trna_rows, "rna_id"),
            "rrna_traceable_rate_locus_id": column_non_empty_rate(rrna_rows, "rrna_locus_id"),
            "lnc_master_overlap": overlap_rate(lnc_rows),
            "trna_master_overlap": overlap_rate(trna_rows),
            "rrna_master_overlap": overlap_rate(rrna_rows),
            "trna_anticodon_non_empty_rate": column_non_empty_rate(trna_rows, "anticodon"),
        },
        "non_empty_rates": {
            "lnc": {
                "rna_id": column_non_empty_rate(lnc_rows, "rna_id"),
                "ensembl_transcript_id": column_non_empty_rate(lnc_rows, "ensembl_transcript_id"),
                "ensembl_gene_id": column_non_empty_rate(lnc_rows, "ensembl_gene_id"),
                "gene_symbol": column_non_empty_rate(lnc_rows, "gene_symbol"),
            },
            "trna": {
                "rna_id": column_non_empty_rate(trna_rows, "rna_id"),
                "trna_label": column_non_empty_rate(trna_rows, "trna_label"),
                "anticodon": column_non_empty_rate(trna_rows, "anticodon"),
                "seq_accession": column_non_empty_rate(trna_rows, "seq_accession"),
            },
            "rrna": {
                "rrna_locus_id": column_non_empty_rate(rrna_rows, "rrna_locus_id"),
                "rna_id": column_non_empty_rate(rrna_rows, "rna_id"),
                "seq_accession": column_non_empty_rate(rrna_rows, "seq_accession"),
                "locus_start": column_non_empty_rate(rrna_rows, "locus_start"),
                "locus_end": column_non_empty_rate(rrna_rows, "locus_end"),
            },
        },
        "distributions": {
            "lnc_mapping_status": dict(lnc_map_dist),
            "trna_anticodon_source": dict(trna_anti_source),
            "trna_locus_type": dict(trna_locus_type),
            "rrna_type": dict(rrna_type_dist),
            "rrna_source_db": dict(rrna_db_dist),
        },
        "scan_stats": {
            "id_mapping": idmap_stats,
            "gtf_lnc": lnc_stats,
        },
    }


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Build RNA type feature pack (lnc/tRNA/rRNA)")
    ap.add_argument("--rna-master", type=Path, default=Path("data/output/rna_master_v1.tsv"))
    ap.add_argument("--gtf", type=Path, default=Path("data/raw/rna/ensembl/Homo_sapiens.GRCh38.115.chr.gtf.gz"))
    ap.add_argument("--id-mapping", type=Path, default=Path("data/raw/rna/rnacentral/id_mapping.tsv.gz"))
    ap.add_argument("--taxon-id", default="9606")
    ap.add_argument("--lnc-output", type=Path, default=Path("data/output/rna_lnc_entries_v1.tsv"))
    ap.add_argument("--trna-output", type=Path, default=Path("data/output/rna_trna_features_v1.tsv"))
    ap.add_argument("--rrna-output", type=Path, default=Path("data/output/rna_rrna_loci_v1.tsv"))
    ap.add_argument("--report", type=Path, default=Path("pipelines/rna_type_features/reports/rna_type_features_v1.metrics.json"))
    ap.add_argument("--source-version", default="RNAcentral:id_mapping.tsv;Ensembl:GRCh38.115.chr.gtf")
    ap.add_argument("--fetch-date", default=FETCH_DATE)
    ap.add_argument("--max-idmap-lines", type=int, default=None)
    ap.add_argument("--max-gtf-lines", type=int, default=None)
    ap.add_argument("--check-inputs", action="store_true")
    return ap.parse_args()


def main() -> int:
    args = parse_args()

    gtf_path = resolve_with_alt(args.gtf)
    idmap_path = resolve_with_alt(args.id_mapping)

    missing = []
    if not args.rna_master.exists():
        missing.append(str(args.rna_master))
    if gtf_path is None:
        missing.append(str(args.gtf))
    if idmap_path is None:
        missing.append(str(args.id_mapping))

    if missing:
        write_json(
            args.report,
            {
                "pipeline": "rna_type_features_v1",
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
                "pipeline": "rna_type_features_v1",
                "status": "inputs_ready",
                "inputs": {
                    "rna_master": str(args.rna_master),
                    "gtf": str(gtf_path),
                    "id_mapping": str(idmap_path),
                },
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            },
        )
        print(f"[OK] inputs ready -> {args.report}")
        return 0

    master_ids = read_master_ids(args.rna_master)
    enst_version_to_urs, enst_base_to_urs, trna_by_urs, rrna_rows, idmap_stats = scan_id_mapping(
        id_mapping_path=idmap_path, taxon_id=args.taxon_id, max_lines=args.max_idmap_lines
    )
    lnc_rows, lnc_stats = build_lnc_rows(
        gtf_path=gtf_path,
        taxon_id=args.taxon_id,
        enst_version_to_urs=enst_version_to_urs,
        enst_base_to_urs=enst_base_to_urs,
        max_lines=args.max_gtf_lines,
    )
    trna_rows = finalize_trna_rows(trna_by_urs)
    rrna_rows.sort(key=lambda r: (r["seq_accession"], int(r["locus_start"]), r["rrna_locus_id"]))

    write_tsv(args.lnc_output, lnc_rows, LNC_COLUMNS, args.source_version, args.fetch_date)
    write_tsv(args.trna_output, trna_rows, TRNA_COLUMNS, args.source_version, args.fetch_date)
    write_tsv(args.rrna_output, rrna_rows, RRNA_COLUMNS, args.source_version, args.fetch_date)

    metrics = build_metrics(
        lnc_rows=lnc_rows,
        trna_rows=trna_rows,
        rrna_rows=rrna_rows,
        master_ids=master_ids,
        idmap_stats=idmap_stats,
        lnc_stats=lnc_stats,
        inputs={
            "rna_master": str(args.rna_master),
            "gtf": str(gtf_path),
            "id_mapping": str(idmap_path),
        },
        outputs={
            "rna_lnc_entries_v1": str(args.lnc_output),
            "rna_trna_features_v1": str(args.trna_output),
            "rna_rrna_loci_v1": str(args.rrna_output),
        },
    )
    write_json(args.report, metrics)

    print(
        "[OK] built type features: "
        f"lnc={len(lnc_rows)} trna={len(trna_rows)} rrna_loci={len(rrna_rows)}"
    )
    print(f"[OK] metrics -> {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
