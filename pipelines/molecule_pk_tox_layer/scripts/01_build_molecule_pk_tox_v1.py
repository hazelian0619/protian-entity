#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from lxml import etree

INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")

FIELDS = [
    "bioavailability",
    "clearance",
    "half_life",
    "volume_of_distribution",
    "ld50",
    "mutagenicity",
]
EVIDENCE_TIERS = {"structured/high", "text_mined/medium", "inferred/low", "none"}
EVIDENCE_RANK = {"none": 0, "inferred/low": 1, "text_mined/medium": 2, "structured/high": 3}


@dataclass
class XrefRec:
    inchikey: str
    drugbank_ids: List[str]
    pubchem_cids: List[str]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def utc_today() -> str:
    return datetime.now(timezone.utc).date().isoformat()


def normalize(x: str) -> str:
    return (x or "").strip()


def split_multi(x: str) -> List[str]:
    t = normalize(x)
    if not t:
        return []
    return [i.strip() for i in t.split(";") if i.strip()]


def valid_inchikey(x: str) -> bool:
    return bool(INCHIKEY_RE.match(normalize(x).upper()))


def normalize_ws(x: str) -> str:
    return re.sub(r"\s+", " ", normalize(x))


def strip_html(x: str) -> str:
    t = re.sub(r"<[^>]+>", "", x or "")
    t = t.replace("&nbsp;", " ").replace("&lt;", "<").replace("&gt;", ">").replace("&amp;", "&")
    return normalize_ws(t)


def first_sentence(x: str, keyword: Optional[str] = None, max_len: int = 320) -> str:
    t = strip_html(x)
    if not t:
        return ""
    sentences = re.split(r"(?<=[\.\!\?])\s+", t)
    if keyword:
        kw = keyword.lower()
        for s in sentences:
            if kw in s.lower():
                return s[:max_len]
    return sentences[0][:max_len]


def write_tsv(path: Path, rows: Iterable[Dict[str, str]], header: Sequence[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    n = 0
    with tmp.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(header), lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in header})
            n += 1
    tmp.replace(path)
    return n


def read_tsv(path: Path, max_rows: Optional[int] = None) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        rows: List[Dict[str, str]] = []
        for i, row in enumerate(r, start=1):
            rows.append(row)
            if max_rows is not None and i >= max_rows:
                break
    return list(r.fieldnames), rows


def join_values(values: Iterable[str]) -> str:
    uniq = sorted({normalize(v) for v in values if normalize(v)})
    return ";".join(uniq)


def detect_drugbank_version(xml_path: Path) -> str:
    try:
        with xml_path.open("r", encoding="utf-8", errors="ignore") as f:
            head = f.read(4096)
    except Exception:
        return "DrugBank:unknown"

    m_root = re.search(r"<drugbank\b[^>]*>", head, flags=re.IGNORECASE | re.DOTALL)
    if not m_root:
        return "DrugBank:unknown"

    root_tag = m_root.group(0)
    m_ver = re.search(r'\bversion="([^"]+)"', root_tag)
    m_exp = re.search(r'\bexported-on="([^"]+)"', root_tag)
    if m_ver and m_exp:
        return f"DrugBank:{m_ver.group(1)}|exported-on:{m_exp.group(1)}"
    if m_ver:
        return f"DrugBank:{m_ver.group(1)}"
    return "DrugBank:unknown"


def choose_col(fieldnames: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    low = {c.lower(): c for c in fieldnames}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    return None


def load_xref(path: Path, max_rows: Optional[int]) -> Tuple[List[XrefRec], Dict[str, Any], Set[str], Set[str]]:
    cols, rows = read_tsv(path, max_rows=max_rows)
    if "inchikey" not in cols:
        raise SystemExit(f"[ERROR] xref missing inchikey: {path}")

    out: List[XrefRec] = []
    bad_inchikey = 0
    all_dbids: Set[str] = set()
    all_cids: Set[str] = set()
    rows_with_dbid = 0
    rows_with_cid = 0

    for row in rows:
        ik = normalize(row.get("inchikey", "")).upper()
        if not valid_inchikey(ik):
            bad_inchikey += 1
            continue

        dbids = [d for d in split_multi(row.get("drugbank_id", "")) if d]
        cids = [c for c in split_multi(row.get("pubchem_cid", "")) if c]
        if dbids:
            rows_with_dbid += 1
            all_dbids.update(dbids)
        if cids:
            rows_with_cid += 1
            all_cids.update(cids)
        out.append(XrefRec(inchikey=ik, drugbank_ids=dbids, pubchem_cids=cids))

    stats = {
        "input_rows": len(rows),
        "valid_rows": len(out),
        "skipped_bad_inchikey_rows": bad_inchikey,
        "rows_with_drugbank_id": rows_with_dbid,
        "rows_with_pubchem_cid": rows_with_cid,
        "unique_drugbank_ids": len(all_dbids),
        "unique_pubchem_cids": len(all_cids),
    }
    return out, stats, all_dbids, all_cids


def get_direct_child_text_by_local(node: etree._Element, local_name: str) -> str:
    for ch in node:
        if etree.QName(ch).localname == local_name:
            return normalize_ws("".join(ch.itertext()))
    return ""


def get_primary_drugbank_id(drug: etree._Element) -> str:
    first = ""
    primary = ""
    for ch in drug:
        if etree.QName(ch).localname != "drugbank-id":
            continue
        text = normalize_ws("".join(ch.itertext()))
        if not text:
            continue
        if not first:
            first = text
        if ch.get("primary", "").lower() == "true":
            primary = text
            break
    return primary or first


def extract_ld50(toxicity_text: str) -> str:
    t = strip_html(toxicity_text)
    if not t:
        return ""
    patterns = [
        r"([^.;\n]{0,100}\bLD\s*50\b[^.;\n]{0,220})",
        r"([^.;\n]{0,100}\blethal dose\b[^.;\n]{0,220})",
    ]
    for p in patterns:
        m = re.search(p, t, flags=re.IGNORECASE)
        if m:
            return normalize_ws(m.group(1))[:320]
    return ""


def infer_mutagenicity(toxicity_text: str) -> str:
    t = strip_html(toxicity_text).lower()
    if not t:
        return ""
    neg_markers = [
        "non-mutagenic",
        "non mutagenic",
        "not mutagenic",
        "no mutagenic",
        "ames test negative",
        "negative ames",
        "not genotoxic",
        "no genotoxicity",
    ]
    pos_markers = [
        "mutagenic",
        "genotoxic",
        "positive ames",
        "ames positive",
    ]
    neg = any(k in t for k in neg_markers)
    pos = any(k in t for k in pos_markers)
    if pos and neg:
        return "mixed"
    if pos:
        return "positive"
    if neg:
        return "negative"
    return ""


def load_drugbank_pk_tox(
    xml_path: Optional[Path], target_dbids: Set[str]
) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Any], str]:
    if xml_path is None or not xml_path.exists():
        return {}, {"available": False, "reason": "drugbank_xml_not_found"}, "DrugBank:unavailable"

    out: Dict[str, Dict[str, str]] = {}
    parsed_top_level = 0
    matched_target = 0
    value_counter = Counter()
    version = detect_drugbank_version(xml_path)

    for _, drug in etree.iterparse(
        str(xml_path),
        events=("end",),
        tag="{*}drug",
        huge_tree=True,
        recover=True,
    ):
        parent = drug.getparent()
        if parent is None or etree.QName(parent).localname != "drugbank":
            drug.clear()
            continue

        parsed_top_level += 1
        dbid = get_primary_drugbank_id(drug)
        if not dbid or (target_dbids and dbid not in target_dbids):
            drug.clear()
            while drug.getprevious() is not None:
                del parent[0]
            continue

        matched_target += 1
        absorption = get_direct_child_text_by_local(drug, "absorption")
        clearance = get_direct_child_text_by_local(drug, "clearance")
        half_life = get_direct_child_text_by_local(drug, "half-life")
        volume = get_direct_child_text_by_local(drug, "volume-of-distribution")
        toxicity = get_direct_child_text_by_local(drug, "toxicity")

        rec = {
            "bioavailability": first_sentence(absorption, keyword="bioavailability"),
            "bioavailability_evidence": "none",
            "bioavailability_source": "none",
            "clearance": first_sentence(clearance),
            "clearance_evidence": "none",
            "clearance_source": "none",
            "half_life": first_sentence(half_life),
            "half_life_evidence": "none",
            "half_life_source": "none",
            "volume_of_distribution": first_sentence(volume),
            "volume_of_distribution_evidence": "none",
            "volume_of_distribution_source": "none",
            "ld50": extract_ld50(toxicity),
            "ld50_evidence": "none",
            "ld50_source": "none",
            "mutagenicity": infer_mutagenicity(toxicity),
            "mutagenicity_evidence": "none",
            "mutagenicity_source": "none",
        }

        if rec["bioavailability"]:
            rec["bioavailability_evidence"] = "structured/high"
            rec["bioavailability_source"] = f"drugbank_xml:{dbid}#absorption"
            value_counter["bioavailability"] += 1
        if rec["clearance"]:
            rec["clearance_evidence"] = "structured/high"
            rec["clearance_source"] = f"drugbank_xml:{dbid}#clearance"
            value_counter["clearance"] += 1
        if rec["half_life"]:
            rec["half_life_evidence"] = "structured/high"
            rec["half_life_source"] = f"drugbank_xml:{dbid}#half-life"
            value_counter["half_life"] += 1
        if rec["volume_of_distribution"]:
            rec["volume_of_distribution_evidence"] = "structured/high"
            rec["volume_of_distribution_source"] = f"drugbank_xml:{dbid}#volume-of-distribution"
            value_counter["volume_of_distribution"] += 1
        if rec["ld50"]:
            rec["ld50_evidence"] = "text_mined/medium"
            rec["ld50_source"] = f"drugbank_xml:{dbid}#toxicity"
            value_counter["ld50"] += 1
        if rec["mutagenicity"]:
            rec["mutagenicity_evidence"] = "inferred/low"
            rec["mutagenicity_source"] = f"drugbank_xml:{dbid}#toxicity"
            value_counter["mutagenicity"] += 1

        out[dbid] = rec

        drug.clear()
        while drug.getprevious() is not None:
            del parent[0]

    report = {
        "available": True,
        "xml_path": str(xml_path),
        "parsed_top_level_drugs": parsed_top_level,
        "matched_target_drugbank_ids": matched_target,
        "target_records_loaded": len(out),
        "field_non_empty_counts": {k: int(value_counter.get(k, 0)) for k in FIELDS},
    }
    return out, report, version


def load_pubchem_tox(path: Optional[Path]) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Any], str]:
    if path is None or not path.exists():
        return {}, {"available": False, "reason": "pubchem_tox_file_not_found"}, "PubChem:unavailable"

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        head = f.readline()
    delimiter = "\t" if "\t" in head else ","

    with path.open("r", encoding="utf-8", newline="", errors="ignore") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        if reader.fieldnames is None:
            return {}, {"available": False, "reason": "pubchem_tox_missing_header"}, "PubChem:unknown"
        cols = list(reader.fieldnames)
        cid_col = choose_col(cols, ["pubchem_cid", "cid", "pubchemcid"])
        if not cid_col:
            return {}, {"available": False, "reason": "pubchem_tox_missing_cid_column"}, "PubChem:unknown"

        col_map = {
            "bioavailability": choose_col(cols, ["bioavailability", "oral_bioavailability"]),
            "clearance": choose_col(cols, ["clearance", "total_clearance"]),
            "half_life": choose_col(cols, ["half_life", "half-life", "t_half"]),
            "volume_of_distribution": choose_col(cols, ["volume_of_distribution", "volume-distribution", "vd"]),
            "ld50": choose_col(cols, ["ld50", "ld_50", "acute_toxicity_ld50"]),
            "mutagenicity": choose_col(cols, ["mutagenicity", "ames", "ames_mutagenicity"]),
        }
        source_ref_col = choose_col(cols, ["source", "source_url", "reference", "url"])
        version_col = choose_col(cols, ["source_version", "version", "dataset_version"])

        rows = list(reader)

    out: Dict[str, Dict[str, str]] = {}
    mapped_rows = 0
    field_non_empty = Counter()
    versions: Set[str] = set()

    for idx, row in enumerate(rows, start=1):
        cid = normalize(row.get(cid_col, ""))
        if not cid:
            continue

        rec: Dict[str, str] = {}
        score = 0
        for field in FIELDS:
            col = col_map.get(field)
            val = normalize(row.get(col, "")) if col else ""
            if val:
                score += 1
                rec[field] = val
                rec[f"{field}_evidence"] = "structured/high"
                ref = normalize(row.get(source_ref_col, "")) if source_ref_col else ""
                if ref:
                    rec[f"{field}_source"] = f"pubchem_tox:{cid}#{field}|{ref}"
                else:
                    rec[f"{field}_source"] = f"pubchem_tox:{cid}#{field}|row={idx}"
                field_non_empty[field] += 1
            else:
                rec[field] = ""
                rec[f"{field}_evidence"] = "none"
                rec[f"{field}_source"] = "none"

        if score == 0:
            continue

        existing = out.get(cid)
        existing_score = sum(1 for f in FIELDS if existing and existing.get(f))
        if existing is None or score > existing_score:
            out[cid] = rec
            mapped_rows += 1

        if version_col:
            v = normalize(row.get(version_col, ""))
            if v:
                versions.add(v)

    version = f"PubChem:{join_values(versions)}" if versions else "PubChem:local_tox_file"
    report = {
        "available": True,
        "file_path": str(path),
        "rows_read": len(rows),
        "rows_with_any_field": mapped_rows,
        "cid_records_loaded": len(out),
        "field_non_empty_counts": {k: int(field_non_empty.get(k, 0)) for k in FIELDS},
        "detected_columns": {
            "cid_col": cid_col,
            "source_ref_col": source_ref_col,
            "version_col": version_col,
            "field_cols": col_map,
        },
    }
    return out, report, version


def pick_best_drugbank_record(dbids: Sequence[str], db_map: Dict[str, Dict[str, str]]) -> Tuple[str, Dict[str, str]]:
    best_id = ""
    best = {}
    best_score = -1
    for dbid in dbids:
        rec = db_map.get(dbid)
        if not rec:
            continue
        score = sum(1 for f in FIELDS if normalize(rec.get(f, "")))
        if score > best_score:
            best_score = score
            best_id = dbid
            best = rec
    return best_id, best


def pick_best_pubchem_record(cids: Sequence[str], cid_map: Dict[str, Dict[str, str]]) -> Tuple[str, Dict[str, str]]:
    best_cid = ""
    best = {}
    best_score = -1
    for cid in cids:
        rec = cid_map.get(cid)
        if not rec:
            continue
        score = sum(1 for f in FIELDS if normalize(rec.get(f, "")))
        if score > best_score:
            best_score = score
            best_cid = cid
            best = rec
    return best_cid, best


def maybe_apply(row: Dict[str, str], field: str, value: str, evidence: str, source: str) -> bool:
    value = normalize(value)
    if not value:
        return False
    if evidence not in EVIDENCE_TIERS:
        return False
    cur_ev = row.get(f"{field}_evidence", "none")
    if EVIDENCE_RANK.get(evidence, 0) >= EVIDENCE_RANK.get(cur_ev, 0):
        row[field] = value
        row[f"{field}_evidence"] = evidence
        row[f"{field}_source"] = normalize(source) or "none"
        return True
    return False


def build_rows(
    xref_rows: Sequence[XrefRec],
    db_map: Dict[str, Dict[str, str]],
    pubchem_map: Dict[str, Dict[str, str]],
    source_version: str,
) -> Tuple[List[Dict[str, str]], Dict[str, Any]]:
    out: List[Dict[str, str]] = []
    field_non_empty_counts = Counter()
    evidence_counts = Counter()
    rows_with_any = 0
    rows_from_db = 0
    rows_from_pubchem = 0

    for rec in xref_rows:
        row: Dict[str, str] = {
            "inchikey": rec.inchikey,
            "drugbank_id": rec.drugbank_ids[0] if rec.drugbank_ids else "",
            "pubchem_cid": rec.pubchem_cids[0] if rec.pubchem_cids else "",
            "fetch_date": utc_today(),
            "source_version": source_version,
        }
        for f in FIELDS:
            row[f] = ""
            row[f"{f}_evidence"] = "none"
            row[f"{f}_source"] = "none"

        selected_dbid, db_pktox = pick_best_drugbank_record(rec.drugbank_ids, db_map)
        if selected_dbid:
            row["drugbank_id"] = selected_dbid
            touched = False
            for field in FIELDS:
                if maybe_apply(
                    row,
                    field,
                    db_pktox.get(field, ""),
                    db_pktox.get(f"{field}_evidence", "none"),
                    db_pktox.get(f"{field}_source", "none"),
                ):
                    touched = True
            if touched:
                rows_from_db += 1

        selected_cid, pubchem_pktox = pick_best_pubchem_record(rec.pubchem_cids, pubchem_map)
        if selected_cid:
            row["pubchem_cid"] = selected_cid
            touched = False
            for field in FIELDS:
                if maybe_apply(
                    row,
                    field,
                    pubchem_pktox.get(field, ""),
                    pubchem_pktox.get(f"{field}_evidence", "none"),
                    pubchem_pktox.get(f"{field}_source", "none"),
                ):
                    touched = True
            if touched:
                rows_from_pubchem += 1

        any_field = False
        for f in FIELDS:
            if normalize(row[f]):
                any_field = True
                field_non_empty_counts[f] += 1
            evidence_counts[row[f"{f}_evidence"]] += 1
        if any_field:
            rows_with_any += 1

        out.append(row)

    metrics = {
        "rows_total": len(out),
        "rows_with_any_pk_tox": rows_with_any,
        "field_non_empty_counts": {k: int(field_non_empty_counts.get(k, 0)) for k in FIELDS},
        "evidence_tier_counts": {k: int(evidence_counts.get(k, 0)) for k in sorted(EVIDENCE_TIERS)},
        "rows_with_drugbank_contrib": rows_from_db,
        "rows_with_pubchem_contrib": rows_from_pubchem,
    }
    return out, metrics


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xref", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--drugbank-xml", type=Path)
    ap.add_argument("--pubchem-tox", type=Path)
    ap.add_argument("--max-rows", type=int, default=None)
    args = ap.parse_args()

    if not args.xref.exists():
        raise SystemExit(f"[ERROR] missing xref input: {args.xref}")

    xref_rows, xref_stats, target_dbids, _ = load_xref(args.xref, max_rows=args.max_rows)

    db_map, db_report, db_version = load_drugbank_pk_tox(args.drugbank_xml, target_dbids=target_dbids)
    pubchem_map, pubchem_report, pubchem_version = load_pubchem_tox(args.pubchem_tox)

    source_versions = ["molecule_xref_core_v2"]
    if db_report.get("available"):
        source_versions.append(db_version)
    if pubchem_report.get("available"):
        source_versions.append(pubchem_version)
    merged_source_version = join_values(source_versions) or "molecule_xref_core_v2"

    rows, build_metrics = build_rows(
        xref_rows=xref_rows,
        db_map=db_map,
        pubchem_map=pubchem_map,
        source_version=merged_source_version,
    )

    header = [
        "inchikey",
        "drugbank_id",
        "pubchem_cid",
        "bioavailability",
        "bioavailability_evidence",
        "bioavailability_source",
        "clearance",
        "clearance_evidence",
        "clearance_source",
        "half_life",
        "half_life_evidence",
        "half_life_source",
        "volume_of_distribution",
        "volume_of_distribution_evidence",
        "volume_of_distribution_source",
        "ld50",
        "ld50_evidence",
        "ld50_source",
        "mutagenicity",
        "mutagenicity_evidence",
        "mutagenicity_source",
        "fetch_date",
        "source_version",
    ]

    written = write_tsv(args.out, rows=rows, header=header)

    report = {
        "name": "molecule_pk_tox_v1.build",
        "created_at": utc_now(),
        "inputs": {
            "xref": str(args.xref),
            "drugbank_xml": str(args.drugbank_xml) if args.drugbank_xml else None,
            "pubchem_tox": str(args.pubchem_tox) if args.pubchem_tox else None,
            "max_rows": args.max_rows,
        },
        "source_status": {
            "drugbank": db_report,
            "pubchem_tox": pubchem_report,
        },
        "source_versions": {
            "drugbank_version": db_version,
            "pubchem_version": pubchem_version,
            "merged_source_version": merged_source_version,
        },
        "xref_stats": xref_stats,
        "metrics": {
            **build_metrics,
            "rows_written": written,
            "blocked_recommended": (build_metrics["rows_with_any_pk_tox"] == 0),
        },
    }

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(
        f"[OK] molecule_pk_tox_v1 build -> {args.out} "
        f"(rows={written}, rows_with_any_pk_tox={build_metrics['rows_with_any_pk_tox']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
