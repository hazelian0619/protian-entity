#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import warnings
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

warnings.filterwarnings("ignore", message="urllib3 .* doesn't match a supported version")
import requests
from lxml import etree

try:
    from requests import RequestsDependencyWarning

    warnings.filterwarnings("ignore", category=RequestsDependencyWarning)
except Exception:
    pass

INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
NUM_RE = re.compile(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")
FLOAT_RE = re.compile(r"^-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?$")

FIELDS = [
    "bioavailability",
    "clearance",
    "half_life",
    "volume_of_distribution",
    "ld50",
    "mutagenicity",
]
NUMERIC_FIELD_SET = {
    "bioavailability",
    "clearance",
    "half_life",
    "volume_of_distribution",
    "ld50",
}

EVIDENCE_TIERS = {"structured/high", "text_mined/medium", "inferred/low", "none"}
EVIDENCE_RANK = {"none": 0, "inferred/low": 1, "text_mined/medium": 2, "structured/high": 3}

CHEMBL_ACTIVITY_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
CHEMBL_STATUS_URL = "https://www.ebi.ac.uk/chembl/api/data/status.json"
CHEMBL_STANDARD_TYPES = "F,CL,T1/2,Vd,LD50,Mutagenicity"
CHEMBL_STD_TO_FIELD = {
    "F": "bioavailability",
    "CL": "clearance",
    "T1/2": "half_life",
    "VD": "volume_of_distribution",
    "Vd": "volume_of_distribution",
    "LD50": "ld50",
    "MUTAGENICITY": "mutagenicity",
    "Mutagenicity": "mutagenicity",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def utc_today() -> str:
    return datetime.now(timezone.utc).date().isoformat()


def norm(x: Any) -> str:
    return str(x or "").strip()


def normalize_ws(x: str) -> str:
    return re.sub(r"\s+", " ", norm(x))


def strip_html(x: str) -> str:
    t = re.sub(r"<[^>]+>", "", x or "")
    t = t.replace("&nbsp;", " ").replace("&lt;", "<").replace("&gt;", ">")
    t = t.replace("&amp;", "&").replace("&plusmn;", "±")
    return normalize_ws(t)


def valid_inchikey(x: str) -> bool:
    return bool(INCHIKEY_RE.match(norm(x).upper()))


def split_multi(x: str) -> List[str]:
    t = norm(x)
    if not t:
        return []
    return [p.strip() for p in t.split(";") if p.strip()]


def join_values(values: Iterable[str]) -> str:
    uniq = sorted({norm(v) for v in values if norm(v)})
    return ";".join(uniq)


def fmt_float(v: float, digits: int = 8) -> str:
    s = f"{v:.{digits}g}"
    return s.replace("e+", "e")


def parse_float(x: str) -> Optional[float]:
    t = norm(x)
    if not t:
        return None
    if not FLOAT_RE.match(t):
        return None
    try:
        return float(t)
    except Exception:
        return None


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


@dataclass
class XrefRec:
    inchikey: str
    chembl_ids: List[str]
    drugbank_ids: List[str]
    pubchem_cids: List[str]


@dataclass
class FieldRec:
    raw: str
    evidence: str
    source_anchor: str
    source_db: str


def detect_drugbank_version(xml_path: Path) -> str:
    try:
        head = xml_path.read_text(encoding="utf-8", errors="ignore")[:4096]
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


def first_sentence(x: str, keyword: Optional[str] = None, max_len: int = 360) -> str:
    t = strip_html(x)
    if not t:
        return ""
    sents = re.split(r"(?<=[\.!?])\s+", t)
    if keyword:
        kw = keyword.lower()
        for s in sents:
            if kw in s.lower():
                return s[:max_len]
    return sents[0][:max_len]


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


def infer_mutagenicity(raw: str) -> str:
    t = strip_html(raw).lower()
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
        "no mutagenicity",
    ]
    pos_markers = [
        "mutagenic",
        "genotoxic",
        "positive ames",
        "ames positive",
        "mutagenicity was measured",
    ]
    neg = any(m in t for m in neg_markers)
    pos = any(m in t for m in pos_markers)
    if pos and neg:
        return "mixed"
    if pos:
        return "positive"
    if neg:
        return "negative"
    return ""


def extract_ld50_sentence(raw: str) -> str:
    t = strip_html(raw)
    if not t:
        return ""
    pats = [
        r"([^.;\n]{0,120}\bLD\s*50\b[^.;\n]{0,260})",
        r"([^.;\n]{0,120}\blethal dose\b[^.;\n]{0,260})",
    ]
    for p in pats:
        m = re.search(p, t, flags=re.IGNORECASE)
        if m:
            return normalize_ws(m.group(1))[:360]
    return ""


def parse_bioavailability_pct(raw: str) -> Optional[float]:
    t = strip_html(raw)
    if not t:
        return None
    m_range = re.search(
        r"([-+]?\d+(?:\.\d+)?)\s*(?:-|to|–|~)\s*([-+]?\d+(?:\.\d+)?)\s*%",
        t,
        flags=re.IGNORECASE,
    )
    if m_range:
        a = float(m_range.group(1))
        b = float(m_range.group(2))
        v = (a + b) / 2.0
        if 0 <= v <= 1000:
            return v
    m_pct = re.search(r"([-+]?\d+(?:\.\d+)?)\s*%", t)
    if m_pct:
        v = float(m_pct.group(1))
        if 0 <= v <= 1000:
            return v
    m_word = re.search(r"([-+]?\d+(?:\.\d+)?)\s*(percent|per\s*cent)\b", t, flags=re.IGNORECASE)
    if m_word:
        v = float(m_word.group(1))
        if 0 <= v <= 1000:
            return v
    if "bioavailability" in t.lower():
        m = NUM_RE.search(t)
        if m:
            v = float(m.group(0))
            if 0 <= v <= 1.0:
                return v * 100.0
    return None


def _clean_unit_token(x: str) -> str:
    t = norm(x).lower()
    t = t.replace("μ", "u").replace("µ", "u").replace("−", "-").replace("·", ".")
    t = t.replace("hours", "h").replace("hour", "h").replace("hrs", "h").replace("hr", "h")
    t = t.replace("minutes", "min").replace("minute", "min").replace("mins", "min")
    t = t.replace("days", "day")
    t = t.replace("weeks", "week")
    t = t.replace(" ", "")
    return t


def _extract_num_unit_pairs(raw: str) -> List[Tuple[float, str]]:
    s = _clean_unit_token(strip_html(raw))
    out: List[Tuple[float, str]] = []
    for m in re.finditer(r"([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", s):
        try:
            val = float(m.group(1))
        except Exception:
            continue
        tail = s[m.end() : m.end() + 28]
        um = re.match(r"([a-z%./()^-]{1,24})", tail)
        if um:
            out.append((val, um.group(1)))
    return out


def parse_half_life_hours(raw: str) -> Optional[float]:
    t = strip_html(raw).lower().replace("μ", "u").replace("µ", "u")
    if not t:
        return None
    m = re.search(
        r"([-+]?\d+(?:\.\d+)?)\s*(?:-|to|–)?\s*([-+]?\d+(?:\.\d+)?)?\s*(hours?|hrs?|hr|h|minutes?|mins?|min|days?|day|d|weeks?|week|w)\b",
        t,
        flags=re.IGNORECASE,
    )
    if m:
        v1 = float(m.group(1))
        v2 = m.group(2)
        if v2 and FLOAT_RE.match(v2):
            v = (v1 + float(v2)) / 2.0
        else:
            v = v1
        unit = m.group(3).lower()
    else:
        m2 = re.search(r"([-+]?\d+(?:\.\d+)?)(h|min|day|week)\b", t)
        if not m2:
            return None
        v = float(m2.group(1))
        unit = m2.group(2).lower()

    if unit.startswith("min"):
        return v / 60.0
    if unit in {"d", "day", "days"}:
        return v * 24.0
    if unit in {"w", "week", "weeks"}:
        return v * 24.0 * 7.0
    return v


def parse_clearance_std(raw: str) -> Optional[Tuple[float, str]]:
    for value, unit in _extract_num_unit_pairs(raw):
        u = unit
        if re.match(r"ml(?:/|\.)min(?:-?1)?(?:/|\.)kg(?:-?1)?", u):
            return value, "mL/min/kg"
        if re.match(r"ul(?:/|\.)min(?:-?1)?(?:/|\.)kg(?:-?1)?", u):
            return value / 1000.0, "mL/min/kg"
        if re.match(r"l(?:/|\.)h(?:-?1)?(?:/|\.)kg(?:-?1)?", u):
            return (value * 1000.0 / 60.0), "mL/min/kg"
        if re.match(r"ml(?:/|\.)h(?:/|\.)kg", u):
            return value / 60.0, "mL/min/kg"

        if re.match(r"ml(?:/|\.)min(?:-?1)?$", u):
            return value, "mL/min"
        if re.match(r"ul(?:/|\.)min(?:-?1)?$", u):
            return value / 1000.0, "mL/min"
        if re.match(r"l(?:/|\.)h(?:-?1)?$", u):
            return value, "L/h"
    return None


def parse_volume_distribution_std(raw: str) -> Optional[Tuple[float, str]]:
    for value, unit in _extract_num_unit_pairs(raw):
        u = unit
        if re.match(r"l(?:/|\.)kg(?:-?1)?$", u):
            return value, "L/kg"
        if re.match(r"ml(?:/|\.)kg(?:-?1)?$", u):
            return value / 1000.0, "L/kg"
        if u.startswith("l/1.73m2"):
            return value, "L/1.73m2"
        if re.match(r"l$", u):
            return value, "L"
    return None


def parse_ld50_std(raw: str) -> Optional[Tuple[float, str]]:
    t = strip_html(raw)
    if not t:
        return None
    scope = t
    m_scope = re.search(r"ld\s*50.{0,180}", t, flags=re.IGNORECASE)
    if m_scope:
        scope = m_scope.group(0)

    s = _clean_unit_token(scope)
    for m in re.finditer(
        r"([-+]?\d+(?:\.\d+)?)(?:\s*(?:-|to|–)\s*[-+]?\d+(?:\.\d+)?)?\s*(mg/kg|mg\.kg-1|g/kg|ug/kg|u?gkg-1)",
        s,
    ):
        v = float(m.group(1))
        u = m.group(2)
        if u.startswith("mg"):
            return v, "mg/kg"
        if u.startswith("g/"):
            return v * 1000.0, "mg/kg"
        if u.startswith("ug") or u.startswith("ugkg"):
            return v / 1000.0, "mg/kg"

    m2 = re.search(
        r"([-+]?\d+(?:\.\d+)?)\s*(mg/kg|mg\.kg-1|g/kg|ug/kg)",
        _clean_unit_token(t),
    )
    if m2:
        v = float(m2.group(1))
        u = m2.group(2)
        if u.startswith("mg"):
            return v, "mg/kg"
        if u.startswith("g/"):
            return v * 1000.0, "mg/kg"
        if u.startswith("ug"):
            return v / 1000.0, "mg/kg"
    return None


def mutagenicity_label_std(raw: str) -> str:
    guess = infer_mutagenicity(raw)
    if guess in {"positive", "negative", "mixed"}:
        return guess
    return "unknown"


def parse_numeric_normalized(field: str, raw: str) -> Tuple[str, str]:
    if not norm(raw):
        return "", ""
    if field == "bioavailability":
        v = parse_bioavailability_pct(raw)
        return (fmt_float(v), "pct") if v is not None else ("", "")
    if field == "clearance":
        p = parse_clearance_std(raw)
        if p is None:
            return "", ""
        return fmt_float(p[0]), p[1]
    if field == "half_life":
        v = parse_half_life_hours(raw)
        return (fmt_float(v), "h") if v is not None else ("", "")
    if field == "volume_of_distribution":
        p = parse_volume_distribution_std(raw)
        if p is None:
            return "", ""
        return fmt_float(p[0]), p[1]
    if field == "ld50":
        p = parse_ld50_std(raw)
        if p is None:
            return "", ""
        return fmt_float(p[0]), p[1]
    return "", ""


def choose_col(fieldnames: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    low = {c.lower(): c for c in fieldnames}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    return None


def load_xref(path: Path, max_rows: Optional[int]) -> Tuple[List[XrefRec], Dict[str, Any], Set[str], Set[str], Set[str]]:
    cols, rows = read_tsv(path, max_rows=max_rows)
    if "inchikey" not in cols:
        raise SystemExit(f"[ERROR] xref missing inchikey: {path}")

    out: List[XrefRec] = []
    bad_inchikey = 0
    all_dbids: Set[str] = set()
    all_cids: Set[str] = set()
    all_chembl: Set[str] = set()

    db_col = choose_col(cols, ["drugbank_id"])
    cid_col = choose_col(cols, ["pubchem_cid", "pubchem_id"])
    chembl_col = choose_col(cols, ["chembl_id", "chembl"])

    rows_with_dbid = 0
    rows_with_cid = 0
    rows_with_chembl = 0

    for row in rows:
        ik = norm(row.get("inchikey", "")).upper()
        if not valid_inchikey(ik):
            bad_inchikey += 1
            continue

        dbids = split_multi(row.get(db_col, "") if db_col else "")
        cids = split_multi(row.get(cid_col, "") if cid_col else "")
        chembl_ids = split_multi(row.get(chembl_col, "") if chembl_col else "")

        if dbids:
            rows_with_dbid += 1
            all_dbids.update(dbids)
        if cids:
            rows_with_cid += 1
            all_cids.update(cids)
        if chembl_ids:
            rows_with_chembl += 1
            all_chembl.update(chembl_ids)

        out.append(XrefRec(inchikey=ik, chembl_ids=chembl_ids, drugbank_ids=dbids, pubchem_cids=cids))

    stats = {
        "input_rows": len(rows),
        "valid_rows": len(out),
        "skipped_bad_inchikey_rows": bad_inchikey,
        "rows_with_drugbank_id": rows_with_dbid,
        "rows_with_pubchem_cid": rows_with_cid,
        "rows_with_chembl_id": rows_with_chembl,
        "unique_drugbank_ids": len(all_dbids),
        "unique_pubchem_cids": len(all_cids),
        "unique_chembl_ids": len(all_chembl),
    }
    return out, stats, all_dbids, all_cids, all_chembl


def resolve_drugbank_xml(cli_path: Optional[Path]) -> Optional[Path]:
    cands: List[Path] = []
    if cli_path:
        cands.append(cli_path)
    cands.extend(
        [
            Path("data/raw/drugbank/drugbank_complete_database.xml"),
            Path("data/raw/drugbank/drugbank_all_full_database.xml"),
            Path("../1218/data/raw/drugbank/drugbank_complete_database_2026-01-04.xml"),
            Path("../1218/data/raw/drugbank/drugbank_all_full_database.xml"),
        ]
    )
    cands.extend(sorted(Path("data/raw/drugbank").glob("*.xml")))
    cands.extend(sorted(Path("../1218/data/raw/drugbank").glob("*.xml")))

    existing = [p for p in cands if p.exists()]
    if not existing:
        return None
    existing = sorted(existing, key=lambda p: p.stat().st_mtime, reverse=True)
    return existing[0]


def load_drugbank_pk_tox(xml_path: Optional[Path], target_dbids: Set[str]) -> Tuple[Dict[str, Dict[str, FieldRec]], Dict[str, Any], str]:
    if xml_path is None or not xml_path.exists():
        return {}, {"available": False, "reason": "drugbank_xml_not_found"}, "DrugBank:unavailable"

    out: Dict[str, Dict[str, FieldRec]] = {}
    parsed_top_level = 0
    matched_target = 0
    value_counter = Counter()
    version = detect_drugbank_version(xml_path)

    for _, drug in etree.iterparse(str(xml_path), events=("end",), tag="{*}drug", huge_tree=True, recover=True):
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

        rec: Dict[str, FieldRec] = {}
        bio = first_sentence(absorption, keyword="bioavailability")
        if bio:
            rec["bioavailability"] = FieldRec(bio, "structured/high", f"drugbank_xml:{dbid}#absorption", "drugbank_xml")
            value_counter["bioavailability"] += 1

        clr = first_sentence(clearance)
        if clr:
            rec["clearance"] = FieldRec(clr, "structured/high", f"drugbank_xml:{dbid}#clearance", "drugbank_xml")
            value_counter["clearance"] += 1

        hl = first_sentence(half_life)
        if hl:
            rec["half_life"] = FieldRec(hl, "structured/high", f"drugbank_xml:{dbid}#half-life", "drugbank_xml")
            value_counter["half_life"] += 1

        vd = first_sentence(volume)
        if vd:
            rec["volume_of_distribution"] = FieldRec(
                vd,
                "structured/high",
                f"drugbank_xml:{dbid}#volume-of-distribution",
                "drugbank_xml",
            )
            value_counter["volume_of_distribution"] += 1

        ld50 = extract_ld50_sentence(toxicity)
        if ld50:
            rec["ld50"] = FieldRec(ld50, "text_mined/medium", f"drugbank_xml:{dbid}#toxicity", "drugbank_xml")
            value_counter["ld50"] += 1

        mut = infer_mutagenicity(toxicity)
        if mut:
            rec["mutagenicity"] = FieldRec(mut, "inferred/low", f"drugbank_xml:{dbid}#toxicity", "drugbank_xml")
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


def _chembl_req_json(session: requests.Session, url: str, params: Dict[str, Any], timeout: int, retries: int = 3) -> Dict[str, Any]:
    last_err: Optional[Exception] = None
    for i in range(retries):
        try:
            r = session.get(url, params=params, timeout=timeout)
            if r.status_code >= 500:
                raise RuntimeError(f"HTTP {r.status_code}: {r.text[:240]}")
            if r.status_code != 200:
                raise RuntimeError(f"HTTP {r.status_code}: {r.text[:240]}")
            return r.json()
        except Exception as e:
            last_err = e
            if i + 1 == retries:
                break
    raise RuntimeError(f"chembl request failed after {retries} tries: {last_err}")


def _chunked(items: Sequence[str], n: int) -> Iterable[List[str]]:
    for i in range(0, len(items), n):
        yield list(items[i : i + n])


def build_chembl_raw(std_type: str, value: str, units: str, relation: str, desc: str) -> str:
    rel = norm(relation)
    val = norm(value)
    u = norm(units)
    if std_type.upper() != "MUTAGENICITY" and not val:
        return ""
    head = normalize_ws(f"{rel}{val} {u}".strip())
    if std_type.upper() == "MUTAGENICITY":
        extra = strip_html(desc)
        if extra:
            return normalize_ws(f"{head}; {extra}" if head else extra)[:360]
    return head[:160]


def score_candidate(field: str, rec: FieldRec, relation: str, year: Optional[int]) -> float:
    score = 0.0
    if field in NUMERIC_FIELD_SET:
        v, u = parse_numeric_normalized(field, rec.raw)
        if v and u:
            score += 10.0
        elif norm(rec.raw):
            score += 3.0
    elif field == "mutagenicity":
        label = mutagenicity_label_std(rec.raw)
        if label != "unknown":
            score += 7.0
        elif norm(rec.raw):
            score += 2.0

    rel = norm(relation)
    if rel in {"", "="}:
        score += 2.0
    else:
        score += 1.0

    if year is not None and year > 1900:
        score += min((year - 1900) / 200.0, 1.0)
    return score


def write_chembl_cache(path: Path, data: Dict[str, Dict[str, FieldRec]], meta: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    header = ["chembl_id", "field", "raw", "evidence", "source_anchor", "source_db"]
    with tmp.open("w", encoding="utf-8", newline="") as f:
        f.write("# meta=" + json.dumps(meta, ensure_ascii=False) + "\n")
        w = csv.DictWriter(f, delimiter="\t", fieldnames=header, lineterminator="\n")
        w.writeheader()
        for chembl_id, fields in sorted(data.items()):
            for field, rec in sorted(fields.items()):
                w.writerow(
                    {
                        "chembl_id": chembl_id,
                        "field": field,
                        "raw": rec.raw,
                        "evidence": rec.evidence,
                        "source_anchor": rec.source_anchor,
                        "source_db": rec.source_db,
                    }
                )
    tmp.replace(path)


def read_chembl_cache(path: Path) -> Tuple[Dict[str, Dict[str, FieldRec]], Dict[str, Any]]:
    if not path.exists():
        return {}, {}
    meta: Dict[str, Any] = {}
    out: Dict[str, Dict[str, FieldRec]] = {}
    with path.open("r", encoding="utf-8", newline="") as f:
        first = f.readline()
        if first.startswith("# meta="):
            try:
                meta = json.loads(first[len("# meta=") :])
            except Exception:
                meta = {}
        else:
            f.seek(0)

        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            return {}, meta
        for row in r:
            chembl_id = norm(row.get("chembl_id", ""))
            field = norm(row.get("field", ""))
            if not chembl_id or not field:
                continue
            out.setdefault(chembl_id, {})[field] = FieldRec(
                raw=norm(row.get("raw", "")),
                evidence=norm(row.get("evidence", "structured/high")) or "structured/high",
                source_anchor=norm(row.get("source_anchor", "none")) or "none",
                source_db=norm(row.get("source_db", "chembl_api")) or "chembl_api",
            )
    return out, meta


def load_chembl_pk_tox(
    target_chembl_ids: Set[str],
    cache_path: Path,
    force_refresh: bool,
    chunk_size: int,
    page_limit: int,
    timeout: int,
) -> Tuple[Dict[str, Dict[str, FieldRec]], Dict[str, Any], str]:
    if not target_chembl_ids:
        return {}, {"available": False, "reason": "no_target_chembl_ids"}, "ChEMBL:unavailable"

    if cache_path.exists() and not force_refresh:
        cached_data, cached_meta = read_chembl_cache(cache_path)
        if cached_data:
            version = norm(cached_meta.get("source_version", "ChEMBL:cache")) or "ChEMBL:cache"
            rep = {
                "available": True,
                "mode": "cache",
                "cache_path": str(cache_path),
                "records_loaded": len(cached_data),
                "field_non_empty_counts": cached_meta.get("field_non_empty_counts", {}),
                "fetch_date": cached_meta.get("fetch_date"),
            }
            return cached_data, rep, version

    session = requests.Session()
    status = _chembl_req_json(session, CHEMBL_STATUS_URL, params={}, timeout=timeout)
    chembl_ver = norm(status.get("chembl_db_version", ""))
    release_date = norm(status.get("chembl_release_date", ""))
    source_version = f"ChEMBL:{chembl_ver}|release_date:{release_date}" if chembl_ver else "ChEMBL:api"

    target_upper = {norm(c).upper() for c in target_chembl_ids if norm(c)}
    target = sorted(target_upper)
    best: Dict[str, Dict[str, Tuple[FieldRec, float]]] = {}
    field_counts = Counter()

    chunks = list(_chunked(target, max(1, chunk_size)))
    request_count = 0
    activities_scanned = 0

    for chunk in chunks:
        base_params = {
            "standard_type__in": CHEMBL_STANDARD_TYPES,
            "molecule_chembl_id__in": ",".join(chunk),
            "limit": page_limit,
        }
        page0 = _chembl_req_json(session, CHEMBL_ACTIVITY_URL, params={**base_params, "offset": 0}, timeout=timeout)
        request_count += 1
        total = int(page0.get("page_meta", {}).get("total_count", 0))

        pages = [page0]
        for offset in range(page_limit, total, page_limit):
            p = _chembl_req_json(session, CHEMBL_ACTIVITY_URL, params={**base_params, "offset": offset}, timeout=timeout)
            pages.append(p)
            request_count += 1

        for page in pages:
            for a in page.get("activities", []):
                activities_scanned += 1
                chembl_id = norm(a.get("molecule_chembl_id", "")).upper()
                if not chembl_id or chembl_id not in target_upper:
                    continue

                std_type = norm(a.get("standard_type", ""))
                field = CHEMBL_STD_TO_FIELD.get(std_type, CHEMBL_STD_TO_FIELD.get(std_type.upper(), ""))
                if field not in FIELDS:
                    continue

                raw = build_chembl_raw(
                    std_type=std_type,
                    value=norm(a.get("standard_value", "")),
                    units=norm(a.get("standard_units", "")),
                    relation=norm(a.get("standard_relation", "")),
                    desc=norm(a.get("assay_description", "")),
                )
                if not raw:
                    continue

                act_id = norm(a.get("activity_id", ""))
                assay_id = norm(a.get("assay_chembl_id", ""))
                source_anchor = f"chembl_api:activity:{act_id}|molecule:{chembl_id}|assay:{assay_id}|standard_type:{std_type}"
                fr = FieldRec(raw=raw, evidence="structured/high", source_anchor=source_anchor, source_db="chembl_api")

                year: Optional[int] = None
                y = a.get("document_year")
                if isinstance(y, int):
                    year = y
                elif norm(y).isdigit():
                    year = int(norm(y))

                score = score_candidate(field, fr, relation=norm(a.get("standard_relation", "")), year=year)

                cur = best.setdefault(chembl_id, {}).get(field)
                if cur is None or score > cur[1]:
                    best[chembl_id][field] = (fr, score)

    out: Dict[str, Dict[str, FieldRec]] = {}
    for chembl_id, fmap in best.items():
        out[chembl_id] = {}
        for field, (rec, _) in fmap.items():
            out[chembl_id][field] = rec
            field_counts[field] += 1

    cache_meta = {
        "created_at": utc_now(),
        "fetch_date": utc_today(),
        "source_version": source_version,
        "records_loaded": len(out),
        "field_non_empty_counts": {k: int(field_counts.get(k, 0)) for k in FIELDS},
        "request_count": request_count,
        "activities_scanned": activities_scanned,
    }
    write_chembl_cache(cache_path, out, cache_meta)

    report = {
        "available": True,
        "mode": "download",
        "request_count": request_count,
        "activities_scanned": activities_scanned,
        "chunk_count": len(chunks),
        "target_chembl_ids": len(target),
        "records_loaded": len(out),
        "field_non_empty_counts": {k: int(field_counts.get(k, 0)) for k in FIELDS},
        "cache_path": str(cache_path),
    }
    return out, report, source_version


def pick_best_map(ids: Sequence[str], mapping: Dict[str, Dict[str, FieldRec]]) -> Tuple[str, Dict[str, FieldRec]]:
    best_id = ""
    best_rec: Dict[str, FieldRec] = {}
    best_score = -1
    for _id in ids:
        m = mapping.get(_id) or mapping.get(norm(_id).upper()) or mapping.get(norm(_id))
        if not m:
            continue
        score = sum(1 for f in FIELDS if f in m and norm(m[f].raw))
        if score > best_score:
            best_score = score
            best_id = _id
            best_rec = m
    return best_id, best_rec


def maybe_apply_field(
    row: Dict[str, str],
    source_db_used: Set[str],
    field: str,
    rec: FieldRec,
    only_if_empty: bool = False,
) -> bool:
    raw = norm(rec.raw)
    if not raw:
        return False
    if rec.evidence not in EVIDENCE_TIERS:
        return False

    if only_if_empty and norm(row.get(field, "")):
        return False

    cur_ev = row.get(f"{field}_evidence", "none")
    if EVIDENCE_RANK.get(rec.evidence, 0) < EVIDENCE_RANK.get(cur_ev, 0):
        return False

    row[field] = raw
    row[f"{field}_evidence"] = rec.evidence
    row[f"{field}_source"] = norm(rec.source_anchor) or "none"
    source_db_used.add(norm(rec.source_db) or "unknown")
    return True


def apply_normalization_columns(row: Dict[str, str]) -> Dict[str, bool]:
    flags: Dict[str, bool] = {}

    b_val, _ = parse_numeric_normalized("bioavailability", row.get("bioavailability", ""))
    row["bioavailability_pct"] = b_val
    flags["bioavailability_pct"] = bool(b_val)

    c_val, c_unit = parse_numeric_normalized("clearance", row.get("clearance", ""))
    row["clearance_value"] = c_val
    row["clearance_unit_std"] = c_unit
    flags["clearance"] = bool(c_val and c_unit)

    h_val, h_unit = parse_numeric_normalized("half_life", row.get("half_life", ""))
    row["half_life_value"] = h_val
    row["half_life_unit_std"] = h_unit
    flags["half_life"] = bool(h_val and h_unit)

    v_val, v_unit = parse_numeric_normalized("volume_of_distribution", row.get("volume_of_distribution", ""))
    row["volume_distribution_value"] = v_val
    row["volume_distribution_unit_std"] = v_unit
    flags["volume_of_distribution"] = bool(v_val and v_unit)

    l_val, l_unit = parse_numeric_normalized("ld50", row.get("ld50", ""))
    row["ld50_value"] = l_val
    row["ld50_unit_std"] = l_unit
    flags["ld50"] = bool(l_val and l_unit)

    row["mutagenicity_label_std"] = mutagenicity_label_std(row.get("mutagenicity", ""))
    flags["mutagenicity"] = row["mutagenicity_label_std"] in {"positive", "negative", "mixed", "unknown"}

    return flags


def build_rows(
    xref_rows: Sequence[XrefRec],
    db_map: Dict[str, Dict[str, FieldRec]],
    chembl_map: Dict[str, Dict[str, FieldRec]],
    source_version: str,
) -> Tuple[List[Dict[str, str]], Dict[str, Any]]:
    out: List[Dict[str, str]] = []

    field_non_empty = Counter()
    evidence_counts = Counter()
    rows_with_any = 0
    rows_with_drugbank_contrib = 0
    rows_with_second_source_contrib = 0

    normalization_raw_counts = Counter()
    normalization_ok_counts = Counter()
    rows_with_any_normalized = 0

    for rec in xref_rows:
        row: Dict[str, str] = {
            "inchikey": rec.inchikey,
            "chembl_id": rec.chembl_ids[0] if rec.chembl_ids else "",
            "drugbank_id": rec.drugbank_ids[0] if rec.drugbank_ids else "",
            "pubchem_cid": rec.pubchem_cids[0] if rec.pubchem_cids else "",
            "fetch_date": utc_today(),
            "source_version": source_version,
            "source": "molecule_xref_core_v2",
        }
        for f in FIELDS:
            row[f] = ""
            row[f"{f}_evidence"] = "none"
            row[f"{f}_source"] = "none"

        source_db_used: Set[str] = set()
        row_db_touched = False
        row_chembl_touched = False

        chosen_dbid, db_fields = pick_best_map(rec.drugbank_ids, db_map)
        if chosen_dbid:
            row["drugbank_id"] = chosen_dbid
            for field in FIELDS:
                fr = db_fields.get(field)
                if fr and maybe_apply_field(row, source_db_used, field, fr, only_if_empty=False):
                    row_db_touched = True

        chosen_chembl, chembl_fields = pick_best_map(rec.chembl_ids, chembl_map)
        if chosen_chembl:
            row["chembl_id"] = chosen_chembl
            for field in FIELDS:
                fr = chembl_fields.get(field)
                if not fr:
                    continue

                # fill missing first
                if maybe_apply_field(row, source_db_used, field, fr, only_if_empty=True):
                    row_chembl_touched = True
                    continue

                # for numeric fields, replace with ChEMBL when current value cannot be normalized but ChEMBL can
                if field in NUMERIC_FIELD_SET and norm(row.get(field, "")):
                    cur_val, cur_unit = parse_numeric_normalized(field, row.get(field, ""))
                    new_val, new_unit = parse_numeric_normalized(field, fr.raw)
                    if (not cur_val or not cur_unit) and (new_val and new_unit):
                        if maybe_apply_field(row, source_db_used, field, fr, only_if_empty=False):
                            row_chembl_touched = True

        if row_db_touched:
            rows_with_drugbank_contrib += 1
        if row_chembl_touched:
            rows_with_second_source_contrib += 1

        norm_flags = apply_normalization_columns(row)

        any_field = False
        any_normalized = False
        for f in FIELDS:
            if norm(row.get(f, "")):
                any_field = True
                field_non_empty[f] += 1

            evidence_counts[row.get(f"{f}_evidence", "none")] += 1

            if f in NUMERIC_FIELD_SET and norm(row.get(f, "")):
                normalization_raw_counts[f] += 1
                key = "bioavailability_pct" if f == "bioavailability" else f
                if norm_flags.get(key):
                    normalization_ok_counts[f] += 1
                    any_normalized = True
            if f == "mutagenicity" and norm(row.get(f, "")):
                normalization_raw_counts[f] += 1
                if norm_flags.get("mutagenicity"):
                    normalization_ok_counts[f] += 1
                    any_normalized = True

        if any_field:
            rows_with_any += 1
        if any_normalized:
            rows_with_any_normalized += 1

        if source_db_used:
            row["source"] = join_values(source_db_used)

        out.append(row)

    normalization_rates = {
        f: (
            (normalization_ok_counts.get(f, 0) / normalization_raw_counts[f])
            if normalization_raw_counts.get(f, 0)
            else 1.0
        )
        for f in ["bioavailability", "clearance", "half_life", "volume_of_distribution", "ld50", "mutagenicity"]
    }

    metrics = {
        "rows_total": len(out),
        "rows_with_any_pk_tox": rows_with_any,
        "rows_with_any_normalized": rows_with_any_normalized,
        "rows_with_drugbank_contrib": rows_with_drugbank_contrib,
        "rows_with_second_source_contrib": rows_with_second_source_contrib,
        "field_non_empty_counts": {k: int(field_non_empty.get(k, 0)) for k in FIELDS},
        "evidence_tier_counts": {k: int(evidence_counts.get(k, 0)) for k in sorted(EVIDENCE_TIERS)},
        "normalization_raw_non_empty_counts": {k: int(normalization_raw_counts.get(k, 0)) for k in ["bioavailability", "clearance", "half_life", "volume_of_distribution", "ld50", "mutagenicity"]},
        "normalization_success_counts": {k: int(normalization_ok_counts.get(k, 0)) for k in ["bioavailability", "clearance", "half_life", "volume_of_distribution", "ld50", "mutagenicity"]},
        "normalization_success_rate_when_raw_present": normalization_rates,
    }
    return out, metrics


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xref", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--drugbank-xml", type=Path, default=None)
    ap.add_argument("--chembl-cache", type=Path, default=Path("pipelines/molecule_pk_tox_v2/.cache/chembl_pk_tox_v2.cache.tsv"))
    ap.add_argument("--force-refresh-chembl", action="store_true")
    ap.add_argument("--chembl-chunk-size", type=int, default=220)
    ap.add_argument("--chembl-page-limit", type=int, default=1000)
    ap.add_argument("--chembl-timeout", type=int, default=60)
    ap.add_argument("--max-rows", type=int, default=None)
    args = ap.parse_args()

    if not args.xref.exists():
        raise SystemExit(f"[ERROR] missing xref input: {args.xref}")

    xref_rows, xref_stats, target_dbids, _, target_chembl_ids = load_xref(args.xref, max_rows=args.max_rows)

    resolved_drugbank_xml = resolve_drugbank_xml(args.drugbank_xml)
    db_map, db_report, db_version = load_drugbank_pk_tox(resolved_drugbank_xml, target_dbids=target_dbids)

    chembl_map, chembl_report, chembl_version = load_chembl_pk_tox(
        target_chembl_ids=target_chembl_ids,
        cache_path=args.chembl_cache,
        force_refresh=args.force_refresh_chembl,
        chunk_size=args.chembl_chunk_size,
        page_limit=args.chembl_page_limit,
        timeout=args.chembl_timeout,
    )

    source_versions = ["molecule_xref_core_v2"]
    if db_report.get("available"):
        source_versions.append(db_version)
    if chembl_report.get("available"):
        source_versions.append(chembl_version)
    merged_source_version = join_values(source_versions) or "molecule_xref_core_v2"

    rows, build_metrics = build_rows(
        xref_rows=xref_rows,
        db_map=db_map,
        chembl_map=chembl_map,
        source_version=merged_source_version,
    )

    header = [
        "inchikey",
        "chembl_id",
        "drugbank_id",
        "pubchem_cid",
        "bioavailability",
        "bioavailability_pct",
        "bioavailability_evidence",
        "bioavailability_source",
        "clearance",
        "clearance_value",
        "clearance_unit_std",
        "clearance_evidence",
        "clearance_source",
        "half_life",
        "half_life_value",
        "half_life_unit_std",
        "half_life_evidence",
        "half_life_source",
        "volume_of_distribution",
        "volume_distribution_value",
        "volume_distribution_unit_std",
        "volume_of_distribution_evidence",
        "volume_of_distribution_source",
        "ld50",
        "ld50_value",
        "ld50_unit_std",
        "ld50_evidence",
        "ld50_source",
        "mutagenicity",
        "mutagenicity_label_std",
        "mutagenicity_evidence",
        "mutagenicity_source",
        "source",
        "source_version",
        "fetch_date",
    ]

    written = write_tsv(args.out, rows=rows, header=header)

    report = {
        "name": "molecule_pk_tox_v2.build",
        "created_at": utc_now(),
        "inputs": {
            "xref": str(args.xref),
            "drugbank_xml": str(resolved_drugbank_xml) if resolved_drugbank_xml else None,
            "chembl_cache": str(args.chembl_cache),
            "force_refresh_chembl": bool(args.force_refresh_chembl),
            "chembl_chunk_size": args.chembl_chunk_size,
            "chembl_page_limit": args.chembl_page_limit,
            "chembl_timeout": args.chembl_timeout,
            "max_rows": args.max_rows,
        },
        "source_status": {
            "drugbank": db_report,
            "chembl_api": chembl_report,
        },
        "source_versions": {
            "drugbank_version": db_version,
            "chembl_version": chembl_version,
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
        f"[OK] molecule_pk_tox_v2 build -> {args.out} "
        f"(rows={written}, rows_with_any_pk_tox={build_metrics['rows_with_any_pk_tox']}, "
        f"rows_with_second_source_contrib={build_metrics['rows_with_second_source_contrib']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
