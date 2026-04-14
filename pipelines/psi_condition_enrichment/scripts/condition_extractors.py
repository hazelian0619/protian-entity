#!/usr/bin/env python3
from __future__ import annotations

import json
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

SOURCE_FIELDS: Sequence[str] = (
    "assay_description",
    "assay_context",
    "activity_comment",
    "data_validity_comment",
)

SOURCE_PRIORITY = {
    "condition_pH": 0,
    "condition_temperature_c": 0,
    "condition_system": 0,
    "condition_context": 0,
    "assay_description": 1,
    "assay_context": 2,
    "activity_comment": 3,
    "data_validity_comment": 4,
}

PH_RANGE_RE = re.compile(
    r"\bp\s*[-/]?\s*h\s*(?:of\s*)?(?:=|:)?\s*(\d{1,2}(?:\.\d{1,2})?)\s*(?:-|–|to|~)\s*(\d{1,2}(?:\.\d{1,2})?)\b",
    re.IGNORECASE,
)
PH_BETWEEN_RE = re.compile(
    r"\bbetween\s+p\s*[-/]?\s*h\s*(\d{1,2}(?:\.\d{1,2})?)\s*(?:and|&)\s*(\d{1,2}(?:\.\d{1,2})?)\b",
    re.IGNORECASE,
)
PH_SINGLE_RE = re.compile(
    r"\bp\s*[-/]?\s*h\s*(?:of\s*)?(?:=|:)?\s*(\d{1,2}(?:\.\d{1,2})?)\b",
    re.IGNORECASE,
)

TEMP_RANGE_RE = re.compile(
    r"(-?\d{1,3}(?:\.\d+)?)\s*(?:-|–|to|~)\s*(-?\d{1,3}(?:\.\d+)?)\s*°?\s*([CFK])\b",
    re.IGNORECASE,
)
TEMP_SINGLE_RE = re.compile(
    r"(-?\d{1,3}(?:\.\d+)?)\s*(?:°\s*)?([CFK])\b",
    re.IGNORECASE,
)
ROOM_TEMP_RE = re.compile(r"\b(room temperature|ambient temperature)\b", re.IGNORECASE)
ICE_TEMP_RE = re.compile(r"\b(on ice|ice-cold|iced)\b", re.IGNORECASE)

INCUBATION_RE = re.compile(r"\bincubat(?:ed|ion)\b|\b\d+\s*(?:h|hr|hrs|hour|hours)\b", re.IGNORECASE)
MEDIA_RE = re.compile(r"\b(DMEM|RPMI|FBS|culture medium|serum[- ]free medium|HBSS|KRH)\b", re.IGNORECASE)

SYSTEM_REGEX = {
    "cell_based": re.compile(r"\b(cell[- ]?based|cell line|cells?)\b", re.IGNORECASE),
    "cell_free": re.compile(r"\b(cell[- ]?free|cell-free)\b", re.IGNORECASE),
    "microsome": re.compile(r"\bmicrosom(?:e|al)\b", re.IGNORECASE),
    "plasma_serum": re.compile(r"\b(plasma|serum)\b", re.IGNORECASE),
    "lysate": re.compile(r"\b(lysate|cytosol|homogenate)\b", re.IGNORECASE),
    "membrane": re.compile(r"\bmembrane\b", re.IGNORECASE),
    "recombinant": re.compile(r"\brecombinant\b", re.IGNORECASE),
    "buffered": re.compile(r"\bbuffer\b|\bPBS\b|\bHEPES\b|\bTRIS\b|\bMOPS\b|\bMES\b", re.IGNORECASE),
    "in_vitro": re.compile(r"\bin\s+vitro\b", re.IGNORECASE),
    "in_vivo": re.compile(r"\bin\s+vivo\b", re.IGNORECASE),
}

BUFFER_RE = re.compile(
    r"\b(Tris(?:-HCl)?|HEPES|PBS|phosphate\s+buffer(?:ed\s+saline)?|MOPS|MES|citrate\s+buffer|acetate\s+buffer|borate\s+buffer|DMEM|RPMI|HBSS|KRH)\b",
    re.IGNORECASE,
)

KNOWN_CELL_LINE_RE = [
    re.compile(r"\b(HEK\-?293|293T|HEK-Blue|CHO|HeLa|Jurkat|A549|U2OS|MCF\-?7|HT\-?29|HCT\-?\d+|MDCK|Vero|COS\-?7|BHK\-?21)\b", re.IGNORECASE),
]
GENERIC_CELL_RE = re.compile(r"\b([A-Za-z][A-Za-z0-9\-]{1,20})\s+cells\b", re.IGNORECASE)
CELL_STOPWORDS = {
    "human",
    "mouse",
    "rat",
    "monkey",
    "primary",
    "normal",
    "fresh",
    "whole",
    "blood",
    "tumor",
    "assay",
    "recombinant",
    "unknown",
}

BUFFER_PH_DEFAULT = {
    "PBS": "7.4",
    "PHOSPHATE_BUFFER": "7.4",
    "PHOSPHATE_BUFFERED_SALINE": "7.4",
    "HEPES": "7.4",
    "TRIS": "7.5",
    "TRIS-HCL": "7.5",
    "DMEM": "7.4",
    "RPMI": "7.4",
    "HBSS": "7.4",
    "KRH": "7.4",
    "MOPS": "7.2",
    "MES": "6.5",
    "CITRATE_BUFFER": "6",
    "ACETATE_BUFFER": "5.5",
    "BORATE_BUFFER": "8.5",
}


@dataclass
class Candidate:
    value: str
    source_field: str
    method: str
    confidence: float
    evidence: str


def non_empty(v: object) -> bool:
    s = str(v or "").strip()
    return s not in {"", "NA", "N/A", "None", "null"}


def _compact_text(v: object) -> str:
    s = str(v or "").strip()
    if not s:
        return ""
    return re.sub(r"\s+", " ", s)


def _fmt_num(x: float) -> str:
    s = f"{x:.2f}".rstrip("0").rstrip(".")
    return "0" if s in {"-0", "-0.0", ""} else s


def _norm_ph(v1: str, v2: str = "") -> str:
    try:
        x = float(v1)
    except Exception:
        return ""
    if x < 0 or x > 14:
        return ""
    if not v2:
        return _fmt_num(x)
    try:
        y = float(v2)
    except Exception:
        return ""
    if y < 0 or y > 14:
        return ""
    lo, hi = sorted([x, y])
    if abs(lo - hi) < 1e-9:
        return _fmt_num(lo)
    return f"{_fmt_num(lo)}-{_fmt_num(hi)}"


def _to_celsius(v: str, unit: str) -> str:
    try:
        x = float(v)
    except Exception:
        return ""
    u = unit.upper()
    c = x if u == "C" else (x - 32.0) * 5.0 / 9.0 if u == "F" else (x - 273.15) if u == "K" else None
    if c is None:
        return ""
    if c < -80 or c > 250:
        return ""
    return _fmt_num(c)


def _norm_temp(v1: str, unit: str, v2: str = "") -> str:
    c1 = _to_celsius(v1, unit)
    if not c1:
        return ""
    if not v2:
        return c1
    c2 = _to_celsius(v2, unit)
    if not c2:
        return ""
    lo, hi = sorted([float(c1), float(c2)])
    if abs(lo - hi) < 1e-9:
        return _fmt_num(lo)
    return f"{_fmt_num(lo)}-{_fmt_num(hi)}"


def _push_candidate(bag: List[Candidate], value: str, source_field: str, method: str, confidence: float, evidence: str) -> None:
    if not non_empty(value):
        return
    bag.append(Candidate(value=value, source_field=source_field, method=method, confidence=round(confidence, 4), evidence=evidence[:120]))


def _normalize_buffer(raw: str) -> str:
    s = _compact_text(raw).upper().replace(" ", "_")
    s = s.replace("PHOSPHATE_BUFFERED_SALINE", "PHOSPHATE_BUFFERED_SALINE")
    if s.startswith("TRIS"):
        return "TRIS"
    return s


def _resolve(cands: List[Candidate]) -> Tuple[str, Candidate | None, bool, List[Candidate]]:
    if not cands:
        return "", None, False, []
    ranked = sorted(
        cands,
        key=lambda c: (
            -float(c.confidence),
            SOURCE_PRIORITY.get(c.source_field, 99),
            c.method,
            c.value,
        ),
    )
    uniq = {c.value for c in ranked}
    return ranked[0].value, ranked[0], len(uniq) > 1, ranked


def extract_condition_bundle(row: Dict[str, str]) -> Dict[str, str]:
    texts: Dict[str, str] = {k: _compact_text(row.get(k, "")) for k in SOURCE_FIELDS}
    joined_text = " | ".join([texts[k] for k in SOURCE_FIELDS if texts[k]])
    lower_joined = joined_text.lower()

    ph_cands: List[Candidate] = []
    temp_cands: List[Candidate] = []

    if non_empty(row.get("condition_pH", "")):
        _push_candidate(ph_cands, _compact_text(row.get("condition_pH", "")), "condition_pH", "v2_existing", 0.99, "v2")
    if non_empty(row.get("condition_temperature_c", "")):
        _push_candidate(
            temp_cands,
            _compact_text(row.get("condition_temperature_c", "")),
            "condition_temperature_c",
            "v2_existing",
            0.99,
            "v2",
        )

    system_categories = set()
    cell_lines = set()
    buffers = set()
    context_source_fields = set()
    rule_hits: List[str] = []

    for src, txt in texts.items():
        if not txt:
            continue
        low = txt.lower()

        if "ph" in low or "p h" in low:
            for m in PH_RANGE_RE.finditer(txt):
                v = _norm_ph(m.group(1), m.group(2))
                _push_candidate(ph_cands, v, src, "ph_explicit_range", 0.95, m.group(0))
            for m in PH_BETWEEN_RE.finditer(txt):
                v = _norm_ph(m.group(1), m.group(2))
                _push_candidate(ph_cands, v, src, "ph_explicit_between", 0.94, m.group(0))
            for m in PH_SINGLE_RE.finditer(txt):
                v = _norm_ph(m.group(1))
                _push_candidate(ph_cands, v, src, "ph_explicit_single", 0.93, m.group(0))

        if any(t in low for t in ["temp", "°", " c", " f", " k", "room temperature", "ambient"]):
            for m in TEMP_RANGE_RE.finditer(txt):
                v = _norm_temp(m.group(1), m.group(3), m.group(2))
                _push_candidate(temp_cands, v, src, "temp_explicit_range", 0.95, m.group(0))
            for m in TEMP_SINGLE_RE.finditer(txt):
                v = _norm_temp(m.group(1), m.group(2))
                _push_candidate(temp_cands, v, src, "temp_explicit_single", 0.93, m.group(0))
            if ROOM_TEMP_RE.search(txt):
                _push_candidate(temp_cands, "25", src, "temp_room_temperature", 0.82, "room temperature")
            if ICE_TEMP_RE.search(txt):
                _push_candidate(temp_cands, "0", src, "temp_on_ice", 0.85, "ice")

        for name, cre in SYSTEM_REGEX.items():
            if cre.search(txt):
                system_categories.add(name)
                context_source_fields.add(src)

        for m in BUFFER_RE.finditer(txt):
            b = _normalize_buffer(m.group(1))
            if b:
                buffers.add(b)
                context_source_fields.add(src)

        for cre in KNOWN_CELL_LINE_RE:
            for m in cre.finditer(txt):
                token = _compact_text(m.group(1)).upper()
                if token:
                    cell_lines.add(token)
                    context_source_fields.add(src)

        for m in GENERIC_CELL_RE.finditer(txt):
            token = _compact_text(m.group(1))
            if not token:
                continue
            tok = token.upper()
            if token.lower() in CELL_STOPWORDS:
                continue
            if len(tok) <= 20:
                cell_lines.add(tok)
                context_source_fields.add(src)

    # pH inference
    if not ph_cands:
        inferred_from_buffers = 0
        for b in sorted(buffers):
            val = BUFFER_PH_DEFAULT.get(b)
            if val:
                _push_candidate(ph_cands, val, "assay_description", "ph_inferred_buffer_default", 0.62, b)
                inferred_from_buffers += 1
        if inferred_from_buffers > 0:
            rule_hits.append("ph_inferred_buffer_default")
        elif (cell_lines or "cell_based" in system_categories) and (INCUBATION_RE.search(joined_text) or MEDIA_RE.search(joined_text)):
            _push_candidate(ph_cands, "7.4", "assay_context", "ph_inferred_cell_culture", 0.48, "cell culture")
            rule_hits.append("ph_inferred_cell_culture")

    # temperature inference
    if not temp_cands:
        if ICE_TEMP_RE.search(joined_text):
            _push_candidate(temp_cands, "0", "activity_comment", "temp_inferred_ice", 0.76, "ice")
            rule_hits.append("temp_inferred_ice")
        elif (cell_lines or "cell_based" in system_categories or MEDIA_RE.search(joined_text)) and (
            INCUBATION_RE.search(joined_text) or MEDIA_RE.search(joined_text)
        ):
            _push_candidate(temp_cands, "37", "assay_context", "temp_inferred_cell_culture", 0.55, "cell culture")
            rule_hits.append("temp_inferred_cell_culture")

    ph_value, ph_best, ph_conflict, ph_ranked = _resolve(ph_cands)
    temp_value, temp_best, temp_conflict, temp_ranked = _resolve(temp_cands)

    if ph_best:
        rule_hits.append(ph_best.method)
    if temp_best:
        rule_hits.append(temp_best.method)

    if cell_lines:
        rule_hits.append("system_cell_line")
    if buffers:
        rule_hits.append("system_buffer")

    # Keep legacy context/system if nothing parsed
    legacy_context = _compact_text(row.get("condition_context", ""))
    legacy_system = _compact_text(row.get("condition_system", ""))

    system_list = sorted(system_categories)
    cell_list = sorted(cell_lines)[:6]
    buffer_list = sorted(buffers)[:6]

    context_parts: List[str] = []
    if ph_value:
        context_parts.append(f"pH={ph_value}")
    if temp_value:
        context_parts.append(f"temperature_c={temp_value}")
    if system_list:
        context_parts.append(f"system={','.join(system_list)}")
    if cell_list:
        context_parts.append(f"cell_line={','.join(cell_list)}")
    if buffer_list:
        context_parts.append(f"buffer={','.join(buffer_list)}")

    condition_system = ";".join(system_list)
    if not condition_system and legacy_system:
        condition_system = legacy_system

    condition_context = "; ".join(context_parts)
    if not condition_context and legacy_context:
        condition_context = legacy_context

    source_fields = set(context_source_fields)
    if ph_best:
        source_fields.add(ph_best.source_field)
    if temp_best:
        source_fields.add(temp_best.source_field)
    if not source_fields and legacy_context:
        source_fields.add("condition_context")

    confidence_parts: List[float] = []
    if ph_best:
        confidence_parts.append(ph_best.confidence)
    if temp_best:
        confidence_parts.append(temp_best.confidence)
    if system_list or cell_list or buffer_list:
        confidence_parts.append(0.72 if (cell_list or buffer_list) else 0.56)
    if not confidence_parts and legacy_context:
        confidence_parts.append(0.4)

    conf = ""
    if confidence_parts:
        conf = _fmt_num(sum(confidence_parts) / len(confidence_parts))

    context_payload: Dict[str, object] = {}
    if ph_value:
        context_payload["pH"] = ph_value
    if temp_value:
        context_payload["temperature_c"] = temp_value
    if system_list:
        context_payload["system"] = system_list
    if cell_list:
        context_payload["cell_line"] = cell_list
    if buffer_list:
        context_payload["buffer"] = buffer_list
    if source_fields:
        context_payload["source_fields"] = sorted(source_fields)
    if legacy_context and legacy_context != condition_context:
        context_payload["legacy_condition_context"] = legacy_context[:220]

    context_json = json.dumps(context_payload, ensure_ascii=False, separators=(",", ":"), sort_keys=True) if context_payload else ""

    conflict_fields = []
    if ph_conflict:
        conflict_fields.append("condition_pH")
    if temp_conflict:
        conflict_fields.append("condition_temperature_c")
    conflict_flag = "true" if conflict_fields else "false"

    candidates_payload = {
        "condition_pH": [c.__dict__ for c in ph_ranked],
        "condition_temperature_c": [c.__dict__ for c in temp_ranked],
        "system": {
            "system": system_list,
            "cell_line": cell_list,
            "buffer": buffer_list,
        },
    }

    source_snapshot = {k: texts.get(k, "")[:280] for k in SOURCE_FIELDS if texts.get(k, "")}

    return {
        "condition_pH": ph_value,
        "condition_temperature_c": temp_value,
        "condition_system": condition_system,
        "condition_context": condition_context,
        "condition_context_json": context_json,
        "condition_extract_confidence": conf,
        "condition_extract_source_field": ";".join(sorted(source_fields)),
        "conflict_flag": conflict_flag,
        "conflict_fields": ";".join(conflict_fields),
        "_audit_candidates_json": json.dumps(candidates_payload, ensure_ascii=False, separators=(",", ":"), sort_keys=True),
        "_source_snapshot_json": json.dumps(source_snapshot, ensure_ascii=False, separators=(",", ":"), sort_keys=True),
        "_rule_hits": ";".join(sorted(set(rule_hits))),
    }
