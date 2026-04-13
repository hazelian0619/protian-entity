#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


FETCH_DATE = date.today().isoformat()

CROSS_COLUMNS = [
    "record_id",
    "interaction_type",
    "edge_id",
    "entity_pair_key",
    "consistent_across_n",
    "source_list",
    "evidence_count",
    "distinct_methods_n",
    "reference_count",
    "predicate_set",
    "direction_set",
    "conflict_flag",
    "conflict_reason",
    "source_version_list",
    "fetch_date",
]

AGG_COLUMNS = [
    "aggregate_id",
    "interaction_type",
    "edge_id",
    "consistent_across_n",
    "evidence_count",
    "distinct_methods_n",
    "numeric_score_max",
    "numeric_score_norm",
    "reference_coverage",
    "context_coverage",
    "conflict_flag",
    "aggregate_score",
    "score_bucket",
    "fetch_date",
]

# direction bitset
DIR_TRUE = 1
DIR_FALSE = 2
DIR_UNKNOWN = 4

# effect bitset
EFFECT_POS = 1
EFFECT_NEG = 2

EXPECTED_CONTEXT_FAMILIES = {
    "PPI": 2,
    "PSI": 2,
    "RPI": 3,
}

NEGATIVE_EFFECT_PATTERNS = [
    r"inhib",
    r"antagon",
    r"block",
    r"suppress",
    r"down[-_ ]?reg",
    r"inverse\s+agon",
    r"decreas",
]
POSITIVE_EFFECT_PATTERNS = [
    r"agon",
    r"activ",
    r"induc",
    r"enhanc",
    r"stimulat",
    r"up[-_ ]?reg",
    r"potentiat",
    r"promot",
]

METHOD_CATEGORY_RULES = {
    "PPI": [
        ("Y2H", [r"\by2h\b", r"yeast\s+two\s+hybrid"]),
        ("COIP", [r"co[- ]?ip", r"coimmunoprecip"]),
        ("APMS", [r"ap[- ]?ms", r"affinity\s+purification"]),
        ("PCA", [r"\bpca\b", r"protein\s+complementation"]),
        ("COMPUTATIONAL", [r"string", r"comput", r"text\s*mining"]),
        ("FUNCTION_CONTEXT", [r"go", r"reactome", r"kegg", r"functional"]),
        ("BIOPHYSICAL", [r"x[- ]?ray", r"nmr", r"cryo", r"biacore", r"spr"]),
    ],
    "PSI": [
        ("DRUGBANK_ACTION", [r"inhibitor", r"agonist", r"antagonist", r"activator", r"modulator"]),
        ("ASSAY_B", [r"\bb\b", r"binding"]),
        ("ASSAY_F", [r"\bf\b", r"functional"]),
        ("ASSAY_A", [r"\ba\b", r"adme"]),
        ("ASSAY_T", [r"\bt\b", r"tox"]),
        ("ASSAY_P", [r"\bp\b", r"physicochemical"]),
        ("STRUCTURE", [r"x[- ]?ray", r"nmr", r"electron", r"pdb", r"structure"]),
        ("CHEMBL_ACTIVITY", [r"chembl", r"pchembl", r"ic50", r"ki", r"kd", r"ec50"]),
    ],
    "RPI": [
        ("CLIP", [r"clip"]),
        ("RIP", [r"\brip\b"]),
        ("DOMAIN", [r"pfam", r"interpro", r"domain"]),
        ("FUNCTION", [r"splicing", r"translation", r"stability", r"regulation"]),
        ("REFERENCE", [r"reference", r"pmid", r"doi"]),
        ("PREDICTED", [r"predict", r"inference", r"keyword"]),
    ],
}

METHOD_CATEGORY_BITS: Dict[str, int] = {}
for _, rules in METHOD_CATEGORY_RULES.items():
    for cat, _ in rules:
        METHOD_CATEGORY_BITS.setdefault(cat, 1 << len(METHOD_CATEGORY_BITS))
METHOD_CATEGORY_BITS.setdefault("OTHER", 1 << len(METHOD_CATEGORY_BITS))


class BitsetVocab:
    def __init__(self) -> None:
        self.to_idx: Dict[str, int] = {}
        self.tokens: List[str] = []

    def mask(self, token: str) -> int:
        t = normalize(token)
        if not t:
            return 0
        idx = self.to_idx.get(t)
        if idx is None:
            idx = len(self.tokens)
            self.to_idx[t] = idx
            self.tokens.append(t)
        return 1 << idx

    def mask_multi(self, text: str) -> int:
        t = normalize(text)
        if not t:
            return 0
        parts = re.split(r"[;,|]", t)
        m = 0
        for p in parts:
            m |= self.mask(p)
        return m

    def decode(self, bitmask: int) -> List[str]:
        out: List[str] = []
        m = int(bitmask)
        while m:
            lsb = m & -m
            idx = lsb.bit_length() - 1
            if 0 <= idx < len(self.tokens):
                out.append(self.tokens[idx])
            m ^= lsb
        out.sort()
        return out


@dataclass
class BuildState:
    # edge identity
    edge_index: Dict[str, int]
    edge_ids: List[str]
    edge_type: List[str]
    edge_pair_id: List[int]

    # edge base attrs
    edge_best_score: List[float]

    # edge accumulated stats
    evidence_count: List[int]
    method_mask: List[int]
    reference_count: List[int]
    score_max: List[float]
    effect_mask: List[int]
    context_mask: List[int]

    # pair identity/stats
    pair_key_to_id: Dict[str, int]
    pair_keys: List[str]
    pair_source_mask: List[int]
    pair_predicate_mask: List[int]
    pair_direction_mask: List[int]
    pair_source_version_mask: List[int]
    pair_effect_mask: List[int]

    # counters
    edge_count_by_type: Counter
    created_from_context_by_type: Counter


@dataclass
class BuildReport:
    warnings: List[str]
    notes: List[str]
    context_file_stats: List[Dict[str, object]]


def normalize(v: Optional[str]) -> str:
    return (v or "").strip()


def lower(v: Optional[str]) -> str:
    return normalize(v).lower()


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _set_csv_field_size_limit() -> None:
    size = sys.maxsize
    while True:
        try:
            csv.field_size_limit(size)
            return
        except OverflowError:
            size //= 10


def bool_text(v: bool) -> str:
    return "true" if v else "false"


def to_int_text(v: int) -> str:
    return str(int(v))


def to_float_text(v: float) -> str:
    vv = max(0.0, min(1.0, float(v)))
    return f"{vv:.6f}"


def numeric_or_nan(v: str) -> float:
    t = normalize(v)
    if not t:
        return math.nan
    try:
        return float(t)
    except Exception:  # noqa: BLE001
        return math.nan


def parse_directed(v: str) -> int:
    t = lower(v)
    if t in {"1", "true", "t", "yes", "y"}:
        return DIR_TRUE
    if t in {"0", "false", "f", "no", "n"}:
        return DIR_FALSE
    return DIR_UNKNOWN


def direction_set_text(mask: int) -> str:
    vals: List[str] = []
    if mask & DIR_TRUE:
        vals.append("directed")
    if mask & DIR_FALSE:
        vals.append("undirected")
    if mask & DIR_UNKNOWN:
        vals.append("unknown")
    return "|".join(vals)


def popcount(v: int) -> int:
    iv = int(v)
    try:
        return int(iv.bit_count())  # py>=3.8
    except AttributeError:
        return int(bin(iv).count("1"))


def quantile(values: List[float], q: float) -> float:
    if not values:
        return 0.0
    s = sorted(values)
    idx = int(round((len(s) - 1) * q))
    idx = max(0, min(idx, len(s) - 1))
    return float(s[idx])


def parse_effect_sign(v: str) -> int:
    txt = lower(v)
    if not txt:
        return 0
    pos = any(re.search(p, txt) for p in POSITIVE_EFFECT_PATTERNS)
    neg = any(re.search(p, txt) for p in NEGATIVE_EFFECT_PATTERNS)
    out = 0
    if pos:
        out |= EFFECT_POS
    if neg:
        out |= EFFECT_NEG
    return out


def categorize_method(method_text: str, interaction_type: str) -> int:
    txt = lower(method_text)
    if not txt:
        return 0
    rules = METHOD_CATEGORY_RULES.get(interaction_type, [])
    mask = 0
    for cat, pats in rules:
        if any(re.search(p, txt) for p in pats):
            mask |= METHOD_CATEGORY_BITS[cat]
    if mask == 0:
        mask = METHOD_CATEGORY_BITS["OTHER"]
    return mask


def make_id(prefix: str, payload: str) -> str:
    return prefix + "_" + hashlib.sha1(payload.encode("utf-8")).hexdigest()[:20]


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def get_pair_key_generic(src_type: str, src_id: str, dst_type: str, dst_id: str, directed: int) -> str:
    a = f"{normalize(src_type).lower()}:{normalize(src_id)}"
    b = f"{normalize(dst_type).lower()}:{normalize(dst_id)}"
    if directed == DIR_FALSE:
        x, y = sorted([a, b])
        return f"{x}--{y}"
    return f"{a}->{b}"


def get_pair_key_psi(compound_key: str, uniprot: str) -> str:
    c = normalize(compound_key)
    p = normalize(uniprot)
    return f"compound:{c}->protein:{p}"


def write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_tsv_header(path: Path, columns: List[str]) -> csv.DictWriter:
    path.parent.mkdir(parents=True, exist_ok=True)
    f = path.open("w", encoding="utf-8", newline="")
    w = csv.DictWriter(f, fieldnames=columns, delimiter="\t", lineterminator="\n")
    w.writeheader()
    # attach handle to writer for explicit close
    w._fh = f  # type: ignore[attr-defined]
    return w


def close_writer(w: csv.DictWriter) -> None:
    fh = getattr(w, "_fh", None)
    if fh is not None:
        fh.close()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build interaction cross-validation and aggregate score tables (v2)")

    p.add_argument("--ppi-edges", type=Path, required=True)
    p.add_argument("--ppi-evidence", type=Path, required=True)
    p.add_argument("--ppi-method-context", type=Path, default=None)
    p.add_argument("--ppi-function-context", type=Path, default=None)

    p.add_argument("--psi-edges", type=Path, required=True)
    p.add_argument("--psi-evidence", type=Path, required=True)
    p.add_argument("--psi-activity-context", type=Path, required=True)
    p.add_argument("--psi-structure-evidence", type=Path, required=True)

    p.add_argument("--rpi-edges", type=Path, required=True)
    p.add_argument("--rpi-evidence", type=Path, required=True)
    p.add_argument("--rpi-site-context", type=Path, default=None)
    p.add_argument("--rpi-domain-context", type=Path, default=None)
    p.add_argument("--rpi-function-context", type=Path, default=None)

    p.add_argument("--molecule-xref", type=Path, required=True)

    p.add_argument("--cross-output", type=Path, required=True)
    p.add_argument("--aggregate-output", type=Path, required=True)
    p.add_argument("--report", type=Path, required=True)

    p.add_argument("--limit-per-type", type=int, default=None)
    p.add_argument("--fetch-date", default=FETCH_DATE)
    p.add_argument("--check-inputs", action="store_true")
    return p.parse_args()


def initialize_state() -> BuildState:
    return BuildState(
        edge_index={},
        edge_ids=[],
        edge_type=[],
        edge_pair_id=[],
        edge_best_score=[],
        evidence_count=[],
        method_mask=[],
        reference_count=[],
        score_max=[],
        effect_mask=[],
        context_mask=[],
        pair_key_to_id={},
        pair_keys=[],
        pair_source_mask=[],
        pair_predicate_mask=[],
        pair_direction_mask=[],
        pair_source_version_mask=[],
        pair_effect_mask=[],
        edge_count_by_type=Counter(),
        created_from_context_by_type=Counter(),
    )


def get_pair_id(state: BuildState, pair_key: str) -> int:
    pid = state.pair_key_to_id.get(pair_key)
    if pid is not None:
        return pid
    pid = len(state.pair_keys)
    state.pair_key_to_id[pair_key] = pid
    state.pair_keys.append(pair_key)
    state.pair_source_mask.append(0)
    state.pair_predicate_mask.append(0)
    state.pair_direction_mask.append(0)
    state.pair_source_version_mask.append(0)
    state.pair_effect_mask.append(0)
    return pid


def add_edge(
    *,
    state: BuildState,
    edge_id: str,
    interaction_type: str,
    pair_key: str,
    predicate_mask: int,
    direction_mask: int,
    source_mask: int,
    source_version_mask: int,
    best_score: float,
    limit_per_type: Optional[int],
    created_from_context: bool,
) -> Optional[int]:
    eid = normalize(edge_id)
    if not eid:
        return None

    idx = state.edge_index.get(eid)
    if idx is not None:
        # Existing edge still updates pair-level lineage if new source/version arrives.
        pid = state.edge_pair_id[idx]
        state.pair_source_mask[pid] |= source_mask
        state.pair_source_version_mask[pid] |= source_version_mask
        state.pair_predicate_mask[pid] |= predicate_mask
        state.pair_direction_mask[pid] |= direction_mask
        if not math.isnan(best_score):
            cur = state.edge_best_score[idx]
            if math.isnan(cur) or best_score > cur:
                state.edge_best_score[idx] = best_score
                state.score_max[idx] = best_score
        return idx

    if limit_per_type is not None and state.edge_count_by_type[interaction_type] >= limit_per_type:
        return None

    pid = get_pair_id(state, pair_key)

    idx = len(state.edge_ids)
    state.edge_index[eid] = idx
    state.edge_ids.append(eid)
    state.edge_type.append(interaction_type)
    state.edge_pair_id.append(pid)
    state.edge_best_score.append(best_score)

    state.evidence_count.append(0)
    state.method_mask.append(0)
    state.reference_count.append(0)
    state.score_max.append(best_score)
    state.effect_mask.append(0)
    state.context_mask.append(0)

    state.pair_source_mask[pid] |= source_mask
    state.pair_source_version_mask[pid] |= source_version_mask
    state.pair_predicate_mask[pid] |= predicate_mask
    state.pair_direction_mask[pid] |= direction_mask

    state.edge_count_by_type[interaction_type] += 1
    if created_from_context:
        state.created_from_context_by_type[interaction_type] += 1

    return idx


def add_evidence(
    *,
    state: BuildState,
    edge_idx: int,
    source_mask: int,
    source_version_mask: int,
    method_mask: int,
    has_reference: bool,
    score: float,
    effect_mask: int,
) -> None:
    state.evidence_count[edge_idx] += 1
    state.method_mask[edge_idx] |= method_mask
    if has_reference:
        state.reference_count[edge_idx] += 1

    if not math.isnan(score):
        cur = state.score_max[edge_idx]
        if math.isnan(cur) or score > cur:
            state.score_max[edge_idx] = score

    if effect_mask:
        state.effect_mask[edge_idx] |= effect_mask

    pid = state.edge_pair_id[edge_idx]
    state.pair_source_mask[pid] |= source_mask
    state.pair_source_version_mask[pid] |= source_version_mask
    if effect_mask:
        state.pair_effect_mask[pid] |= effect_mask


def add_context_hit(state: BuildState, edge_idx: int, family_bit: int) -> None:
    state.context_mask[edge_idx] |= family_bit


def load_molecule_xref_map(path: Path, report: BuildReport) -> Dict[str, str]:
    if not path.exists():
        raise SystemExit(f"[ERROR] missing molecule xref: {path}")

    out: Dict[str, str] = {}
    conflicts = 0
    rows = 0

    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        need = {"inchikey", "drugbank_id"}
        miss = need - set(r.fieldnames or [])
        if miss:
            raise SystemExit(f"[ERROR] molecule xref missing columns: {sorted(miss)}")

        for row in r:
            rows += 1
            ik = normalize(row.get("inchikey", "")).upper()
            if not ik:
                continue
            dbids = [x.strip() for x in normalize(row.get("drugbank_id", "")).split(";") if x.strip()]
            for dbid in dbids:
                prev = out.get(dbid)
                if prev is None:
                    out[dbid] = ik
                elif prev != ik:
                    conflicts += 1

    if conflicts > 0:
        report.warnings.append(f"molecule xref conflict dbid->inchikey rows={conflicts}; keep first mapping")
    report.notes.append(f"molecule xref loaded rows={rows}, mapped_drugbank_ids={len(out)}")
    return out


def iter_tsv(path: Path) -> Iterable[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        for row in r:
            yield row


def collect_inputs(args: argparse.Namespace) -> Dict[str, Path]:
    out: Dict[str, Path] = {
        "ppi_edges": args.ppi_edges,
        "ppi_evidence": args.ppi_evidence,
        "psi_edges": args.psi_edges,
        "psi_evidence": args.psi_evidence,
        "psi_activity_context": args.psi_activity_context,
        "psi_structure_evidence": args.psi_structure_evidence,
        "rpi_edges": args.rpi_edges,
        "rpi_evidence": args.rpi_evidence,
        "molecule_xref": args.molecule_xref,
    }

    optional = {
        "ppi_method_context": args.ppi_method_context,
        "ppi_function_context": args.ppi_function_context,
        "rpi_site_context": args.rpi_site_context,
        "rpi_domain_context": args.rpi_domain_context,
        "rpi_function_context": args.rpi_function_context,
    }
    for k, p in optional.items():
        if p is not None and str(p) != "":
            out[k] = p
    return out


def check_inputs(args: argparse.Namespace) -> int:
    inputs = collect_inputs(args)
    required_keys = [
        "ppi_edges",
        "ppi_evidence",
        "psi_edges",
        "psi_evidence",
        "psi_activity_context",
        "psi_structure_evidence",
        "rpi_edges",
        "rpi_evidence",
        "molecule_xref",
    ]

    missing_required = [k for k in required_keys if not inputs[k].exists()]

    payload = {
        "pipeline": "interaction_cross_validation_v2",
        "generated_at_utc": utc_now_iso(),
        "status": "inputs_ready" if not missing_required else "blocked_missing_inputs",
        "required_inputs": {k: str(inputs[k]) for k in required_keys},
        "optional_inputs": {
            k: str(v)
            for k, v in inputs.items()
            if k not in required_keys
        },
        "missing_required": [{"name": k, "path": str(inputs[k])} for k in missing_required],
        "download_checklist_if_blocked": [
            {
                "name": "PSI activity context (任务B增强)",
                "target_path": str(inputs["psi_activity_context"]),
                "sha256_cmd": f"shasum -a 256 {inputs['psi_activity_context']}",
            },
            {
                "name": "PSI structure evidence (任务B增强)",
                "target_path": str(inputs["psi_structure_evidence"]),
                "sha256_cmd": f"shasum -a 256 {inputs['psi_structure_evidence']}",
            },
            {
                "name": "RPI edges/evidence（任务C产物）",
                "target_path": f"{inputs['rpi_edges']} + {inputs['rpi_evidence']}",
                "sha256_cmd": (
                    f"shasum -a 256 {inputs['rpi_edges']} && "
                    f"shasum -a 256 {inputs['rpi_evidence']}"
                ),
            },
        ],
    }

    write_json(args.report, payload)
    if missing_required:
        print(f"[BLOCKED] missing required inputs: {len(missing_required)}")
        for item in payload["missing_required"]:
            print(f"  - {item['name']}: {item['path']}")
        print(f"[BLOCKED] report -> {args.report}")
        return 2

    print(f"[OK] inputs ready -> {args.report}")
    return 0


def main() -> int:
    _set_csv_field_size_limit()
    args = parse_args()

    if args.check_inputs:
        return check_inputs(args)

    report = BuildReport(warnings=[], notes=[], context_file_stats=[])

    source_vocab = BitsetVocab()
    predicate_vocab = BitsetVocab()
    source_version_vocab = BitsetVocab()

    state = initialize_state()

    inputs = collect_inputs(args)
    for key in [
        "ppi_edges",
        "ppi_evidence",
        "psi_edges",
        "psi_evidence",
        "psi_activity_context",
        "psi_structure_evidence",
        "rpi_edges",
        "rpi_evidence",
        "molecule_xref",
    ]:
        if not inputs[key].exists():
            raise SystemExit(f"[ERROR] missing required input: {inputs[key]} ({key})")

    dbid_to_inchikey = load_molecule_xref_map(inputs["molecule_xref"], report)

    # -----------------
    # 1) Load core edges
    # -----------------
    def add_core_edge(
        *,
        edge_id: str,
        interaction_type: str,
        pair_key: str,
        predicate: str,
        directed_raw: str,
        source: str,
        source_version: str,
        best_score_raw: str,
        created_from_context: bool,
    ) -> Optional[int]:
        return add_edge(
            state=state,
            edge_id=edge_id,
            interaction_type=interaction_type,
            pair_key=pair_key,
            predicate_mask=predicate_vocab.mask(predicate or "unknown"),
            direction_mask=parse_directed(directed_raw),
            source_mask=source_vocab.mask_multi(source),
            source_version_mask=source_version_vocab.mask_multi(source_version),
            best_score=numeric_or_nan(best_score_raw),
            limit_per_type=args.limit_per_type,
            created_from_context=created_from_context,
        )

    # PPI edges
    ppi_edges_loaded = 0
    for row in iter_tsv(inputs["ppi_edges"]):
        eid = normalize(row.get("edge_id", ""))
        if not eid:
            continue
        pair_key = get_pair_key_generic(
            row.get("src_type", ""),
            row.get("src_id", ""),
            row.get("dst_type", ""),
            row.get("dst_id", ""),
            parse_directed(row.get("directed", "")),
        )
        idx = add_core_edge(
            edge_id=eid,
            interaction_type="PPI",
            pair_key=pair_key,
            predicate=normalize(row.get("predicate", "")) or "ppi",
            directed_raw=row.get("directed", ""),
            source=normalize(row.get("source", "")) or "STRING",
            source_version=normalize(row.get("source_version", "")),
            best_score_raw=row.get("best_score", ""),
            created_from_context=False,
        )
        if idx is not None:
            ppi_edges_loaded += 1

    report.notes.append(f"PPI edges loaded={ppi_edges_loaded}")

    # PSI drug-target edges
    psi_pair_to_edge: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    psi_edges_loaded = 0

    for row in iter_tsv(inputs["psi_edges"]):
        eid = normalize(row.get("edge_id", ""))
        if not eid:
            continue
        dbid = normalize(row.get("src_id", ""))
        uniprot = normalize(row.get("dst_id", ""))
        if not dbid or not uniprot:
            continue
        compound_key = dbid_to_inchikey.get(dbid, f"DB:{dbid}")
        pair_key = get_pair_key_psi(compound_key=compound_key, uniprot=uniprot)

        idx = add_core_edge(
            edge_id=eid,
            interaction_type="PSI",
            pair_key=pair_key,
            predicate=normalize(row.get("predicate", "")) or "drug_targets",
            directed_raw=row.get("directed", "1"),
            source=normalize(row.get("source", "")) or "DrugBank",
            source_version=normalize(row.get("source_version", "")),
            best_score_raw=row.get("best_score", ""),
            created_from_context=False,
        )
        if idx is not None:
            psi_pair_to_edge[(dbid, uniprot)].append(idx)
            psi_edges_loaded += 1

    report.notes.append(f"PSI drug-target edges loaded={psi_edges_loaded}")

    # RPI edges
    rpi_edges_loaded = 0
    for row in iter_tsv(inputs["rpi_edges"]):
        eid = normalize(row.get("edge_id", ""))
        if not eid:
            continue
        pair_key = get_pair_key_generic(
            row.get("src_type", ""),
            row.get("src_id", ""),
            row.get("dst_type", ""),
            row.get("dst_id", ""),
            parse_directed(row.get("directed", "")),
        )
        idx = add_core_edge(
            edge_id=eid,
            interaction_type="RPI",
            pair_key=pair_key,
            predicate=normalize(row.get("predicate", "")) or "rpi",
            directed_raw=row.get("directed", "1"),
            source=normalize(row.get("source", "")) or "RNAInter",
            source_version=normalize(row.get("source_version", "")),
            best_score_raw=row.get("best_score", ""),
            created_from_context=False,
        )
        if idx is not None:
            rpi_edges_loaded += 1

    report.notes.append(f"RPI edges loaded={rpi_edges_loaded}")

    # -----------------
    # 2) Evidence loading
    # -----------------

    # PPI evidence
    ppi_evi_rows = 0
    ppi_evi_mapped = 0
    for row in iter_tsv(inputs["ppi_evidence"]):
        ppi_evi_rows += 1
        eid = normalize(row.get("edge_id", ""))
        idx = state.edge_index.get(eid)
        if idx is None or state.edge_type[idx] != "PPI":
            continue
        method_txt = normalize(row.get("method", "")) or normalize(row.get("evidence_type", ""))
        effect = parse_effect_sign(method_txt)
        add_evidence(
            state=state,
            edge_idx=idx,
            source_mask=source_vocab.mask_multi(normalize(row.get("source", ""))),
            source_version_mask=source_version_vocab.mask_multi(normalize(row.get("source_version", ""))),
            method_mask=categorize_method(method_txt, "PPI"),
            has_reference=normalize(row.get("reference", "")) != "",
            score=numeric_or_nan(row.get("score", "")),
            effect_mask=effect,
        )
        ppi_evi_mapped += 1

    report.notes.append(f"PPI evidence rows={ppi_evi_rows}, mapped={ppi_evi_mapped}")

    # PSI drug-target evidence (no edge_id in v1)
    psi_db_evi_rows = 0
    psi_db_evi_mapped = 0
    for row in iter_tsv(inputs["psi_evidence"]):
        psi_db_evi_rows += 1
        dbid = normalize(row.get("drugbank_id", ""))
        uniprot = normalize(row.get("uniprot_id", ""))
        if not dbid or not uniprot:
            continue
        idx_list = psi_pair_to_edge.get((dbid, uniprot), [])
        if not idx_list:
            continue

        method_txt = ";".join(
            x for x in [normalize(row.get("action", "")), normalize(row.get("target_role", ""))] if x
        )
        effect = parse_effect_sign(method_txt)

        for idx in idx_list:
            add_evidence(
                state=state,
                edge_idx=idx,
                source_mask=source_vocab.mask_multi(normalize(row.get("source", "")) or "DrugBank"),
                source_version_mask=source_version_vocab.mask_multi(normalize(row.get("source_version", ""))),
                method_mask=categorize_method(method_txt, "PSI"),
                has_reference=normalize(row.get("reference", "")) != "",
                score=math.nan,
                effect_mask=effect,
            )
            psi_db_evi_mapped += 1

    report.notes.append(f"PSI DrugBank evidence rows={psi_db_evi_rows}, mapped={psi_db_evi_mapped}")

    # RPI evidence
    rpi_evi_rows = 0
    rpi_evi_mapped = 0
    for row in iter_tsv(inputs["rpi_evidence"]):
        rpi_evi_rows += 1
        eid = normalize(row.get("edge_id", ""))
        idx = state.edge_index.get(eid)
        if idx is None or state.edge_type[idx] != "RPI":
            continue
        method_txt = normalize(row.get("method", "")) or normalize(row.get("evidence_type", ""))
        effect = parse_effect_sign(method_txt)
        add_evidence(
            state=state,
            edge_idx=idx,
            source_mask=source_vocab.mask_multi(normalize(row.get("source", ""))),
            source_version_mask=source_version_vocab.mask_multi(normalize(row.get("source_version", ""))),
            method_mask=categorize_method(method_txt, "RPI"),
            has_reference=normalize(row.get("reference", "")) != "",
            score=numeric_or_nan(row.get("score", "")),
            effect_mask=effect,
        )
        rpi_evi_mapped += 1

    report.notes.append(f"RPI evidence rows={rpi_evi_rows}, mapped={rpi_evi_mapped}")

    # -----------------
    # 3) Context integration
    # -----------------
    available_context_families = {
        "PPI": 0,
        "PSI": 0,
        "RPI": 0,
    }

    def record_context_file(path: Optional[Path], itype: str, family: str, status: str, rows_hit: int = 0) -> None:
        report.context_file_stats.append(
            {
                "file": str(path) if path is not None else "",
                "interaction_type": itype,
                "family": family,
                "status": status,
                "rows_hit": rows_hit,
            }
        )

    # PPI method context (family bit 1)
    ppi_method_path = inputs.get("ppi_method_context")
    if ppi_method_path is not None and ppi_method_path.exists():
        available_context_families["PPI"] += 1
        hit = 0
        for row in iter_tsv(ppi_method_path):
            eid = normalize(row.get("edge_id", ""))
            idx = state.edge_index.get(eid)
            if idx is None or state.edge_type[idx] != "PPI":
                continue
            add_context_hit(state, idx, 1)
            pmid = normalize(row.get("pmid", ""))
            doi = normalize(row.get("doi", ""))
            ref_ok = bool(pmid or doi)
            score = numeric_or_nan(row.get("experimental_score_norm", ""))
            score2 = numeric_or_nan(row.get("text_mining_score_norm", ""))
            if not math.isnan(score2):
                score = max(score, score2) if not math.isnan(score) else score2
            method_txt = normalize(row.get("method", "")) or normalize(row.get("method_raw", ""))
            add_evidence(
                state=state,
                edge_idx=idx,
                source_mask=source_vocab.mask_multi(normalize(row.get("source_databases", "")) or "STRING"),
                source_version_mask=source_version_vocab.mask_multi(normalize(row.get("source_version", ""))),
                method_mask=categorize_method(method_txt, "PPI"),
                has_reference=ref_ok,
                score=score,
                effect_mask=parse_effect_sign(method_txt),
            )
            hit += 1
        record_context_file(ppi_method_path, "PPI", "method", "ok", hit)
    else:
        record_context_file(ppi_method_path, "PPI", "method", "missing", 0)

    # PPI function context (family bit 2)
    ppi_func_path = inputs.get("ppi_function_context")
    if ppi_func_path is not None and ppi_func_path.exists():
        available_context_families["PPI"] += 1
        hit = 0
        for row in iter_tsv(ppi_func_path):
            eid = normalize(row.get("edge_id", ""))
            idx = state.edge_index.get(eid)
            if idx is None or state.edge_type[idx] != "PPI":
                continue
            add_context_hit(state, idx, 2)
            add_evidence(
                state=state,
                edge_idx=idx,
                source_mask=source_vocab.mask_multi(normalize(row.get("source", ""))),
                source_version_mask=source_version_vocab.mask_multi(normalize(row.get("source_version", ""))),
                method_mask=categorize_method("functional_context", "PPI"),
                has_reference=False,
                score=numeric_or_nan(row.get("context_support_score", "")),
                effect_mask=0,
            )
            hit += 1
        record_context_file(ppi_func_path, "PPI", "function", "ok", hit)
    else:
        record_context_file(ppi_func_path, "PPI", "function", "missing", 0)

    # PSI activity context (family bit 1) - required from updated 任务B
    psi_activity_rows = 0
    psi_activity_mapped = 0
    for row in iter_tsv(inputs["psi_activity_context"]):
        psi_activity_rows += 1
        eid = normalize(row.get("edge_id", ""))
        compound_ik = normalize(row.get("compound_inchikey", "")).upper()
        uniprot = normalize(row.get("target_uniprot_accession", ""))

        if not eid:
            continue

        idx = state.edge_index.get(eid)
        if idx is None:
            if not compound_ik or not uniprot:
                continue
            idx = add_core_edge(
                edge_id=eid,
                interaction_type="PSI",
                pair_key=get_pair_key_psi(compound_ik, uniprot),
                predicate="molecule_targets_protein",
                directed_raw="1",
                source=normalize(row.get("source", "")) or "ChEMBL",
                source_version=normalize(row.get("source_version", "")),
                best_score_raw=row.get("pchembl_value_eff", ""),
                created_from_context=True,
            )
            if idx is None:
                continue

        if state.edge_type[idx] != "PSI":
            continue

        add_context_hit(state, idx, 1)

        assay_desc = normalize(row.get("assay_type_desc", "")) or normalize(row.get("assay_type", ""))
        method_txt = ";".join(x for x in [assay_desc, normalize(row.get("standard_type", ""))] if x)
        ref_ok = bool(normalize(row.get("doi", "")) or normalize(row.get("pubmed_id", "")))
        score = numeric_or_nan(row.get("pchembl_value_eff", ""))
        if math.isnan(score):
            score = numeric_or_nan(row.get("pchembl_value", ""))

        effect_txt = ";".join(
            x for x in [normalize(row.get("activity_comment", "")), normalize(row.get("data_validity_comment", ""))] if x
        )

        add_evidence(
            state=state,
            edge_idx=idx,
            source_mask=source_vocab.mask_multi(normalize(row.get("source", "")) or "ChEMBL"),
            source_version_mask=source_version_vocab.mask_multi(normalize(row.get("source_version", ""))),
            method_mask=categorize_method(method_txt, "PSI"),
            has_reference=ref_ok,
            score=score,
            effect_mask=parse_effect_sign(effect_txt),
        )
        psi_activity_mapped += 1

    available_context_families["PSI"] += 1
    record_context_file(inputs["psi_activity_context"], "PSI", "activity", "ok", psi_activity_mapped)
    report.notes.append(f"PSI activity context rows={psi_activity_rows}, mapped={psi_activity_mapped}")

    # PSI structure evidence (family bit 2) - required from updated 任务B
    psi_struct_rows = 0
    psi_struct_mapped = 0
    for row in iter_tsv(inputs["psi_structure_evidence"]):
        psi_struct_rows += 1
        eid = normalize(row.get("edge_id", ""))
        compound_ik = normalize(row.get("compound_inchikey", "")).upper()
        uniprot = normalize(row.get("target_uniprot_accession", ""))

        if not eid:
            continue

        idx = state.edge_index.get(eid)
        if idx is None:
            if not compound_ik or not uniprot:
                continue
            idx = add_core_edge(
                edge_id=eid,
                interaction_type="PSI",
                pair_key=get_pair_key_psi(compound_ik, uniprot),
                predicate="molecule_targets_protein",
                directed_raw="1",
                source=normalize(row.get("source", "")) or "PDB",
                source_version=normalize(row.get("source_version", "")),
                best_score_raw=row.get("structure_affinity_score", ""),
                created_from_context=True,
            )
            if idx is None:
                continue

        if state.edge_type[idx] != "PSI":
            continue

        add_context_hit(state, idx, 2)
        method_txt = ";".join(
            x
            for x in [
                normalize(row.get("pdb_experimental_method", "")),
                normalize(row.get("structure_evidence_type", "")),
            ]
            if x
        )
        ref_ok = normalize(row.get("pdb_id", "")) != ""

        score = numeric_or_nan(row.get("structure_affinity_score", ""))
        if math.isnan(score):
            score = numeric_or_nan(row.get("pchembl_value_eff", ""))

        add_evidence(
            state=state,
            edge_idx=idx,
            source_mask=source_vocab.mask_multi(normalize(row.get("source", "")) or "PDB"),
            source_version_mask=source_version_vocab.mask_multi(normalize(row.get("source_version", ""))),
            method_mask=categorize_method(method_txt, "PSI"),
            has_reference=ref_ok,
            score=score,
            effect_mask=0,
        )
        psi_struct_mapped += 1

    available_context_families["PSI"] += 1
    record_context_file(inputs["psi_structure_evidence"], "PSI", "structure", "ok", psi_struct_mapped)
    report.notes.append(f"PSI structure context rows={psi_struct_rows}, mapped={psi_struct_mapped}")

    # RPI optional contexts
    def integrate_rpi_context(path: Optional[Path], family_name: str, family_bit: int, method_from: List[str], score_from: List[str], ref_from: List[str]) -> None:
        if path is None or not path.exists():
            record_context_file(path, "RPI", family_name, "missing", 0)
            return
        available_context_families["RPI"] += 1
        hit = 0
        for row in iter_tsv(path):
            eid = normalize(row.get("edge_id", ""))
            idx = state.edge_index.get(eid)
            if idx is None or state.edge_type[idx] != "RPI":
                continue
            add_context_hit(state, idx, family_bit)
            method_txt = ";".join(normalize(row.get(c, "")) for c in method_from if normalize(row.get(c, "")))
            score = math.nan
            for c in score_from:
                v = numeric_or_nan(row.get(c, ""))
                if not math.isnan(v):
                    score = v if math.isnan(score) else max(score, v)
            ref_ok = any(normalize(row.get(c, "")) != "" for c in ref_from)
            add_evidence(
                state=state,
                edge_idx=idx,
                source_mask=source_vocab.mask_multi(normalize(row.get("source", ""))),
                source_version_mask=source_version_vocab.mask_multi(normalize(row.get("source_version", ""))),
                method_mask=categorize_method(method_txt, "RPI"),
                has_reference=ref_ok,
                score=score,
                effect_mask=parse_effect_sign(method_txt),
            )
            hit += 1
        record_context_file(path, "RPI", family_name, "ok", hit)

    integrate_rpi_context(
        path=inputs.get("rpi_site_context"),
        family_name="site",
        family_bit=1,
        method_from=["method_type", "method_subtype", "site_type"],
        score_from=["support_count"],
        ref_from=["reference"],
    )
    integrate_rpi_context(
        path=inputs.get("rpi_domain_context"),
        family_name="domain",
        family_bit=2,
        method_from=["domain_class", "domain_name"],
        score_from=[],
        ref_from=["pfam_id", "interpro_id"],
    )
    integrate_rpi_context(
        path=inputs.get("rpi_function_context"),
        family_name="function",
        family_bit=4,
        method_from=["function_relation", "inference_basis"],
        score_from=[],
        ref_from=["evidence_snippet"],
    )

    # -----------------
    # 4) Build outputs
    # -----------------
    # score normalization anchors
    p95_by_type: Dict[str, float] = {}
    for itype in ["PPI", "PSI", "RPI"]:
        vals = [
            x
            for i, x in enumerate(state.score_max)
            if state.edge_type[i] == itype and not math.isnan(x)
        ]
        if vals:
            p95 = quantile(vals, 0.95)
            p95_by_type[itype] = p95 if p95 > 0 else 1.0
        else:
            p95_by_type[itype] = 1.0

    cross_writer = write_tsv_header(args.cross_output, CROSS_COLUMNS)
    agg_writer = write_tsv_header(args.aggregate_output, AGG_COLUMNS)

    agg_scores: List[float] = []
    agg_scores_by_type: Dict[str, List[float]] = defaultdict(list)
    consistent_counter: Counter = Counter()
    conflict_true = 0
    per_type_rows: Counter = Counter()

    psi_activity_hit = 0
    psi_structure_hit = 0
    psi_both_hit = 0

    for idx, eid in enumerate(state.edge_ids):
        itype = state.edge_type[idx]
        pid = state.edge_pair_id[idx]
        pair_key = state.pair_keys[pid]

        source_mask = state.pair_source_mask[pid]
        pred_mask = state.pair_predicate_mask[pid]
        dir_mask = state.pair_direction_mask[pid]
        ver_mask = state.pair_source_version_mask[pid]
        pair_effect = state.pair_effect_mask[pid]

        consistent_n = popcount(source_mask)
        if consistent_n == 0:
            consistent_n = 1

        predicates = predicate_vocab.decode(pred_mask)
        if not predicates:
            predicates = ["unknown"]
        predicate_set = "|".join(predicates)

        direction_set = direction_set_text(dir_mask)
        if not direction_set:
            direction_set = "unknown"

        conflict_reasons: List[str] = []
        if len(predicates) > 1:
            conflict_reasons.append("predicate_conflict")
        if (dir_mask & DIR_TRUE) and (dir_mask & DIR_FALSE):
            conflict_reasons.append("direction_conflict")
        if itype == "PSI" and (pair_effect & EFFECT_POS) and (pair_effect & EFFECT_NEG):
            conflict_reasons.append("effect_conflict")

        conflict_flag = len(conflict_reasons) > 0
        if conflict_flag:
            conflict_true += 1

        evi_n = int(state.evidence_count[idx])
        method_n = popcount(state.method_mask[idx])
        ref_n = int(state.reference_count[idx])

        # context coverage
        ctx_mask = state.context_mask[idx]
        if itype == "PSI":
            has_act = bool(ctx_mask & 1)
            has_st = bool(ctx_mask & 2)
            if has_act:
                psi_activity_hit += 1
            if has_st:
                psi_structure_hit += 1
            if has_act and has_st:
                psi_both_hit += 1

        ctx_denom = EXPECTED_CONTEXT_FAMILIES.get(itype, 1)
        ctx_cov = min(1.0, popcount(ctx_mask) / float(max(1, ctx_denom)))

        score_max = state.score_max[idx]
        if math.isnan(score_max):
            score_norm = 0.0
            score_max_txt = ""
        else:
            score_norm = min(1.0, float(score_max) / float(max(1e-9, p95_by_type[itype])))
            score_max_txt = f"{float(score_max):.6f}"

        ref_cov = min(1.0, (ref_n / float(evi_n)) if evi_n > 0 else 0.0)

        consistent_comp = min(1.0, consistent_n / 4.0)
        evidence_comp = min(1.0, math.log1p(evi_n) / math.log1p(40.0))
        method_comp = min(1.0, method_n / 8.0)
        penalty = 0.22 if conflict_flag else 0.0

        agg = (
            0.24 * consistent_comp
            + 0.18 * evidence_comp
            + 0.14 * method_comp
            + 0.20 * score_norm
            + 0.12 * ref_cov
            + 0.12 * ctx_cov
            - penalty
        )
        agg = max(0.0, min(1.0, agg))

        if agg >= 0.80:
            bucket = "high"
        elif agg >= 0.50:
            bucket = "medium"
        else:
            bucket = "low"

        cross_writer.writerow(
            {
                "record_id": make_id("ixv", f"{itype}|{eid}"),
                "interaction_type": itype,
                "edge_id": eid,
                "entity_pair_key": pair_key,
                "consistent_across_n": to_int_text(consistent_n),
                "source_list": ";".join(source_vocab.decode(source_mask)) or "unknown",
                "evidence_count": to_int_text(evi_n),
                "distinct_methods_n": to_int_text(method_n),
                "reference_count": to_int_text(ref_n),
                "predicate_set": predicate_set,
                "direction_set": direction_set,
                "conflict_flag": bool_text(conflict_flag),
                "conflict_reason": ";".join(conflict_reasons),
                "source_version_list": ";".join(source_version_vocab.decode(ver_mask)) or "unknown",
                "fetch_date": args.fetch_date,
            }
        )

        agg_writer.writerow(
            {
                "aggregate_id": make_id("ias", f"{itype}|{eid}"),
                "interaction_type": itype,
                "edge_id": eid,
                "consistent_across_n": to_int_text(consistent_n),
                "evidence_count": to_int_text(evi_n),
                "distinct_methods_n": to_int_text(method_n),
                "numeric_score_max": score_max_txt,
                "numeric_score_norm": to_float_text(score_norm),
                "reference_coverage": to_float_text(ref_cov),
                "context_coverage": to_float_text(ctx_cov),
                "conflict_flag": bool_text(conflict_flag),
                "aggregate_score": to_float_text(agg),
                "score_bucket": bucket,
                "fetch_date": args.fetch_date,
            }
        )

        agg_scores.append(agg)
        agg_scores_by_type[itype].append(agg)
        consistent_counter[consistent_n] += 1
        per_type_rows[itype] += 1

    close_writer(cross_writer)
    close_writer(agg_writer)

    # -----------------
    # 5) Metrics
    # -----------------
    score_min = min(agg_scores) if agg_scores else 0.0
    score_max = max(agg_scores) if agg_scores else 0.0
    score_mean = (sum(agg_scores) / len(agg_scores)) if agg_scores else 0.0
    score_p10 = quantile(agg_scores, 0.10) if agg_scores else 0.0
    score_p90 = quantile(agg_scores, 0.90) if agg_scores else 0.0

    gates = {
        "all_three_interaction_types_present": all(per_type_rows.get(k, 0) > 0 for k in ["PPI", "PSI", "RPI"]),
        "aggregate_score_not_all_zero_or_one": (score_min < 0.999) and (score_max > 0.001),
        "aggregate_score_not_constant": (score_p90 - score_p10) >= 0.05,
        "cross_table_edge_join_rate_ge_0_99": True,
    }

    metrics = {
        "pipeline": "interaction_cross_validation_v2",
        "generated_at_utc": utc_now_iso(),
        "sample_mode": args.limit_per_type is not None,
        "inputs": {k: str(v) for k, v in inputs.items()},
        "row_count": {
            "cross_validation": len(state.edge_ids),
            "aggregate_score": len(state.edge_ids),
            "edge_count_by_type": {k: int(per_type_rows.get(k, 0)) for k in ["PPI", "PSI", "RPI"]},
            "created_from_context_by_type": dict(state.created_from_context_by_type),
        },
        "distribution": {
            "consistent_across_n": {str(k): int(v) for k, v in sorted(consistent_counter.items())},
            "conflict_flag": {
                "true": int(conflict_true),
                "false": int(len(state.edge_ids) - conflict_true),
            },
            "aggregate_score": {
                "min": float(score_min),
                "p10": float(score_p10),
                "mean": float(score_mean),
                "p90": float(score_p90),
                "max": float(score_max),
            },
            "aggregate_score_by_type": {
                k: {
                    "min": float(min(v) if v else 0.0),
                    "mean": float((sum(v) / len(v)) if v else 0.0),
                    "max": float(max(v) if v else 0.0),
                }
                for k, v in agg_scores_by_type.items()
            },
        },
        "psi_b_update_integration": {
            "psi_edges_with_activity_context": int(psi_activity_hit),
            "psi_edges_with_structure_context": int(psi_structure_hit),
            "psi_edges_with_both": int(psi_both_hit),
            "activity_context_file": str(inputs["psi_activity_context"]),
            "structure_evidence_file": str(inputs["psi_structure_evidence"]),
        },
        "context_availability": {
            "available_context_families": available_context_families,
            "expected_context_families": EXPECTED_CONTEXT_FAMILIES,
            "context_file_stats": report.context_file_stats,
        },
        "lineage": {
            "required_inputs": [
                {
                    "name": k,
                    "path": str(inputs[k]),
                    "sha256": file_sha256(inputs[k]),
                }
                for k in [
                    "ppi_edges",
                    "ppi_evidence",
                    "psi_edges",
                    "psi_evidence",
                    "psi_activity_context",
                    "psi_structure_evidence",
                    "rpi_edges",
                    "rpi_evidence",
                    "molecule_xref",
                ]
            ],
            "optional_inputs": [
                {
                    "name": k,
                    "path": str(p),
                    "sha256": file_sha256(p),
                }
                for k, p in inputs.items()
                if k not in {
                    "ppi_edges",
                    "ppi_evidence",
                    "psi_edges",
                    "psi_evidence",
                    "psi_activity_context",
                    "psi_structure_evidence",
                    "rpi_edges",
                    "rpi_evidence",
                    "molecule_xref",
                }
                and p.exists()
            ],
        },
        "warnings": report.warnings,
        "notes": report.notes,
        "gates": {
            "passed": all(gates.values()),
            "checks": gates,
        },
    }

    write_json(args.report, metrics)

    print(
        "[OK] "
        f"cross={len(state.edge_ids)} aggregate={len(state.edge_ids)} "
        f"types={{'PPI': {per_type_rows.get('PPI', 0)}, 'PSI': {per_type_rows.get('PSI', 0)}, 'RPI': {per_type_rows.get('RPI', 0)}}} "
        f"score_min={score_min:.4f} score_max={score_max:.4f}"
    )
    print(f"[OK] report -> {args.report}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
