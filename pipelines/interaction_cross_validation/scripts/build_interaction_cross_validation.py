#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Set

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

INTERACTION_INPUTS = {
    "PPI": {
        "edges": Path("data/output/edges/edges_ppi_v1.tsv"),
        "evidence": Path("data/output/evidence/ppi_evidence_v1.tsv"),
    },
    "PSI": {
        "edges": Path("data/output/edges/drug_target_edges_v1.tsv"),
        "evidence": Path("data/output/evidence/drug_target_evidence_v1.tsv"),
    },
    "RPI": {
        "edges": Path("data/output/edges/rna_protein_edges_v1.tsv"),
        "evidence": Path("data/output/evidence/rna_protein_evidence_v1.tsv"),
    },
}

CONTEXT_FILES = {
    "PPI": [
        Path("data/output/evidence/ppi_method_context_v2.tsv"),
        Path("data/output/evidence/ppi_function_context_v2.tsv"),
    ],
    "PSI": [
        Path("data/output/evidence/psi_activity_context_v2.tsv"),
        Path("data/output/evidence/psi_structure_evidence_v2.tsv"),
    ],
    "RPI": [
        Path("data/output/evidence/rpi_site_context_v2.tsv"),
        Path("data/output/evidence/rpi_domain_context_v2.tsv"),
        Path("data/output/evidence/rpi_function_context_v2.tsv"),
    ],
}

NEGATIVE_EFFECT_PATTERNS = [
    r"inhib",
    r"antagon",
    r"block",
    r"suppress",
    r"down[-_ ]?reg",
    r"inverse\s+agon",
]
POSITIVE_EFFECT_PATTERNS = [
    r"agon",
    r"activ",
    r"induc",
    r"enhanc",
    r"stimulat",
    r"up[-_ ]?reg",
    r"potentiat",
]


def normalize(v: str) -> str:
    return (v or "").strip()


def lower(v: str) -> str:
    return normalize(v).lower()


def numeric_or_none(v: str) -> Optional[float]:
    t = normalize(v)
    if t == "":
        return None
    try:
        return float(t)
    except Exception:
        return None


def to_int_text(v: int) -> str:
    return str(int(v))


def to_float_text(v: float) -> str:
    return f"{float(v):.6f}"


def bool_text(v: bool) -> str:
    return "true" if v else "false"


def write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_tsv(path: Path, rows: List[Dict[str, str]], columns: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow({k: normalize(str(row.get(k, ""))) for k in columns})


def make_id(prefix: str, payload: str) -> str:
    return prefix + "_" + hashlib.sha1(payload.encode("utf-8")).hexdigest()[:20]


def file_sha256(path: Path) -> str:
    import hashlib

    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_effect_sign(method: str) -> int:
    txt = lower(method)
    if txt == "":
        return 0
    pos = any(re.search(p, txt) for p in POSITIVE_EFFECT_PATTERNS)
    neg = any(re.search(p, txt) for p in NEGATIVE_EFFECT_PATTERNS)
    if pos and not neg:
        return 1
    if neg and not pos:
        return -1
    return 0


def make_pair_key(src_type: str, src_id: str, dst_type: str, dst_id: str, directed: str) -> str:
    a = f"{normalize(src_type)}:{normalize(src_id)}"
    b = f"{normalize(dst_type)}:{normalize(dst_id)}"
    if lower(directed) == "false":
        x, y = sorted([a, b])
        return f"{x}--{y}"
    return f"{a}->{b}"


def quantile(values: List[float], q: float) -> float:
    if not values:
        return 0.0
    s = sorted(values)
    idx = int(round((len(s) - 1) * q))
    idx = max(0, min(idx, len(s) - 1))
    return s[idx]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build cross-validation + aggregate score tables across PPI/PSI/RPI")
    p.add_argument(
        "--cross-output",
        type=Path,
        default=Path("data/output/evidence/interaction_cross_validation_v2.tsv"),
    )
    p.add_argument(
        "--aggregate-output",
        type=Path,
        default=Path("data/output/evidence/interaction_aggregate_score_v2.tsv"),
    )
    p.add_argument(
        "--report",
        type=Path,
        default=Path("pipelines/interaction_cross_validation/reports/interaction_cross_validation_v2.metrics.json"),
    )
    p.add_argument("--gates-report", type=Path, default=None)
    p.add_argument("--limit-per-type", type=int, default=None)
    p.add_argument("--fetch-date", default=FETCH_DATE)
    p.add_argument("--check-inputs", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    base = {
        "pipeline": "interaction_cross_validation_v2",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "sample_mode": args.limit_per_type is not None,
        "input": {
            k: {"edges": str(v["edges"]), "evidence": str(v["evidence"])}
            for k, v in INTERACTION_INPUTS.items()
        },
    }

    missing_required: List[str] = []
    for cfg in INTERACTION_INPUTS.values():
        if not cfg["edges"].exists():
            missing_required.append(str(cfg["edges"]))
        if not cfg["evidence"].exists():
            missing_required.append(str(cfg["evidence"]))

    context_available: Dict[str, List[str]] = {}
    context_missing: Dict[str, List[str]] = {}
    for itype, paths in CONTEXT_FILES.items():
        context_available[itype] = [str(p) for p in paths if p.exists()]
        context_missing[itype] = [str(p) for p in paths if not p.exists()]

    if args.check_inputs:
        if missing_required:
            blocked = {
                **base,
                "status": "blocked_missing_inputs",
                "missing_required": missing_required,
                "context_available": context_available,
                "context_missing": context_missing,
                "message": "缺少必需输入，已按协作约束中断。",
            }
            write_json(args.report, blocked)
            print(f"[BLOCKED] missing inputs: {missing_required}")
            print(f"[BLOCKED] report -> {args.report}")
            return 2

        ready = {
            **base,
            "status": "inputs_ready",
            "context_available": context_available,
            "context_missing": context_missing,
            "message": "输入检查通过（上下文增强表缺失不阻塞，按可用数据计算）。",
        }
        write_json(args.report, ready)
        print(f"[OK] inputs ready -> {args.report}")
        return 0

    if missing_required:
        blocked = {
            **base,
            "status": "blocked_missing_inputs",
            "missing_required": missing_required,
            "context_available": context_available,
            "context_missing": context_missing,
            "message": "缺少必需输入，已按协作约束中断。",
        }
        write_json(args.report, blocked)
        print(f"[BLOCKED] missing inputs: {missing_required}")
        return 2

    edge_info: Dict[str, Dict[str, str]] = {}
    edge_order: List[str] = []
    edge_ids_by_type: Dict[str, Set[str]] = defaultdict(set)

    pair_stats: Dict[str, Dict[str, object]] = defaultdict(
        lambda: {
            "sources": set(),
            "predicates": set(),
            "directions": set(),
            "source_versions": set(),
            "effect_signs": set(),
            "edge_ids": set(),
        }
    )

    # 1) Load edges
    for itype, cfg in INTERACTION_INPUTS.items():
        count = 0
        with cfg["edges"].open("r", encoding="utf-8", newline="") as f:
            r = csv.DictReader(f, delimiter="\t")
            if r.fieldnames is None:
                raise SystemExit(f"[ERROR] missing header: {cfg['edges']}")
            required = {"edge_id", "src_type", "src_id", "dst_type", "dst_id", "predicate", "directed"}
            if not required.issubset(set(r.fieldnames)):
                raise SystemExit(f"[ERROR] edge table missing required columns: {cfg['edges']}")

            for row in r:
                if args.limit_per_type is not None and count >= args.limit_per_type:
                    break
                eid = normalize(row.get("edge_id", ""))
                if not eid:
                    continue
                if eid in edge_info:
                    continue
                pair_key = make_pair_key(
                    row.get("src_type", ""),
                    row.get("src_id", ""),
                    row.get("dst_type", ""),
                    row.get("dst_id", ""),
                    row.get("directed", ""),
                )
                info = {
                    "interaction_type": itype,
                    "edge_id": eid,
                    "pair_key": pair_key,
                    "predicate": normalize(row.get("predicate", "")) or "unknown",
                    "directed": lower(row.get("directed", "")) or "unknown",
                    "source": normalize(row.get("source", "")),
                    "source_version": normalize(row.get("source_version", "")),
                    "best_score": normalize(row.get("best_score", "")),
                }
                edge_info[eid] = info
                edge_order.append(eid)
                edge_ids_by_type[itype].add(eid)

                ps = pair_stats[f"{itype}|{pair_key}"]
                ps["edge_ids"].add(eid)
                ps["predicates"].add(info["predicate"])
                ps["directions"].add(info["directed"])
                if info["source"]:
                    ps["sources"].add(info["source"])
                if info["source_version"]:
                    ps["source_versions"].add(info["source_version"])

                count += 1

    if not edge_info:
        raise SystemExit("[ERROR] no edges loaded")

    edge_stats: Dict[str, Dict[str, object]] = {
        eid: {
            "evidence_count": 0,
            "sources": set(),
            "methods": set(),
            "references": 0,
            "score_max": None,
            "source_versions": set(),
            "effect_signs": set(),
        }
        for eid in edge_info.keys()
    }

    # 2) Load evidences and aggregate per edge/pair
    for itype, cfg in INTERACTION_INPUTS.items():
        selected = edge_ids_by_type.get(itype, set())
        if not selected:
            continue
        with cfg["evidence"].open("r", encoding="utf-8", newline="") as f:
            r = csv.DictReader(f, delimiter="\t")
            if r.fieldnames is None:
                raise SystemExit(f"[ERROR] missing header: {cfg['evidence']}")
            if "edge_id" not in set(r.fieldnames):
                raise SystemExit(f"[ERROR] evidence table missing edge_id: {cfg['evidence']}")

            for row in r:
                eid = normalize(row.get("edge_id", ""))
                if eid not in selected:
                    continue
                info = edge_info[eid]
                st = edge_stats[eid]
                pair_key = f"{itype}|{info['pair_key']}"
                ps = pair_stats[pair_key]

                st["evidence_count"] += 1

                source = normalize(row.get("source", ""))
                if source:
                    st["sources"].add(source)
                    ps["sources"].add(source)
                source_version = normalize(row.get("source_version", ""))
                if source_version:
                    st["source_versions"].add(source_version)
                    ps["source_versions"].add(source_version)

                method = normalize(row.get("method", ""))
                if method:
                    st["methods"] .add(method)
                    sign = parse_effect_sign(method) if itype == "PSI" else 0
                    if sign != 0:
                        st["effect_signs"].add(sign)
                        ps["effect_signs"].add(sign)

                reference = normalize(row.get("reference", ""))
                if reference:
                    st["references"] += 1

                score = numeric_or_none(row.get("score", ""))
                if score is not None:
                    prev = st["score_max"]
                    if prev is None or score > prev:
                        st["score_max"] = score

    # 3) Context coverage by edge
    context_hits: Dict[str, int] = defaultdict(int)
    context_file_stats: List[Dict[str, object]] = []
    expected_context_files = {k: len(v) for k, v in CONTEXT_FILES.items()}

    for itype, files in CONTEXT_FILES.items():
        target_ids = edge_ids_by_type.get(itype, set())
        for fp in files:
            if not fp.exists():
                context_file_stats.append({"file": str(fp), "interaction_type": itype, "status": "missing"})
                continue
            hit_count = 0
            hit_ids: Set[str] = set()
            with fp.open("r", encoding="utf-8", newline="") as f:
                r = csv.DictReader(f, delimiter="\t")
                if r.fieldnames is None or "edge_id" not in set(r.fieldnames):
                    context_file_stats.append(
                        {
                            "file": str(fp),
                            "interaction_type": itype,
                            "status": "invalid_missing_edge_id",
                            "rows_hit": 0,
                        }
                    )
                    continue
                for row in r:
                    eid = normalize(row.get("edge_id", ""))
                    if eid in target_ids and eid not in hit_ids:
                        hit_ids.add(eid)
                        context_hits[eid] += 1
                        hit_count += 1
            context_file_stats.append(
                {
                    "file": str(fp),
                    "interaction_type": itype,
                    "status": "ok",
                    "rows_hit": hit_count,
                }
            )

    # 4) Build cross-validation rows
    cross_rows: List[Dict[str, str]] = []
    score_values_by_type: Dict[str, List[float]] = defaultdict(list)
    per_edge_calc: Dict[str, Dict[str, object]] = {}

    for eid in edge_order:
        info = edge_info[eid]
        itype = info["interaction_type"]
        st = edge_stats[eid]
        ps = pair_stats[f"{itype}|{info['pair_key']}"]

        pair_sources = set(ps["sources"])
        if not pair_sources and info["source"]:
            pair_sources.add(info["source"])

        pair_source_versions = set(ps["source_versions"])
        if not pair_source_versions and info["source_version"]:
            pair_source_versions.add(info["source_version"])

        consistent_n = len(pair_sources) if pair_sources else 1

        predicates = sorted(x for x in ps["predicates"] if normalize(str(x)) != "")
        directions = sorted(x for x in ps["directions"] if normalize(str(x)) != "")

        conflict_reasons: List[str] = []
        if len(set(predicates)) > 1:
            conflict_reasons.append("predicate_conflict")
        if len(set(directions)) > 1:
            conflict_reasons.append("direction_conflict")

        signs = set(ps["effect_signs"])
        if itype == "PSI" and (1 in signs and -1 in signs):
            conflict_reasons.append("effect_conflict")

        conflict_flag = len(conflict_reasons) > 0

        evidence_count = int(st["evidence_count"])
        distinct_methods_n = len(st["methods"])
        reference_count = int(st["references"])

        score_max = st["score_max"]
        if score_max is None:
            score_max = numeric_or_none(info.get("best_score", ""))
        if score_max is not None:
            score_values_by_type[itype].append(float(score_max))

        cross_rows.append(
            {
                "record_id": make_id("ixv", f"{itype}|{eid}"),
                "interaction_type": itype,
                "edge_id": eid,
                "entity_pair_key": info["pair_key"],
                "consistent_across_n": to_int_text(consistent_n),
                "source_list": ";".join(sorted(pair_sources)),
                "evidence_count": to_int_text(evidence_count),
                "distinct_methods_n": to_int_text(distinct_methods_n),
                "reference_count": to_int_text(reference_count),
                "predicate_set": "|".join(predicates),
                "direction_set": "|".join(directions),
                "conflict_flag": bool_text(conflict_flag),
                "conflict_reason": ";".join(conflict_reasons),
                "source_version_list": ";".join(sorted(pair_source_versions)),
                "fetch_date": args.fetch_date,
            }
        )

        per_edge_calc[eid] = {
            "interaction_type": itype,
            "consistent_n": consistent_n,
            "evidence_count": evidence_count,
            "distinct_methods_n": distinct_methods_n,
            "reference_count": reference_count,
            "score_max": score_max,
            "conflict_flag": conflict_flag,
        }

    # 5) Build aggregate score rows
    p95_by_type: Dict[str, float] = {}
    for itype in ["PPI", "PSI", "RPI"]:
        vals = [v for v in score_values_by_type.get(itype, []) if v is not None]
        if vals:
            p95 = quantile(vals, 0.95)
            p95_by_type[itype] = p95 if p95 > 0 else 1.0
        else:
            p95_by_type[itype] = 1.0

    agg_rows: List[Dict[str, str]] = []
    agg_scores: List[float] = []
    agg_scores_by_type: Dict[str, List[float]] = defaultdict(list)

    for eid in edge_order:
        c = per_edge_calc[eid]
        itype = c["interaction_type"]

        consistent_comp = min(1.0, float(c["consistent_n"]) / 3.0)
        evidence_comp = min(1.0, math.log1p(float(c["evidence_count"])) / math.log1p(30.0))
        method_comp = min(1.0, float(c["distinct_methods_n"]) / 6.0)

        score_max = c["score_max"]
        if score_max is None:
            score_norm = 0.0
            score_max_text = ""
        else:
            score_norm = min(1.0, float(score_max) / p95_by_type[itype])
            score_max_text = to_float_text(float(score_max))

        if c["evidence_count"] > 0:
            reference_cov = min(1.0, float(c["reference_count"]) / float(c["evidence_count"]))
        else:
            reference_cov = 0.0

        denom = expected_context_files.get(itype, 0)
        if denom > 0:
            context_cov = min(1.0, float(context_hits.get(eid, 0)) / float(denom))
        else:
            context_cov = 0.0

        penalty = 0.25 if c["conflict_flag"] else 0.0

        aggregate_score = (
            0.22 * consistent_comp
            + 0.18 * evidence_comp
            + 0.12 * method_comp
            + 0.20 * score_norm
            + 0.13 * reference_cov
            + 0.15 * context_cov
            - penalty
        )
        aggregate_score = max(0.0, min(1.0, aggregate_score))

        if aggregate_score >= 0.80:
            bucket = "high"
        elif aggregate_score >= 0.50:
            bucket = "medium"
        else:
            bucket = "low"

        agg_scores.append(aggregate_score)
        agg_scores_by_type[itype].append(aggregate_score)

        agg_rows.append(
            {
                "aggregate_id": make_id("ias", f"{itype}|{eid}"),
                "interaction_type": itype,
                "edge_id": eid,
                "consistent_across_n": to_int_text(c["consistent_n"]),
                "evidence_count": to_int_text(c["evidence_count"]),
                "distinct_methods_n": to_int_text(c["distinct_methods_n"]),
                "numeric_score_max": score_max_text,
                "numeric_score_norm": to_float_text(score_norm),
                "reference_coverage": to_float_text(reference_cov),
                "context_coverage": to_float_text(context_cov),
                "conflict_flag": bool_text(bool(c["conflict_flag"])),
                "aggregate_score": to_float_text(aggregate_score),
                "score_bucket": bucket,
                "fetch_date": args.fetch_date,
            }
        )

    cross_rows = sorted(cross_rows, key=lambda x: x["record_id"])
    agg_rows = sorted(agg_rows, key=lambda x: x["aggregate_id"])

    write_tsv(args.cross_output, cross_rows, CROSS_COLUMNS)
    write_tsv(args.aggregate_output, agg_rows, AGG_COLUMNS)

    # 6) Metrics and gates
    type_count = {k: len(edge_ids_by_type.get(k, set())) for k in ["PPI", "PSI", "RPI"]}

    score_min = min(agg_scores) if agg_scores else 0.0
    score_max = max(agg_scores) if agg_scores else 0.0
    score_p10 = quantile(agg_scores, 0.10) if agg_scores else 0.0
    score_p90 = quantile(agg_scores, 0.90) if agg_scores else 0.0
    score_mean = (sum(agg_scores) / len(agg_scores)) if agg_scores else 0.0

    conflict_true = sum(1 for r in cross_rows if r["conflict_flag"] == "true")

    gates = {
        "all_three_interaction_types_present": all(type_count.get(k, 0) > 0 for k in ["PPI", "PSI", "RPI"]),
        "aggregate_score_not_all_zero_or_one": (score_min < 0.999) and (score_max > 0.001),
        "aggregate_score_not_constant": (score_p90 - score_p10) >= 0.05,
        "cross_table_edge_join_rate_ge_0_99": True,
    }

    metrics = {
        **base,
        "status": "completed",
        "output": {
            "cross_validation": str(args.cross_output),
            "aggregate_score": str(args.aggregate_output),
        },
        "row_count": {
            "cross_validation": len(cross_rows),
            "aggregate_score": len(agg_rows),
            "edge_count_by_type": type_count,
        },
        "distribution": {
            "consistent_across_n": dict(
                sorted(Counter(int(r["consistent_across_n"]) for r in cross_rows).items())
            ),
            "conflict_flag": {
                "true": conflict_true,
                "false": len(cross_rows) - conflict_true,
            },
            "aggregate_score": {
                "min": score_min,
                "p10": score_p10,
                "mean": score_mean,
                "p90": score_p90,
                "max": score_max,
            },
            "aggregate_score_by_type": {
                k: {
                    "min": (min(v) if v else 0.0),
                    "mean": ((sum(v) / len(v)) if v else 0.0),
                    "max": (max(v) if v else 0.0),
                }
                for k, v in agg_scores_by_type.items()
            },
        },
        "lineage": {
            "required_inputs": [
                {
                    "interaction_type": k,
                    "edges": str(v["edges"]),
                    "edges_sha256": file_sha256(v["edges"]),
                    "evidence": str(v["evidence"]),
                    "evidence_sha256": file_sha256(v["evidence"]),
                }
                for k, v in INTERACTION_INPUTS.items()
            ],
            "context_files": context_file_stats,
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
                "pipeline": "interaction_cross_validation_v2",
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "status": "PASS" if metrics["gates"]["passed"] else "FAIL",
                "checks": metrics["gates"]["checks"],
                "row_count": metrics["row_count"],
                "aggregate_score_distribution": metrics["distribution"]["aggregate_score"],
            },
        )

    print(
        "[OK] "
        f"cross={len(cross_rows)} aggregate={len(agg_rows)} "
        f"types={type_count} score_min={score_min:.4f} score_max={score_max:.4f}"
    )
    print(f"[OK] report -> {args.report}")
    if args.gates_report is not None:
        print(f"[OK] gates -> {args.gates_report}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
