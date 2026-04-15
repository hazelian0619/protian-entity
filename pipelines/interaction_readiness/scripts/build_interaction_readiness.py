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
from typing import Any, Dict, Iterable, List, Optional, Tuple

DATE_TAG = "2026-04-12"
STATUS_DONE = "已完成"
STATUS_PARTIAL = "部分"
STATUS_MISSING = "缺失"


@dataclass
class ArtifactSpec:
    key: str
    role: str
    required: bool
    patterns: List[str]


INTERACTION_SPECS: Dict[str, Dict[str, Any]] = {
    "PPI": {
        "description": "Protein-Protein Interaction",
        "artifacts": [
            ArtifactSpec(
                "edges",
                "table",
                True,
                [
                    "data/output/edges/edges_ppi_v1.tsv",
                    "pipelines/edges_ppi/data/output/edges/edges_ppi_v1.tsv",
                ],
            ),
            ArtifactSpec(
                "evidence",
                "table",
                True,
                [
                    "data/output/evidence/ppi_evidence_v1.tsv",
                    "pipelines/edges_ppi/data/output/evidence/ppi_evidence_v1.tsv",
                ],
            ),
            ArtifactSpec("contract_edges", "contract", True, ["pipelines/edges_ppi/contracts/edges_ppi_v1.json"]),
            ArtifactSpec("contract_evidence", "contract", True, ["pipelines/edges_ppi/contracts/ppi_evidence_v1.json"]),
            ArtifactSpec("report_validation_edges", "report", False, ["pipelines/edges_ppi/reports/edges_ppi_v1.validation.json"]),
            ArtifactSpec("report_validation_evidence", "report", False, ["pipelines/edges_ppi/reports/ppi_evidence_v1.validation.json"]),
            ArtifactSpec("report_qa", "report", False, ["pipelines/edges_ppi/reports/edges_ppi_v1.qa.json"]),
        ],
        "expected_columns": {
            "edges": [
                "edge_id",
                "src_type",
                "src_id",
                "dst_type",
                "dst_id",
                "predicate",
                "directed",
                "source",
                "source_version",
                "fetch_date",
            ],
            "evidence": [
                "evidence_id",
                "edge_id",
                "evidence_type",
                "method",
                "score",
                "reference",
                "source",
                "source_version",
                "fetch_date",
            ],
        },
        "pipeline_hints": ["pipelines/edges_ppi/README.md", "pipelines/edges_ppi/run.sh"],
    },
    "PSI": {
        "description": "Protein-Small molecule Interaction",
        "artifacts": [
            ArtifactSpec(
                "drug_target_edges",
                "table",
                True,
                [
                    "data/output/edges/drug_target_edges_v1.tsv",
                    "pipelines/drugbank/data/output/edges/drug_target_edges_v1.tsv",
                ],
            ),
            ArtifactSpec(
                "drug_target_evidence",
                "table",
                True,
                [
                    "data/output/evidence/drug_target_evidence_v1.tsv",
                    "pipelines/drugbank/data/output/evidence/drug_target_evidence_v1.tsv",
                ],
            ),
            ArtifactSpec(
                "molecules_m3_psi",
                "table",
                True,
                [
                    "data/output/molecules/*psi*.tsv",
                    "data/output/molecules/*m3*psi*.tsv",
                    "data/output/molecules_m3*psi*.tsv",
                ],
            ),
            ArtifactSpec(
                "drug_xref_molecules",
                "table",
                False,
                [
                    "data/output/drugbank/drug_xref_molecules_v1.tsv",
                    "pipelines/drugbank/data/output/drugbank/drug_xref_molecules_v1.tsv",
                ],
            ),
            ArtifactSpec("contract_edges", "contract", False, ["pipelines/drugbank/contracts/drug_target_edges_v1.json"]),
            ArtifactSpec("contract_evidence", "contract", False, ["pipelines/drugbank/contracts/drug_target_evidence_v1.json"]),
            ArtifactSpec("report_validation_edges", "report", False, ["pipelines/drugbank/reports/drug_target_edges_v1.validation.json"]),
            ArtifactSpec("report_validation_evidence", "report", False, ["pipelines/drugbank/reports/drug_target_evidence_v1.validation.json"]),
        ],
        "expected_columns": {
            "drug_target_edges": [
                "edge_id",
                "src_type",
                "src_id",
                "dst_type",
                "dst_id",
                "predicate",
                "directed",
                "source",
                "source_version",
                "fetch_date",
            ],
            "drug_target_evidence": [
                "evidence_id",
                "edge_id",
                "evidence_type",
                "method",
                "score",
                "reference",
                "source",
                "source_version",
                "fetch_date",
            ],
        },
        "pipeline_hints": ["pipelines/drugbank/README.md", "pipelines/drugbank/run.sh", "pipelines/molecules/README.md", "pipelines/molecules/run.sh"],
    },
    "RPI": {
        "description": "RNA-Protein Interaction",
        "artifacts": [
            ArtifactSpec("edges", "table", True, ["data/output/edges/rna_protein_edges_v1.tsv"]),
            ArtifactSpec("evidence", "table", True, ["data/output/evidence/rna_protein_evidence_v1.tsv"]),
            ArtifactSpec("contract_edges", "contract", True, ["pipelines/rna_rpi/contracts/rna_protein_edges_v1.json"]),
            ArtifactSpec("contract_evidence", "contract", True, ["pipelines/rna_rpi/contracts/rna_protein_evidence_v1.json"]),
            ArtifactSpec("report_metrics", "report", True, ["pipelines/rna_rpi/reports/rna_rpi_v1.metrics.json"]),
            ArtifactSpec("report_gates", "report", True, ["pipelines/rna_rpi/reports/rna_rpi_v1.gates.json"]),
            ArtifactSpec("report_manifest", "report", True, ["pipelines/rna_rpi/reports/rna_rpi_v1.manifest.json"]),
            ArtifactSpec("report_validation_edges", "report", True, ["pipelines/rna_rpi/reports/rna_protein_edges_v1.validation.json"]),
            ArtifactSpec("report_validation_evidence", "report", True, ["pipelines/rna_rpi/reports/rna_protein_evidence_v1.validation.json"]),
        ],
        "expected_columns": {
            "edges": [
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
            ],
            "evidence": [
                "evidence_id",
                "edge_id",
                "evidence_type",
                "method",
                "score",
                "reference",
                "source",
                "source_version",
                "fetch_date",
            ],
        },
        "pipeline_hints": ["pipelines/rna_rpi/README.md", "pipelines/rna_rpi/run.sh"],
    },
}

SOURCE_MAINTAINER = {
    "starbase": "Sun Yat-sen University (ENCORI/starBase)",
    "rnainter": "RNAInter Team",
    "npinter": "NPInter Team",
    "string": "STRING Consortium",
    "biogrid": "BioGRID Team",
    "intact": "IntAct/EBI",
    "drugbank": "DrugBank (OMx Personal Health Analytics)",
    "chembl": "ChEMBL (EMBL-EBI)",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build interaction L2 local readiness report")
    p.add_argument("--repo-root", type=Path, default=Path("."))
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--doc-out", type=Path, required=True)
    p.add_argument("--mode", choices=["sample", "full"], default="full")
    p.add_argument("--sample-rows", type=int, default=200)
    p.add_argument("--offrepo-root", action="append", default=[])
    return p.parse_args()


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def rel(p: Path, root: Path) -> str:
    try:
        return str(p.resolve().relative_to(root.resolve()))
    except Exception:
        return str(p)


def find_first(root: Path, patterns: Iterable[str]) -> Optional[Path]:
    for patt in patterns:
        hits = sorted(root.glob(patt))
        for h in hits:
            if h.is_file():
                return h
    return None


def find_all(root: Path, patterns: Iterable[str]) -> List[Path]:
    out: List[Path] = []
    for patt in patterns:
        for h in sorted(root.glob(patt)):
            if h.is_file() and h not in out:
                out.append(h)
    return out


def detect_delim(path: Path) -> str:
    head = path.read_text(encoding="utf-8", errors="replace")[:4096]
    return "\t" if head.count("\t") >= head.count(",") else ","


def is_numeric(v: str) -> bool:
    try:
        float(v)
        return True
    except Exception:
        return False


def infer_col_type(values: List[str]) -> str:
    vals = [x for x in values if x != ""]
    if not vals:
        return "string"
    lv = [x.lower() for x in vals]
    if all(x in {"true", "false", "0", "1", "yes", "no"} for x in lv):
        return "bool"
    if all(re.fullmatch(r"[-+]?\d+", x) for x in vals):
        return "int"
    if all(is_numeric(x) for x in vals):
        return "float"
    if sum(1 for x in vals if re.fullmatch(r"\d{4}-\d{2}-\d{2}", x)) / len(vals) >= 0.8:
        return "date"
    return "string"


def summarize_table(path: Path, mode: str, sample_rows: int) -> Dict[str, Any]:
    delim = detect_delim(path)
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        cols = reader.fieldnames or []
        seen: List[Dict[str, str]] = []
        values: Dict[str, List[str]] = {c: [] for c in cols}
        row_count = 0
        non_empty: Counter[str] = Counter()
        for row in reader:
            row_count += 1
            if len(seen) < 1:
                seen.append({c: (row.get(c, "") or "") for c in cols})
            take_for_profile = row_count <= sample_rows if mode == "sample" else row_count <= max(sample_rows, 1000)
            if take_for_profile:
                for c in cols:
                    v = (row.get(c, "") or "").strip()
                    values[c].append(v)
                    if v != "":
                        non_empty[c] += 1
            if mode == "sample" and row_count >= sample_rows:
                break

    inferred = {c: infer_col_type(values[c]) for c in cols}
    profile_rows = min(row_count, sample_rows if mode == "sample" else max(sample_rows, 1000))
    non_empty_rate = {c: (non_empty[c] / profile_rows if profile_rows > 0 else 0.0) for c in cols}

    return {
        "path": str(path),
        "delimiter": "tsv" if delim == "\t" else "csv",
        "columns": cols,
        "row_count": row_count,
        "example": seen[0] if seen else {},
        "inferred_types": inferred,
        "non_empty_rate": non_empty_rate,
    }


def load_json_safe(path: Optional[Path]) -> Dict[str, Any]:
    if path is None or not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def combine_status(values: List[str]) -> str:
    if not values:
        return STATUS_MISSING
    if all(v == STATUS_DONE for v in values):
        return STATUS_DONE
    if all(v == STATUS_MISSING for v in values):
        return STATUS_MISSING
    return STATUS_PARTIAL


def status_from_ratio(ok: int, total: int) -> str:
    if total == 0 or ok == 0:
        return STATUS_MISSING
    if ok == total:
        return STATUS_DONE
    return STATUS_PARTIAL


def assess_interaction(
    name: str,
    conf: Dict[str, Any],
    repo_root: Path,
    mode: str,
    sample_rows: int,
    offrepo_roots: List[Path],
) -> Dict[str, Any]:
    artifacts = conf["artifacts"]
    expected_cols = conf.get("expected_columns", {})

    found: Dict[str, Optional[Path]] = {}
    offrepo_hits: Dict[str, List[str]] = {}

    for spec in artifacts:
        hit = find_first(repo_root, spec.patterns)
        found[spec.key] = hit
        if hit is None and offrepo_roots:
            candidates: List[str] = []
            for r in offrepo_roots:
                for c in find_all(r, spec.patterns):
                    candidates.append(str(c))
            if candidates:
                offrepo_hits[spec.key] = candidates

    table_profiles: Dict[str, Any] = {}
    for spec in artifacts:
        if spec.role == "table" and found.get(spec.key) is not None:
            table_profiles[spec.key] = summarize_table(found[spec.key], mode=mode, sample_rows=sample_rows)

    # 组件
    required_components = [s for s in artifacts if s.required and s.role == "table"]
    completed_components = sum(1 for s in required_components if found.get(s.key) is not None)
    component_status = status_from_ratio(completed_components, len(required_components))

    # 属性字段
    schema_hits = 0
    schema_total = 0
    missing_columns_by_component: Dict[str, List[str]] = {}
    for key, cols in expected_cols.items():
        schema_total += 1
        profile = table_profiles.get(key)
        if not profile:
            missing_columns_by_component[key] = cols
            continue
        present = set(profile["columns"])
        miss = [c for c in cols if c not in present]
        missing_columns_by_component[key] = miss
        if not miss:
            schema_hits += 1
    field_status = status_from_ratio(schema_hits, schema_total)

    # 数据类型
    dtype_status = STATUS_MISSING
    if table_profiles:
        bad = 0
        total = 0
        for _, p in table_profiles.items():
            for _, t in p["inferred_types"].items():
                total += 1
                if t not in {"string", "int", "float", "bool", "date"}:
                    bad += 1
        dtype_status = STATUS_DONE if total > 0 and bad == 0 else STATUS_PARTIAL

    # 获取方式
    hints = conf.get("pipeline_hints", [])
    hint_hits = sum(1 for h in hints if (repo_root / h).exists())
    acquisition_status = status_from_ratio(hint_hits, len(hints)) if hints else STATUS_MISSING

    # 存储形式
    storage_status = STATUS_DONE if table_profiles else STATUS_MISSING

    # 示例
    example_status = STATUS_DONE if any(p.get("example") for p in table_profiles.values()) else STATUS_MISSING

    # 指标/用途
    metrics_report = None
    for k in ["report_metrics", "report_qa", "report_gates"]:
        if found.get(k) is not None:
            metrics_report = found[k]
            break
    metric_status = STATUS_DONE if metrics_report else STATUS_MISSING
    usage_status = metric_status

    # 层次结构: evidence.edge_id 能否 join edges.edge_id
    hierarchy_status = STATUS_MISSING
    hierarchy = {}
    edge_keys = ["edges", "drug_target_edges"]
    evidence_keys = ["evidence", "drug_target_evidence"]
    edge_profile = next((table_profiles[k] for k in edge_keys if k in table_profiles), None)
    ev_profile = next((table_profiles[k] for k in evidence_keys if k in table_profiles), None)
    if edge_profile and ev_profile and "edge_id" in edge_profile["columns"] and "edge_id" in ev_profile["columns"]:
        edge_ids = set()
        with Path(edge_profile["path"]).open("r", encoding="utf-8", errors="replace", newline="") as f:
            r = csv.DictReader(f, delimiter="\t" if edge_profile["delimiter"] == "tsv" else ",")
            for i, row in enumerate(r, start=1):
                edge_ids.add((row.get("edge_id", "") or "").strip())
                if mode == "sample" and i >= sample_rows:
                    break
        total = 0
        hit = 0
        with Path(ev_profile["path"]).open("r", encoding="utf-8", errors="replace", newline="") as f:
            r = csv.DictReader(f, delimiter="\t" if ev_profile["delimiter"] == "tsv" else ",")
            for i, row in enumerate(r, start=1):
                eid = (row.get("edge_id", "") or "").strip()
                if eid:
                    total += 1
                    if eid in edge_ids:
                        hit += 1
                if mode == "sample" and i >= sample_rows:
                    break
        join_rate = hit / total if total else 0.0
        hierarchy = {"evidence_edge_join_rate": join_rate, "checked_rows": total}
        hierarchy_status = STATUS_DONE if join_rate >= 0.99 else (STATUS_PARTIAL if join_rate > 0 else STATUS_MISSING)

    # 语义关系
    semantic_status = STATUS_MISSING
    semantic_required = {
        "edge": ["predicate", "directed"],
        "evidence": ["evidence_type"],
    }
    sem_hit = 0
    sem_total = 0
    if edge_profile:
        sem_total += len(semantic_required["edge"])
        sem_hit += sum(1 for c in semantic_required["edge"] if c in edge_profile["columns"])
    if ev_profile:
        sem_total += len(semantic_required["evidence"])
        sem_hit += sum(1 for c in semantic_required["evidence"] if c in ev_profile["columns"])
    if sem_total > 0:
        semantic_status = STATUS_DONE if sem_hit == sem_total else STATUS_PARTIAL

    # 数据来源 + 维护机构 + 更新频率
    source_status = STATUS_MISSING
    maintainer_status = STATUS_MISSING
    update_freq_status = STATUS_MISSING
    source_values: Counter[str] = Counter()
    fetch_dates: Counter[str] = Counter()
    for prof in table_profiles.values():
        cols = set(prof["columns"])
        if {"source", "source_version", "fetch_date"}.issubset(cols):
            source_status = STATUS_DONE
        path = Path(prof["path"])
        with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
            r = csv.DictReader(f, delimiter="\t" if prof["delimiter"] == "tsv" else ",")
            for i, row in enumerate(r, start=1):
                src = (row.get("source", "") or "").strip()
                if src:
                    source_values[src.lower()] += 1
                dt = (row.get("fetch_date", "") or "").strip()
                if dt:
                    fetch_dates[dt] += 1
                if mode == "sample" and i >= sample_rows:
                    break
    if source_values:
        mapped = 0
        for s in source_values:
            if any(k in s for k in SOURCE_MAINTAINER):
                mapped += 1
        maintainer_status = STATUS_DONE if mapped > 0 else STATUS_PARTIAL
    if fetch_dates:
        update_freq_status = STATUS_PARTIAL if len(fetch_dates) == 1 else STATUS_DONE

    # 当前状态 / 证据路径
    dimension_status = {
        "组件": component_status,
        "属性字段": field_status,
        "数据类型": dtype_status,
        "获取方式": acquisition_status,
        "存储形式": storage_status,
        "示例": example_status,
        "用途/指标": usage_status,
        "指标": metric_status,
        "层次结构": hierarchy_status,
        "语义关系": semantic_status,
        "数据来源": source_status,
        "维护机构": maintainer_status,
        "更新频率": update_freq_status,
        "证据路径": STATUS_DONE if any(found.values()) else STATUS_MISSING,
    }
    current_status = combine_status(list(dimension_status.values()))
    dimension_status["当前状态"] = current_status

    evidence_paths = {
        "in_repo": {k: rel(v, repo_root) for k, v in found.items() if v is not None},
        "off_repo_candidates": offrepo_hits,
    }

    source_preview = [{"source": k, "count": v, "maintainer": _map_maintainer(k)} for k, v in source_values.most_common(5)]

    return {
        "interaction": name,
        "description": conf["description"],
        "status": current_status,
        "dimension_status": dimension_status,
        "artifacts": {
            "required_tables": [s.key for s in required_components],
            "found_in_repo": {k: rel(v, repo_root) for k, v in found.items() if v is not None},
            "missing_required": [s.key for s in artifacts if s.required and found.get(s.key) is None],
            "off_repo_candidates": offrepo_hits,
        },
        "schema": {
            "expected_columns": expected_cols,
            "missing_columns_by_component": missing_columns_by_component,
        },
        "table_profiles": table_profiles,
        "metrics": {
            "component_completion": {
                "done": completed_components,
                "total": len(required_components),
                "ratio": (completed_components / len(required_components)) if required_components else 0.0,
            },
            "schema_completion": {
                "done": schema_hits,
                "total": schema_total,
                "ratio": (schema_hits / schema_total) if schema_total else 0.0,
            },
            **hierarchy,
            "source_preview": source_preview,
            "fetch_date_preview": [{"fetch_date": k, "count": v} for k, v in fetch_dates.most_common(5)],
        },
        "evidence_path": evidence_paths,
    }


def _map_maintainer(source_text: str) -> str:
    s = source_text.lower()
    for k, v in SOURCE_MAINTAINER.items():
        if k in s:
            return v
    return "未知"


def build_top_gaps(interactions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    gaps: List[Dict[str, Any]] = []

    for item in interactions:
        it = item["interaction"]
        missing_required = item["artifacts"]["missing_required"]
        for key in missing_required:
            sev = "阻塞级"
            effort = "L"
            if "contract" in key:
                effort = "S"
            elif "molecules_m3" in key:
                effort = "M"
            gaps.append(
                {
                    "interaction": it,
                    "gap": f"缺少必需产物: {key}",
                    "severity": sev,
                    "effort": effort,
                    "evidence": item["evidence_path"],
                    "suggestion": "补齐对应 pipeline 产物并运行 contract+validation。",
                }
            )

        dim = item["dimension_status"]
        for dkey in ["属性字段", "层次结构", "语义关系", "维护机构", "更新频率"]:
            if dim.get(dkey) in {STATUS_MISSING, STATUS_PARTIAL}:
                sev = "阻塞级" if dkey in {"属性字段", "层次结构", "语义关系"} else "可延期"
                effort = "M" if dkey in {"属性字段", "层次结构"} else "S"
                gaps.append(
                    {
                        "interaction": it,
                        "gap": f"{dkey}为{dim.get(dkey)}",
                        "severity": sev,
                        "effort": effort,
                        "evidence": item["metrics"],
                        "suggestion": "按统一 L2 字段规范补齐并在本地复跑 QA。",
                    }
                )

    def score(g: Dict[str, Any]) -> Tuple[int, int]:
        s = 0 if g["severity"] == "阻塞级" else 1
        e = {"S": 0, "M": 1, "L": 2}.get(g["effort"], 1)
        return (s, e)

    gaps_sorted = sorted(gaps, key=score)
    dedup: List[Dict[str, Any]] = []
    seen = set()
    for g in gaps_sorted:
        k = (g["interaction"], g["gap"])
        if k in seen:
            continue
        seen.add(k)
        dedup.append(g)
    blockers = [g for g in dedup if g["severity"] == "阻塞级"]
    deferables = [g for g in dedup if g["severity"] == "可延期"]

    picked: List[Dict[str, Any]] = []
    if deferables:
        # 保留可延期缺口可见性：优先放入最多 8 个阻塞项，给可延期至少 2 个位置
        picked.extend(blockers[:8])
        picked.extend(deferables[: max(0, 10 - len(picked))])
        if len(picked) < 10:
            picked.extend(blockers[8 : 8 + (10 - len(picked))])
    else:
        picked.extend(blockers[:10])
    return picked[:10]


def build_preupload_gates(interactions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    by_name = {x["interaction"]: x for x in interactions}

    gates = []

    for name in ["PPI", "PSI", "RPI"]:
        i = by_name.get(name, {})
        miss = i.get("artifacts", {}).get("missing_required", ["<missing interaction>"])
        passed = len(miss) == 0
        gates.append(
            {
                "gate_id": f"{name.lower()}_required_artifacts_ready",
                "description": f"{name} 必需产物（edges/evidence/核心衍生）齐全",
                "passed": passed,
                "evidence": {"missing_required": miss},
            }
        )

    rpi = by_name.get("RPI", {})
    rpi_metrics = rpi.get("metrics", {})
    join_rate = rpi_metrics.get("evidence_edge_join_rate", 0.0)
    gates.append(
        {
            "gate_id": "rpi_edge_evidence_hierarchy_ge_0_99",
            "description": "RPI evidence.edge_id 与 edges.edge_id 连接率 >= 0.99",
            "passed": join_rate >= 0.99,
            "evidence": {"join_rate": join_rate},
        }
    )

    all_schema_done = all(x["dimension_status"].get("属性字段") == STATUS_DONE for x in interactions)
    gates.append(
        {
            "gate_id": "interaction_schema_uniformity",
            "description": "三类互作均满足统一 L2 字段口径",
            "passed": all_schema_done,
            "evidence": {k: v["dimension_status"].get("属性字段") for k, v in by_name.items()},
        }
    )

    has_blockers = any(g["severity"] == "阻塞级" for g in build_top_gaps(interactions))
    gates.append(
        {
            "gate_id": "zero_blocking_gaps",
            "description": "Top gaps 中阻塞级缺口为 0",
            "passed": not has_blockers,
            "evidence": {"has_blocking_gaps": has_blockers},
        }
    )

    return gates


def render_markdown(report: Dict[str, Any]) -> str:
    lines: List[str] = []
    lines.append(f"# INTERACTION L2 Local Readiness ({DATE_TAG})")
    lines.append("")
    lines.append(f"- 生成时间(UTC): {report['generated_at_utc']}")
    lines.append(f"- 扫描模式: {report['mode']}")
    lines.append(f"- 总体结论: **{report['overall_status']}**")
    lines.append("")

    lines.append("## 1) 三类互作同口径状态表")
    lines.append("")
    header = [
        "互作",
        "组件",
        "属性字段",
        "数据类型",
        "获取方式",
        "存储形式",
        "示例",
        "用途/指标",
        "指标",
        "层次结构",
        "语义关系",
        "数据来源",
        "维护机构",
        "更新频率",
        "证据路径",
        "当前状态",
    ]
    lines.append("| " + " | ".join(header) + " |")
    lines.append("|" + "---|" * len(header))
    for item in report["interactions"]:
        d = item["dimension_status"]
        row = [
            item["interaction"],
            d.get("组件", ""),
            d.get("属性字段", ""),
            d.get("数据类型", ""),
            d.get("获取方式", ""),
            d.get("存储形式", ""),
            d.get("示例", ""),
            d.get("用途/指标", ""),
            d.get("指标", ""),
            d.get("层次结构", ""),
            d.get("语义关系", ""),
            d.get("数据来源", ""),
            d.get("维护机构", ""),
            d.get("更新频率", ""),
            d.get("证据路径", ""),
            d.get("当前状态", ""),
        ]
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")

    lines.append("## 2) 关键指标摘要")
    lines.append("")
    for item in report["interactions"]:
        m = item.get("metrics", {})
        cc = m.get("component_completion", {})
        sc = m.get("schema_completion", {})
        lines.append(f"### {item['interaction']}")
        lines.append(f"- 组件完成度: {cc.get('done', 0)}/{cc.get('total', 0)} ({cc.get('ratio', 0.0):.2%})")
        lines.append(f"- 字段完成度: {sc.get('done', 0)}/{sc.get('total', 0)} ({sc.get('ratio', 0.0):.2%})")
        if "evidence_edge_join_rate" in m:
            lines.append(f"- evidence-edge join_rate: {m['evidence_edge_join_rate']:.4f}")
        if m.get("source_preview"):
            src = ", ".join(f"{x['source']}({x['count']})" for x in m["source_preview"][:3])
            lines.append(f"- 来源预览: {src}")
        lines.append(f"- 证据路径: {json.dumps(item['evidence_path']['in_repo'], ensure_ascii=False)}")
        lines.append("")

    lines.append("## 3) Top10 缺口（阻塞级/可延期 + 工作量）")
    lines.append("")
    lines.append("| # | 互作 | 缺口 | 分级 | 工作量 | 建议 |")
    lines.append("|---|---|---|---|---|---|")
    for i, g in enumerate(report["top10_gaps"], start=1):
        lines.append(
            f"| {i} | {g['interaction']} | {g['gap']} | {g['severity']} | {g['effort']} | {g['suggestion']} |"
        )
    lines.append("")

    lines.append("## 4) 本地通过后再上传：前置门槛")
    lines.append("")
    lines.append("| Gate | 说明 | 结果 |")
    lines.append("|---|---|---|")
    for g in report["preupload_gates"]:
        lines.append(f"| {g['gate_id']} | {g['description']} | {'PASS' if g['passed'] else 'FAIL'} |")
    lines.append("")

    lines.append("## 5) 收口建议")
    lines.append("")
    if report.get("overall_status") == "PASS":
        lines.append("1. 已完成三类互作本地收口，建议冻结当前 manifest 与 gates 作为上传前基线。")
        lines.append("2. 维持统一 L2 字段口径（method/reference/score/source_version/fetch_date），后续增量按同合同校验。")
        lines.append("3. 当前仅剩可延期项（更新频率分层），可在下一轮引入多时间窗快照后提升为“已完成”。")
    else:
        lines.append("1. 先补齐 PPI/PSI 的本地 edges+evidence+contract+validation 产物，再做统一上传。")
        lines.append("2. 保持三类互作统一 L2 字段（method/reference/score/source_version/fetch_date）。")
        lines.append("3. 所有缺口关闭后，重新运行 interaction_readiness 全量核查并锁定 manifest。")
    lines.append("")

    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    offrepo_roots = [Path(x).resolve() for x in args.offrepo_root]

    interactions = [
        assess_interaction(
            name=k,
            conf=v,
            repo_root=repo_root,
            mode=args.mode,
            sample_rows=args.sample_rows,
            offrepo_roots=offrepo_roots,
        )
        for k, v in INTERACTION_SPECS.items()
    ]

    top10 = build_top_gaps(interactions)
    gates = build_preupload_gates(interactions)

    blocking = sum(1 for g in top10 if g["severity"] == "阻塞级")
    overall_status = "PASS" if blocking == 0 and all(g["passed"] for g in gates) else "FAIL"

    report = {
        "report_id": f"interaction_l2_readiness_{DATE_TAG}",
        "generated_at_utc": now_utc(),
        "mode": args.mode,
        "repo_root": str(repo_root),
        "dimensions_required": [
            "组件",
            "属性字段",
            "数据类型",
            "获取方式",
            "存储形式",
            "示例",
            "用途/指标",
            "指标",
            "层次结构",
            "语义关系",
            "数据来源",
            "维护机构",
            "更新频率",
            "当前状态",
            "证据路径",
        ],
        "interactions": interactions,
        "top10_gaps": top10,
        "preupload_gates": gates,
        "summary": {
            "interaction_status": {x["interaction"]: x["status"] for x in interactions},
            "blocking_gap_count": blocking,
            "gate_pass_count": sum(1 for g in gates if g["passed"]),
            "gate_total": len(gates),
        },
        "overall_status": overall_status,
    }

    write_json(args.report, report)

    md = render_markdown(report)
    args.doc_out.parent.mkdir(parents=True, exist_ok=True)
    args.doc_out.write_text(md, encoding="utf-8")

    print(
        "[OK] interaction_readiness "
        f"mode={args.mode} overall={overall_status} blocking_gaps={blocking} "
        f"report={args.report} doc={args.doc_out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
