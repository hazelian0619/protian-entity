#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def tsv_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8", newline="") as f:
        return max(0, sum(1 for _ in f) - 1)


def read_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build local release bundle summary for interaction datasets")
    p.add_argument("--release-tag", default="2026-04-13")
    p.add_argument(
        "--summary-out",
        type=Path,
        default=Path("pipelines/interaction_release_local/reports/interaction_release_local_2026-04-13.summary.json"),
    )
    p.add_argument(
        "--gates-out",
        type=Path,
        default=Path("pipelines/interaction_release_local/reports/interaction_release_local_2026-04-13.gates.json"),
    )
    p.add_argument(
        "--checklist-out",
        type=Path,
        default=Path("../_coord/互作补洗_本地收口Checklist_2026-04-13.md"),
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()

    required_tables = {
        "PPI_edges_v1": Path("data/output/edges/edges_ppi_v1.tsv"),
        "PPI_evidence_v1": Path("data/output/evidence/ppi_evidence_v1.tsv"),
        "PPI_method_context_v2": Path("data/output/evidence/ppi_method_context_v2.tsv"),
        "PPI_function_context_v2": Path("data/output/evidence/ppi_function_context_v2.tsv"),
        "PSI_edges_v1": Path("data/output/edges/drug_target_edges_v1.tsv"),
        "PSI_evidence_v1": Path("data/output/evidence/drug_target_evidence_v1.tsv"),
        "PSI_activity_context_v2": Path("data/output/evidence/psi_activity_context_v2.tsv"),
        "PSI_structure_evidence_v2": Path("data/output/evidence/psi_structure_evidence_v2.tsv"),
        "RPI_edges_v1": Path("data/output/edges/rna_protein_edges_v1.tsv"),
        "RPI_evidence_v1": Path("data/output/evidence/rna_protein_evidence_v1.tsv"),
        "RPI_site_context_v2": Path("data/output/evidence/rpi_site_context_v2.tsv"),
        "RPI_domain_context_v2": Path("data/output/evidence/rpi_domain_context_v2.tsv"),
        "RPI_function_context_v2": Path("data/output/evidence/rpi_function_context_v2.tsv"),
        "interaction_cross_validation_v2": Path("data/output/evidence/interaction_cross_validation_v2.tsv"),
        "interaction_aggregate_score_v2": Path("data/output/evidence/interaction_aggregate_score_v2.tsv"),
        "interaction_ontology_mapping_v2": Path("data/output/evidence/interaction_ontology_mapping_v2.tsv"),
    }

    qa_reports = {
        "A_ppi_semantic": Path("pipelines/ppi_semantic_enrichment/reports/ppi_semantic_enrichment_v2.qa.json"),
        "B_psi_activity_structure": Path("pipelines/psi_activity_structure_enrichment/reports/psi_activity_structure_enrichment_v2.qa.json"),
        "C_rpi_site_domain": Path("pipelines/rpi_site_domain_enrichment/reports/rpi_site_domain_enrichment_v2.gates.json"),
        "D_interaction_cross_validation": Path("pipelines/interaction_cross_validation/reports/interaction_cross_validation_v2.qa.json"),
        "E_interaction_ontology": Path("pipelines/interaction_ontology_mapping/reports/interaction_ontology_mapping_v2.gates.json"),
    }

    table_stats: Dict[str, Dict] = {}
    missing_tables: List[str] = []

    for name, path in required_tables.items():
        if not path.exists():
            missing_tables.append(name)
            table_stats[name] = {"path": str(path), "exists": False}
            continue
        table_stats[name] = {
            "path": str(path),
            "exists": True,
            "rows": tsv_rows(path),
            "sha256": sha256(path),
        }

    qa_stats: Dict[str, Dict] = {}
    missing_qa: List[str] = []

    for name, path in qa_reports.items():
        if not path.exists():
            missing_qa.append(name)
            qa_stats[name] = {"path": str(path), "exists": False, "passed": False}
            continue
        obj = read_json(path)
        passed = bool(obj.get("passed", False) if "passed" in obj else obj.get("status") == "PASS")
        qa_stats[name] = {
            "path": str(path),
            "exists": True,
            "passed": passed,
        }

    qa_all_pass = all(v.get("passed", False) for v in qa_stats.values()) and len(qa_stats) > 0
    all_tables_present = len(missing_tables) == 0

    checks = {
        "required_tables_present": all_tables_present,
        "all_qa_reports_present": len(missing_qa) == 0,
        "all_qa_pass": qa_all_pass,
    }

    status = "PASS" if all(checks.values()) else "FAIL"

    summary = {
        "pipeline": "interaction_release_local",
        "release_tag": args.release_tag,
        "generated_at_utc": now_iso(),
        "status": status,
        "checks": checks,
        "missing_tables": missing_tables,
        "missing_qa_reports": missing_qa,
        "tables": table_stats,
        "qa_reports": qa_stats,
        "known_boundaries": [
            "本次提审聚焦已实现数据集，未覆盖新增研发缺口（如复合物成员/taxonomy 扩展）。",
            "跨仓依赖已收敛到 1218 本地产物路径。",
        ],
    }

    gates = {
        "pipeline": "interaction_release_local",
        "release_tag": args.release_tag,
        "generated_at_utc": now_iso(),
        "status": status,
        "checks": checks,
        "table_count": len(required_tables),
        "qa_report_count": len(qa_reports),
    }

    checklist_lines = [
        f"# 互作补洗本地收口 Checklist（{args.release_tag}）",
        "",
        f"- 生成时间（UTC）：{summary['generated_at_utc']}",
        f"- 收口状态：**{status}**",
        "",
        "## 1) 必需结果表",
    ]
    for name, v in table_stats.items():
        mark = "✅" if v.get("exists") else "❌"
        extra = f"rows={v.get('rows')}" if v.get("exists") else "MISSING"
        checklist_lines.append(f"- {mark} `{name}` -> `{v['path']}` ({extra})")

    checklist_lines += ["", "## 2) QA/Gates", ""]
    for name, v in qa_stats.items():
        mark = "✅" if v.get("passed") else "❌"
        checklist_lines.append(f"- {mark} `{name}` -> `{v['path']}`")

    checklist_lines += [
        "",
        "## 3) 总结",
        "",
        f"- required_tables_present: `{checks['required_tables_present']}`",
        f"- all_qa_reports_present: `{checks['all_qa_reports_present']}`",
        f"- all_qa_pass: `{checks['all_qa_pass']}`",
    ]

    write_json(args.summary_out, summary)
    write_json(args.gates_out, gates)
    write_text(args.checklist_out, "\n".join(checklist_lines) + "\n")

    print(f"[{status}] summary -> {args.summary_out}")
    print(f"[{status}] gates   -> {args.gates_out}")
    print(f"[{status}] checklist -> {args.checklist_out}")

    return 0 if status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
