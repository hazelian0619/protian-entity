#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


def write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def row_count_tsv(path: Path) -> int:
    out = subprocess.check_output(["wc", "-l", str(path)]).decode("utf-8").strip().split()[0]
    return max(0, int(out) - 1)


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build release package summary/gates for interaction dataset")
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--gates", type=Path, required=True)
    p.add_argument("--checklist", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()

    tables = {
        "ppi_method_context_v2": Path("data/output/evidence/ppi_method_context_v2.tsv"),
        "ppi_function_context_v2": Path("data/output/evidence/ppi_function_context_v2.tsv"),
        "psi_activity_context_v2": Path("data/output/evidence/psi_activity_context_v2.tsv"),
        "psi_structure_evidence_v2": Path("data/output/evidence/psi_structure_evidence_v2.tsv"),
        "rpi_site_context_v2": Path("data/output/evidence/rpi_site_context_v2.tsv"),
        "rpi_domain_context_v2": Path("data/output/evidence/rpi_domain_context_v2.tsv"),
        "rpi_function_context_v2": Path("data/output/evidence/rpi_function_context_v2.tsv"),
        "interaction_cross_validation_v2": Path("data/output/evidence/interaction_cross_validation_v2.tsv"),
        "interaction_aggregate_score_v2": Path("data/output/evidence/interaction_aggregate_score_v2.tsv"),
        "interaction_ontology_mapping_v2": Path("data/output/evidence/interaction_ontology_mapping_v2.tsv"),
    }

    qa_reports = {
        "ppi_A": Path("pipelines/ppi_semantic_enrichment/reports/ppi_semantic_enrichment_v2.qa.json"),
        "psi_B": Path("pipelines/psi_activity_structure_enrichment/reports/psi_activity_structure_enrichment_v2.qa.json"),
        "cross_D": Path("pipelines/interaction_cross_validation/reports/interaction_cross_validation_v2.qa.json"),
        "rpi_C_upstream": Path("pipelines/interaction_release_local/reports/upstream_rpi_site_domain_enrichment_v2.gates.json"),
        "ontology_E_upstream": Path("pipelines/interaction_release_local/reports/upstream_interaction_ontology_mapping_v2.gates.json"),
    }

    missing_tables = [k for k, p in tables.items() if not p.exists()]
    missing_reports = [k for k, p in qa_reports.items() if not p.exists()]

    if missing_tables or missing_reports:
        payload = {
            "pipeline": "interaction_release_local",
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "status": "blocked_missing_inputs",
            "missing_tables": missing_tables,
            "missing_reports": missing_reports,
        }
        write_json(args.report, payload)
        write_json(
            args.gates,
            {
                "pipeline": "interaction_release_local",
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "status": "FAIL",
                "checks": {
                    "all_required_tables_exist": False,
                    "all_required_reports_exist": False,
                },
            },
        )
        print(f"[BLOCKED] missing_tables={len(missing_tables)} missing_reports={len(missing_reports)}")
        return 2

    table_stats: Dict[str, Dict[str, object]] = {}
    for name, path in tables.items():
        table_stats[name] = {
            "path": str(path),
            "rows": row_count_tsv(path),
            "size_bytes": path.stat().st_size,
            "modified_time_utc": datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat(),
        }

    ppi_qa = load_json(qa_reports["ppi_A"])
    psi_qa = load_json(qa_reports["psi_B"])
    cross_qa = load_json(qa_reports["cross_D"])
    c_gates = load_json(qa_reports["rpi_C_upstream"])
    e_gates = load_json(qa_reports["ontology_E_upstream"])

    checks = {
        "all_required_tables_exist": True,
        "all_required_reports_exist": True,
        "ppi_A_pass": bool(ppi_qa.get("passed", False)),
        "psi_B_pass": bool(psi_qa.get("passed", False)),
        "cross_D_pass": bool(cross_qa.get("passed", False)),
        "rpi_C_pass": str(c_gates.get("status", "")).upper() == "PASS",
        "ontology_E_pass": str(e_gates.get("status", "")).upper() == "PASS",
    }

    status = "PASS" if all(checks.values()) else "FAIL"

    coverage_snapshot = {
        "ppi_method_rows": table_stats["ppi_method_context_v2"]["rows"],
        "psi_activity_rows": table_stats["psi_activity_context_v2"]["rows"],
        "rpi_site_rows": table_stats["rpi_site_context_v2"]["rows"],
        "cross_rows": table_stats["interaction_cross_validation_v2"]["rows"],
        "ontology_rows": table_stats["interaction_ontology_mapping_v2"]["rows"],
    }

    report_payload = {
        "pipeline": "interaction_release_local",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "checks": checks,
        "table_stats": table_stats,
        "coverage_snapshot": coverage_snapshot,
        "qa_summary": {
            "ppi_A": ppi_qa.get("metrics", {}),
            "psi_B": psi_qa.get("metrics", {}),
            "cross_D": cross_qa.get("metrics", {}),
            "rpi_C_gates": c_gates.get("checks", c_gates.get("gates", {})),
            "ontology_E_gates": e_gates.get("checks", e_gates.get("gates", {})),
        },
    }
    write_json(args.report, report_payload)

    gates_payload = {
        "pipeline": "interaction_release_local",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "checks": checks,
        "coverage_snapshot": coverage_snapshot,
    }
    write_json(args.gates, gates_payload)

    checklist_lines = [
        "# 互作补洗本地收口 Checklist（2026-04-13）",
        "",
        f"- 总体状态：**{status}**",
        "",
        "## A/B/C/D/E 已实现结果表",
    ]
    for name in [
        "ppi_method_context_v2",
        "ppi_function_context_v2",
        "psi_activity_context_v2",
        "psi_structure_evidence_v2",
        "rpi_site_context_v2",
        "rpi_domain_context_v2",
        "rpi_function_context_v2",
        "interaction_cross_validation_v2",
        "interaction_aggregate_score_v2",
        "interaction_ontology_mapping_v2",
    ]:
        s = table_stats[name]
        checklist_lines.append(f"- [x] `{s['path']}`（rows={s['rows']}）")

    checklist_lines += [
        "",
        "## 质量门禁",
        f"- [x] A(PPI语义增强) QA PASS: {checks['ppi_A_pass']}",
        f"- [x] B(PSI活性结构增强) QA PASS: {checks['psi_B_pass']}",
        f"- [x] C(RPI位点结构域功能) Gates PASS: {checks['rpi_C_pass']}",
        f"- [x] D(跨库一致性与聚合评分) QA PASS: {checks['cross_D_pass']}",
        f"- [x] E(本体标准化) Gates PASS: {checks['ontology_E_pass']}",
        "",
        "## 说明",
        "- 本清单用于提审前本地收口；暂不处理新增缺口，只聚焦已实现成果并线。",
    ]
    args.checklist.parent.mkdir(parents=True, exist_ok=True)
    args.checklist.write_text("\n".join(checklist_lines) + "\n", encoding="utf-8")

    print(f"[{status}] release package -> {args.report}")
    print(f"[{status}] gates -> {args.gates}")
    print(f"[OK] checklist -> {args.checklist}")
    return 0 if status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
