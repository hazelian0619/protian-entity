#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


def write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Sync interaction V2 artifacts from another repo into local 1218 workspace")
    p.add_argument("--source-root", type=Path, default=Path("../protian-entity"))
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--check-only", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    src_root = args.source_root.resolve()

    required_local = [
        Path("data/output/evidence/ppi_method_context_v2.tsv"),
        Path("data/output/evidence/ppi_function_context_v2.tsv"),
        Path("data/output/evidence/psi_activity_context_v2.tsv"),
        Path("data/output/evidence/psi_structure_evidence_v2.tsv"),
        Path("data/output/evidence/interaction_cross_validation_v2.tsv"),
        Path("data/output/evidence/interaction_aggregate_score_v2.tsv"),
    ]

    sync_pairs = [
        (
            src_root / "data/output/evidence/rpi_site_context_v2.tsv",
            Path("data/output/evidence/rpi_site_context_v2.tsv"),
            "RPI site context v2",
        ),
        (
            src_root / "data/output/evidence/rpi_domain_context_v2.tsv",
            Path("data/output/evidence/rpi_domain_context_v2.tsv"),
            "RPI domain context v2",
        ),
        (
            src_root / "data/output/evidence/rpi_function_context_v2.tsv",
            Path("data/output/evidence/rpi_function_context_v2.tsv"),
            "RPI function context v2",
        ),
        (
            src_root / "data/output/evidence/interaction_ontology_mapping_v2.tsv",
            Path("data/output/evidence/interaction_ontology_mapping_v2.tsv"),
            "Interaction ontology mapping v2",
        ),
        (
            src_root / "pipelines/rpi_site_domain_enrichment/reports/rpi_site_domain_enrichment_v2.gates.json",
            Path("pipelines/interaction_release_local/reports/upstream_rpi_site_domain_enrichment_v2.gates.json"),
            "Upstream RPI gates report",
        ),
        (
            src_root / "pipelines/interaction_ontology_mapping/reports/interaction_ontology_mapping_v2.gates.json",
            Path("pipelines/interaction_release_local/reports/upstream_interaction_ontology_mapping_v2.gates.json"),
            "Upstream ontology gates report",
        ),
    ]

    missing_local = [str(p) for p in required_local if not p.exists()]
    missing_upstream = [str(src) for src, _, _ in sync_pairs if not src.exists()]

    if missing_local or missing_upstream:
        payload = {
            "pipeline": "interaction_release_local",
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "status": "blocked_missing_inputs",
            "missing_local": missing_local,
            "missing_upstream": missing_upstream,
            "source_root": str(src_root),
            "message": "缺少提审必需输入，已中断。",
        }
        write_json(args.report, payload)
        print(f"[BLOCKED] local_missing={len(missing_local)} upstream_missing={len(missing_upstream)}")
        return 2

    copied: List[Dict[str, str]] = []
    skipped: List[Dict[str, str]] = []

    for src, dst, desc in sync_pairs:
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists() and dst.stat().st_size == src.stat().st_size:
            skipped.append({"description": desc, "src": str(src), "dst": str(dst), "reason": "same_size_exists"})
            continue
        if not args.check_only:
            shutil.copy2(src, dst)
        copied.append({"description": desc, "src": str(src), "dst": str(dst)})

    payload = {
        "pipeline": "interaction_release_local",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "ready" if args.check_only else "synced",
        "source_root": str(src_root),
        "copied_count": len(copied),
        "skipped_count": len(skipped),
        "copied": copied,
        "skipped": skipped,
    }
    write_json(args.report, payload)
    print(f"[OK] sync complete: copied={len(copied)} skipped={len(skipped)} -> {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
