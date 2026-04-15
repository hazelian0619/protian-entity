#!/usr/bin/env python3
"""Check consistency between manifests, local tables, and validation reports."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import glob
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional


def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def count_tsv_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        try:
            next(reader)
        except StopIteration:
            return 0
        return sum(1 for _ in reader)


def build_validation_row_map() -> Dict[str, int]:
    mapping: Dict[str, int] = {}
    for report in glob.glob("pipelines/**/reports/*.validation.json", recursive=True):
        bn = os.path.basename(report)
        if ".smoke." in bn or ".sample." in bn:
            continue
        try:
            d = load_json(Path(report))
        except Exception:
            continue
        table = d.get("table")
        row_count = d.get("row_count")
        if table and isinstance(row_count, int):
            mapping[table] = row_count
    return mapping


def tables_from_manifest(manifest_path: Path) -> List[Dict[str, Any]]:
    d = load_json(manifest_path)
    out: List[Dict[str, Any]] = []

    # rna-like: tables[] with source_path + rows
    for t in d.get("tables", []) if isinstance(d.get("tables"), list) else []:
        src = t.get("source_path") or t.get("path")
        rows = t.get("rows") or t.get("row_count")
        name = t.get("name") or Path(str(src)).name if src else t.get("name")
        if src:
            out.append({"name": name, "path": src, "manifest_rows": rows, "source": "tables"})

    # interaction-like: items[] with raw_path + raw_rows
    for t in d.get("items", []) if isinstance(d.get("items"), list) else []:
        src = t.get("raw_path")
        rows = t.get("raw_rows")
        name = t.get("table") or Path(str(src)).name if src else t.get("table")
        if src:
            out.append({"name": name, "path": src, "manifest_rows": rows, "source": "items"})

    # generic artifacts
    for t in d.get("artifacts", []) if isinstance(d.get("artifacts"), list) else []:
        src = t.get("path")
        rows = t.get("rows") or t.get("row_count")
        name = t.get("name") or Path(str(src)).name if src else t.get("name")
        if src and any(str(src).endswith(ext) for ext in [".tsv", ".csv"]):
            out.append({"name": name, "path": src, "manifest_rows": rows, "source": "artifacts"})

    # de-duplicate by path
    dedup: Dict[str, Dict[str, Any]] = {}
    for item in out:
        dedup[item["path"]] = item
    return list(dedup.values())


def check_release_assets_mode(
    product: str,
    latest: Dict[str, Any],
    manifest_path: Path,
    issues: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """Validate release-asset metadata without requiring local raw tables."""
    d = load_json(manifest_path)
    summary: Dict[str, Any] = {
        "product": product,
        "manifest": str(manifest_path),
        "mode": "release_assets",
        "assets_checked": 0,
        "asset_issues": 0,
    }

    # interaction-style: items[].release_assets[]
    if isinstance(d.get("items"), list):
        for item in d["items"]:
            table = item.get("table")
            release_assets = item.get("release_assets", [])
            release_available = item.get("release_available", None)
            if release_available is False:
                issues.append(
                    {
                        "product": product,
                        "severity": "error",
                        "type": "release_asset_missing",
                        "table": table,
                        "manifest": str(manifest_path),
                    }
                )
                summary["asset_issues"] += 1
            if not release_assets:
                issues.append(
                    {
                        "product": product,
                        "severity": "error",
                        "type": "release_asset_empty",
                        "table": table,
                        "manifest": str(manifest_path),
                    }
                )
                summary["asset_issues"] += 1
                continue
            for asset in release_assets:
                summary["assets_checked"] += 1
                if not asset.get("name") or not asset.get("sha256"):
                    issues.append(
                        {
                            "product": product,
                            "severity": "error",
                            "type": "release_asset_missing_metadata",
                            "table": table,
                            "asset": asset,
                            "manifest": str(manifest_path),
                        }
                    )
                    summary["asset_issues"] += 1

    # molecule-style: assets[]
    if isinstance(d.get("assets"), list):
        for asset in d["assets"]:
            summary["assets_checked"] += 1
            if not asset.get("name") or not asset.get("sha256"):
                issues.append(
                    {
                        "product": product,
                        "severity": "error",
                        "type": "release_asset_missing_metadata",
                        "asset": asset,
                        "manifest": str(manifest_path),
                    }
                )
                summary["asset_issues"] += 1

    # optional delivery report validation (interaction)
    delivery_report = latest.get("delivery_report_path")
    if delivery_report:
        dp = Path(delivery_report)
        if not dp.exists():
            issues.append(
                {
                    "product": product,
                    "severity": "error",
                    "type": "delivery_report_missing",
                    "path": delivery_report,
                }
            )
            summary["asset_issues"] += 1
        else:
            dv = load_json(dp)
            status = str(dv.get("status", "")).upper()
            if status not in {"PASS", "OK"}:
                issues.append(
                    {
                        "product": product,
                        "severity": "error",
                        "type": "delivery_report_not_pass",
                        "path": delivery_report,
                        "status": dv.get("status"),
                    }
                )
                summary["asset_issues"] += 1

    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Check release consistency across manifests and validation reports")
    parser.add_argument("--index", default="release/index.json", help="Unified release index path")
    parser.add_argument("--out", default="release/consistency_report.json", help="Output report path")
    parser.add_argument("--allow-issues", action="store_true", help="Return exit code 0 even when issues are found")
    args = parser.parse_args()

    index_path = Path(args.index)
    if not index_path.exists():
        raise SystemExit(f"[ERROR] missing index: {index_path}")

    index = load_json(index_path)
    validation_rows = build_validation_row_map()

    checks: List[Dict[str, Any]] = []
    issues: List[Dict[str, Any]] = []

    for p in index.get("products", []):
        product = p.get("product")
        latest = p.get("latest", {})
        consistency_mode = latest.get("consistency_mode", "strict_local")
        manifest_rel = latest.get("manifest_path")
        if not manifest_rel:
            checks.append({"product": product, "status": "skip", "reason": "no manifest_path"})
            continue

        manifest_path = Path(manifest_rel)
        if not manifest_path.exists():
            issues.append(
                {
                    "product": product,
                    "severity": "error",
                    "type": "manifest_missing",
                    "manifest_path": manifest_rel,
                }
            )
            continue

        if consistency_mode == "release_assets":
            checks.append(check_release_assets_mode(product, latest, manifest_path, issues))
            continue

        tables = tables_from_manifest(manifest_path)
        product_summary = {
            "product": product,
            "manifest": manifest_rel,
            "mode": "strict_local",
            "tables_checked": 0,
            "table_issues": 0,
        }

        for t in tables:
            table_path = Path(t["path"])
            product_summary["tables_checked"] += 1

            if not table_path.exists():
                issues.append(
                    {
                        "product": product,
                        "severity": "error",
                        "type": "table_missing",
                        "table": t["name"],
                        "path": t["path"],
                    }
                )
                product_summary["table_issues"] += 1
                continue

            actual_rows = count_tsv_rows(table_path)
            manifest_rows = t.get("manifest_rows")
            validation_row_count: Optional[int] = validation_rows.get(t["path"])

            if isinstance(manifest_rows, int) and actual_rows != manifest_rows:
                issues.append(
                    {
                        "product": product,
                        "severity": "error",
                        "type": "manifest_row_mismatch",
                        "table": t["name"],
                        "path": t["path"],
                        "manifest_rows": manifest_rows,
                        "actual_rows": actual_rows,
                    }
                )
                product_summary["table_issues"] += 1

            if isinstance(validation_row_count, int) and actual_rows != validation_row_count:
                issues.append(
                    {
                        "product": product,
                        "severity": "error",
                        "type": "validation_row_mismatch",
                        "table": t["name"],
                        "path": t["path"],
                        "validation_rows": validation_row_count,
                        "actual_rows": actual_rows,
                    }
                )
                product_summary["table_issues"] += 1

        checks.append(product_summary)

    report = {
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "index": str(index_path),
        "summary": {
            "products": len(index.get("products", [])),
            "checks": len(checks),
            "issues": len(issues),
            "status": "PASS" if len(issues) == 0 else "FAIL",
        },
        "checks": checks,
        "issues": issues,
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"[{'PASS' if len(issues)==0 else 'FAIL'}] issues={len(issues)} report={out_path}")
    if issues:
        for issue in issues[:20]:
            print(" -", issue)

    if issues and not args.allow_issues:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
