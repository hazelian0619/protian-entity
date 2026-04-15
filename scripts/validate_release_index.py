#!/usr/bin/env python3
"""Validate release/index.json structure and referenced paths.

Primary checks:
1) Schema validation (jsonschema, if available)
2) Path integrity for current/manifest/checksum/quality reports
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List


def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def basic_structure_errors(index: Dict[str, Any]) -> List[str]:
    errors: List[str] = []
    if not isinstance(index, dict):
        return ["index root must be an object"]

    for k in ["schema_version", "generated_at_utc", "products"]:
        if k not in index:
            errors.append(f"missing root key: {k}")

    products = index.get("products")
    if not isinstance(products, list):
        errors.append("products must be an array")
        return errors

    for i, product in enumerate(products):
        if not isinstance(product, dict):
            errors.append(f"products[{i}] must be an object")
            continue
        if "product" not in product:
            errors.append(f"products[{i}] missing key: product")
        latest = product.get("latest")
        if not isinstance(latest, dict):
            errors.append(f"products[{i}].latest must be an object")
            continue
        for req in ["version", "distribution_mode"]:
            if req not in latest:
                errors.append(f"products[{i}].latest missing key: {req}")

    return errors


def path_integrity_errors(
    index: Dict[str, Any],
    repo_root: Path,
    strict_github_release_paths: bool = False,
) -> List[str]:
    errors: List[str] = []
    for p in index.get("products", []):
        pname = p.get("product", "<unknown>")

        current_path = p.get("current_path")
        if current_path and not (repo_root / current_path).exists():
            errors.append(f"{pname}: current_path missing: {current_path}")

        latest = p.get("latest", {}) if isinstance(p.get("latest"), dict) else {}
        dist_mode = latest.get("distribution_mode")
        is_github_release = dist_mode == "github_release"

        # For release-based products, local bundle files may not be checked into git.
        # Enforce release metadata fields first; local file paths are optional unless strict mode is enabled.
        if is_github_release:
            if not latest.get("release_tag"):
                errors.append(f"{pname}: missing release_tag for github_release")
            if not latest.get("release_url"):
                errors.append(f"{pname}: missing release_url for github_release")

        for key in ["manifest_path", "checksum_path", "delivery_report_path"]:
            rel = latest.get(key)
            if not rel:
                continue
            if is_github_release and not strict_github_release_paths and key in {"manifest_path", "checksum_path"}:
                continue
            if not (repo_root / rel).exists():
                errors.append(f"{pname}: {key} missing: {rel}")

        for art in latest.get("artifacts", []) if isinstance(latest.get("artifacts"), list) else []:
            rel = art.get("path")
            if rel and not (repo_root / rel).exists():
                errors.append(f"{pname}: artifact path missing: {rel}")

        for report in latest.get("quality_reports", []) if isinstance(latest.get("quality_reports"), list) else []:
            if not report:
                continue
            if is_github_release and not strict_github_release_paths:
                continue
            if not (repo_root / report).exists():
                errors.append(f"{pname}: quality report missing: {report}")

    return errors


def schema_errors(index: Dict[str, Any], schema_path: Path) -> List[str]:
    if not schema_path.exists():
        return [f"schema missing: {schema_path}"]

    try:
        import jsonschema  # type: ignore
    except Exception:
        return [
            "jsonschema is not installed; install via `python3 -m pip install jsonschema` "
            "or rely on basic structure checks only"
        ]

    schema = load_json(schema_path)
    validator = jsonschema.Draft202012Validator(schema)
    errs = sorted(validator.iter_errors(index), key=lambda e: list(e.path))

    out: List[str] = []
    for e in errs:
        path = ".".join(str(x) for x in e.path)
        loc = path if path else "<root>"
        out.append(f"schema error at {loc}: {e.message}")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate release index and referenced paths")
    parser.add_argument("--index", default="release/index.json", help="release index path")
    parser.add_argument("--schema", default="release/schema/index.schema.json", help="json schema path")
    parser.add_argument("--repo-root", default=".", help="repository root for path checks")
    parser.add_argument("--allow-missing-jsonschema", action="store_true", help="do not fail when jsonschema package is missing")
    parser.add_argument(
        "--strict-github-release-paths",
        action="store_true",
        help="require manifest/checksum/quality report local paths to exist even for github_release products",
    )
    args = parser.parse_args()

    index_path = Path(args.index)
    repo_root = Path(args.repo_root)

    if not index_path.exists():
        raise SystemExit(f"[ERROR] missing index: {index_path}")

    index = load_json(index_path)

    errors: List[str] = []
    errors.extend(basic_structure_errors(index))

    schema_errs = schema_errors(index, Path(args.schema))
    if schema_errs:
        missing_pkg_only = (
            len(schema_errs) == 1
            and schema_errs[0].startswith("jsonschema is not installed")
        )
        if not (missing_pkg_only and args.allow_missing_jsonschema):
            errors.extend(schema_errs)

    errors.extend(
        path_integrity_errors(
            index,
            repo_root,
            strict_github_release_paths=args.strict_github_release_paths,
        )
    )

    if errors:
        print(f"[FAIL] release index validation errors: {len(errors)}")
        for e in errors[:50]:
            print(" -", e)
        return 1

    print("[PASS] release index is valid and referenced paths exist")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
