#!/usr/bin/env python3
"""Build a unified machine-readable release index from products/*/current.json."""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path
from typing import Any, Dict, List


def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Build unified release index")
    parser.add_argument("--products-root", default="products", help="Products root directory")
    parser.add_argument("--out", default="release/index.json", help="Output index path")
    args = parser.parse_args()

    products_root = Path(args.products_root)
    out_path = Path(args.out)

    products: List[Dict[str, Any]] = []
    for current_path in sorted(products_root.glob("*/current.json")):
        product_dir = current_path.parent
        current = load_json(current_path)
        product_name = current.get("product", product_dir.name)

        product_meta_path = product_dir / "product.json"
        product_meta: Dict[str, Any] = {}
        if product_meta_path.exists():
            product_meta = load_json(product_meta_path)

        latest = current.get("latest", {})
        manifest_path = latest.get("manifest_path")
        manifest_exists = bool(manifest_path and Path(manifest_path).exists())

        products.append(
            {
                "product": product_name,
                "display_name": current.get("display_name", product_name),
                "status": current.get("status", "unknown"),
                "metadata": product_meta,
                "current_path": str(current_path),
                "latest": latest,
                "checks": {
                    "manifest_exists": manifest_exists,
                    "quality_report_count": len(latest.get("quality_reports", [])),
                },
            }
        )

    index: Dict[str, Any] = {
        "schema_version": "1.0.0",
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "products": products,
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(index, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"[OK] wrote {out_path} with {len(products)} products")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
