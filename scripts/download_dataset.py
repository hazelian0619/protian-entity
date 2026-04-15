#!/usr/bin/env python3
"""Download datasets by product/version using release/index.json.

Supports:
- GitHub Release assets download
- Repository snapshot direct URLs (for non-release products)
- Checksum verification from manifest and/or SHA256SUMS
- Interaction chunk merge (*.part.000 ...)
- Optional decompression (.gz, .zst)
"""

from __future__ import annotations

import argparse
import fnmatch
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


GITHUB_API = "https://api.github.com"


def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download_file(url: str, out_path: Path, timeout: int = 60) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "bio-entity-kg-foundation-downloader/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp, out_path.open("wb") as f:
        shutil.copyfileobj(resp, f)


def github_release_assets(repo: str, tag: str) -> List[Dict[str, Any]]:
    url = f"{GITHUB_API}/repos/{repo}/releases/tags/{tag}"
    req = urllib.request.Request(url, headers={"User-Agent": "bio-entity-kg-foundation-downloader/1.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        payload = json.loads(resp.read().decode("utf-8"))
    return payload.get("assets", [])


def manifest_checksum_map(manifest_path: Optional[Path]) -> Dict[str, str]:
    checksums: Dict[str, str] = {}
    if not manifest_path or not manifest_path.exists():
        return checksums

    data = load_json(manifest_path)

    # RNA manifest
    for key in ["tables", "reports"]:
        items = data.get(key)
        if isinstance(items, list):
            for item in items:
                if isinstance(item, dict) and item.get("name") and item.get("sha256"):
                    checksums[item["name"]] = str(item["sha256"])

    # Molecule manifest
    items = data.get("assets")
    if isinstance(items, list):
        for item in items:
            if isinstance(item, dict) and item.get("name") and item.get("sha256"):
                checksums[item["name"]] = str(item["sha256"])

    # Interaction release assets manifest
    items = data.get("items")
    if isinstance(items, list):
        for item in items:
            if not isinstance(item, dict):
                continue
            for asset in item.get("release_assets", []):
                if isinstance(asset, dict) and asset.get("name") and asset.get("sha256"):
                    checksums[asset["name"]] = str(asset["sha256"])

    return checksums


def parse_sha256sums(path: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not path.exists():
        return out
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            out[parts[1].lstrip("*")] = parts[0]
    return out


def merge_chunk_files(directory: Path) -> List[Path]:
    pattern = re.compile(r"^(?P<base>.+)\.part\.(?P<num>\d+)$")
    groups: Dict[str, List[Tuple[int, Path]]] = {}

    for p in directory.iterdir():
        if not p.is_file():
            continue
        m = pattern.match(p.name)
        if not m:
            continue
        base = m.group("base")
        num = int(m.group("num"))
        groups.setdefault(base, []).append((num, p))

    merged: List[Path] = []
    for base, parts in groups.items():
        parts_sorted = sorted(parts, key=lambda x: x[0])
        out_path = directory / base
        with out_path.open("wb") as out_f:
            for _, part_path in parts_sorted:
                with part_path.open("rb") as in_f:
                    shutil.copyfileobj(in_f, out_f)
        merged.append(out_path)
    return merged


def decompress_files(directory: Path) -> List[Path]:
    created: List[Path] = []
    zstd_exists = shutil.which("zstd") is not None

    for p in sorted(directory.iterdir()):
        if not p.is_file():
            continue

        if p.name.endswith(".gz"):
            out = directory / p.name[:-3]
            with gzip.open(p, "rb") as src, out.open("wb") as dst:
                shutil.copyfileobj(src, dst)
            created.append(out)

        elif p.name.endswith(".zst"):
            out = directory / p.name[:-4]
            if zstd_exists:
                subprocess.run(["zstd", "-d", "-f", str(p), "-o", str(out)], check=True)
                created.append(out)
            else:
                print(f"[WARN] zstd not found; skip decompress: {p.name}")

    return created


def filter_assets(assets: List[Dict[str, Any]], patterns: List[str]) -> List[Dict[str, Any]]:
    if not patterns:
        return assets
    out = []
    for a in assets:
        name = a.get("name", "")
        if any(fnmatch.fnmatch(name, pat) for pat in patterns):
            out.append(a)
    return out


def load_index_product(index_path: Path, product: str) -> Dict[str, Any]:
    idx = load_json(index_path)
    for p in idx.get("products", []):
        if p.get("product") == product:
            return p
    raise SystemExit(f"[ERROR] product not found in index: {product}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Download bio-entity KG dataset assets")
    parser.add_argument("--product", required=True, choices=["protein", "rna", "molecule", "interaction"])
    parser.add_argument("--version", default="latest", help="Version to fetch (currently supports latest)")
    parser.add_argument("--repo", default="hazelian0619/protian-entity", help="GitHub repo owner/name")
    parser.add_argument("--index", default="release/index.json", help="Unified release index")
    parser.add_argument("--dest", default="downloads", help="Destination directory")
    parser.add_argument("--asset-glob", action="append", default=[], help="Filter assets by glob pattern")
    parser.add_argument("--skip-verify", action="store_true", help="Skip sha256 verification")
    parser.add_argument("--merge-chunks", action="store_true", help="Merge chunked assets (*.part.000...) after download")
    parser.add_argument("--decompress", action="store_true", help="Decompress .gz/.zst after download")
    args = parser.parse_args()

    if args.version != "latest":
        raise SystemExit("[ERROR] only --version latest is currently supported")

    product_meta = load_index_product(Path(args.index), args.product)
    latest = product_meta.get("latest", {})
    version = latest.get("version", "latest")

    out_dir = Path(args.dest) / args.product / version
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = Path(latest.get("manifest_path")) if latest.get("manifest_path") else None
    manifest_sha = manifest_checksum_map(manifest_path)

    downloaded: List[Path] = []

    if latest.get("distribution_mode") == "github_release":
        tag = latest.get("release_tag")
        if not tag:
            raise SystemExit(f"[ERROR] product {args.product} missing release_tag")

        assets = github_release_assets(args.repo, tag)
        assets = filter_assets(assets, args.asset_glob)
        if not assets:
            print(f"[WARN] no assets matched for {args.product}:{tag}")

        print(f"[INFO] downloading {len(assets)} assets from release {tag}")
        for asset in assets:
            name = asset["name"]
            url = asset["browser_download_url"]
            target = out_dir / name
            print(f"[DL] {name}")
            download_file(url, target)
            downloaded.append(target)

    else:
        artifacts = latest.get("artifacts", [])
        if args.asset_glob:
            artifacts = [
                a for a in artifacts
                if any(fnmatch.fnmatch(str(a.get("name", "")), pat) for pat in args.asset_glob)
            ]
        if not artifacts:
            raise SystemExit(f"[ERROR] product {args.product} has no artifacts in current.json")

        print(f"[INFO] downloading {len(artifacts)} repository artifacts")
        for a in artifacts:
            name = a["name"]
            url = a["download_url"]
            target = out_dir / name
            print(f"[DL] {name}")
            download_file(url, target)
            downloaded.append(target)

    # optional checksum verification
    if not args.skip_verify:
        sha_failures = []

        # 1) manifest checksums (preferred)
        expected = dict(manifest_sha)

        # 2) SHA256SUMS from downloaded assets (if present)
        sha_sums_file = out_dir / "SHA256SUMS.txt"
        expected.update(parse_sha256sums(sha_sums_file))

        if expected:
            for f in downloaded:
                exp = expected.get(f.name)
                if not exp:
                    continue
                got = sha256_file(f)
                if got.lower() != exp.lower():
                    sha_failures.append((f.name, exp, got))

        if sha_failures:
            print("[FAIL] checksum verification failed:")
            for name, exp, got in sha_failures[:20]:
                print(f"  - {name}: expected={exp} got={got}")
            return 2
        print("[OK] checksum verification passed (where checksum available)")

    if args.merge_chunks:
        merged = merge_chunk_files(out_dir)
        if merged:
            print(f"[OK] merged {len(merged)} chunked files")
            for p in merged:
                print("  -", p.name)

    if args.decompress:
        created = decompress_files(out_dir)
        if created:
            print(f"[OK] decompressed {len(created)} files")
            for p in created:
                print("  -", p.name)

    print(f"[DONE] downloaded assets to {out_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except urllib.error.HTTPError as e:
        print(f"[ERROR] HTTP {e.code}: {e.reason}", file=sys.stderr)
        raise
