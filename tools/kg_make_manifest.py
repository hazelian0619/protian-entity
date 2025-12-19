#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def git_commit(repo_root: Path) -> str | None:
    try:
        out = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(repo_root))
        return out.decode("utf-8").strip()
    except Exception:
        return None


def find_repo_root(start: Path) -> Path | None:
    for parent in [start] + list(start.parents):
        if (parent / ".git").exists():
            return parent
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-version", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("files", nargs="+", type=Path)
    args = ap.parse_args()

    now = datetime.now(timezone.utc).isoformat()
    repo_root = find_repo_root(Path.cwd())
    commit = git_commit(repo_root) if repo_root else None

    artifacts: List[Dict[str, Any]] = []
    for f in args.files:
        st = f.stat()
        artifacts.append(
            {
                "path": str(f),
                "size_bytes": st.st_size,
                "sha256": sha256_file(f),
                "modified_time": datetime.fromtimestamp(st.st_mtime, tz=timezone.utc).isoformat(),
            }
        )

    manifest: Dict[str, Any] = {
        "data_version": args.data_version,
        "created_at": now,
        "git_commit": commit,
        "artifacts": artifacts,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"[OK] manifest -> {args.out} (artifacts={len(artifacts)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
