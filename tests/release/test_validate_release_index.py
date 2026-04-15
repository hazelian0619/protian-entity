from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "validate_release_index.py"


def test_validate_release_index_pass(tmp_path: Path) -> None:
    (tmp_path / "products" / "alpha").mkdir(parents=True)
    (tmp_path / "products" / "alpha" / "current.json").write_text("{}", encoding="utf-8")

    (tmp_path / "reports").mkdir(parents=True)
    (tmp_path / "reports" / "a.validation.json").write_text("{}", encoding="utf-8")

    (tmp_path / "manifest.json").write_text(json.dumps({"tables": []}), encoding="utf-8")

    index = {
        "schema_version": "1.0.0",
        "generated_at_utc": "2026-01-01T00:00:00Z",
        "products": [
            {
                "product": "alpha",
                "current_path": "products/alpha/current.json",
                "latest": {
                    "version": "v1",
                    "distribution_mode": "repository_snapshot",
                    "manifest_path": "manifest.json",
                    "quality_reports": ["reports/a.validation.json"],
                },
            }
        ],
    }

    (tmp_path / "release").mkdir(parents=True)
    (tmp_path / "release" / "index.json").write_text(json.dumps(index), encoding="utf-8")

    cmd = [
        sys.executable,
        str(SCRIPT),
        "--index",
        "release/index.json",
        "--schema",
        str(REPO_ROOT / "release" / "schema" / "index.schema.json"),
        "--repo-root",
        ".",
        "--allow-missing-jsonschema",
    ]
    proc = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
