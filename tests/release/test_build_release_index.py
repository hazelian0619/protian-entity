from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "build_release_index.py"


def test_build_release_index(tmp_path: Path) -> None:
    products = tmp_path / "products"
    products.mkdir(parents=True)

    # product A with manifest present
    pa = products / "alpha"
    pa.mkdir()
    (pa / "product.json").write_text(json.dumps({"product": "alpha", "owner": "team-a"}), encoding="utf-8")
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "alpha_manifest.json").write_text("{}", encoding="utf-8")
    (pa / "current.json").write_text(
        json.dumps(
            {
                "product": "alpha",
                "display_name": "Alpha",
                "status": "stable",
                "latest": {
                    "version": "v1",
                    "distribution_mode": "repository_snapshot",
                    "manifest_path": "data/alpha_manifest.json",
                    "quality_reports": ["x.json", "y.json"],
                },
            }
        ),
        encoding="utf-8",
    )

    # product B without product.json
    pb = products / "beta"
    pb.mkdir()
    (pb / "current.json").write_text(
        json.dumps(
            {
                "product": "beta",
                "latest": {
                    "version": "v2",
                    "distribution_mode": "github_release",
                    "manifest_path": "missing/manifest.json",
                },
            }
        ),
        encoding="utf-8",
    )

    out = tmp_path / "release" / "index.json"
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--products-root",
        str(products),
        "--out",
        str(out),
    ]
    subprocess.run(cmd, cwd=tmp_path, check=True)

    idx = json.loads(out.read_text(encoding="utf-8"))
    assert idx["schema_version"] == "1.0.0"
    assert len(idx["products"]) == 2

    alpha = next(p for p in idx["products"] if p["product"] == "alpha")
    beta = next(p for p in idx["products"] if p["product"] == "beta")

    assert alpha["checks"]["manifest_exists"] is True
    assert alpha["checks"]["quality_report_count"] == 2
    assert alpha["metadata"]["owner"] == "team-a"

    assert beta["checks"]["manifest_exists"] is False
