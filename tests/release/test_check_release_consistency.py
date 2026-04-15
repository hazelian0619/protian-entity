from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "check_release_consistency.py"


def _write_tsv(path: Path, rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write("c1\tc2\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


def test_strict_local_pass(tmp_path: Path) -> None:
    table = tmp_path / "data" / "table.tsv"
    _write_tsv(table, [["a", "1"], ["b", "2"]])

    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps({"tables": [{"name": "table", "source_path": "data/table.tsv", "rows": 2}]}),
        encoding="utf-8",
    )

    val_report = tmp_path / "pipelines" / "x" / "reports" / "table.validation.json"
    val_report.parent.mkdir(parents=True, exist_ok=True)
    val_report.write_text(
        json.dumps({"table": "data/table.tsv", "row_count": 2, "passed": True}),
        encoding="utf-8",
    )

    index = tmp_path / "release" / "index.json"
    index.parent.mkdir(parents=True, exist_ok=True)
    index.write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "generated_at_utc": "2026-01-01T00:00:00Z",
                "products": [
                    {
                        "product": "alpha",
                        "latest": {
                            "version": "v1",
                            "distribution_mode": "repository_snapshot",
                            "manifest_path": "manifest.json",
                        },
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    out = tmp_path / "release" / "consistency_report.json"
    cmd = [sys.executable, str(SCRIPT), "--index", "release/index.json", "--out", "release/consistency_report.json"]
    subprocess.run(cmd, cwd=tmp_path, check=True)

    report = json.loads(out.read_text(encoding="utf-8"))
    assert report["summary"]["status"] == "PASS"
    assert report["summary"]["issues"] == 0


def test_release_assets_fail_on_missing_sha(tmp_path: Path) -> None:
    manifest = tmp_path / "interaction_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "items": [
                    {
                        "table": "ppi_edges",
                        "release_available": True,
                        "release_assets": [{"name": "ppi_edges.tsv.zst"}],
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    index = tmp_path / "release" / "index.json"
    index.parent.mkdir(parents=True, exist_ok=True)
    index.write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "generated_at_utc": "2026-01-01T00:00:00Z",
                "products": [
                    {
                        "product": "interaction",
                        "latest": {
                            "version": "v1",
                            "distribution_mode": "github_release",
                            "consistency_mode": "release_assets",
                            "manifest_path": "interaction_manifest.json",
                        },
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    cmd = [sys.executable, str(SCRIPT), "--index", "release/index.json", "--out", "release/out.json"]
    proc = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True)

    assert proc.returncode == 1
    report = json.loads((tmp_path / "release" / "out.json").read_text(encoding="utf-8"))
    assert report["summary"]["status"] == "FAIL"
    assert report["summary"]["issues"] >= 1
