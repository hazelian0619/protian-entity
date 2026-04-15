from __future__ import annotations

import importlib.util
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO_ROOT / "scripts" / "download_dataset.py"


spec = importlib.util.spec_from_file_location("download_dataset", MODULE_PATH)
assert spec and spec.loader
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_parse_sha256sums(tmp_path: Path) -> None:
    p = tmp_path / "SHA256SUMS.txt"
    p.write_text(
        "aaabbbccc111  file_a.tsv.zst\n"
        "ddd222 *file_b.tsv.gz\n",
        encoding="utf-8",
    )
    got = mod.parse_sha256sums(p)
    assert got["file_a.tsv.zst"] == "aaabbbccc111"
    assert got["file_b.tsv.gz"] == "ddd222"


def test_merge_chunk_files(tmp_path: Path) -> None:
    d = tmp_path / "dl"
    d.mkdir()
    (d / "x.tsv.zst.part.000").write_bytes(b"abc")
    (d / "x.tsv.zst.part.001").write_bytes(b"def")
    (d / "x.tsv.zst.part.002").write_bytes(b"ghi")

    merged = mod.merge_chunk_files(d)
    assert len(merged) == 1
    out = d / "x.tsv.zst"
    assert out.exists()
    assert out.read_bytes() == b"abcdefghi"


def test_manifest_checksum_map(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "tables": [{"name": "a.tsv", "sha256": "111"}],
                "assets": [{"name": "b.tsv.zst", "sha256": "222"}],
                "items": [
                    {
                        "release_assets": [
                            {"name": "c.tsv.zst", "sha256": "333"},
                        ]
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    got = mod.manifest_checksum_map(manifest)
    assert got["a.tsv"] == "111"
    assert got["b.tsv.zst"] == "222"
    assert got["c.tsv.zst"] == "333"
