#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

RFAM_RE = re.compile(r"^RF\d{5}$")
FLOAT_RE = re.compile(r"^-?\d+(?:\.\d+)?$")
FETCH_DATE = date.today().isoformat()

COV_COLUMNS = ["rfam_id", "cm_file", "ga_threshold", "source_version", "fetch_date"]
PRED_COLUMNS = [
    "rna_id",
    "model_tool",
    "structure_file",
    "confidence",
    "energy",
    "source_version",
    "fetch_date",
]


def open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def detect_delimiter(path: Path) -> str:
    with open_maybe_gzip(path) as f:
        head = f.read(8192)
    return "\t" if head.count("\t") >= head.count(",") else ","


def write_tsv(path: Path, rows: List[Dict[str, str]], columns: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: str(row.get(c, "")).strip() for c in columns})


def write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def normalize_numeric(value: str) -> str:
    v = str(value or "").strip()
    if v == "":
        return ""
    if FLOAT_RE.match(v):
        return v
    try:
        return f"{float(v):.4f}".rstrip("0").rstrip(".")
    except Exception:
        return ""


def to_repo_relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(Path.cwd().resolve()))
    except Exception:
        return str(path)


def parse_family_row(line: str) -> Optional[Tuple[str, str]]:
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 9:
        return None
    rfam_id = cols[0].strip()
    if not RFAM_RE.match(rfam_id):
        return None

    # Rfam family.txt 常见为 [6,7,8] = gathering/trusted/noise cutoff
    ga = ""
    for idx in (6, 7, 8):
        if idx < len(cols):
            ga_candidate = normalize_numeric(cols[idx])
            if ga_candidate:
                ga = ga_candidate
                break
    return rfam_id, ga


def detect_cm_source(cm_bundle: Path, cm_dir: Path) -> Dict[str, object]:
    if cm_bundle.exists() and cm_bundle.is_file():
        return {"mode": "bundle", "bundle": cm_bundle, "cm_map": {}}

    cm_map: Dict[str, Path] = {}
    if cm_dir.exists() and cm_dir.is_dir():
        for p in sorted(cm_dir.rglob("*")):
            if not p.is_file():
                continue
            m = re.search(r"(RF\d{5})", p.name)
            if not m:
                continue
            low_name = p.name.lower()
            if not (low_name.endswith(".cm") or low_name.endswith(".cm.gz")):
                continue
            rfam_id = m.group(1)
            if rfam_id not in cm_map:
                cm_map[rfam_id] = p
        if cm_map:
            return {"mode": "per_family", "bundle": None, "cm_map": cm_map}

    return {"mode": "missing", "bundle": None, "cm_map": {}}


def build_covariance_rows(
    family_path: Path,
    cm_info: Dict[str, object],
    source_version: str,
    fetch_date: str,
    limit: Optional[int],
) -> Tuple[List[Dict[str, str]], Dict[str, int]]:
    rows: List[Dict[str, str]] = []
    seen: set[str] = set()
    missing_cm = 0
    parsed_family_rows = 0

    mode = str(cm_info.get("mode"))
    bundle = cm_info.get("bundle")
    cm_map = cm_info.get("cm_map") or {}

    with open_maybe_gzip(family_path) as f:
        for line in f:
            parsed = parse_family_row(line)
            if not parsed:
                continue
            parsed_family_rows += 1
            rfam_id, ga = parsed
            if rfam_id in seen:
                continue
            seen.add(rfam_id)

            cm_ref = ""
            if mode == "bundle" and isinstance(bundle, Path):
                cm_ref = f"{to_repo_relative(bundle)}#{rfam_id}"
            elif mode == "per_family":
                cm_path = cm_map.get(rfam_id)
                if isinstance(cm_path, Path):
                    cm_ref = to_repo_relative(cm_path)

            if not cm_ref:
                missing_cm += 1
                continue

            rows.append(
                {
                    "rfam_id": rfam_id,
                    "cm_file": cm_ref,
                    "ga_threshold": ga,
                    "source_version": source_version,
                    "fetch_date": fetch_date,
                }
            )

            if limit is not None and len(rows) >= limit:
                break

    stats = {
        "parsed_family_rows": parsed_family_rows,
        "unique_rfam_ids": len(seen),
        "rows_emitted": len(rows),
        "rows_missing_cm": missing_cm,
    }
    return rows, stats


def infer_model_tool(path: Path) -> str:
    low_parts = [x.lower() for x in path.parts]
    low_name = path.name.lower()
    if any("rnafold" in p for p in low_parts) or "rnafold" in low_name:
        return "RNAfold"
    if any("rhofold" in p for p in low_parts) or "rhofold" in low_name or "rho-fold" in low_name:
        return "RhoFold"
    if path.suffix.lower() in {".dbn", ".ct", ".bpseq"}:
        return "RNAfold"
    if path.suffix.lower() in {".pdb", ".cif", ".mmcif"}:
        return "RhoFold"
    return "Unknown"


def infer_rna_id(path: Path) -> str:
    name = path.stem
    name = re.sub(r"\.(dbn|ct|bpseq|pdb|cif|mmcif)$", "", name, flags=re.IGNORECASE)
    name = re.sub(r"(?:[_\-.](?:rnafold|rhofold|rhofold))$", "", name, flags=re.IGNORECASE)
    return name.strip("_-. ")


def parse_energy_from_dbn(path: Path) -> str:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as f:
            lines = [next(f).strip() for _ in range(5)]
    except Exception:
        return ""

    joined = " ".join(lines)
    # RNAfold 常见: ".... (-3.40)"
    m = re.search(r"\(([+-]?\d+(?:\.\d+)?)\)", joined)
    if m:
        return normalize_numeric(m.group(1))
    return ""


def load_sidecar_meta(structure_path: Path) -> Dict[str, str]:
    candidates = [
        structure_path.with_suffix(structure_path.suffix + ".json"),
        structure_path.with_suffix(".json"),
        structure_path.with_name(structure_path.name + ".meta.json"),
    ]
    for c in candidates:
        if not c.exists() or not c.is_file():
            continue
        try:
            obj = json.loads(c.read_text(encoding="utf-8"))
            if not isinstance(obj, dict):
                continue
            return {k: "" if v is None else str(v) for k, v in obj.items()}
        except Exception:
            continue
    return {}


def parse_predicted_manifest(path: Path, default_source_version: str, fetch_date: str) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    delim = detect_delimiter(path)
    with open_maybe_gzip(path) as f:
        reader = csv.DictReader(f, delimiter=delim)
        if reader.fieldnames is None:
            return rows
        for row in reader:
            rna_id = str(row.get("rna_id", "")).strip()
            model_tool = str(row.get("model_tool", "")).strip()
            structure_file = str(row.get("structure_file", "")).strip()
            if not (rna_id and model_tool and structure_file):
                continue
            struct_path = Path(structure_file)
            if not struct_path.is_absolute():
                struct_path = (path.parent / struct_path).resolve()
            rows.append(
                {
                    "rna_id": rna_id,
                    "model_tool": model_tool,
                    "structure_file": to_repo_relative(struct_path),
                    "confidence": normalize_numeric(str(row.get("confidence", ""))),
                    "energy": normalize_numeric(str(row.get("energy", ""))),
                    "source_version": str(row.get("source_version", "")).strip() or default_source_version,
                    "fetch_date": str(row.get("fetch_date", "")).strip() or fetch_date,
                }
            )
    return rows


def discover_predicted_structures(
    root: Path,
    default_source_version: str,
    fetch_date: str,
    limit: Optional[int],
) -> List[Dict[str, str]]:
    if not root.exists() or not root.is_dir():
        return []

    allowed = {".dbn", ".ct", ".bpseq", ".pdb", ".cif", ".mmcif"}
    rows: List[Dict[str, str]] = []

    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        if p.suffix.lower() not in allowed:
            continue

        meta = load_sidecar_meta(p)

        rna_id = (meta.get("rna_id") or "").strip() or infer_rna_id(p)
        if not rna_id:
            continue

        model_tool = (meta.get("model_tool") or "").strip() or infer_model_tool(p)
        confidence = normalize_numeric(meta.get("confidence", ""))
        energy = normalize_numeric(meta.get("energy", ""))
        if energy == "" and p.suffix.lower() == ".dbn":
            energy = parse_energy_from_dbn(p)

        rows.append(
            {
                "rna_id": rna_id,
                "model_tool": model_tool,
                "structure_file": to_repo_relative(p),
                "confidence": confidence,
                "energy": energy,
                "source_version": (meta.get("source_version") or "").strip() or default_source_version,
                "fetch_date": (meta.get("fetch_date") or "").strip() or fetch_date,
            }
        )

        if limit is not None and len(rows) >= limit:
            break

    return rows


def detect_predicted_rows(
    predicted_root: Path,
    predicted_manifest: Path,
    source_version: str,
    fetch_date: str,
    limit: Optional[int],
) -> Tuple[List[Dict[str, str]], str]:
    if predicted_manifest.exists() and predicted_manifest.is_file():
        rows = parse_predicted_manifest(predicted_manifest, source_version, fetch_date)
        if limit is not None:
            rows = rows[:limit]
        if rows:
            return rows, "manifest"

    rows = discover_predicted_structures(predicted_root, source_version, fetch_date, limit)
    if rows:
        return rows, "scan"

    return [], "missing"


def make_download_checklist(args: argparse.Namespace) -> List[Dict[str, str]]:
    return [
        {
            "file": "Rfam family metadata",
            "required_path": str(args.rfam_family),
            "suggested_source": "https://ftp.ebi.ac.uk/pub/databases/Rfam/current_database_files/family.txt.gz",
            "sha256_command": f"sha256sum {args.rfam_family}",
        },
        {
            "file": "Rfam covariance model bundle",
            "required_path": str(args.rfam_cm_bundle),
            "suggested_source": "https://ftp.ebi.ac.uk/pub/databases/Rfam/current_database_files/Rfam.cm.gz",
            "sha256_command": f"sha256sum {args.rfam_cm_bundle}",
        },
        {
            "file": "Rfam per-family CM directory（可替代 bundle）",
            "required_path": str(args.rfam_cm_dir),
            "suggested_source": "本地拆分的 RFxxxxx.cm / RFxxxxx.cm.gz 文件目录",
            "sha256_command": f"find {args.rfam_cm_dir} -type f -name '*.cm*' -print0 | xargs -0 sha256sum",
        },
        {
            "file": "预测结构目录",
            "required_path": str(args.predicted_root),
            "suggested_source": "放置 RNAfold/RhoFold 结构产物（.dbn/.ct/.bpseq/.pdb/.cif）",
            "sha256_command": f"find {args.predicted_root} -type f -print0 | xargs -0 sha256sum",
        },
        {
            "file": "预测结构登记清单（可选，优先）",
            "required_path": str(args.predicted_manifest),
            "suggested_source": "本地生成 manifest.tsv，至少包含 rna_id,model_tool,structure_file",
            "sha256_command": f"sha256sum {args.predicted_manifest}",
        },
    ]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build RNA covariance model index + predicted structure register")
    p.add_argument("--rfam-family", type=Path, default=Path("data/raw/rna/rfam/family.txt.gz"))
    p.add_argument("--rfam-cm-bundle", type=Path, default=Path("data/raw/rna/rfam/Rfam.cm.gz"))
    p.add_argument("--rfam-cm-dir", type=Path, default=Path("data/raw/rna/rfam/cm"))
    p.add_argument("--predicted-root", type=Path, default=Path("data/raw/rna/predicted_structures"))
    p.add_argument(
        "--predicted-manifest",
        type=Path,
        default=Path("data/raw/rna/predicted_structures/manifest.tsv"),
    )
    p.add_argument("--cov-output", type=Path, default=Path("data/output/rna_covariance_models_index_v1.tsv"))
    p.add_argument("--pred-output", type=Path, default=Path("data/output/rna_predicted_structures_v1.tsv"))
    p.add_argument(
        "--report",
        type=Path,
        default=Path("pipelines/rna_structure_models/reports/rna_structure_models_v1.metrics.json"),
    )
    p.add_argument("--source-version-rfam", default="Rfam:current")
    p.add_argument("--source-version-pred", default="RNAfold/RhoFold:local")
    p.add_argument("--fetch-date", default=FETCH_DATE)
    p.add_argument("--limit-cov", type=int, default=None)
    p.add_argument("--limit-pred", type=int, default=None)
    p.add_argument("--check-inputs", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    cm_info = detect_cm_source(args.rfam_cm_bundle, args.rfam_cm_dir)
    pred_rows_probe, pred_source = detect_predicted_rows(
        args.predicted_root,
        args.predicted_manifest,
        args.source_version_pred,
        args.fetch_date,
        limit=1,
    )

    base = {
        "pipeline": "rna_structure_models_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": {
            "rfam_family": str(args.rfam_family),
            "rfam_cm_bundle": str(args.rfam_cm_bundle),
            "rfam_cm_dir": str(args.rfam_cm_dir),
            "predicted_root": str(args.predicted_root),
            "predicted_manifest": str(args.predicted_manifest),
        },
    }

    missing: List[str] = []
    if not args.rfam_family.exists():
        missing.append(str(args.rfam_family))
    if str(cm_info.get("mode")) == "missing":
        missing.append(f"CM source missing: {args.rfam_cm_bundle} OR {args.rfam_cm_dir}")
    if pred_source == "missing" or len(pred_rows_probe) == 0:
        missing.append(f"predicted structures missing: {args.predicted_manifest} OR {args.predicted_root}")

    if missing:
        blocked = {
            **base,
            "status": "blocked_missing_inputs",
            "missing_required": missing,
            "cm_mode": cm_info.get("mode"),
            "predicted_detection": pred_source,
            "download_checklist": make_download_checklist(args),
            "message": "缺输入数据，已按约束停止。",
        }
        write_json(args.report, blocked)
        print(f"[BLOCKED] missing inputs -> {missing}")
        print(f"[BLOCKED] report -> {args.report}")
        return 2

    if args.check_inputs:
        ready = {
            **base,
            "status": "inputs_ready",
            "cm_mode": cm_info.get("mode"),
            "predicted_detection": pred_source,
            "message": "输入检查通过。",
        }
        write_json(args.report, ready)
        print(f"[OK] inputs ready -> {args.report}")
        return 0

    cov_rows, cov_stats = build_covariance_rows(
        family_path=args.rfam_family,
        cm_info=cm_info,
        source_version=args.source_version_rfam,
        fetch_date=args.fetch_date,
        limit=args.limit_cov,
    )
    pred_rows, pred_source_runtime = detect_predicted_rows(
        args.predicted_root,
        args.predicted_manifest,
        args.source_version_pred,
        args.fetch_date,
        limit=args.limit_pred,
    )

    # 去重：预测结构按 (rna_id, model_tool, structure_file)
    pred_seen = set()
    pred_rows_unique: List[Dict[str, str]] = []
    for row in pred_rows:
        key = (row["rna_id"], row["model_tool"], row["structure_file"])
        if key in pred_seen:
            continue
        pred_seen.add(key)
        pred_rows_unique.append(row)

    write_tsv(args.cov_output, cov_rows, COV_COLUMNS)
    write_tsv(args.pred_output, pred_rows_unique, PRED_COLUMNS)

    metrics = {
        **base,
        "status": "completed",
        "output": {
            "covariance_table": str(args.cov_output),
            "predicted_table": str(args.pred_output),
        },
        "covariance": {
            **cov_stats,
            "cm_mode": cm_info.get("mode"),
        },
        "predicted": {
            "rows_emitted": len(pred_rows_unique),
            "raw_rows_before_dedup": len(pred_rows),
            "detection_mode": pred_source_runtime,
            "tools": sorted({r["model_tool"] for r in pred_rows_unique}),
        },
        "acceptance": {
            "covariance_rows_gt_0": len(cov_rows) > 0,
            "predicted_rows_gt_0": len(pred_rows_unique) > 0,
        },
    }
    write_json(args.report, metrics)

    passed = metrics["acceptance"]["covariance_rows_gt_0"] and metrics["acceptance"]["predicted_rows_gt_0"]
    print(
        f"[{'PASS' if passed else 'WARN'}] covariance_rows={len(cov_rows)} "
        f"predicted_rows={len(pred_rows_unique)} cm_mode={cm_info.get('mode')} pred_mode={pred_source_runtime}"
    )
    print(f"[OK] covariance -> {args.cov_output}")
    print(f"[OK] predicted -> {args.pred_output}")
    print(f"[OK] report -> {args.report}")
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
