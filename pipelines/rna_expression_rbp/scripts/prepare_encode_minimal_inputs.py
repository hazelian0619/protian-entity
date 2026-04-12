#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean
from typing import Dict, Iterable, List, Tuple


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Prepare minimal ENCODE-derived inputs for RNA expression/RBP pipeline.")
    p.add_argument("--rna-master", type=Path, default=Path("data/output/rna_master_v1.tsv"))
    p.add_argument("--expression-quant", type=Path, default=Path(".tmp/ENCFF937GQE.tsv"))
    p.add_argument("--rbp-peak-bed", type=Path, default=Path(".tmp/ENCFF593RED.bed.gz"))
    p.add_argument(
        "--gtf",
        type=Path,
        default=Path("data/raw/rna/ensembl/Homo_sapiens.GRCh38.115.chr.gtf.gz"),
    )
    p.add_argument("--expression-out", type=Path, default=Path("data/raw/rna/encode/expression_evidence.tsv.gz"))
    p.add_argument("--rbp-out", type=Path, default=Path("data/raw/rna/encode/rbp_sites.tsv.gz"))
    p.add_argument(
        "--report",
        type=Path,
        default=Path("pipelines/rna_expression_rbp/reports/rna_expression_rbp_v1.prep.metrics.json"),
    )
    p.add_argument("--bin-size", type=int, default=1_000_000)
    return p.parse_args()


def ensure_exists(path: Path) -> None:
    if not path.exists():
        raise SystemExit(f"[ERROR] missing input: {path}")


def parse_attrs(attr: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for m in re.finditer(r'([A-Za-z0-9_]+)\s+"([^"]+)"', attr):
        out[m.group(1)] = m.group(2)
    return out


def normalize_chr(chrom: str) -> str:
    c = (chrom or "").strip()
    if c == "":
        return c
    if c.startswith("chr"):
        if c == "chrMT":
            return "chrM"
        return c
    if c == "MT":
        return "chrM"
    return f"chr{c}"


def evidence_from_value(v: float) -> str:
    if v >= 100:
        return "high"
    if v >= 10:
        return "medium"
    return "low"


def evidence_from_peaks(n: int) -> str:
    if n >= 10:
        return "high"
    if n >= 3:
        return "medium"
    return "low"


def load_master_maps(rna_master: Path) -> Tuple[set[str], Dict[str, List[str]]]:
    master_ids: set[str] = set()
    gene_to_rna: Dict[str, List[str]] = defaultdict(list)
    with rna_master.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None or "rna_id" not in set(reader.fieldnames):
            raise SystemExit("[ERROR] rna_master 缺少 rna_id 列")
        for row in reader:
            rid = (row.get("rna_id") or "").strip()
            if not rid:
                continue
            master_ids.add(rid)
            if (row.get("rna_type") or "").strip().lower() != "mrna":
                continue
            g = (row.get("ensembl_gene_id") or "").strip()
            if not g.startswith("ENSG"):
                continue
            g0 = g.split(".", 1)[0]
            gene_to_rna[g0].append(rid)
    return master_ids, gene_to_rna


def load_gene_quant(path: Path) -> Dict[str, float]:
    out: Dict[str, float] = {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            gid = cols[0].split(".", 1)[0]
            if not gid.startswith("ENSG"):
                continue
            nums: List[float] = []
            for x in cols[1:]:
                try:
                    nums.append(float(x))
                except Exception:
                    continue
            if not nums:
                continue
            out[gid] = float(mean(nums))
    return out


def write_expression_input(
    out_path: Path,
    gene_quant: Dict[str, float],
    gene_to_rna: Dict[str, List[str]],
) -> Dict[str, int]:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "rna_id",
        "biosample_tissue",
        "biosample_cell_type",
        "expression_value",
        "expression_unit",
        "assay",
        "sample_id",
        "experiment_id",
        "source_dataset",
        "source_file",
        "evidence_level",
    ]
    rows = 0
    mapped_genes = 0
    with gzip.open(out_path, "wt", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        w.writeheader()
        for gene_id, value in gene_quant.items():
            targets = gene_to_rna.get(gene_id, [])
            if not targets:
                continue
            mapped_genes += 1
            lvl = evidence_from_value(value)
            val = f"{value:.6f}".rstrip("0").rstrip(".")
            for rid in targets:
                w.writerow(
                    {
                        "rna_id": rid,
                        "biosample_tissue": "blood",
                        "biosample_cell_type": "CD4-positive, alpha-beta T cell",
                        "expression_value": val,
                        "expression_unit": "counts_mean",
                        "assay": "total RNA-seq",
                        "sample_id": "ENCFF937GQE",
                        "experiment_id": "ENCSR140GZF",
                        "source_dataset": "ENCODE",
                        "source_file": "ENCFF937GQE.tsv",
                        "evidence_level": lvl,
                    }
                )
                rows += 1
    return {"expression_rows": rows, "matched_genes": mapped_genes, "quant_genes": len(gene_quant)}


def build_tx_bins(gtf: Path, master_ids: set[str], bin_size: int):
    bins: Dict[str, Dict[int, List[Tuple[int, int, str]]]] = defaultdict(lambda: defaultdict(list))
    tx_count = 0
    with gzip.open(gtf, "rt", encoding="utf-8", errors="replace", newline="") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            attrs = parse_attrs(cols[8])
            tid = attrs.get("transcript_id", "")
            if not tid.startswith("ENST"):
                continue
            rid = f"{tid.split('.', 1)[0]}_9606"
            if rid not in master_ids:
                continue
            chrom = normalize_chr(cols[0])
            try:
                start = int(cols[3]) - 1
                end = int(cols[4])
            except Exception:
                continue
            if end <= start:
                continue
            b0 = start // bin_size
            b1 = (end - 1) // bin_size
            for b in range(b0, b1 + 1):
                bins[chrom][b].append((start, end, rid))
            tx_count += 1
    return bins, tx_count


def write_rbp_input(
    out_path: Path,
    bed_path: Path,
    bins,
    bin_size: int,
) -> Dict[str, int]:
    agg: Dict[str, Dict[str, float]] = {}
    peaks_total = 0
    overlaps_total = 0
    with gzip.open(bed_path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            chrom = normalize_chr(cols[0])
            try:
                start = int(cols[1])
                end = int(cols[2])
            except Exception:
                continue
            if end <= start:
                continue
            peaks_total += 1

            try:
                signal = float(cols[6]) if len(cols) > 6 else 0.0
            except Exception:
                signal = 0.0
            try:
                qvalue = float(cols[8]) if len(cols) > 8 else 0.0
            except Exception:
                qvalue = 0.0

            c_bins = bins.get(chrom)
            if not c_bins:
                continue
            b0 = start // bin_size
            b1 = (end - 1) // bin_size
            seen_rids: set[str] = set()
            for b in range(b0, b1 + 1):
                for tx_start, tx_end, rid in c_bins.get(b, []):
                    if rid in seen_rids:
                        continue
                    if start < tx_end and end > tx_start:
                        seen_rids.add(rid)
                        overlaps_total += 1
                        cur = agg.get(rid)
                        if cur is None:
                            agg[rid] = {"peak_count": 1.0, "max_score": signal, "min_q": qvalue}
                        else:
                            cur["peak_count"] += 1.0
                            cur["max_score"] = max(cur["max_score"], signal)
                            cur["min_q"] = min(cur["min_q"], qvalue)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "rna_id",
        "rbp_symbol",
        "biosample_cell_type",
        "assay",
        "peak_count",
        "score",
        "evalue",
        "sample_id",
        "experiment_id",
        "source_dataset",
        "source_file",
        "evidence_level",
    ]
    rows = 0
    with gzip.open(out_path, "wt", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        w.writeheader()
        for rid in sorted(agg.keys()):
            peak_count = int(agg[rid]["peak_count"])
            lvl = evidence_from_peaks(peak_count)
            score = f"{agg[rid]['max_score']:.6f}".rstrip("0").rstrip(".")
            evalue = f"{agg[rid]['min_q']:.6f}".rstrip("0").rstrip(".")
            w.writerow(
                {
                    "rna_id": rid,
                    "rbp_symbol": "TARDBP",
                    "biosample_cell_type": "K562",
                    "assay": "eCLIP",
                    "peak_count": str(peak_count),
                    "score": score,
                    "evalue": evalue,
                    "sample_id": "ENCFF593RED",
                    "experiment_id": "ENCSR720BJU",
                    "source_dataset": "ENCODE",
                    "source_file": "ENCFF593RED.bed.gz",
                    "evidence_level": lvl,
                }
            )
            rows += 1

    return {
        "rbp_rows": rows,
        "unique_rna_with_peaks": len(agg),
        "peaks_total": peaks_total,
        "peak_transcript_overlaps": overlaps_total,
    }


def write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    ensure_exists(args.rna_master)
    ensure_exists(args.expression_quant)
    ensure_exists(args.rbp_peak_bed)
    ensure_exists(args.gtf)

    master_ids, gene_to_rna = load_master_maps(args.rna_master)
    gene_quant = load_gene_quant(args.expression_quant)
    expr_stats = write_expression_input(args.expression_out, gene_quant, gene_to_rna)

    bins, tx_with_coords = build_tx_bins(args.gtf, master_ids, args.bin_size)
    rbp_stats = write_rbp_input(args.rbp_out, args.rbp_peak_bed, bins, args.bin_size)

    report = {
        "pipeline": "rna_expression_rbp_v1_preparation",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": {
            "rna_master": str(args.rna_master),
            "expression_quant": str(args.expression_quant),
            "rbp_peak_bed": str(args.rbp_peak_bed),
            "gtf": str(args.gtf),
        },
        "output": {
            "expression_input": str(args.expression_out),
            "rbp_input": str(args.rbp_out),
        },
        "stats": {
            "master_ids": len(master_ids),
            "master_genes_for_mrna": len(gene_to_rna),
            "transcripts_with_gtf_coords": tx_with_coords,
            **expr_stats,
            **rbp_stats,
        },
        "notes": [
            "Expression input is ENCODE RNA-seq gene quant (ENCSR140GZF / ENCFF937GQE) expanded to transcript-level via ensembl_gene_id join.",
            "RBP input is ENCODE eCLIP peak overlap (ENCSR720BJU / ENCFF593RED) intersected with Ensembl transcript intervals (GRCh38.115).",
        ],
    }
    write_json(args.report, report)
    print(f"[OK] expression input -> {args.expression_out}")
    print(f"[OK] rbp input -> {args.rbp_out}")
    print(f"[OK] report -> {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

