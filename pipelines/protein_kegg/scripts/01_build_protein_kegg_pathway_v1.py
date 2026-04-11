#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import ssl
import time
import urllib.error
import urllib.request
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple


KEGG_INFO_URL = "https://rest.kegg.jp/info/pathway"
KEGG_CONV_URL = "https://rest.kegg.jp/conv/hsa/ncbi-geneid"
KEGG_LINK_URL = "https://rest.kegg.jp/link/pathway/hsa"
KEGG_LIST_URL = "https://rest.kegg.jp/list/pathway/hsa"


@dataclass
class BuildStats:
    master_rows: int = 0
    master_uniprot_non_empty: int = 0
    master_rows_with_parseable_ncbi: int = 0
    distinct_ncbi_ids: int = 0
    ncbi_ids_mapped_to_kegg_gene: int = 0
    ncbi_ids_without_kegg_gene: int = 0
    kegg_genes_with_pathway: int = 0
    mapped_uniprot_count: int = 0
    output_rows: int = 0
    output_unique_pathways: int = 0
    output_unique_kegg_genes: int = 0
    pathway_name_fallback_count: int = 0


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def utc_today_date() -> str:
    return datetime.now(timezone.utc).date().isoformat()


def _atomic_write_tsv(path: Path, header: List[str], rows: Iterable[List[str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    n = 0
    with tmp.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for row in rows:
            w.writerow(row)
            n += 1
    tmp.replace(path)
    return n


def _fetch_text_with_cache(
    url: str,
    cache_path: Path,
    timeout_sec: int,
    retries: int,
    force_refresh: bool,
) -> str:
    if cache_path.exists() and not force_refresh:
        return cache_path.read_text(encoding="utf-8")

    cache_path.parent.mkdir(parents=True, exist_ok=True)

    last_err: Optional[BaseException] = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout_sec) as resp:
                data = resp.read().decode("utf-8")
            cache_path.write_text(data, encoding="utf-8")
            return data
        except urllib.error.URLError as e:
            reason = getattr(e, "reason", None)
            # Some local Python envs fail to resolve CA bundles. Fallback to
            # explicit insecure TLS context so pipeline can still proceed.
            if isinstance(reason, ssl.SSLCertVerificationError):
                try:
                    insecure_ctx = ssl._create_unverified_context()
                    with urllib.request.urlopen(url, timeout=timeout_sec, context=insecure_ctx) as resp:
                        data = resp.read().decode("utf-8")
                    cache_path.write_text(data, encoding="utf-8")
                    return data
                except Exception as insecure_err:
                    last_err = insecure_err
                if attempt < retries:
                    time.sleep(min(10, 1.5 * attempt))
                continue

            last_err = e
            if attempt < retries:
                time.sleep(min(10, 1.5 * attempt))
        except (TimeoutError, ConnectionError, OSError) as e:
            last_err = e
            if attempt < retries:
                time.sleep(min(10, 1.5 * attempt))

    raise RuntimeError(
        "KEGG API fetch failed after retries. "
        "Manual download/authorization may be required. "
        f"URL={url}; cache_target={cache_path}; last_error={last_err!r}. "
        "Please manually place file at this cache_target and verify with: "
        f"sha256sum {cache_path}"
    )


def _parse_pathway_source_version(info_text: str) -> str:
    for line in info_text.splitlines():
        m = re.search(r"Release\s+(.+)$", line)
        if m:
            return f"KEGG pathway {m.group(1).strip()}"
    return "KEGG pathway (version_unknown)"


def _parse_ncbi_to_kegg_gene(conv_text: str) -> Dict[str, Set[str]]:
    out: Dict[str, Set[str]] = defaultdict(set)
    for raw in conv_text.splitlines():
        if not raw.strip():
            continue
        parts = raw.split("\t")
        if len(parts) != 2:
            continue
        left, right = parts
        if left.startswith("ncbi-geneid:") and right.startswith("hsa:"):
            ncbi_id = left.split(":", 1)[1].strip()
            kegg_gene = right.strip()
            if ncbi_id and kegg_gene:
                out[ncbi_id].add(kegg_gene)
    return out


def _parse_kegg_gene_to_pathways(link_text: str) -> Dict[str, Set[str]]:
    out: Dict[str, Set[str]] = defaultdict(set)
    for raw in link_text.splitlines():
        if not raw.strip():
            continue
        parts = raw.split("\t")
        if len(parts) != 2:
            continue
        left, right = parts
        if left.startswith("hsa:") and right.startswith("path:hsa"):
            pathway_id = right.split(":", 1)[1].strip()
            if pathway_id:
                out[left.strip()].add(pathway_id)
    return out


def _parse_pathway_names(list_text: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for raw in list_text.splitlines():
        if not raw.strip():
            continue
        parts = raw.split("\t", 1)
        if len(parts) != 2:
            continue
        pid, pname = parts
        pid = pid.strip()
        pname = pname.strip()
        if pid:
            out[pid] = pname
    return out


def _normalize_ncbi_gene_id(v: str) -> Optional[str]:
    x = (v or "").strip()
    if x in {"", "NA", "N/A", "None", "null"}:
        return None

    # Most rows use forms like "12345" or "12345.0"
    m = re.fullmatch(r"(\d+)(?:\.0+)?", x)
    if m:
        return m.group(1)

    # Fallback: if a decimal-like representation exists, truncate to integer token.
    m2 = re.fullmatch(r"(\d+)\.\d+", x)
    if m2:
        return m2.group(1)

    return None


def _load_master_ncbi_index(master_path: Path, max_rows: Optional[int]) -> Tuple[Dict[str, Set[str]], BuildStats]:
    stats = BuildStats()
    ncbi_to_uniprots: Dict[str, Set[str]] = defaultdict(set)

    with master_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"uniprot_id", "ncbi_gene_id"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] protein master missing columns: {sorted(missing)}")

        for i, row in enumerate(reader, start=1):
            if max_rows is not None and i > max_rows:
                break

            stats.master_rows += 1
            uniprot = (row.get("uniprot_id") or "").strip()
            if uniprot:
                stats.master_uniprot_non_empty += 1

            ncbi_norm = _normalize_ncbi_gene_id(row.get("ncbi_gene_id", ""))
            if uniprot and ncbi_norm:
                ncbi_to_uniprots[ncbi_norm].add(uniprot)
                stats.master_rows_with_parseable_ncbi += 1

    stats.distinct_ncbi_ids = len(ncbi_to_uniprots)
    return ncbi_to_uniprots, stats


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-master", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--cache-dir", type=Path, required=True)
    ap.add_argument("--max-rows", type=int, default=None)
    ap.add_argument("--timeout-sec", type=int, default=60)
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--force-refresh", action="store_true")
    args = ap.parse_args()

    if not args.input_master.exists():
        raise SystemExit(f"[ERROR] missing input master table: {args.input_master}")

    mode = "sample" if args.max_rows is not None else "full"

    info_text = _fetch_text_with_cache(
        KEGG_INFO_URL,
        args.cache_dir / "kegg_info_pathway.txt",
        timeout_sec=args.timeout_sec,
        retries=args.retries,
        force_refresh=args.force_refresh,
    )
    conv_text = _fetch_text_with_cache(
        KEGG_CONV_URL,
        args.cache_dir / "kegg_conv_hsa_ncbi-geneid.tsv",
        timeout_sec=args.timeout_sec,
        retries=args.retries,
        force_refresh=args.force_refresh,
    )
    link_text = _fetch_text_with_cache(
        KEGG_LINK_URL,
        args.cache_dir / "kegg_link_pathway_hsa.tsv",
        timeout_sec=args.timeout_sec,
        retries=args.retries,
        force_refresh=args.force_refresh,
    )
    list_text = _fetch_text_with_cache(
        KEGG_LIST_URL,
        args.cache_dir / "kegg_list_pathway_hsa.tsv",
        timeout_sec=args.timeout_sec,
        retries=args.retries,
        force_refresh=args.force_refresh,
    )

    source_version = _parse_pathway_source_version(info_text)
    fetch_date = utc_today_date()

    ncbi_to_kegg_genes = _parse_ncbi_to_kegg_gene(conv_text)
    kegg_gene_to_pathways = _parse_kegg_gene_to_pathways(link_text)
    pathway_names = _parse_pathway_names(list_text)

    ncbi_to_uniprots, stats = _load_master_ncbi_index(args.input_master, args.max_rows)

    rows_set: Set[Tuple[str, str, str, str, str, str]] = set()
    mapped_uniprots: Set[str] = set()
    mapped_ncbi_ids = 0
    missing_ncbi_ids = 0
    kegg_genes_with_pathway: Set[str] = set()
    fallback_name_count = 0

    for ncbi_id, uniprots in ncbi_to_uniprots.items():
        kegg_genes = ncbi_to_kegg_genes.get(ncbi_id)
        if not kegg_genes:
            missing_ncbi_ids += 1
            continue

        mapped_ncbi_ids += 1
        for kegg_gene_id in sorted(kegg_genes):
            pathways = kegg_gene_to_pathways.get(kegg_gene_id)
            if not pathways:
                continue

            kegg_genes_with_pathway.add(kegg_gene_id)
            for pathway_id in sorted(pathways):
                pathway_name = pathway_names.get(pathway_id)
                if not pathway_name:
                    pathway_name = pathway_id
                    fallback_name_count += 1

                for uniprot_id in sorted(uniprots):
                    rows_set.add(
                        (
                            uniprot_id,
                            kegg_gene_id,
                            pathway_id,
                            pathway_name,
                            source_version,
                            fetch_date,
                        )
                    )
                    mapped_uniprots.add(uniprot_id)

    rows_sorted = sorted(rows_set, key=lambda x: (x[0], x[2], x[1]))

    header = [
        "uniprot_id",
        "kegg_gene_id",
        "kegg_pathway_id",
        "pathway_name",
        "source_version",
        "fetch_date",
    ]

    n_rows = _atomic_write_tsv(args.out, header, rows_sorted)

    stats.ncbi_ids_mapped_to_kegg_gene = mapped_ncbi_ids
    stats.ncbi_ids_without_kegg_gene = missing_ncbi_ids
    stats.kegg_genes_with_pathway = len(kegg_genes_with_pathway)
    stats.mapped_uniprot_count = len(mapped_uniprots)
    stats.output_rows = n_rows
    stats.output_unique_pathways = len({row[2] for row in rows_sorted})
    stats.output_unique_kegg_genes = len({row[1] for row in rows_sorted})
    stats.pathway_name_fallback_count = fallback_name_count

    mapped_rate = (
        (stats.mapped_uniprot_count / stats.master_uniprot_non_empty)
        if stats.master_uniprot_non_empty
        else 0.0
    )

    report = {
        "name": "protein_kegg_pathway_v1.build",
        "created_at": utc_now_iso(),
        "mode": mode,
        "inputs": {
            "master": str(args.input_master),
            "max_rows": args.max_rows,
        },
        "kegg_resources": {
            "info_url": KEGG_INFO_URL,
            "conv_url": KEGG_CONV_URL,
            "link_url": KEGG_LINK_URL,
            "list_url": KEGG_LIST_URL,
            "cache_dir": str(args.cache_dir),
            "source_version": source_version,
            "fetch_date": fetch_date,
        },
        "metrics": {
            "master_rows": stats.master_rows,
            "master_uniprot_non_empty": stats.master_uniprot_non_empty,
            "master_rows_with_parseable_ncbi": stats.master_rows_with_parseable_ncbi,
            "distinct_ncbi_ids": stats.distinct_ncbi_ids,
            "ncbi_ids_mapped_to_kegg_gene": stats.ncbi_ids_mapped_to_kegg_gene,
            "ncbi_ids_without_kegg_gene": stats.ncbi_ids_without_kegg_gene,
            "kegg_genes_with_pathway": stats.kegg_genes_with_pathway,
            "mapped_uniprot_count": stats.mapped_uniprot_count,
            "mapped_uniprot_rate": mapped_rate,
            "output_rows": stats.output_rows,
            "output_unique_pathways": stats.output_unique_pathways,
            "output_unique_kegg_genes": stats.output_unique_kegg_genes,
            "pathway_name_fallback_count": stats.pathway_name_fallback_count,
        },
        "notes": {
            "ncbi_gene_id_normalization": "accepts integer and *.0 forms; malformed values are ignored",
            "source": "KEGG REST API",
        },
    }

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"[OK] mode={mode} out={args.out} rows={n_rows}")
    print(f"[OK] report={args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
