#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import time
import warnings
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

import requests
from requests.exceptions import RequestException, RequestsDependencyWarning

warnings.filterwarnings("ignore", category=RequestsDependencyWarning)

INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
ZINC_ID_RE = re.compile(r"^ZINC\d+$")

DEFAULT_BASE_URL = "http://files.docking.org/catalogs/"
DEFAULT_TIERS = ["10", "20", "30", "40", "50"]
USER_AGENT = "Mozilla/5.0 (compatible; protian-entity-zinc-direct-xref/1.0)"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def utc_today() -> str:
    return datetime.now(timezone.utc).date().isoformat()


def normalize(x: str) -> str:
    return (x or "").strip()


def valid_inchikey(x: str) -> bool:
    return bool(INCHIKEY_RE.match(normalize(x).upper()))


def valid_zinc_id(x: str) -> bool:
    return bool(ZINC_ID_RE.match(normalize(x).upper()))


def join_values(values: Iterable[str]) -> str:
    return ";".join(sorted({normalize(v) for v in values if normalize(v)}))


def write_tsv(path: Path, rows: Iterable[Dict[str, str]], header: Sequence[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    n = 0
    with tmp.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(header), lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in header})
            n += 1
    tmp.replace(path)
    return n


def read_tsv(path: Path, max_rows: Optional[int] = None) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"[ERROR] missing header: {path}")
        rows: List[Dict[str, str]] = []
        for i, row in enumerate(r, start=1):
            rows.append(row)
            if max_rows is not None and i >= max_rows:
                break
    return list(r.fieldnames), rows


def load_xref_inchikey(path: Path, max_rows: Optional[int]) -> Tuple[Set[str], Dict[str, Any]]:
    cols, rows = read_tsv(path, max_rows=max_rows)
    if "inchikey" not in cols:
        raise SystemExit(f"[ERROR] xref missing inchikey: {path}")

    out: Set[str] = set()
    bad = 0
    for r in rows:
        ik = normalize(r.get("inchikey", "")).upper()
        if not valid_inchikey(ik):
            bad += 1
            continue
        out.add(ik)

    return out, {
        "input_rows": len(rows),
        "valid_inchikey_rows": len(out),
        "skipped_bad_inchikey_rows": bad,
    }


def session_get_text(session: requests.Session, url: str, timeout: int, retries: int) -> str:
    last_err: Optional[str] = None
    for attempt in range(retries + 1):
        try:
            resp = session.get(url, timeout=timeout)
            if resp.status_code == 200:
                return resp.text
            last_err = f"HTTP {resp.status_code}"
            if resp.status_code in {403, 404}:
                break
        except Exception as e:
            last_err = f"{type(e).__name__}: {e}"
        time.sleep(min(2.5, 0.4 * (attempt + 1)))
    raise RuntimeError(last_err or "unknown_http_error")


def parse_hrefs(index_html: str) -> List[str]:
    return re.findall(r'href="([^"]+)"', index_html)


def list_vendor_dirs(session: requests.Session, base_url: str, tier: str, timeout: int, retries: int) -> List[str]:
    url = f"{base_url}{tier}/"
    html = session_get_text(session, url, timeout=timeout, retries=retries)
    out: List[str] = []
    seen: Set[str] = set()
    for href in parse_hrefs(html):
        if href.startswith("?"):
            continue
        if href.startswith("/"):
            continue
        if not href.endswith("/"):
            continue
        name = href.strip("/").strip()
        if not name:
            continue
        if name in {".", "..", "list", "tranches", "source", "errors", "filtered", "unique"}:
            continue
        if name not in seen:
            seen.add(name)
            out.append(name)
    return sorted(out)


def iter_response_lines(resp: requests.Response, gzipped: bool) -> Iterator[str]:
    if gzipped:
        gz = gzip.GzipFile(fileobj=resp.raw)
        for raw in gz:
            if not raw:
                continue
            yield raw.decode("utf-8", errors="replace").rstrip("\n")
    else:
        for raw in resp.iter_lines(decode_unicode=False):
            if raw is None:
                continue
            if isinstance(raw, bytes):
                yield raw.decode("utf-8", errors="replace").rstrip("\n")
            else:
                yield str(raw).rstrip("\n")


def fetch_info_matches(
    session: requests.Session,
    info_url: str,
    *,
    tier: str,
    vendor: str,
    xref_inchikeys: Set[str],
    timeout: int,
    retries: int,
) -> Tuple[List[Dict[str, str]], Dict[str, Any], str]:
    """Return (matched_rows, metrics, status). status in {ok, not_found, failed}."""
    gzipped = info_url.endswith(".gz")
    metrics = Counter()
    last_err = ""

    for attempt in range(retries + 1):
        matched: List[Dict[str, str]] = []
        try:
            with session.get(info_url, stream=True, timeout=timeout) as resp:
                if resp.status_code == 404:
                    return [], {"status_code": 404}, "not_found"
                if resp.status_code != 200:
                    last_err = f"HTTP {resp.status_code}"
                    raise RuntimeError(last_err)

                for line in iter_response_lines(resp, gzipped=gzipped):
                    metrics["lines_total"] += 1
                    s = line.strip()
                    if not s or s.startswith("#"):
                        continue

                    parts = s.split("\t")
                    if len(parts) < 3:
                        parts = s.split()
                    if len(parts) < 3:
                        metrics["lines_bad_short"] += 1
                        continue

                    supplier_code = normalize(parts[0])
                    zinc_id = normalize(parts[1]).upper()
                    inchikey = normalize(parts[2]).upper()
                    tranche = normalize(parts[3]) if len(parts) >= 4 else ""
                    availability_tag = normalize(" ".join(parts[4:])) if len(parts) >= 5 else ""

                    if not valid_inchikey(inchikey):
                        metrics["lines_bad_inchikey"] += 1
                        continue
                    if not valid_zinc_id(zinc_id):
                        metrics["lines_bad_zinc_id"] += 1
                        continue

                    if inchikey not in xref_inchikeys:
                        metrics["lines_not_in_xref"] += 1
                        continue

                    matched.append(
                        {
                            "inchikey": inchikey,
                            "zinc_id": zinc_id,
                            "supplier_code": supplier_code,
                            "supplier_name": vendor,
                            "catalog_tier": tier,
                            "zinc_tranche": tranche,
                            "availability_tag": availability_tag,
                            "_info_url": info_url,
                        }
                    )
                    metrics["matched_rows"] += 1

            return matched, {k: int(v) for k, v in metrics.items()}, "ok"

        except Exception as e:
            last_err = f"{type(e).__name__}: {e}"
            metrics["fetch_errors"] += 1
            if attempt < retries:
                time.sleep(min(3.0, 0.5 * (attempt + 1)))

    return [], {**{k: int(v) for k, v in metrics.items()}, "error": last_err}, "failed"


def infer_tier_label(tier: str) -> str:
    labels = {
        "10": "tier10",
        "20": "tier20",
        "30": "tier30",
        "40": "tier40",
        "50": "tier50",
    }
    return labels.get(tier, f"tier{tier}")


def infer_purchase_label(tier: str, availability_tag: str) -> str:
    tag = normalize(availability_tag).lower()
    if any(k in tag for k in ["in-stock", "for sale", "forsale", "purchas"]):
        return "purchasable"
    if tag.startswith("in-"):
        return tag.split()[0]
    if tier in {"40", "50"}:
        return "likely_purchasable"
    return "unknown"


def fetch_source_smiles(
    session: requests.Session,
    base_url: str,
    vendor: str,
    needed_codes: Set[str],
    timeout: int,
    retries: int,
) -> Tuple[Dict[str, str], Dict[str, Any]]:
    if not needed_codes:
        return {}, {"status": "skip_empty_needed_codes"}

    candidates = [
        (f"{base_url}source/{vendor}.src.txt", False),
        (f"{base_url}source/{vendor}.src.txt.gz", True),
    ]

    metrics = Counter()
    metrics["needed_codes"] = len(needed_codes)
    last_err = ""

    for src_url, force_gz in candidates:
        gzipped = force_gz or src_url.endswith(".gz")
        for attempt in range(retries + 1):
            found: Dict[str, str] = {}
            try:
                with session.get(src_url, stream=True, timeout=timeout) as resp:
                    if resp.status_code == 404:
                        metrics["source_404"] += 1
                        break
                    if resp.status_code != 200:
                        last_err = f"HTTP {resp.status_code}"
                        raise RuntimeError(last_err)

                    for line in iter_response_lines(resp, gzipped=gzipped):
                        metrics["source_lines_total"] += 1
                        s = line.strip()
                        if not s or s.startswith("#"):
                            continue
                        parts = s.split()
                        if len(parts) < 2:
                            metrics["source_lines_bad_short"] += 1
                            continue

                        smiles = normalize(parts[0])
                        code = normalize(parts[1])
                        if not smiles or not code:
                            continue
                        if code in needed_codes and code not in found:
                            found[code] = smiles
                            metrics["source_codes_found"] += 1
                            if len(found) >= len(needed_codes):
                                break

                metrics["source_url_used"] += 1
                return found, {
                    **{k: int(v) for k, v in metrics.items()},
                    "status": "ok",
                    "source_url": src_url,
                    "source_error": "",
                }
            except Exception as e:
                last_err = f"{type(e).__name__}: {e}"
                metrics["source_fetch_errors"] += 1
                if attempt < retries:
                    time.sleep(min(2.5, 0.4 * (attempt + 1)))

    return {}, {
        **{k: int(v) for k, v in metrics.items()},
        "status": "failed",
        "source_url": "",
        "source_error": last_err,
    }


def non_empty_rate(rows: List[Dict[str, str]], col: str) -> Dict[str, Any]:
    n = sum(1 for r in rows if normalize(r.get(col, "")) != "")
    total = len(rows)
    return {"non_empty": n, "total": total, "rate": (n / total) if total else 0.0}


def maybe_load_baseline_rate(path: Optional[Path], fallback: float = 0.0183) -> Tuple[float, str]:
    if path is None or not path.exists():
        return fallback, "fallback_constant"
    try:
        d = json.loads(path.read_text(encoding="utf-8"))
        rate = float(d.get("row_metrics", {}).get("zinc_3d_rate", fallback))
        return rate, f"from:{path}"
    except Exception:
        return fallback, "fallback_constant_parse_error"


def write_manual_download_report(path: Path, base_url: str) -> None:
    report = {
        "name": "molecule_zinc_xref_v1.manual_download",
        "reason": "failed to auto fetch ZINC catalog info mapping files",
        "checklist": [
            {
                "name": "ZINC catalogs root",
                "url": base_url,
                "hint": "Use browser or curl with User-Agent to avoid 403",
            },
            {
                "name": "Tier examples",
                "urls": [
                    f"{base_url}40/",
                    f"{base_url}50/",
                    f"{base_url}source/",
                ],
            },
            {
                "name": "Example manual download command",
                "cmd": "curl -A 'Mozilla/5.0' -O http://files.docking.org/catalogs/50/enamine/enamine.info.txt.gz",
            },
        ],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xref", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--coverage-report", type=Path, required=True)
    ap.add_argument("--conflict-report", type=Path, required=True)
    ap.add_argument("--manual-report", type=Path, required=True)
    ap.add_argument("--zinc-base-url", default=DEFAULT_BASE_URL)
    ap.add_argument("--tiers", nargs="*", default=DEFAULT_TIERS)
    ap.add_argument("--max-rows", type=int, default=None)
    ap.add_argument("--max-vendors", type=int, default=None)
    ap.add_argument("--http-timeout", type=int, default=35)
    ap.add_argument("--http-retries", type=int, default=2)
    ap.add_argument("--baseline-coverage-report", type=Path, default=Path("pipelines/molecule_3d_registry/reports/molecule_3d_registry_v1.coverage.json"))
    args = ap.parse_args()

    if not args.xref.exists():
        raise SystemExit(f"[ERROR] missing xref input: {args.xref}")

    created_at = utc_now()
    fetch_date = utc_today()

    xref_inchikeys, xref_stats = load_xref_inchikey(args.xref, args.max_rows)

    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT, "Accept": "text/html,application/octet-stream,*/*"})

    scan_stats = Counter()
    per_vendor_scan: Dict[str, Any] = {}

    all_matched: List[Dict[str, str]] = []
    discovered_vendor_count = 0

    tiers = [normalize(t).strip("/") for t in args.tiers if normalize(t)]

    for tier in tiers:
        try:
            vendors = list_vendor_dirs(session, args.zinc_base_url, tier=tier, timeout=args.http_timeout, retries=args.http_retries)
        except Exception as e:
            per_vendor_scan[f"tier_{tier}"] = {"status": "failed", "error": f"{type(e).__name__}: {e}"}
            continue

        if args.max_vendors is not None and args.max_vendors >= 0:
            vendors = vendors[: args.max_vendors]

        discovered_vendor_count += len(vendors)
        per_vendor_scan[f"tier_{tier}"] = {"status": "ok", "vendors": len(vendors)}

        for vendor in vendors:
            scan_stats["vendor_attempted"] += 1
            info_candidates = [
                f"{args.zinc_base_url}{tier}/{vendor}/{vendor}.info.txt.gz",
                f"{args.zinc_base_url}{tier}/{vendor}/{vendor}.info.txt",
            ]

            matched_rows: List[Dict[str, str]] = []
            used_url = ""
            status = "not_found"
            metrics: Dict[str, Any] = {}

            for info_url in info_candidates:
                got, met, st = fetch_info_matches(
                    session,
                    info_url,
                    tier=tier,
                    vendor=vendor,
                    xref_inchikeys=xref_inchikeys,
                    timeout=args.http_timeout,
                    retries=args.http_retries,
                )
                if st == "ok":
                    matched_rows = got
                    metrics = met
                    used_url = info_url
                    status = "ok"
                    break
                if st == "failed":
                    metrics = met
                    status = "failed"
                else:
                    metrics = met

            if status == "ok":
                scan_stats["vendor_info_ok"] += 1
                if matched_rows:
                    scan_stats["vendor_with_match"] += 1
                    all_matched.extend(matched_rows)
                    scan_stats["matched_rows_total"] += len(matched_rows)
            elif status == "failed":
                scan_stats["vendor_info_failed"] += 1
            else:
                scan_stats["vendor_info_not_found"] += 1

            per_vendor_scan[f"{tier}/{vendor}"] = {
                "status": status,
                "info_url": used_url,
                "metrics": metrics,
            }

    if discovered_vendor_count == 0 or int(scan_stats.get("vendor_info_ok", 0)) == 0:
        write_manual_download_report(args.manual_report, base_url=args.zinc_base_url)
        print(f"[ERROR] unable to fetch any ZINC info files; manual checklist -> {args.manual_report}")
        return 2

    # Dedupe + collect needed supplier codes for smiles mapping
    dedup_key: Set[Tuple[str, str, str, str, str]] = set()
    dedup_rows: List[Dict[str, str]] = []
    needed_codes_by_vendor: Dict[str, Set[str]] = defaultdict(set)

    for r in all_matched:
        key = (
            r["inchikey"],
            r["zinc_id"],
            r["supplier_name"],
            r["supplier_code"],
            r["catalog_tier"],
        )
        if key in dedup_key:
            scan_stats["dedup_rows_dropped"] += 1
            continue
        dedup_key.add(key)
        dedup_rows.append(r)
        code = normalize(r.get("supplier_code", ""))
        if code:
            needed_codes_by_vendor[r["supplier_name"]].add(code)

    # Load smiles from source/<vendor>.src.txt where possible
    smiles_by_vendor_code: Dict[Tuple[str, str], str] = {}
    source_fetch_report: Dict[str, Any] = {}

    for vendor, codes in sorted(needed_codes_by_vendor.items()):
        smiles_map, rep = fetch_source_smiles(
            session,
            base_url=args.zinc_base_url,
            vendor=vendor,
            needed_codes=codes,
            timeout=args.http_timeout,
            retries=args.http_retries,
        )
        source_fetch_report[vendor] = rep
        for code, smiles in smiles_map.items():
            smiles_by_vendor_code[(vendor, code)] = smiles

    # Build output rows
    out_rows: List[Dict[str, str]] = []
    xref_version = "v2" if "_v2" in args.xref.name else "v1"

    for r in dedup_rows:
        vendor = r["supplier_name"]
        code = normalize(r.get("supplier_code", ""))
        smiles = normalize(smiles_by_vendor_code.get((vendor, code), ""))

        tier = normalize(r.get("catalog_tier", ""))
        tier_label = infer_tier_label(tier)
        availability_tag = normalize(r.get("availability_tag", ""))
        purchase_label = infer_purchase_label(tier, availability_tag)

        source_url = normalize(r.get("_info_url", ""))
        source_version = f"ZINC:catalogs/{tier}/{vendor}/{vendor}.info;ZINC:catalogs/source/{vendor}.src"

        out_rows.append(
            {
                "inchikey": r["inchikey"],
                "zinc_id": r["zinc_id"],
                "supplier_code": code,
                "supplier_name": vendor,
                "smiles": smiles,
                "catalog_tier": tier,
                "catalog_tier_label": tier_label,
                "zinc_tranche": normalize(r.get("zinc_tranche", "")),
                "availability_tag": availability_tag,
                "purchase_label": purchase_label,
                "mapping_method": "zinc_catalog_info_inchikey_exact",
                "xref_version": xref_version,
                "source": source_url or "http://files.docking.org/catalogs/",
                "source_version": source_version,
                "fetch_date": fetch_date,
            }
        )

    header = [
        "inchikey",
        "zinc_id",
        "supplier_code",
        "supplier_name",
        "smiles",
        "catalog_tier",
        "catalog_tier_label",
        "zinc_tranche",
        "availability_tag",
        "purchase_label",
        "mapping_method",
        "xref_version",
        "source",
        "source_version",
        "fetch_date",
    ]

    rows_written = write_tsv(args.out, out_rows, header)

    # coverage + conflict
    mapped_inchikey = {r["inchikey"] for r in out_rows}
    coverage_rate = (len(mapped_inchikey) / len(xref_inchikeys)) if xref_inchikeys else 0.0
    baseline_rate, baseline_source = maybe_load_baseline_rate(args.baseline_coverage_report)

    inchikey_to_zinc: Dict[str, Set[str]] = defaultdict(set)
    inchikey_to_vendor: Dict[str, Set[str]] = defaultdict(set)
    zinc_to_inchikey: Dict[str, Set[str]] = defaultdict(set)
    for r in out_rows:
        ik = r["inchikey"]
        zid = r["zinc_id"]
        v = r["supplier_name"]
        inchikey_to_zinc[ik].add(zid)
        inchikey_to_vendor[ik].add(v)
        zinc_to_inchikey[zid].add(ik)

    multi_zinc_inchikey = {ik: len(v) for ik, v in inchikey_to_zinc.items() if len(v) > 1}
    multi_vendor_inchikey = {ik: len(v) for ik, v in inchikey_to_vendor.items() if len(v) > 1}
    multi_inchikey_zinc = {z: len(v) for z, v in zinc_to_inchikey.items() if len(v) > 1}

    conflict_report = {
        "name": "molecule_zinc_xref_v1.conflict_audit",
        "created_at": created_at,
        "metrics": {
            "rows": rows_written,
            "inchikey_with_multi_zinc": len(multi_zinc_inchikey),
            "inchikey_with_multi_vendor": len(multi_vendor_inchikey),
            "zinc_with_multi_inchikey": len(multi_inchikey_zinc),
        },
        "samples": {
            "inchikey_with_multi_zinc": [
                {
                    "inchikey": ik,
                    "zinc_count": len(inchikey_to_zinc[ik]),
                    "zinc_ids": sorted(inchikey_to_zinc[ik])[:10],
                    "vendors": sorted(inchikey_to_vendor[ik])[:10],
                }
                for ik in sorted(multi_zinc_inchikey.keys())[:100]
            ],
            "zinc_with_multi_inchikey": [
                {
                    "zinc_id": zid,
                    "inchikey_count": len(zinc_to_inchikey[zid]),
                    "inchikey_sample": sorted(zinc_to_inchikey[zid])[:10],
                }
                for zid in sorted(multi_inchikey_zinc.keys())[:100]
            ],
        },
    }

    args.conflict_report.parent.mkdir(parents=True, exist_ok=True)
    args.conflict_report.write_text(json.dumps(conflict_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    coverage_report = {
        "name": "molecule_zinc_xref_v1.coverage",
        "created_at": created_at,
        "xref_rows": len(xref_inchikeys),
        "rows": rows_written,
        "mapped_inchikey": len(mapped_inchikey),
        "coverage_rate": coverage_rate,
        "baseline_rate": baseline_rate,
        "baseline_source": baseline_source,
        "coverage_delta_vs_baseline": coverage_rate - baseline_rate,
        "coverage_fold_vs_baseline": (coverage_rate / baseline_rate) if baseline_rate > 0 else None,
        "non_empty_rates": {
            "source": non_empty_rate(out_rows, "source"),
            "source_version": non_empty_rate(out_rows, "source_version"),
            "fetch_date": non_empty_rate(out_rows, "fetch_date"),
            "smiles": non_empty_rate(out_rows, "smiles"),
        },
        "catalog_tier_distribution": dict(Counter(r["catalog_tier"] for r in out_rows)),
        "purchase_label_distribution": dict(Counter(r["purchase_label"] for r in out_rows)),
    }
    args.coverage_report.parent.mkdir(parents=True, exist_ok=True)
    args.coverage_report.write_text(json.dumps(coverage_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    build_report = {
        "name": "molecule_zinc_xref_v1.build",
        "created_at": created_at,
        "sample_mode": args.max_rows is not None or args.max_vendors is not None,
        "max_rows": args.max_rows,
        "max_vendors": args.max_vendors,
        "inputs": {
            "xref": str(args.xref),
            "zinc_base_url": args.zinc_base_url,
            "tiers": tiers,
        },
        "output": str(args.out),
        "metrics": {
            "rows_written": rows_written,
            **xref_stats,
            **{k: int(v) for k, v in scan_stats.items()},
            "discovered_vendor_count": discovered_vendor_count,
            "smiles_joined_rows": int(sum(1 for r in out_rows if normalize(r.get("smiles", "")) != "")),
            "coverage_rate": coverage_rate,
            "mapped_inchikey": len(mapped_inchikey),
            "baseline_rate": baseline_rate,
            "coverage_delta_vs_baseline": coverage_rate - baseline_rate,
        },
        "vendor_scan": per_vendor_scan,
        "source_fetch": source_fetch_report,
        "coverage_report": str(args.coverage_report),
        "conflict_report": str(args.conflict_report),
    }

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(build_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"[OK] build -> {args.out} (rows={rows_written})")
    print(
        "[OK] coverage mapped_inchikey={} / {} = {:.4f} (baseline {:.4f})".format(
            len(mapped_inchikey),
            len(xref_inchikeys),
            coverage_rate,
            baseline_rate,
        )
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
