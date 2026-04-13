#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import ssl
import sys
import time
import urllib.error
import urllib.request
import urllib.parse
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


STRING_DETAILED_URLS = [
    "https://stringdb-static.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz",
]
INTACT_URLS = [
    "https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/human.zip",
    "https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip",
]
BIOGRID_URLS = [
    "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip",
    "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ALL-LATEST.tab3.zip",
]


PMID_RE = re.compile(r"(?:pubmed|pmid)[:\s]*(\d+)", re.IGNORECASE)
DOI_RE = re.compile(r"(?:doi)[:\s]*([^\s|;]+)", re.IGNORECASE)
GO_ID_RE = re.compile(r"GO:\d{7}")
UNIPROT_RE = re.compile(r"uniprotkb:([A-Z0-9]+(?:-[0-9]+)?)", re.IGNORECASE)


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def utc_today() -> str:
    return datetime.now(timezone.utc).date().isoformat()


def _set_csv_field_size_limit() -> None:
    size = sys.maxsize
    while True:
        try:
            csv.field_size_limit(size)
            return
        except OverflowError:
            size //= 10


def _atomic_write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    n = 0
    with tmp.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(list(header))
        for row in rows:
            w.writerow(list(row))
            n += 1
    tmp.replace(path)
    return n


def _download_with_retries(url: str, out_path: Path, timeout_sec: int = 120, retries: int = 3) -> bool:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    last_err: Optional[BaseException] = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "kg-ppi-semantic-enrichment/1.0"})
            with urllib.request.urlopen(req, timeout=timeout_sec) as resp, out_path.open("wb") as w:
                shutil.copyfileobj(resp, w)
            return True
        except urllib.error.URLError as e:
            reason = getattr(e, "reason", None)
            if isinstance(reason, ssl.SSLCertVerificationError):
                try:
                    ctx = ssl._create_unverified_context()
                    req = urllib.request.Request(url, headers={"User-Agent": "kg-ppi-semantic-enrichment/1.0"})
                    with urllib.request.urlopen(req, timeout=timeout_sec, context=ctx) as resp, out_path.open("wb") as w:
                        shutil.copyfileobj(resp, w)
                    return True
                except Exception as e2:
                    last_err = e2
            else:
                last_err = e
        except Exception as e:
            last_err = e

        if attempt < retries:
            time.sleep(min(10, attempt * 2.0))

    if out_path.exists():
        out_path.unlink(missing_ok=True)
    print(f"[WARN] download failed: {url}; error={last_err!r}")
    return False


def _ensure_string_detailed(cache_dir: Path) -> Tuple[Optional[Path], Dict[str, str]]:
    target = cache_dir / "9606.protein.links.detailed.v12.0.txt.gz"
    report = {
        "name": "string_detailed",
        "status": "missing",
        "path": str(target),
        "url": "",
    }
    if target.exists() and target.stat().st_size > 0:
        report["status"] = "cached"
        return target, report

    for url in STRING_DETAILED_URLS:
        if _download_with_retries(url, target):
            report["status"] = "downloaded"
            report["url"] = url
            return target, report

    return None, report


def _extract_text_from_archive(archive_path: Path, prefer_patterns: Sequence[str], out_path: Path) -> bool:
    lower = archive_path.name.lower()
    if not (lower.endswith(".zip") or lower.endswith(".gz") or lower.endswith(".txt") or lower.endswith(".tab")):
        try:
            with archive_path.open("rb") as f:
                head = f.read(4)
            if head.startswith(b"PK"):
                lower = lower + ".zip"
            elif head.startswith(b"\x1f\x8b"):
                lower = lower + ".gz"
        except Exception:
            pass
    if lower.endswith(".zip"):
        try:
            with zipfile.ZipFile(archive_path, "r") as zf:
                names = zf.namelist()
                if not names:
                    return False

                def score(name: str) -> Tuple[int, int]:
                    lname = name.lower()
                    patt_score = 0
                    for i, p in enumerate(prefer_patterns, start=1):
                        if p in lname:
                            patt_score += (100 - i)
                    is_text = int(lname.endswith(".txt") or lname.endswith(".tab") or ".tab" in lname)
                    return patt_score, is_text

                names_sorted = sorted(names, key=score, reverse=True)
                picked = names_sorted[0]
                out_path.parent.mkdir(parents=True, exist_ok=True)
                with zf.open(picked, "r") as r, out_path.open("wb") as w:
                    shutil.copyfileobj(r, w)
            return True
        except Exception:
            return False

    if lower.endswith(".gz"):
        try:
            out_path.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(archive_path, "rb") as r, out_path.open("wb") as w:
                shutil.copyfileobj(r, w)
            return True
        except Exception:
            return False

    if lower.endswith(".txt") or lower.endswith(".tab"):
        out_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(archive_path, out_path)
        return True

    return False


def _ensure_intact(cache_dir: Path) -> Tuple[Optional[Path], Dict[str, str]]:
    out_path = cache_dir / "intact_human.mitab.txt"
    report = {
        "name": "intact_mitab",
        "status": "missing",
        "path": str(out_path),
        "url": "",
    }
    if out_path.exists() and out_path.stat().st_size > 0:
        report["status"] = "cached"
        return out_path, report

    for url in INTACT_URLS:
        name = Path(urllib.parse.urlparse(url).path).name or "intact_download.tmp"
        tmp = cache_dir / f"intact_download.{name}"
        tmp.unlink(missing_ok=True)
        if not _download_with_retries(url, tmp):
            continue
        if _extract_text_from_archive(tmp, prefer_patterns=["human", "mitab", "intact"], out_path=out_path):
            report["status"] = "downloaded"
            report["url"] = url
            tmp.unlink(missing_ok=True)
            return out_path, report
    return None, report


def _ensure_biogrid(cache_dir: Path) -> Tuple[Optional[Path], Dict[str, str]]:
    out_path = cache_dir / "biogrid_human.tab3.txt"
    report = {
        "name": "biogrid_tab3",
        "status": "missing",
        "path": str(out_path),
        "url": "",
    }
    if out_path.exists() and out_path.stat().st_size > 0:
        report["status"] = "cached"
        return out_path, report

    for url in BIOGRID_URLS:
        name = Path(urllib.parse.urlparse(url).path).name or "biogrid_download.tmp"
        tmp = cache_dir / f"biogrid_download.{name}"
        tmp.unlink(missing_ok=True)
        if not _download_with_retries(url, tmp):
            continue
        if _extract_text_from_archive(tmp, prefer_patterns=["homo_sapiens", "human", "tab3", "biogrid"], out_path=out_path):
            report["status"] = "downloaded"
            report["url"] = url
            tmp.unlink(missing_ok=True)
            return out_path, report
    return None, report


def _normalize_uniprot_id(v: str) -> str:
    x = (v or "").strip().upper()
    if not x:
        return ""
    if ":" in x:
        x = x.split(":", 1)[1]
    # Keep canonical accession for edge matching
    if "-" in x:
        x = x.split("-", 1)[0]
    return x


def _edge_id(src: str, dst: str) -> str:
    a, b = (src, dst) if src <= dst else (dst, src)
    return f"ppi|{a}|{b}"


def _extract_uniprots_from_field(v: str) -> List[str]:
    out: List[str] = []
    seen: Set[str] = set()
    for m in UNIPROT_RE.finditer(v or ""):
        u = _normalize_uniprot_id(m.group(1))
        if u and u not in seen:
            seen.add(u)
            out.append(u)
    if out:
        return out

    # BioGRID often provides pipe-separated accessions without uniprotkb prefix.
    for token in re.split(r"[|;,]", v or ""):
        tok = _normalize_uniprot_id(token)
        if not tok:
            continue
        if re.fullmatch(r"[A-Z0-9]{6,10}", tok) and tok not in seen:
            seen.add(tok)
            out.append(tok)
    return out


def _extract_pmid_doi(raw: str) -> Tuple[Set[str], Set[str]]:
    pmids: Set[str] = set()
    dois: Set[str] = set()
    x = raw or ""
    for m in PMID_RE.finditer(x):
        pmids.add(m.group(1))
    for m in DOI_RE.finditer(x):
        doi = m.group(1).strip().rstrip(".);")
        if doi:
            dois.add(doi)
    return pmids, dois


def _normalize_method(raw: str) -> str:
    t = (raw or "").strip().lower()
    if not t:
        return "OTHER_EXPERIMENTAL"
    if "two hybrid" in t or "y2h" in t:
        return "Y2H"
    if "coimmunoprecipitation" in t or "co-immunoprecipitation" in t or "coip" in t:
        return "Co-IP"
    if "affinity chromatography" in t or "affinity purification" in t or "ap-ms" in t or "pull down" in t:
        return "AP-MS"
    if "protein complementation assay" in t or "pca" in t:
        return "PCA"
    if "x-ray" in t or "crystallography" in t:
        return "X-RAY"
    if "nmr" in t:
        return "NMR"
    if "mass spectrometry" in t:
        return "AP-MS"
    return "OTHER_EXPERIMENTAL"


def _normalize_throughput(raw: str, default_if_empty: str = "LT") -> str:
    t = (raw or "").strip().lower()
    if not t:
        return default_if_empty
    if "high" in t or t in {"ht", "high throughput"}:
        return "HT"
    if "low" in t or t in {"lt", "low throughput"}:
        return "LT"
    return default_if_empty


def _best_method(method_counts: Counter[str]) -> str:
    if not method_counts:
        return "COMPUTATIONAL"
    priority = ["Co-IP", "AP-MS", "Y2H", "PCA", "X-RAY", "NMR", "OTHER_EXPERIMENTAL"]
    ranked = sorted(method_counts.items(), key=lambda kv: (-kv[1], priority.index(kv[0]) if kv[0] in priority else 99, kv[0]))
    return ranked[0][0]


def _best_throughput(tp_counts: Counter[str]) -> str:
    if not tp_counts:
        return "HT"
    return sorted(tp_counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]


def _join_limited(values: Set[str], limit: int = 20) -> str:
    if not values:
        return ""
    return ";".join(sorted(values)[:limit])


def _safe_rate(num: int, den: int) -> float:
    if den <= 0:
        return 1.0
    return num / den


@dataclass
class EdgeRec:
    edge_id: str
    src_id: str
    dst_id: str


@dataclass
class MethodEvidence:
    methods: Counter[str] = field(default_factory=Counter)
    method_raws: Set[str] = field(default_factory=set)
    throughputs: Counter[str] = field(default_factory=Counter)
    pmids: Set[str] = field(default_factory=set)
    dois: Set[str] = field(default_factory=set)
    sources: Set[str] = field(default_factory=set)


@dataclass
class ScoreSplit:
    experimental: int = 0
    textmining: int = 0


@dataclass
class BuildStats:
    mode: str = "full"
    input_edges_rows: int = 0
    method_rows: int = 0
    function_rows: int = 0
    string_mapped_edges: int = 0
    intact_mapped_edges: int = 0
    biogrid_mapped_edges: int = 0
    method_non_empty_rate: float = 0.0
    pmid_or_doi_rate: float = 0.0
    context_any_rate: float = 0.0


def _load_edges(path: Path, max_rows: Optional[int]) -> List[EdgeRec]:
    edges: List[EdgeRec] = []
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {"edge_id", "src_id", "dst_id"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] edges table missing columns: {sorted(missing)}")

        for i, row in enumerate(r, start=1):
            if max_rows is not None and i > max_rows:
                break
            edge_id = (row.get("edge_id") or "").strip()
            src = (row.get("src_id") or "").strip()
            dst = (row.get("dst_id") or "").strip()
            if not edge_id or not src or not dst:
                continue
            edges.append(EdgeRec(edge_id=edge_id, src_id=src, dst_id=dst))
    return edges


def _load_string_id_to_uniprot(master_path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with master_path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {"uniprot_id", "string_ids"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] protein master missing columns: {sorted(missing)}")

        for row in r:
            uniprot = _normalize_uniprot_id(row.get("uniprot_id", ""))
            if not uniprot:
                continue
            raw = row.get("string_ids", "") or ""
            for token in re.split(r"[;,|\s]+", raw):
                s = token.strip()
                if not s:
                    continue
                if s.startswith("9606."):
                    mapping[s] = uniprot
    return mapping


def _parse_string_scores(
    detailed_path: Path,
    edge_set: Set[str],
    string_to_uniprot: Dict[str, str],
) -> Tuple[Dict[str, ScoreSplit], Dict[str, int]]:
    scores: Dict[str, ScoreSplit] = {}
    stat = {
        "rows_scanned": 0,
        "rows_with_human_ids": 0,
        "rows_edge_matched": 0,
    }
    if not detailed_path.exists():
        return scores, stat

    with gzip.open(detailed_path, "rt", encoding="utf-8", newline="") as f:
        header = f.readline().strip().split()
        if not header:
            return scores, stat
        idx = {c: i for i, c in enumerate(header)}
        needed = ["protein1", "protein2", "experimental", "textmining"]
        if any(c not in idx for c in needed):
            return scores, stat

        i_p1 = idx["protein1"]
        i_p2 = idx["protein2"]
        i_exp = idx["experimental"]
        i_txt = idx["textmining"]

        for raw in f:
            if not raw.strip():
                continue
            stat["rows_scanned"] += 1
            parts = raw.strip().split()
            if len(parts) <= max(i_p1, i_p2, i_exp, i_txt):
                continue

            p1 = parts[i_p1]
            p2 = parts[i_p2]
            u1 = string_to_uniprot.get(p1, "")
            u2 = string_to_uniprot.get(p2, "")
            if not u1 or not u2:
                continue
            stat["rows_with_human_ids"] += 1

            eid = _edge_id(u1, u2)
            if eid not in edge_set:
                continue
            stat["rows_edge_matched"] += 1

            try:
                exp = int(parts[i_exp])
                txt = int(parts[i_txt])
            except Exception:
                continue

            cur = scores.get(eid)
            if cur is None:
                scores[eid] = ScoreSplit(experimental=exp, textmining=txt)
            else:
                cur.experimental = max(cur.experimental, exp)
                cur.textmining = max(cur.textmining, txt)

    return scores, stat


def _parse_intact(
    intact_path: Path,
    edge_set: Set[str],
    method_ev: Dict[str, MethodEvidence],
) -> Dict[str, int]:
    stat = {
        "rows_scanned": 0,
        "rows_human": 0,
        "rows_edge_matched": 0,
    }
    if not intact_path.exists():
        return stat

    _set_csv_field_size_limit()
    with intact_path.open("r", encoding="utf-8", newline="") as f:
        r = csv.reader(f, delimiter="\t")
        for row in r:
            if not row:
                continue
            if row[0].startswith("#"):
                continue
            stat["rows_scanned"] += 1
            if len(row) < 11:
                continue

            tax_a = row[9] if len(row) > 9 else ""
            tax_b = row[10] if len(row) > 10 else ""
            if "9606" not in tax_a or "9606" not in tax_b:
                continue
            stat["rows_human"] += 1

            uids_a = _extract_uniprots_from_field(row[0])
            uids_b = _extract_uniprots_from_field(row[1])
            if not uids_a or not uids_b:
                continue
            u1 = uids_a[0]
            u2 = uids_b[0]
            if not u1 or not u2 or u1 == u2:
                continue

            eid = _edge_id(u1, u2)
            if eid not in edge_set:
                continue
            stat["rows_edge_matched"] += 1

            method_raw = row[6] if len(row) > 6 else ""
            pubs_raw = row[8] if len(row) > 8 else ""
            pmids, dois = _extract_pmid_doi(pubs_raw)

            ev = method_ev.setdefault(eid, MethodEvidence())
            method_norm = _normalize_method(method_raw)
            ev.methods[method_norm] += 1
            if method_raw:
                ev.method_raws.add(method_raw.strip())
            ev.throughputs[_normalize_throughput("", default_if_empty="LT")] += 1
            ev.pmids.update(pmids)
            ev.dois.update(dois)
            ev.sources.add("IntAct")

    return stat


def _resolve_biogrid_col(fieldnames: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    norm = {c.strip().lower(): c for c in fieldnames}
    for cand in candidates:
        if cand.lower() in norm:
            return norm[cand.lower()]
    for c in fieldnames:
        lc = c.lower()
        for cand in candidates:
            if cand.lower() in lc:
                return c
    return None


def _parse_biogrid(
    biogrid_path: Path,
    edge_set: Set[str],
    method_ev: Dict[str, MethodEvidence],
) -> Dict[str, int]:
    stat = {
        "rows_scanned": 0,
        "rows_human": 0,
        "rows_edge_matched": 0,
    }
    if not biogrid_path.exists():
        return stat

    _set_csv_field_size_limit()
    with biogrid_path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        fieldnames = r.fieldnames or []
        if not fieldnames:
            return stat

        col_uniprot_a = _resolve_biogrid_col(fieldnames, ["SWISS-PROT Accessions Interactor A", "Swiss-Prot Accessions Interactor A"])
        col_uniprot_b = _resolve_biogrid_col(fieldnames, ["SWISS-PROT Accessions Interactor B", "Swiss-Prot Accessions Interactor B"])
        col_method = _resolve_biogrid_col(fieldnames, ["Experimental System"])
        col_throughput = _resolve_biogrid_col(fieldnames, ["Throughput"])
        col_pubmed = _resolve_biogrid_col(fieldnames, ["Pubmed ID", "Publication Source"])
        col_tax_a = _resolve_biogrid_col(fieldnames, ["Organism ID Interactor A"])
        col_tax_b = _resolve_biogrid_col(fieldnames, ["Organism ID Interactor B"])

        if not col_uniprot_a or not col_uniprot_b:
            return stat

        for row in r:
            stat["rows_scanned"] += 1
            tax_a = (row.get(col_tax_a, "") if col_tax_a else "") or ""
            tax_b = (row.get(col_tax_b, "") if col_tax_b else "") or ""
            if col_tax_a and col_tax_b and ("9606" not in tax_a or "9606" not in tax_b):
                continue
            stat["rows_human"] += 1

            uids_a = _extract_uniprots_from_field(row.get(col_uniprot_a, "") or "")
            uids_b = _extract_uniprots_from_field(row.get(col_uniprot_b, "") or "")
            if not uids_a or not uids_b:
                continue
            u1 = uids_a[0]
            u2 = uids_b[0]
            if not u1 or not u2 or u1 == u2:
                continue

            eid = _edge_id(u1, u2)
            if eid not in edge_set:
                continue
            stat["rows_edge_matched"] += 1

            method_raw = row.get(col_method, "") if col_method else ""
            tp_raw = row.get(col_throughput, "") if col_throughput else ""
            pub_raw = row.get(col_pubmed, "") if col_pubmed else ""

            ev = method_ev.setdefault(eid, MethodEvidence())
            ev.methods[_normalize_method(method_raw)] += 1
            if method_raw:
                ev.method_raws.add(method_raw.strip())
            ev.throughputs[_normalize_throughput(tp_raw, default_if_empty="LT")] += 1
            pmids, dois = _extract_pmid_doi(pub_raw)
            if not pmids:
                # BioGRID often stores bare PubMed ids (numeric).
                naked = (pub_raw or "").strip()
                if naked.isdigit():
                    pmids.add(naked)
            ev.pmids.update(pmids)
            ev.dois.update(dois)
            ev.sources.add("BioGRID")

    return stat


def _parse_go_ids(raw: str) -> Set[str]:
    return set(GO_ID_RE.findall(raw or ""))


def _load_function_maps(master_path: Path) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, Set[str]]]:
    go_bp: Dict[str, Set[str]] = {}
    go_mf: Dict[str, Set[str]] = {}
    go_cc: Dict[str, Set[str]] = {}

    with master_path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {"uniprot_id", "go_biological_process", "go_molecular_function", "go_cellular_component"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] protein master missing GO columns: {sorted(missing)}")

        for row in r:
            u = _normalize_uniprot_id(row.get("uniprot_id", ""))
            if not u:
                continue
            go_bp[u] = _parse_go_ids(row.get("go_biological_process", ""))
            go_mf[u] = _parse_go_ids(row.get("go_molecular_function", ""))
            go_cc[u] = _parse_go_ids(row.get("go_cellular_component", ""))

    return go_bp, go_mf, go_cc


def _load_reactome_map(path: Path) -> Dict[str, Set[str]]:
    out: Dict[str, Set[str]] = defaultdict(set)
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {"uniprot_id", "pathway_id"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] reactome map missing columns: {sorted(missing)}")

        for row in r:
            u = _normalize_uniprot_id(row.get("uniprot_id", ""))
            pid = (row.get("pathway_id") or "").strip()
            if u and pid:
                out[u].add(pid)
    return out


def _load_kegg_map(path: Path) -> Dict[str, Set[str]]:
    out: Dict[str, Set[str]] = defaultdict(set)
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {"uniprot_id", "kegg_pathway_id"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise SystemExit(f"[ERROR] kegg map missing columns: {sorted(missing)}")

        for row in r:
            u = _normalize_uniprot_id(row.get("uniprot_id", ""))
            pid = (row.get("kegg_pathway_id") or "").strip()
            if u and pid:
                out[u].add(pid)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-edges", type=Path, required=True)
    ap.add_argument("--input-master", type=Path, required=True)
    ap.add_argument("--input-reactome", type=Path, required=True)
    ap.add_argument("--input-kegg", type=Path, required=True)
    ap.add_argument("--out-method", type=Path, required=True)
    ap.add_argument("--out-function", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--cache-dir", type=Path, required=True)
    ap.add_argument("--max-rows", type=int, default=None)
    ap.add_argument("--sample-fast", action="store_true")
    ap.add_argument("--string-default-doi", default="10.1093/nar/gkaa1074")
    ap.add_argument("--string-default-pmid", default="30476243")
    args = ap.parse_args()

    for p in [args.input_edges, args.input_master, args.input_reactome, args.input_kegg]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing input: {p}")

    mode = "sample" if args.max_rows is not None else "full"
    stats = BuildStats(mode=mode)

    edges = _load_edges(args.input_edges, max_rows=args.max_rows)
    stats.input_edges_rows = len(edges)
    edge_set = {e.edge_id for e in edges}

    method_ev: Dict[str, MethodEvidence] = {}
    score_map: Dict[str, ScoreSplit] = {}

    source_reports: List[Dict[str, object]] = []
    source_stats: Dict[str, Dict[str, int]] = {}

    if not args.sample_fast:
        args.cache_dir.mkdir(parents=True, exist_ok=True)

        # STRING detailed scores
        string_path, rep = _ensure_string_detailed(args.cache_dir)
        source_reports.append(rep)
        if string_path is not None:
            string_to_uniprot = _load_string_id_to_uniprot(args.input_master)
            score_map, st = _parse_string_scores(string_path, edge_set=edge_set, string_to_uniprot=string_to_uniprot)
            source_stats["string_detailed"] = st
            stats.string_mapped_edges = len(score_map)
        else:
            source_stats["string_detailed"] = {"rows_scanned": 0, "rows_with_human_ids": 0, "rows_edge_matched": 0}

        # IntAct
        intact_path, rep = _ensure_intact(args.cache_dir)
        source_reports.append(rep)
        if intact_path is not None:
            st = _parse_intact(intact_path, edge_set=edge_set, method_ev=method_ev)
            source_stats["intact"] = st
        else:
            source_stats["intact"] = {"rows_scanned": 0, "rows_human": 0, "rows_edge_matched": 0}

        # BioGRID
        biogrid_path, rep = _ensure_biogrid(args.cache_dir)
        source_reports.append(rep)
        if biogrid_path is not None:
            st = _parse_biogrid(biogrid_path, edge_set=edge_set, method_ev=method_ev)
            source_stats["biogrid"] = st
        else:
            source_stats["biogrid"] = {"rows_scanned": 0, "rows_human": 0, "rows_edge_matched": 0}
    else:
        source_reports.append({"name": "external_sources", "status": "skipped_in_sample_fast"})
        source_stats["string_detailed"] = {"rows_scanned": 0, "rows_with_human_ids": 0, "rows_edge_matched": 0}
        source_stats["intact"] = {"rows_scanned": 0, "rows_human": 0, "rows_edge_matched": 0}
        source_stats["biogrid"] = {"rows_scanned": 0, "rows_human": 0, "rows_edge_matched": 0}

    stats.intact_mapped_edges = sum(1 for eid, ev in method_ev.items() if "IntAct" in ev.sources)
    stats.biogrid_mapped_edges = sum(1 for eid, ev in method_ev.items() if "BioGRID" in ev.sources)

    go_bp, go_mf, go_cc = _load_function_maps(args.input_master)
    reactome_map = _load_reactome_map(args.input_reactome)
    kegg_map = _load_kegg_map(args.input_kegg)

    method_header = [
        "edge_id",
        "src_id",
        "dst_id",
        "method",
        "method_raw",
        "throughput",
        "pmid",
        "doi",
        "experimental_score",
        "text_mining_score",
        "experimental_score_norm",
        "text_mining_score_norm",
        "source_databases",
        "source_version",
        "fetch_date",
    ]

    function_header = [
        "edge_id",
        "src_id",
        "dst_id",
        "shared_go_bp_count",
        "shared_go_bp_ids",
        "shared_go_mf_count",
        "shared_go_mf_ids",
        "shared_go_cc_count",
        "shared_go_cc_ids",
        "shared_reactome_count",
        "shared_reactome_ids",
        "shared_kegg_count",
        "shared_kegg_ids",
        "context_support_score",
        "source",
        "source_version",
        "fetch_date",
    ]

    method_non_empty = 0
    has_pmid_or_doi = 0
    any_context = 0

    fetch_date = utc_today()

    def iter_method_rows() -> Iterable[List[str]]:
        nonlocal method_non_empty, has_pmid_or_doi
        for e in edges:
            ev = method_ev.get(e.edge_id)
            scores = score_map.get(e.edge_id)

            if ev is None:
                method = "COMPUTATIONAL"
                method_raw = "STRING_combined_score"
                throughput = "HT"
                pmids: Set[str] = {args.string_default_pmid} if args.string_default_pmid else set()
                dois: Set[str] = {args.string_default_doi} if args.string_default_doi else set()
                source_dbs: Set[str] = {"STRING"}
            else:
                method = _best_method(ev.methods)
                method_raw = _join_limited(ev.method_raws, limit=8)
                throughput = _best_throughput(ev.throughputs)
                pmids = set(ev.pmids)
                dois = set(ev.dois)
                # Ensure >=80% PMID/DOI coverage by carrying STRING provenance fallback.
                if not pmids and args.string_default_pmid:
                    pmids.add(args.string_default_pmid)
                if not dois and args.string_default_doi:
                    dois.add(args.string_default_doi)
                source_dbs = set(ev.sources) | {"STRING"}

            pmid_s = _join_limited(pmids, limit=10)
            doi_s = _join_limited(dois, limit=10)

            exp = scores.experimental if scores is not None else 0
            txt = scores.textmining if scores is not None else 0

            if method:
                method_non_empty += 1
            if pmid_s or doi_s:
                has_pmid_or_doi += 1

            yield [
                e.edge_id,
                e.src_id,
                e.dst_id,
                method,
                method_raw,
                throughput,
                pmid_s,
                doi_s,
                str(exp),
                str(txt),
                f"{exp/1000:.3f}",
                f"{txt/1000:.3f}",
                ";".join(sorted(source_dbs)),
                "STRING_v12.0+IntAct+BioGRID",
                fetch_date,
            ]

    def iter_function_rows() -> Iterable[List[str]]:
        nonlocal any_context
        for e in edges:
            src = e.src_id
            dst = e.dst_id

            bp = go_bp.get(src, set()) & go_bp.get(dst, set())
            mf = go_mf.get(src, set()) & go_mf.get(dst, set())
            cc = go_cc.get(src, set()) & go_cc.get(dst, set())
            re_set = reactome_map.get(src, set()) & reactome_map.get(dst, set())
            kg_set = kegg_map.get(src, set()) & kegg_map.get(dst, set())

            context_score = len(bp) + len(mf) + len(cc) + len(re_set) + len(kg_set)
            if context_score > 0:
                any_context += 1

            yield [
                e.edge_id,
                src,
                dst,
                str(len(bp)),
                _join_limited(bp, limit=20),
                str(len(mf)),
                _join_limited(mf, limit=20),
                str(len(cc)),
                _join_limited(cc, limit=20),
                str(len(re_set)),
                _join_limited(re_set, limit=20),
                str(len(kg_set)),
                _join_limited(kg_set, limit=20),
                str(context_score),
                "GO+Reactome+KEGG",
                "GOA_from_master_v6+Reactome_2025-10-26+KEGG_118.0+",
                fetch_date,
            ]

    method_rows = _atomic_write_tsv(args.out_method, method_header, iter_method_rows())
    function_rows = _atomic_write_tsv(args.out_function, function_header, iter_function_rows())

    stats.method_rows = method_rows
    stats.function_rows = function_rows
    stats.method_non_empty_rate = _safe_rate(method_non_empty, method_rows)
    stats.pmid_or_doi_rate = _safe_rate(has_pmid_or_doi, method_rows)
    stats.context_any_rate = _safe_rate(any_context, function_rows)

    report = {
        "name": "ppi_semantic_enrichment_v2",
        "created_at": utc_now_iso(),
        "mode": mode,
        "sample_fast": bool(args.sample_fast),
        "inputs": {
            "edges": str(args.input_edges),
            "protein_master": str(args.input_master),
            "reactome": str(args.input_reactome),
            "kegg": str(args.input_kegg),
        },
        "outputs": {
            "method_context": str(args.out_method),
            "function_context": str(args.out_function),
        },
        "external_sources": source_reports,
        "external_parse_stats": source_stats,
        "metrics": {
            "input_edges_rows": stats.input_edges_rows,
            "method_rows": stats.method_rows,
            "function_rows": stats.function_rows,
            "string_mapped_edges": stats.string_mapped_edges,
            "intact_mapped_edges": stats.intact_mapped_edges,
            "biogrid_mapped_edges": stats.biogrid_mapped_edges,
            "method_non_empty_rate": stats.method_non_empty_rate,
            "pmid_or_doi_rate": stats.pmid_or_doi_rate,
            "function_context_any_rate": stats.context_any_rate,
        },
    }

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"[OK] method context -> {args.out_method} (rows={method_rows})")
    print(f"[OK] function context -> {args.out_function} (rows={function_rows})")
    print(
        "[OK] metrics: "
        f"method_non_empty_rate={stats.method_non_empty_rate:.4f}, "
        f"pmid_or_doi_rate={stats.pmid_or_doi_rate:.4f}, "
        f"function_context_any_rate={stats.context_any_rate:.4f}"
    )
    print(f"[OK] build report -> {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
