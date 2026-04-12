#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sqlite3
import tarfile
import urllib.error
import urllib.request
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
UNIPROT_RE = re.compile(r"^[A-Z0-9]{6,10}(?:-\d+)?$")
NUM_RE = re.compile(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")
VALUE_UNIT_RE = re.compile(
    r"^(IC50|Ki|Kd|EC50)\s*([<>~=]{0,2})\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*([fpnum]?M|M)$",
    re.IGNORECASE,
)

ALLOWED_STANDARD = {"IC50", "Ki", "Kd", "EC50"}
ALLOWED_RELATION = {"=", "<", ">", "<=", ">=", "~"}
UNIT_TO_NM = {
    "M": 1e9,
    "MM": 1e6,
    "UM": 1e3,
    "NM": 1.0,
    "PM": 1e-3,
    "FM": 1e-6,
}

EVIDENCE_HEADER = [
    "evidence_id",
    "compound_inchikey",
    "compound_chembl_id",
    "target_uniprot_accession",
    "target_chembl_id",
    "source_db",
    "source_record_id",
    "assay_type",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_unit",
    "normalized_nM",
    "confidence",
    "confidence_score",
    "reference_doi",
    "reference_pubmed_id",
    "reference_pdb_id",
    "source_raw",
    "source_version",
    "fetch_date",
]

EDGE_HEADER = [
    "edge_id",
    "compound_inchikey",
    "target_uniprot_accession",
    "standard_type",
    "best_normalized_nM",
    "mean_normalized_nM",
    "max_normalized_nM",
    "best_confidence",
    "best_confidence_score",
    "evidence_count",
    "doc_count",
    "source_count",
    "source_dbs",
    "source_version",
    "fetch_date",
]

CONFLICT_HEADER = [
    "compound_inchikey",
    "target_uniprot_accession",
    "standard_type",
    "source_dbs",
    "source_count",
    "evidence_count",
    "min_nM",
    "max_nM",
    "fold_change",
    "conflict_flag",
    "conflict_level",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def utc_today() -> str:
    return datetime.now(timezone.utc).date().isoformat()


def norm(x: object) -> str:
    return str(x or "").strip()


def fmt_float(x: float, digits: int = 8) -> str:
    s = f"{x:.{digits}g}"
    if "e" in s:
        return s.replace("e+", "e")
    return s


def valid_inchikey(x: str) -> bool:
    return bool(INCHIKEY_RE.match(norm(x).upper()))


def normalize_inchikey(x: object) -> str:
    return norm(x).upper()


def normalize_uniprot(x: object) -> str:
    t = norm(x).upper()
    if not t:
        return ""
    # split mixed tokens, keep first valid
    tokens = []
    for token in re.split(r"[,\s;/|]+", t):
        v = token.strip().upper()
        if not v:
            continue
        tokens.append(v)
        if UNIPROT_RE.match(v):
            return v
    return tokens[0] if tokens else ""


def normalize_relation(raw: object) -> str:
    s = norm(raw)
    if not s:
        return "="
    mapping = {
        "==": "=",
        "~=": "~",
        "≈": "~",
        "=<": "<=",
        "=>": ">=",
    }
    s = mapping.get(s, s)
    if s in ALLOWED_RELATION:
        return s
    if "<" in s and "=" in s:
        return "<="
    if ">" in s and "=" in s:
        return ">="
    if "<" in s:
        return "<"
    if ">" in s:
        return ">"
    return "~"


def parse_float(raw: object) -> Optional[float]:
    t = norm(raw)
    if not t:
        return None
    try:
        return float(t)
    except Exception:
        return None


def parse_value_and_relation(raw: str, default_unit: str = "nM") -> Optional[Tuple[str, float, str]]:
    t = norm(raw)
    if not t:
        return None
    s = t.replace("μ", "u").replace("µ", "u").replace(" ", "")
    relation = "="
    if s.startswith("<="):
        relation = "<="
        s = s[2:]
    elif s.startswith(">="):
        relation = ">="
        s = s[2:]
    elif s.startswith("<"):
        relation = "<"
        s = s[1:]
    elif s.startswith(">"):
        relation = ">"
        s = s[1:]
    elif s.startswith("="):
        relation = "="
        s = s[1:]
    elif s.startswith("~"):
        relation = "~"
        s = s[1:]

    m = NUM_RE.search(s)
    if not m:
        return None
    value = float(m.group(0))
    tail = s[m.end():].upper()
    unit = default_unit.upper()
    if tail:
        if tail in UNIT_TO_NM:
            unit = tail
    return relation, value, unit


def parse_pdbbind_affinity(raw: str) -> Optional[Tuple[str, str, float, str]]:
    s = norm(raw).replace("μ", "u").replace("µ", "u")
    if not s:
        return None
    m = VALUE_UNIT_RE.match(s.replace(" ", ""))
    if not m:
        return None
    std = m.group(1)
    rel = normalize_relation(m.group(2) or "=")
    value = float(m.group(3))
    unit = m.group(4).upper()
    std_norm = {"KI": "Ki", "KD": "Kd", "IC50": "IC50", "EC50": "EC50"}[std.upper()]
    return std_norm, rel, value, unit


def to_nM(value: float, unit: str) -> Optional[float]:
    u = norm(unit).upper()
    if u not in UNIT_TO_NM:
        return None
    out = value * UNIT_TO_NM[u]
    if out <= 0:
        return None
    return out


def score_to_label(score: float) -> str:
    if score >= 0.85:
        return "high"
    if score >= 0.70:
        return "medium"
    return "low"


def sha1_token(parts: Sequence[str]) -> str:
    return hashlib.sha1("|".join(parts).encode("utf-8")).hexdigest()


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def month_candidates(back_months: int = 8) -> List[str]:
    now = datetime.now(timezone.utc)
    out: List[str] = []
    for i in range(back_months + 1):
        d = (now.replace(day=1) - timedelta(days=31 * i)).replace(day=1)
        out.append(d.strftime("%Y%m"))
    # keep order, de-dup in case month roll weirdness
    uniq: List[str] = []
    seen = set()
    for x in out:
        if x not in seen:
            uniq.append(x)
            seen.add(x)
    return uniq


def download_if_missing(path: Path, urls: Sequence[str], timeout: int = 60) -> Tuple[Optional[Path], Dict[str, object]]:
    report: Dict[str, object] = {
        "target_path": str(path),
        "status": "missing",
        "attempted_urls": [],
        "error": "",
    }
    if path.exists() and path.stat().st_size > 0:
        report["status"] = "ready_local"
        report["size_bytes"] = path.stat().st_size
        return path, report

    ensure_parent(path)
    for url in urls:
        report["attempted_urls"].append(url)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                data = r.read()
            if not data:
                continue
            path.write_bytes(data)
            report["status"] = "downloaded"
            report["download_url"] = url
            report["size_bytes"] = path.stat().st_size
            return path, report
        except urllib.error.HTTPError as e:
            report["error"] = f"{e.code} {e.reason}"
            continue
        except Exception as e:  # noqa: BLE001
            report["error"] = str(e)
            continue
    report["status"] = "unavailable"
    return None, report


def choose_bindingdb_urls() -> List[str]:
    base = "https://www.bindingdb.org/rwd/bind/downloads"
    urls: List[str] = []
    for mm in month_candidates(10):
        urls.append(f"{base}/BindingDB_PubChem_{mm}_tsv.zip")
    for mm in month_candidates(10):
        urls.append(f"{base}/BindingDB_All_{mm}_tsv.zip")
    return urls


def parse_bindingdb_source_version(path: Path) -> str:
    name = path.name
    if name.endswith(".zip"):
        name = name[:-4]
    return name or "BindingDB_unknown"


def iter_bindingdb_rows(zip_path: Path) -> Tuple[Iterator[Dict[str, str]], List[str]]:
    zf = zipfile.ZipFile(zip_path)
    names = [n for n in zf.namelist() if n.lower().endswith(".tsv")]
    if not names:
        zf.close()
        raise SystemExit(f"[ERROR] no tsv inside zip: {zip_path}")
    target = names[0]
    raw = zf.open(target)
    text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")
    reader = csv.DictReader(text, delimiter="\t")
    if reader.fieldnames is None:
        zf.close()
        raise SystemExit(f"[ERROR] invalid bindingdb header: {zip_path}")
    fieldnames = list(reader.fieldnames)

    def _gen() -> Iterator[Dict[str, str]]:
        try:
            for row in reader:
                yield row
        finally:
            text.close()
            raw.close()
            zf.close()

    return _gen(), fieldnames


def parse_pdb_ids(raw: str) -> List[str]:
    out: List[str] = []
    for token in re.split(r"[,\s;/|]+", norm(raw).upper()):
        if re.match(r"^[0-9][A-Z0-9]{3}$", token):
            out.append(token)
    return sorted(set(out))


def parse_pdbbind_tar(path: Path) -> Tuple[Dict[str, str], List[Tuple[str, str, float, str, str]]]:
    with tarfile.open(path, "r:gz") as tf:
        names = set(tf.getnames())
        affinity_file = None
        name_file = None
        for n in names:
            if n.endswith("INDEX_general_PL.2020"):
                affinity_file = n
            if n.endswith("INDEX_general_PL_name.2020"):
                name_file = n
        if not affinity_file or not name_file:
            raise SystemExit(f"[ERROR] pdbbind index files missing in {path}")

        pdb_to_uniprot: Dict[str, str] = {}
        f_name = tf.extractfile(name_file)
        assert f_name is not None
        for raw_line in f_name.read().decode("utf-8", errors="replace").splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            pdb = parts[0].upper()
            uniprot = normalize_uniprot(parts[2])
            if re.match(r"^[0-9][A-Z0-9]{3}$", pdb) and uniprot:
                pdb_to_uniprot[pdb] = uniprot

        affinity_rows: List[Tuple[str, str, float, str, str]] = []
        f_aff = tf.extractfile(affinity_file)
        assert f_aff is not None
        for raw_line in f_aff.read().decode("utf-8", errors="replace").splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            lhs = line.split("//", 1)[0].strip()
            parts = lhs.split()
            if len(parts) < 4:
                continue
            pdb = parts[0].upper()
            affinity_token = "".join(parts[3:])
            parsed = parse_pdbbind_affinity(affinity_token)
            if not parsed:
                continue
            std, rel, value, unit = parsed
            affinity_rows.append((pdb, std, value, unit, rel))

    return pdb_to_uniprot, affinity_rows


def make_temp_db(path: Path) -> sqlite3.Connection:
    ensure_parent(path)
    if path.exists():
        path.unlink()
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute(
        """
        CREATE TABLE evidence_norm (
            evidence_id TEXT PRIMARY KEY,
            compound_inchikey TEXT NOT NULL,
            target_uniprot_accession TEXT NOT NULL,
            standard_type TEXT NOT NULL,
            normalized_nM REAL NOT NULL,
            source_db TEXT NOT NULL,
            confidence_score REAL NOT NULL,
            confidence TEXT NOT NULL,
            reference_key TEXT
        )
        """
    )
    conn.commit()
    return conn


def chembl_query(max_rows: Optional[int]) -> str:
    base = """
        SELECT
            activity_id,
            assay_id,
            doc_id,
            compound_chembl_id,
            compound_inchikey,
            target_chembl_id,
            target_uniprot_accession,
            assay_type,
            standard_type,
            standard_relation,
            standard_value,
            standard_units,
            standard_value_nM,
            doi,
            pubmed_id,
            source,
            evidence_score_v1
        FROM psi_evidence_v1
        WHERE standard_type IN ('IC50','Ki','Kd','EC50')
          AND compound_inchikey IS NOT NULL
          AND target_uniprot_accession IS NOT NULL
          AND TRIM(compound_inchikey) != ''
          AND TRIM(target_uniprot_accession) != ''
          AND standard_value_nM IS NOT NULL
    """
    if max_rows is not None and max_rows > 0:
        base += f"\n LIMIT {int(max_rows)}"
    return base


def count_chembl_baseline(conn: sqlite3.Connection) -> Dict[str, int]:
    q_rows = """
        SELECT COUNT(*)
        FROM psi_evidence_v1
        WHERE standard_type IN ('IC50','Ki','Kd','EC50')
          AND compound_inchikey IS NOT NULL
          AND target_uniprot_accession IS NOT NULL
          AND TRIM(compound_inchikey) != ''
          AND TRIM(target_uniprot_accession) != ''
          AND standard_value_nM IS NOT NULL
          AND CAST(standard_value_nM AS REAL) > 0
    """
    q_edges = """
        SELECT COUNT(*)
        FROM (
            SELECT compound_inchikey, target_uniprot_accession, standard_type
            FROM psi_evidence_v1
            WHERE standard_type IN ('IC50','Ki','Kd','EC50')
              AND compound_inchikey IS NOT NULL
              AND target_uniprot_accession IS NOT NULL
              AND TRIM(compound_inchikey) != ''
              AND TRIM(target_uniprot_accession) != ''
              AND standard_value_nM IS NOT NULL
              AND CAST(standard_value_nM AS REAL) > 0
            GROUP BY compound_inchikey, target_uniprot_accession, standard_type
        )
    """
    return {
        "chembl_baseline_rows": int(conn.execute(q_rows).fetchone()[0]),
        "chembl_baseline_edges": int(conn.execute(q_edges).fetchone()[0]),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chembl-m3-db", type=Path, required=True)
    ap.add_argument("--out-evidence", type=Path, required=True)
    ap.add_argument("--out-edges", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--conflict-audit-tsv", type=Path, required=True)
    ap.add_argument("--conflict-audit-json", type=Path, required=True)
    ap.add_argument("--manual-download-report", type=Path, required=True)
    ap.add_argument("--bindingdb-zip", type=Path, required=True)
    ap.add_argument("--pdbbind-index-tar", type=Path, required=True)
    ap.add_argument("--cache-dir", type=Path, required=True)
    ap.add_argument("--fetch-date", default=utc_today())
    ap.add_argument("--max-chembl-rows", type=int, default=0)
    ap.add_argument("--no-auto-download", action="store_true")
    args = ap.parse_args()

    for p in [
        args.out_evidence,
        args.out_edges,
        args.report,
        args.conflict_audit_tsv,
        args.conflict_audit_json,
        args.manual_download_report,
    ]:
        ensure_parent(p)
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    if not args.chembl_m3_db.exists():
        raise SystemExit(f"[ERROR] missing required input: {args.chembl_m3_db}")

    download_notes: Dict[str, object] = {}
    manual_items: List[Dict[str, str]] = []

    bindingdb_path: Optional[Path] = args.bindingdb_zip if args.bindingdb_zip.exists() else None
    if bindingdb_path is None and not args.no_auto_download:
        resolved, note = download_if_missing(args.bindingdb_zip, choose_bindingdb_urls())
        bindingdb_path = resolved
        download_notes["bindingdb"] = note
    elif bindingdb_path is None:
        download_notes["bindingdb"] = {
            "status": "skipped_no_auto_download",
            "target_path": str(args.bindingdb_zip),
        }
    else:
        download_notes["bindingdb"] = {
            "status": "ready_local",
            "target_path": str(bindingdb_path),
            "size_bytes": bindingdb_path.stat().st_size,
        }

    if bindingdb_path is None:
        manual_items.append(
            {
                "name": "BindingDB PubChem TSV ZIP",
                "expected_path": str(args.bindingdb_zip),
                "example_url": "https://www.bindingdb.org/rwd/bind/downloads/BindingDB_PubChem_202604_tsv.zip",
                "sha256_cmd": f"sha256sum {args.bindingdb_zip}",
            }
        )

    pdbbind_path: Optional[Path] = args.pdbbind_index_tar if args.pdbbind_index_tar.exists() else None
    if pdbbind_path is None and not args.no_auto_download:
        resolved, note = download_if_missing(
            args.pdbbind_index_tar,
            [
                "https://static.pdbbind-plus.org.cn/v2020-renew_website/PDBbind_v2020_plain_text_index.tar.gz",
            ],
        )
        pdbbind_path = resolved
        download_notes["pdbbind"] = note
    elif pdbbind_path is None:
        download_notes["pdbbind"] = {
            "status": "skipped_no_auto_download",
            "target_path": str(args.pdbbind_index_tar),
        }
    else:
        download_notes["pdbbind"] = {
            "status": "ready_local",
            "target_path": str(pdbbind_path),
            "size_bytes": pdbbind_path.stat().st_size,
        }

    if pdbbind_path is None:
        manual_items.append(
            {
                "name": "PDBbind v2020 plain text index",
                "expected_path": str(args.pdbbind_index_tar),
                "example_url": "https://static.pdbbind-plus.org.cn/v2020-renew_website/PDBbind_v2020_plain_text_index.tar.gz",
                "sha256_cmd": f"sha256sum {args.pdbbind_index_tar}",
            }
        )

    args.manual_download_report.write_text(
        json.dumps(
            {
                "name": "molecule_activity_fusion_v2.manual_download",
                "created_at": utc_now(),
                "optional_missing_count": len(manual_items),
                "download_notes": download_notes,
                "download_checklist": manual_items,
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    tmp_db = args.cache_dir / "molecule_activity_fusion_tmp.sqlite"
    agg_conn = make_temp_db(tmp_db)
    insert_sql = """
        INSERT OR REPLACE INTO evidence_norm(
            evidence_id, compound_inchikey, target_uniprot_accession, standard_type,
            normalized_nM, source_db, confidence_score, confidence, reference_key
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    insert_batch: List[Tuple[str, str, str, str, float, str, float, str, str]] = []

    source_counts = Counter()
    source_rows = Counter()
    skipped = Counter()
    notes: List[str] = []
    source_versions = set()
    pdb_to_inchikey_counts: Dict[str, Counter[str]] = defaultdict(Counter)

    max_chembl_rows = args.max_chembl_rows if args.max_chembl_rows and args.max_chembl_rows > 0 else None

    with args.out_evidence.open("w", encoding="utf-8", newline="") as f_ev:
        w_ev = csv.DictWriter(f_ev, delimiter="\t", fieldnames=EVIDENCE_HEADER, lineterminator="\n")
        w_ev.writeheader()

        # -------- ChEMBL --------
        with sqlite3.connect(args.chembl_m3_db) as cconn:
            baseline = count_chembl_baseline(cconn)
            cur = cconn.execute(chembl_query(max_chembl_rows))
            for row in cur:
                (
                    activity_id,
                    assay_id,
                    doc_id,
                    compound_chembl_id,
                    compound_inchikey,
                    target_chembl_id,
                    target_uniprot_accession,
                    assay_type,
                    standard_type,
                    standard_relation,
                    standard_value,
                    standard_units,
                    standard_value_nM,
                    doi,
                    pubmed_id,
                    source,
                    evidence_score_v1,
                ) = row
                ik = normalize_inchikey(compound_inchikey)
                uniprot = norm(target_uniprot_accession).upper()
                std = norm(standard_type)
                if not valid_inchikey(ik) or not uniprot or std not in ALLOWED_STANDARD:
                    skipped["chembl_invalid_key_or_type"] += 1
                    continue

                nm = parse_float(standard_value_nM)
                if nm is None or nm <= 0:
                    skipped["chembl_missing_nm"] += 1
                    continue

                conf_score = parse_float(evidence_score_v1)
                if conf_score is None:
                    conf_score = 0.70
                conf_score = max(0.0, min(conf_score, 0.99))
                conf_label = score_to_label(conf_score)

                relation = normalize_relation(standard_relation)
                std_value = parse_float(standard_value)
                if std_value is None:
                    std_value = nm
                std_unit = norm(standard_units) or "nM"
                if std_unit.lower() == "nm":
                    std_unit = "nM"

                evidence_id = sha1_token(["chembl", str(activity_id), ik, uniprot, std, fmt_float(nm)])
                out = {
                    "evidence_id": evidence_id,
                    "compound_inchikey": ik,
                    "compound_chembl_id": norm(compound_chembl_id).upper(),
                    "target_uniprot_accession": uniprot,
                    "target_chembl_id": norm(target_chembl_id).upper(),
                    "source_db": "chembl",
                    "source_record_id": str(activity_id),
                    "assay_type": norm(assay_type).upper() or "UNK",
                    "standard_type": std,
                    "standard_relation": relation,
                    "standard_value": fmt_float(std_value),
                    "standard_unit": std_unit,
                    "normalized_nM": fmt_float(nm),
                    "confidence": conf_label,
                    "confidence_score": f"{conf_score:.3f}",
                    "reference_doi": norm(doi),
                    "reference_pubmed_id": norm(pubmed_id),
                    "reference_pdb_id": "",
                    "source_raw": norm(source) or "chembl_36",
                    "source_version": norm(source) or "chembl_36",
                    "fetch_date": args.fetch_date,
                }
                w_ev.writerow(out)

                ref_key = norm(doi) or norm(pubmed_id) or str(doc_id or "")
                insert_batch.append(
                    (
                        evidence_id,
                        ik,
                        uniprot,
                        std,
                        float(nm),
                        "chembl",
                        conf_score,
                        conf_label,
                        ref_key,
                    )
                )
                if len(insert_batch) >= 12000:
                    agg_conn.executemany(insert_sql, insert_batch)
                    agg_conn.commit()
                    insert_batch.clear()

                source_counts["chembl"] += 1
                source_rows["chembl"] += 1
                source_versions.add(norm(source) or "chembl_36")

        # -------- BindingDB --------
        bindingdb_stats = Counter()
        bindingdb_version = ""
        if bindingdb_path and bindingdb_path.exists():
            bindingdb_version = parse_bindingdb_source_version(bindingdb_path)
            source_versions.add(bindingdb_version)
            b_iter, b_fields = iter_bindingdb_rows(bindingdb_path)
            uniprot_cols = [c for c in b_fields if "Primary ID of Target Chain" in c]
            if "UniProt (SwissProt) Primary ID of Target Chain 1" in b_fields:
                uniprot_cols.insert(0, uniprot_cols.pop(uniprot_cols.index("UniProt (SwissProt) Primary ID of Target Chain 1")))
            if "UniProt (TrEMBL) Primary ID of Target Chain 1" in b_fields and "UniProt (TrEMBL) Primary ID of Target Chain 1" in uniprot_cols:
                uniprot_cols.insert(1, uniprot_cols.pop(uniprot_cols.index("UniProt (TrEMBL) Primary ID of Target Chain 1")))

            seen_bindingdb_keys = set()
            affinity_cols = [
                ("Ki (nM)", "Ki"),
                ("IC50 (nM)", "IC50"),
                ("Kd (nM)", "Kd"),
                ("EC50 (nM)", "EC50"),
            ]
            for row in b_iter:
                bindingdb_stats["rows_total"] += 1
                ik = normalize_inchikey(row.get("Ligand InChI Key", ""))
                if not valid_inchikey(ik):
                    bindingdb_stats["skip_invalid_inchikey"] += 1
                    continue
                uniprot = ""
                for c in uniprot_cols:
                    uniprot = normalize_uniprot(row.get(c, ""))
                    if uniprot:
                        break
                if not uniprot:
                    bindingdb_stats["skip_missing_uniprot"] += 1
                    continue

                reactant_id = norm(row.get("BindingDB Reactant_set_id", ""))
                doi = norm(row.get("Article DOI", "")) or norm(row.get("BindingDB Entry DOI", ""))
                pmid = norm(row.get("PMID", ""))
                pdb_ids = parse_pdb_ids(row.get("PDB ID(s) for Ligand-Target Complex", ""))
                for pdb_id in pdb_ids:
                    pdb_to_inchikey_counts[pdb_id][ik] += 1

                for col, std in affinity_cols:
                    raw = norm(row.get(col, ""))
                    if not raw:
                        continue
                    parsed = parse_value_and_relation(raw, default_unit="nM")
                    if not parsed:
                        bindingdb_stats["skip_unparseable_affinity"] += 1
                        continue
                    rel, value, unit = parsed
                    nm = to_nM(value, unit)
                    if nm is None:
                        bindingdb_stats["skip_invalid_nm"] += 1
                        continue

                    dedup_key = (ik, uniprot, std, rel, round(nm, 8), doi, pmid, reactant_id)
                    if dedup_key in seen_bindingdb_keys:
                        bindingdb_stats["skip_duplicate_record"] += 1
                        continue
                    seen_bindingdb_keys.add(dedup_key)

                    conf_score = 0.72
                    if doi or pmid:
                        conf_score += 0.15
                    if rel == "=":
                        conf_score += 0.05
                    conf_score = min(conf_score, 0.97)
                    conf_label = score_to_label(conf_score)

                    source_record_id = (
                        f"{reactant_id}:{std}:{bindingdb_stats['rows_total']}"
                        if reactant_id
                        else f"{std}:{ik}:{uniprot}:{bindingdb_stats['rows_total']}"
                    )
                    evidence_id = sha1_token(
                        ["bindingdb", source_record_id, ik, uniprot, std, fmt_float(nm), rel]
                    )
                    out = {
                        "evidence_id": evidence_id,
                        "compound_inchikey": ik,
                        "compound_chembl_id": norm(row.get("ChEMBL ID of Ligand", "")).upper(),
                        "target_uniprot_accession": uniprot,
                        "target_chembl_id": "",
                        "source_db": "bindingdb",
                        "source_record_id": source_record_id,
                        "assay_type": "B",
                        "standard_type": std,
                        "standard_relation": rel,
                        "standard_value": fmt_float(value),
                        "standard_unit": "nM" if unit == "NM" else unit,
                        "normalized_nM": fmt_float(nm),
                        "confidence": conf_label,
                        "confidence_score": f"{conf_score:.3f}",
                        "reference_doi": doi,
                        "reference_pubmed_id": pmid,
                        "reference_pdb_id": ";".join(pdb_ids),
                        "source_raw": norm(row.get("Curation/DataSource", "")) or "BindingDB",
                        "source_version": bindingdb_version,
                        "fetch_date": args.fetch_date,
                    }
                    w_ev.writerow(out)

                    ref_key = doi or pmid or source_record_id
                    insert_batch.append(
                        (
                            evidence_id,
                            ik,
                            uniprot,
                            std,
                            float(nm),
                            "bindingdb",
                            conf_score,
                            conf_label,
                            ref_key,
                        )
                    )
                    if len(insert_batch) >= 12000:
                        agg_conn.executemany(insert_sql, insert_batch)
                        agg_conn.commit()
                        insert_batch.clear()

                    source_counts["bindingdb"] += 1
                    source_rows["bindingdb"] += 1
                    bindingdb_stats["evidence_rows"] += 1
        else:
            notes.append("BindingDB unavailable: skipped")
            bindingdb_stats = Counter({"rows_total": 0, "evidence_rows": 0})

        # -------- PDBbind --------
        pdbbind_stats = Counter()
        if pdbbind_path and pdbbind_path.exists():
            try:
                pdb_to_uniprot, affinity_rows = parse_pdbbind_tar(pdbbind_path)
                source_versions.add("PDBbind_v2020")
                pdbbind_stats["affinity_rows_parsed"] = len(affinity_rows)
                pdbbind_stats["pdb_uniprot_map_size"] = len(pdb_to_uniprot)
                seen_pdbbind_keys = set()
                for pdb_id, std, value, unit, rel in affinity_rows:
                    pdbbind_stats["rows_total"] += 1
                    if pdb_id not in pdb_to_uniprot:
                        pdbbind_stats["skip_missing_uniprot"] += 1
                        continue
                    if pdb_id not in pdb_to_inchikey_counts:
                        pdbbind_stats["skip_missing_inchikey_mapping"] += 1
                        continue
                    uniprot = pdb_to_uniprot[pdb_id]
                    ik_counter = pdb_to_inchikey_counts[pdb_id]
                    if not ik_counter:
                        pdbbind_stats["skip_empty_inchikey_mapping"] += 1
                        continue
                    ik, top_cnt = ik_counter.most_common(1)[0]
                    total_cnt = sum(ik_counter.values())
                    ambiguous = len(ik_counter) > 1 and (top_cnt / total_cnt) < 0.7

                    nm = to_nM(value, unit)
                    if nm is None:
                        pdbbind_stats["skip_invalid_nm"] += 1
                        continue

                    dedup_key = (pdb_id, ik, uniprot, std, rel, round(nm, 8))
                    if dedup_key in seen_pdbbind_keys:
                        pdbbind_stats["skip_duplicate_record"] += 1
                        continue
                    seen_pdbbind_keys.add(dedup_key)

                    conf_score = 0.78
                    if top_cnt >= 2:
                        conf_score += 0.05
                    if ambiguous:
                        conf_score -= 0.10
                    conf_score = max(0.55, min(conf_score, 0.95))
                    conf_label = score_to_label(conf_score)

                    source_record_id = f"{pdb_id}:{std}"
                    evidence_id = sha1_token(
                        ["pdbbind", source_record_id, ik, uniprot, std, fmt_float(nm), rel]
                    )
                    out = {
                        "evidence_id": evidence_id,
                        "compound_inchikey": ik,
                        "compound_chembl_id": "",
                        "target_uniprot_accession": uniprot,
                        "target_chembl_id": "",
                        "source_db": "pdbbind",
                        "source_record_id": source_record_id,
                        "assay_type": "B",
                        "standard_type": std,
                        "standard_relation": rel,
                        "standard_value": fmt_float(value),
                        "standard_unit": unit if unit != "NM" else "nM",
                        "normalized_nM": fmt_float(nm),
                        "confidence": conf_label,
                        "confidence_score": f"{conf_score:.3f}",
                        "reference_doi": "",
                        "reference_pubmed_id": "",
                        "reference_pdb_id": pdb_id,
                        "source_raw": "PDBbind_general_PL_index + BindingDB_PDB_ligand_mapping",
                        "source_version": "PDBbind_v2020",
                        "fetch_date": args.fetch_date,
                    }
                    w_ev.writerow(out)

                    ref_key = pdb_id
                    insert_batch.append(
                        (
                            evidence_id,
                            ik,
                            uniprot,
                            std,
                            float(nm),
                            "pdbbind",
                            conf_score,
                            conf_label,
                            ref_key,
                        )
                    )
                    if len(insert_batch) >= 12000:
                        agg_conn.executemany(insert_sql, insert_batch)
                        agg_conn.commit()
                        insert_batch.clear()

                    source_counts["pdbbind"] += 1
                    source_rows["pdbbind"] += 1
                    pdbbind_stats["evidence_rows"] += 1
            except Exception as e:  # noqa: BLE001
                notes.append(f"PDBbind parse failed: {e}")
        else:
            notes.append("PDBbind unavailable: skipped")
            pdbbind_stats = Counter({"rows_total": 0, "evidence_rows": 0})

    if insert_batch:
        agg_conn.executemany(insert_sql, insert_batch)
        agg_conn.commit()
        insert_batch.clear()

    # build aggregation indexes
    agg_conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_ev_group ON evidence_norm(compound_inchikey, target_uniprot_accession, standard_type)"
    )
    agg_conn.execute("CREATE INDEX IF NOT EXISTS idx_ev_source ON evidence_norm(source_db)")
    agg_conn.commit()

    # -------- edges + conflict audit --------
    edge_rows = 0
    conflict_groups = 0
    conflict_flagged = 0
    multi_source_groups = 0
    source_version_joined = ";".join(sorted(source_versions))

    q_agg = """
        SELECT
            compound_inchikey,
            target_uniprot_accession,
            standard_type,
            COUNT(*) AS evidence_count,
            COUNT(DISTINCT source_db) AS source_count,
            GROUP_CONCAT(DISTINCT source_db) AS source_dbs,
            MIN(normalized_nM) AS min_nM,
            AVG(normalized_nM) AS mean_nM,
            MAX(normalized_nM) AS max_nM,
            MAX(confidence_score) AS max_confidence,
            COUNT(DISTINCT CASE WHEN reference_key IS NOT NULL AND reference_key != '' THEN reference_key END) AS doc_count
        FROM evidence_norm
        GROUP BY compound_inchikey, target_uniprot_accession, standard_type
    """

    with args.out_edges.open("w", encoding="utf-8", newline="") as f_edge, args.conflict_audit_tsv.open(
        "w", encoding="utf-8", newline=""
    ) as f_conf:
        w_edge = csv.DictWriter(f_edge, delimiter="\t", fieldnames=EDGE_HEADER, lineterminator="\n")
        w_conf = csv.DictWriter(f_conf, delimiter="\t", fieldnames=CONFLICT_HEADER, lineterminator="\n")
        w_edge.writeheader()
        w_conf.writeheader()

        for (
            ik,
            uniprot,
            std,
            n_ev,
            n_src,
            srcs_raw,
            min_nm,
            mean_nm,
            max_nm,
            max_conf,
            doc_count,
        ) in agg_conn.execute(q_agg):
            edge_id = sha1_token([ik, uniprot, std])
            src_tokens = sorted({s.strip() for s in norm(srcs_raw).split(",") if s.strip()})
            src_join = ";".join(src_tokens)
            max_conf_f = float(max_conf or 0.0)

            w_edge.writerow(
                {
                    "edge_id": edge_id,
                    "compound_inchikey": ik,
                    "target_uniprot_accession": uniprot,
                    "standard_type": std,
                    "best_normalized_nM": fmt_float(float(min_nm)),
                    "mean_normalized_nM": fmt_float(float(mean_nm)),
                    "max_normalized_nM": fmt_float(float(max_nm)),
                    "best_confidence": score_to_label(max_conf_f),
                    "best_confidence_score": f"{max_conf_f:.3f}",
                    "evidence_count": str(int(n_ev)),
                    "doc_count": str(int(doc_count or 0)),
                    "source_count": str(int(n_src)),
                    "source_dbs": src_join,
                    "source_version": source_version_joined,
                    "fetch_date": args.fetch_date,
                }
            )
            edge_rows += 1

            if int(n_src) > 1:
                multi_source_groups += 1
                min_f = float(min_nm)
                max_f = float(max_nm)
                fold = (max_f / min_f) if min_f > 0 else 0.0
                flag = fold >= 10.0
                if flag:
                    conflict_flagged += 1
                level = "none"
                if fold >= 100:
                    level = "high"
                elif fold >= 10:
                    level = "medium"
                elif fold >= 3:
                    level = "low"

                w_conf.writerow(
                    {
                        "compound_inchikey": ik,
                        "target_uniprot_accession": uniprot,
                        "standard_type": std,
                        "source_dbs": src_join,
                        "source_count": str(int(n_src)),
                        "evidence_count": str(int(n_ev)),
                        "min_nM": fmt_float(min_f),
                        "max_nM": fmt_float(max_f),
                        "fold_change": fmt_float(fold),
                        "conflict_flag": "1" if flag else "0",
                        "conflict_level": level,
                    }
                )
                conflict_groups += 1

    total_evidence_rows = int(
        agg_conn.execute("SELECT COUNT(*) FROM evidence_norm").fetchone()[0]
    )
    chembl_rows_v2 = int(
        agg_conn.execute("SELECT COUNT(*) FROM evidence_norm WHERE source_db='chembl'").fetchone()[0]
    )
    chembl_edges_v2 = int(
        agg_conn.execute(
            """
            SELECT COUNT(*)
            FROM (
                SELECT compound_inchikey, target_uniprot_accession, standard_type
                FROM evidence_norm
                WHERE source_db='chembl'
                GROUP BY compound_inchikey, target_uniprot_accession, standard_type
            )
            """
        ).fetchone()[0]
    )
    agg_conn.close()

    conflict_json = {
        "name": "molecule_activity_conflict_audit_v2",
        "created_at": utc_now(),
        "table": str(args.conflict_audit_tsv),
        "summary": {
            "multi_source_groups": multi_source_groups,
            "conflict_groups_in_audit_table": conflict_groups,
            "flagged_conflicts_fold_ge_10": conflict_flagged,
            "flagged_conflict_rate": (conflict_flagged / conflict_groups) if conflict_groups else 0.0,
        },
    }
    args.conflict_audit_json.write_text(json.dumps(conflict_json, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    build_report = {
        "name": "molecule_activity_fusion_v2.build",
        "created_at": utc_now(),
        "inputs": {
            "chembl_m3_db": str(args.chembl_m3_db),
            "bindingdb_zip": str(bindingdb_path) if bindingdb_path else "",
            "pdbbind_index_tar": str(pdbbind_path) if pdbbind_path else "",
        },
        "download_notes": download_notes,
        "manual_download_report": str(args.manual_download_report),
        "max_chembl_rows": int(args.max_chembl_rows or 0),
        "metrics": {
            "total_evidence_rows": total_evidence_rows,
            "total_edge_rows": edge_rows,
            "source_rows": dict(source_rows),
            "chembl_rows_in_v2": chembl_rows_v2,
            "chembl_edges_in_v2": chembl_edges_v2,
            "multi_source_groups": multi_source_groups,
            "conflict_groups_in_audit_table": conflict_groups,
            "flagged_conflicts_fold_ge_10": conflict_flagged,
            **baseline,
        },
        "source_stats": {
            "bindingdb": dict(bindingdb_stats),
            "pdbbind": dict(pdbbind_stats),
        },
        "skipped": dict(skipped),
        "source_versions": sorted(source_versions),
        "notes": notes,
        "outputs": {
            "evidence": str(args.out_evidence),
            "edges": str(args.out_edges),
            "conflict_audit_tsv": str(args.conflict_audit_tsv),
            "conflict_audit_json": str(args.conflict_audit_json),
        },
    }
    args.report.write_text(json.dumps(build_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(
        f"[DONE] evidence_rows={total_evidence_rows} edges={edge_rows} "
        f"sources={dict(source_rows)} conflicts={conflict_flagged}/{conflict_groups}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
