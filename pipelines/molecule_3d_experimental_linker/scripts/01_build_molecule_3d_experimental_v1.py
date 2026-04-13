#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import requests

RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"
USER_AGENT = "kg-molecule-3d-experimental-linker/1.0"

INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
PDB_RE = re.compile(r"^[0-9A-Z]{4}$")

GRAPHQL_QUERY = """
query($ids:[String!]!) {
  entries(entry_ids:$ids) {
    rcsb_id
    nonpolymer_entities {
      rcsb_nonpolymer_entity_container_identifiers {
        entity_id
        nonpolymer_comp_id
      }
      nonpolymer_comp {
        chem_comp { id }
        rcsb_chem_comp_descriptor {
          InChIKey
        }
      }
    }
  }
}
""".strip()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def utc_today() -> str:
    return datetime.now(timezone.utc).date().isoformat()


def normalize(x: str) -> str:
    return (x or "").strip()


def valid_inchikey(x: str) -> bool:
    return bool(INCHIKEY_RE.match(normalize(x).upper()))


def valid_pdb_id(x: str) -> bool:
    return bool(PDB_RE.match(normalize(x).upper()))


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


@dataclass
class TargetRec:
    pdb_id: str
    uniprot_id: str
    experimental_method: str
    resolution: str
    release_date: str


@dataclass
class LigandRec:
    pdb_id: str
    entity_id: str
    ligand_ccd_id: str
    inchikey: str


def load_xref_inchikeys(path: Path, max_rows: Optional[int]) -> Tuple[Set[str], Dict[str, Any]]:
    cols, rows = read_tsv(path, max_rows=max_rows)
    if "inchikey" not in cols:
        raise SystemExit(f"[ERROR] xref table missing inchikey: {path}")

    out: Set[str] = set()
    bad_ik = 0
    for row in rows:
        ik = normalize(row.get("inchikey", "")).upper()
        if not valid_inchikey(ik):
            bad_ik += 1
            continue
        out.add(ik)

    stats = {
        "input_rows": len(rows),
        "valid_inchikey_rows": len(out),
        "skipped_bad_inchikey_rows": bad_ik,
    }
    return out, stats


def load_protein_pdb_targets(path: Path, max_pdb: Optional[int]) -> Tuple[Dict[str, List[TargetRec]], Dict[str, Any]]:
    cols, rows = read_tsv(path, max_rows=None)
    required = {"pdb_id", "uniprot_id", "experimental_method", "resolution", "release_date"}
    missing = required - set(cols)
    if missing:
        raise SystemExit(f"[ERROR] protein PDB table missing columns: {sorted(missing)}")

    by_pdb: Dict[str, List[TargetRec]] = defaultdict(list)
    target_seen: Set[Tuple[str, str]] = set()

    metrics = Counter()

    for row in rows:
        metrics["input_rows"] += 1
        pdb_id = normalize(row.get("pdb_id", "")).upper()
        uniprot_id = normalize(row.get("uniprot_id", ""))

        if not valid_pdb_id(pdb_id):
            metrics["skip_bad_pdb_id"] += 1
            continue
        if not uniprot_id:
            metrics["skip_empty_uniprot"] += 1
            continue

        key = (pdb_id, uniprot_id)
        if key in target_seen:
            metrics["skip_duplicate_target_pair"] += 1
            continue
        target_seen.add(key)

        by_pdb[pdb_id].append(
            TargetRec(
                pdb_id=pdb_id,
                uniprot_id=uniprot_id,
                experimental_method=normalize(row.get("experimental_method", "")),
                resolution=normalize(row.get("resolution", "")),
                release_date=normalize(row.get("release_date", "")),
            )
        )

    ordered_pdb = sorted(by_pdb.keys())
    if max_pdb is not None and max_pdb >= 0:
        keep = set(ordered_pdb[:max_pdb])
        by_pdb = {p: by_pdb[p] for p in ordered_pdb if p in keep}
        metrics["max_pdb_applied"] = int(max_pdb)

    metrics["unique_pdb"] = len(by_pdb)
    metrics["unique_target_pairs"] = sum(len(v) for v in by_pdb.values())

    return by_pdb, {k: int(v) for k, v in metrics.items()}


def _fetch_ligands_batch(session: requests.Session, pdb_ids: Sequence[str], timeout: int) -> Dict[str, List[LigandRec]]:
    payload = {"query": GRAPHQL_QUERY, "variables": {"ids": list(pdb_ids)}}
    resp = session.post(RCSB_GRAPHQL_URL, json=payload, timeout=timeout)
    if resp.status_code != 200:
        raise RuntimeError(f"HTTP {resp.status_code}")

    data = resp.json()
    if data.get("errors"):
        raise RuntimeError(f"graphql_errors={data.get('errors')}")

    out: Dict[str, List[LigandRec]] = defaultdict(list)

    for entry in (data.get("data") or {}).get("entries") or []:
        pdb_id = normalize(str((entry or {}).get("rcsb_id") or "")).upper()
        if not valid_pdb_id(pdb_id):
            continue

        for np in (entry or {}).get("nonpolymer_entities") or []:
            ids = (np or {}).get("rcsb_nonpolymer_entity_container_identifiers") or {}
            comp = (np or {}).get("nonpolymer_comp") or {}
            descriptor = comp.get("rcsb_chem_comp_descriptor") or {}
            chem_comp = comp.get("chem_comp") or {}

            inchikey = normalize(str(descriptor.get("InChIKey") or "")).upper()
            if not valid_inchikey(inchikey):
                continue

            entity_id = normalize(str(ids.get("entity_id") or ""))
            ligand_ccd_id = normalize(str(ids.get("nonpolymer_comp_id") or chem_comp.get("id") or "")).upper()

            out[pdb_id].append(
                LigandRec(
                    pdb_id=pdb_id,
                    entity_id=entity_id,
                    ligand_ccd_id=ligand_ccd_id,
                    inchikey=inchikey,
                )
            )

    return out


def fetch_pdb_ligand_map(
    pdb_ids: Sequence[str],
    *,
    batch_size: int,
    timeout: int,
    sleep_seconds: float,
    retries: int,
) -> Tuple[Dict[str, List[LigandRec]], Dict[str, Any]]:
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})

    ordered = sorted({p for p in pdb_ids if valid_pdb_id(p)})

    out: Dict[str, List[LigandRec]] = defaultdict(list)
    failed_singletons: Set[str] = set()
    metrics = Counter()

    def fetch_recursive(chunk: List[str]) -> None:
        if not chunk:
            return

        for attempt in range(retries + 1):
            metrics["graphql_http_calls"] += 1
            try:
                got = _fetch_ligands_batch(session, chunk, timeout=timeout)
                for pdb_id, rows in got.items():
                    out[pdb_id].extend(rows)
                metrics["graphql_http_success_calls"] += 1
                return
            except Exception:
                metrics["graphql_http_failed_calls"] += 1
                if attempt < retries:
                    time.sleep(0.35 * (attempt + 1))
                else:
                    break

        if len(chunk) == 1:
            failed_singletons.add(chunk[0])
            return

        mid = len(chunk) // 2
        fetch_recursive(chunk[:mid])
        fetch_recursive(chunk[mid:])

    for i in range(0, len(ordered), batch_size):
        chunk = ordered[i : i + batch_size]
        fetch_recursive(chunk)
        if sleep_seconds > 0:
            time.sleep(sleep_seconds)

    # per-entry ligand dedupe (same ligand can appear multiple times in API payload)
    deduped: Dict[str, List[LigandRec]] = {}
    for pdb_id, rows in out.items():
        seen: Set[Tuple[str, str, str]] = set()
        keep: List[LigandRec] = []
        for r in rows:
            key = (r.entity_id, r.ligand_ccd_id, r.inchikey)
            if key in seen:
                continue
            seen.add(key)
            keep.append(r)
        deduped[pdb_id] = keep

    report = {
        "requested_pdb_ids": len(ordered),
        "resolved_pdb_ids": len(deduped),
        "failed_pdb_ids": len(failed_singletons),
        "failed_pdb_id_sample": sorted(list(failed_singletons))[:100],
        "pdb_with_ligands": sum(1 for _, rows in deduped.items() if rows),
        "ligand_rows_total": int(sum(len(v) for v in deduped.values())),
        **{k: int(v) for k, v in metrics.items()},
    }
    return deduped, report


def non_empty_rate(rows: List[Dict[str, str]], col: str) -> Dict[str, Any]:
    n = sum(1 for r in rows if normalize(r.get(col, "")) != "")
    total = len(rows)
    return {
        "non_empty": n,
        "total": total,
        "rate": (n / total) if total else 0.0,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xref", type=Path, required=True)
    ap.add_argument("--protein-pdb", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--coverage-report", type=Path, required=True)
    ap.add_argument("--max-rows", type=int, default=None)
    ap.add_argument("--max-pdb", type=int, default=None)
    ap.add_argument("--batch-size", type=int, default=120)
    ap.add_argument("--http-timeout", type=int, default=30)
    ap.add_argument("--http-retries", type=int, default=2)
    ap.add_argument("--http-sleep", type=float, default=0.08)
    args = ap.parse_args()

    if not args.xref.exists():
        raise SystemExit(f"[ERROR] missing xref input: {args.xref}")
    if not args.protein_pdb.exists():
        raise SystemExit(f"[ERROR] missing protein pdb input: {args.protein_pdb}")

    created_at = utc_now()
    fetch_date = utc_today()

    xref_iks, xref_stats = load_xref_inchikeys(args.xref, args.max_rows)
    pdb_targets, protein_stats = load_protein_pdb_targets(args.protein_pdb, args.max_pdb)

    ligand_by_pdb, graphql_report = fetch_pdb_ligand_map(
        sorted(pdb_targets.keys()),
        batch_size=max(1, args.batch_size),
        timeout=max(10, args.http_timeout),
        retries=max(0, args.http_retries),
        sleep_seconds=max(0.0, args.http_sleep),
    )

    header = [
        "inchikey",
        "uniprot_id",
        "pdb_id",
        "ligand_ccd_id",
        "pdb_entity_id",
        "experimental_method",
        "resolution",
        "release_date",
        "target_mapping_source",
        "ligand_mapping_source",
        "source_type",
        "evidence_scope",
        "xref_version",
        "fetch_date",
        "source_version",
    ]

    out_rows: List[Dict[str, str]] = []
    row_seen: Set[Tuple[str, str, str, str, str]] = set()
    counts = Counter()

    for pdb_id, targets in pdb_targets.items():
        ligands = ligand_by_pdb.get(pdb_id, [])
        if not ligands:
            counts["pdb_without_ligand"] += 1
            continue

        for lig in ligands:
            counts["ligand_candidates_total"] += 1
            if lig.inchikey not in xref_iks:
                counts["ligand_not_in_xref"] += 1
                continue

            for t in targets:
                # keep quality gate in build stage to make contract & QA deterministic
                if not t.experimental_method:
                    counts["skip_empty_experimental_method"] += 1
                    continue
                if not t.resolution:
                    counts["skip_empty_resolution"] += 1
                    continue

                key = (lig.inchikey, t.uniprot_id, pdb_id, lig.ligand_ccd_id, lig.entity_id)
                if key in row_seen:
                    counts["skip_duplicate_output_row"] += 1
                    continue
                row_seen.add(key)

                out_rows.append(
                    {
                        "inchikey": lig.inchikey,
                        "uniprot_id": t.uniprot_id,
                        "pdb_id": pdb_id,
                        "ligand_ccd_id": lig.ligand_ccd_id,
                        "pdb_entity_id": lig.entity_id,
                        "experimental_method": t.experimental_method,
                        "resolution": t.resolution,
                        "release_date": t.release_date,
                        "target_mapping_source": "RCSB_DATA_API:pdb_structures_v1",
                        "ligand_mapping_source": "RCSB_GRAPHQL:entries.nonpolymer_entities",
                        "source_type": "experimental",
                        "evidence_scope": "pdb_ligand_target_complex",
                        "xref_version": "v2" if "_v2" in args.xref.name else "v1",
                        "fetch_date": fetch_date,
                        "source_version": "RCSB:GraphQL_nonpolymer_entities;RCSB:pdb_structures_v1",
                    }
                )
                counts["rows_written_candidate"] += 1

    rows_written = write_tsv(args.out, out_rows, header)

    source_type_distribution = Counter(r["source_type"] for r in out_rows)

    coverage = {
        "name": "molecule_3d_experimental_v1.coverage",
        "created_at": created_at,
        "rows": rows_written,
        "unique_inchikey": len({r["inchikey"] for r in out_rows}),
        "unique_pdb_id": len({r["pdb_id"] for r in out_rows}),
        "unique_uniprot_id": len({r["uniprot_id"] for r in out_rows}),
        "non_empty_rates": {
            "pdb_id": non_empty_rate(out_rows, "pdb_id"),
            "experimental_method": non_empty_rate(out_rows, "experimental_method"),
            "resolution": non_empty_rate(out_rows, "resolution"),
            "target_mapping_source": non_empty_rate(out_rows, "target_mapping_source"),
            "ligand_mapping_source": non_empty_rate(out_rows, "ligand_mapping_source"),
            "source_type": non_empty_rate(out_rows, "source_type"),
            "source_version": non_empty_rate(out_rows, "source_version"),
        },
        "backlink": {
            "in_xref_rows": sum(1 for r in out_rows if r["inchikey"] in xref_iks),
            "total_rows": rows_written,
            "rate": (sum(1 for r in out_rows if r["inchikey"] in xref_iks) / rows_written) if rows_written else 0.0,
        },
        "source_type_distribution": dict(source_type_distribution),
    }

    args.coverage_report.parent.mkdir(parents=True, exist_ok=True)
    args.coverage_report.write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    build_report = {
        "name": "molecule_3d_experimental_v1.build",
        "created_at": created_at,
        "sample_mode": args.max_rows is not None or args.max_pdb is not None,
        "max_rows": args.max_rows,
        "max_pdb": args.max_pdb,
        "inputs": {
            "xref": str(args.xref),
            "protein_pdb": str(args.protein_pdb),
        },
        "output": str(args.out),
        "metrics": {
            "rows_written": rows_written,
            **xref_stats,
            **protein_stats,
            **{k: int(v) for k, v in counts.items()},
            "source_type_distribution": dict(source_type_distribution),
        },
        "rcsb_graphql": graphql_report,
        "coverage_report": str(args.coverage_report),
    }

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(build_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"[OK] build -> {args.out} (rows={rows_written})")
    print(
        "[OK] mapped experimental rows={} unique_inchikey={} unique_pdb={}".format(
            rows_written,
            coverage["unique_inchikey"],
            coverage["unique_pdb_id"],
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
