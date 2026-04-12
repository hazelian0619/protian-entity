#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
CHEMBL_RE = re.compile(r"^CHEMBL\d+$")
CHEBI_RE = re.compile(r"CHEBI:\d+")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def utc_today() -> str:
    return datetime.now(timezone.utc).date().isoformat()


def normalize(x: str) -> str:
    return (x or "").strip()


def split_multi(x: str) -> List[str]:
    t = normalize(x)
    if not t:
        return []
    return [s.strip() for s in t.split(";") if s.strip()]


def valid_inchikey(x: str) -> bool:
    return bool(INCHIKEY_RE.match(normalize(x).upper()))


def valid_chembl(x: str) -> bool:
    return bool(CHEMBL_RE.match(normalize(x).upper()))


def join_values(values: Iterable[str]) -> str:
    uniq = sorted({normalize(v) for v in values if normalize(v)})
    return ";".join(uniq)


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
class XrefRec:
    inchikey: str
    chembl_ids: Set[str]
    drugbank_ids: Set[str]


def load_xref(path: Path, max_rows: Optional[int]) -> Tuple[List[XrefRec], Dict[str, Any]]:
    cols, rows = read_tsv(path, max_rows=max_rows)
    if "inchikey" not in cols:
        raise SystemExit(f"[ERROR] xref missing inchikey: {path}")

    out: List[XrefRec] = []
    bad_ik = 0
    for row in rows:
        ik = normalize(row.get("inchikey", "")).upper()
        if not valid_inchikey(ik):
            bad_ik += 1
            continue

        chembl_ids = {c.upper() for c in split_multi(row.get("chembl_id", "")) if valid_chembl(c.upper())}
        drugbank_ids = {d for d in split_multi(row.get("drugbank_id", "")) if d}

        out.append(XrefRec(inchikey=ik, chembl_ids=chembl_ids, drugbank_ids=drugbank_ids))

    stats = {
        "input_rows": len(rows),
        "valid_rows": len(out),
        "skipped_bad_inchikey_rows": bad_ik,
        "rows_with_chembl": sum(1 for r in out if r.chembl_ids),
        "rows_with_drugbank": sum(1 for r in out if r.drugbank_ids),
    }
    return out, stats


def load_drugbank_groups(drug_master_path: Path) -> Tuple[Dict[str, Set[str]], Dict[str, Any]]:
    if not drug_master_path.exists():
        raise FileNotFoundError(str(drug_master_path))

    _, rows = read_tsv(drug_master_path)
    dbid_to_groups: Dict[str, Set[str]] = defaultdict(set)
    versions: Set[str] = set()

    for row in rows:
        dbid = normalize(row.get("drugbank_id", ""))
        if not dbid:
            continue
        for g in split_multi(row.get("groups", "")):
            dbid_to_groups[dbid].add(g.lower())

        sv = normalize(row.get("source_version", ""))
        if sv:
            versions.add(sv)

    report = {
        "drugbank_ids": len(dbid_to_groups),
        "version": join_values(f"DrugBank:{v}" for v in versions) or "DrugBank:unknown",
    }
    return dbid_to_groups, report


def load_zinc_from_registry(registry_path: Optional[Path]) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Any]]:
    if registry_path is None or not registry_path.exists():
        return {}, {"available": False, "warning": "molecule_3d_registry_v1.tsv not found"}

    _, rows = read_tsv(registry_path)
    out: Dict[str, Dict[str, str]] = {}
    source_versions: Set[str] = set()

    for row in rows:
        ik = normalize(row.get("inchikey", "")).upper()
        if not valid_inchikey(ik):
            continue
        out[ik] = {
            "zinc_id": normalize(row.get("zinc_id", "")),
            "zinc_3d_available": normalize(row.get("zinc_3d_available", "0")),
            "source_version": normalize(row.get("source_version", "")),
        }
        if out[ik]["source_version"]:
            source_versions.add(out[ik]["source_version"])

    report = {
        "available": True,
        "rows": len(rows),
        "mapped_inchikey": len(out),
        "source_version": join_values(source_versions) or "molecule_3d_registry_v1:unknown",
    }
    return out, report


def parse_chebi_obo(
    obo_path: Path,
    target_inchikeys: Set[str],
) -> Tuple[Dict[str, Set[str]], Dict[str, str], Dict[str, List[str]], Dict[str, Any]]:
    if not obo_path.exists():
        raise FileNotFoundError(str(obo_path))

    inchikey_to_chebi: Dict[str, Set[str]] = defaultdict(set)
    name_by_id: Dict[str, str] = {}
    parents_by_id: Dict[str, List[str]] = {}

    data_version = "unknown"

    cur_id = ""
    cur_name = ""
    cur_parents: List[str] = []
    cur_ik = ""
    cur_obsolete = False
    in_term = False

    counts = Counter()

    def flush_term() -> None:
        nonlocal cur_id, cur_name, cur_parents, cur_ik, cur_obsolete
        if not cur_id:
            return
        counts["term_blocks"] += 1

        if not cur_obsolete:
            name_by_id[cur_id] = cur_name
            parents_by_id[cur_id] = sorted(set(cur_parents))
            if cur_ik and cur_ik in target_inchikeys:
                inchikey_to_chebi[cur_ik].add(cur_id)
                counts["target_inchikey_mapped_terms"] += 1

        cur_id = ""
        cur_name = ""
        cur_parents = []
        cur_ik = ""
        cur_obsolete = False

    with obo_path.open("r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.rstrip("\n")

            if line.startswith("data-version:"):
                data_version = line.split(":", 1)[1].strip()

            if line == "[Term]":
                flush_term()
                in_term = True
                continue

            if line.startswith("[") and line != "[Term]":
                flush_term()
                in_term = False
                continue

            if not in_term:
                continue

            if line.startswith("id:"):
                cur_id = line.split(":", 1)[1].strip()
            elif line.startswith("name:"):
                cur_name = line.split(":", 1)[1].strip()
            elif line.startswith("is_a:"):
                parent = line.split("is_a:", 1)[1].strip().split(" ! ", 1)[0].strip()
                if parent.startswith("CHEBI:"):
                    cur_parents.append(parent)
            elif line.startswith("is_obsolete:") and line.split(":", 1)[1].strip().lower() == "true":
                cur_obsolete = True
            elif line.startswith("property_value:") and "inchi_key_string" in line:
                # property_value: chemrof:inchi_key_string "XGEA..." xsd:string
                parts = line.split('"')
                if len(parts) >= 2:
                    ik = parts[1].strip().upper()
                    if valid_inchikey(ik):
                        cur_ik = ik

    flush_term()

    counts["target_inchikey_hits"] = len(inchikey_to_chebi)
    counts["terms_loaded"] = len(name_by_id)
    counts["parent_nodes"] = len(parents_by_id)

    report = {
        "data_version": data_version,
        **{k: int(v) for k, v in counts.items()},
    }
    return inchikey_to_chebi, name_by_id, parents_by_id, report


def shortest_path_to_root(term_id: str, parents_by_id: Dict[str, List[str]], memo: Dict[str, List[str]], visiting: Set[str]) -> List[str]:
    if term_id in memo:
        return memo[term_id]
    if term_id in visiting:
        # cycle guard
        return [term_id]

    visiting.add(term_id)
    parents = parents_by_id.get(term_id, [])

    if not parents:
        path = [term_id]
    else:
        candidates: List[List[str]] = []
        for p in parents:
            p_path = shortest_path_to_root(p, parents_by_id, memo, visiting)
            candidates.append([term_id] + p_path)
        # shortest path first; tie-break lexicographically for determinism
        candidates.sort(key=lambda x: (len(x), ">".join(x)))
        path = candidates[0]

    visiting.remove(term_id)
    memo[term_id] = path
    return path


def load_atc_map_from_chembl(chembl_db: Path) -> Tuple[Dict[str, Dict[str, Set[str]]], Dict[str, Any]]:
    if not chembl_db.exists():
        raise FileNotFoundError(str(chembl_db))

    map_by_chembl: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))

    with sqlite3.connect(chembl_db) as conn:
        release_row = conn.execute(
            "SELECT chembl_release FROM chembl_release ORDER BY chembl_release_id DESC LIMIT 1"
        ).fetchone()
        chembl_release = release_row[0] if release_row else "CHEMBL_UNKNOWN"

        q = """
            SELECT UPPER(TRIM(md.chembl_id)) AS chembl_id,
                   mac.level5,
                   ac.level1,
                   ac.level1_description,
                   ac.level2,
                   ac.level2_description
            FROM molecule_dictionary md
            JOIN molecule_atc_classification mac ON md.molregno = mac.molregno
            LEFT JOIN atc_classification ac ON ac.level5 = mac.level5
            WHERE md.chembl_id IS NOT NULL
        """

        for chembl_id, level5, level1, level1_desc, level2, level2_desc in conn.execute(q):
            c = normalize(chembl_id).upper()
            if not valid_chembl(c):
                continue
            if normalize(level5):
                map_by_chembl[c]["atc_code"].add(normalize(level5))
            if normalize(level1):
                map_by_chembl[c]["atc_level1_code"].add(normalize(level1))
            if normalize(level1_desc):
                map_by_chembl[c]["atc_level1_desc"].add(normalize(level1_desc))
            if normalize(level2):
                map_by_chembl[c]["atc_level2_code"].add(normalize(level2))
            if normalize(level2_desc):
                map_by_chembl[c]["atc_level2_desc"].add(normalize(level2_desc))

    report = {
        "chembl_release": chembl_release,
        "mapped_chembl_ids": len(map_by_chembl),
        "atc_pairs": int(sum(len(v.get("atc_code", set())) for v in map_by_chembl.values())),
    }
    return map_by_chembl, report


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xref", type=Path, required=True)
    ap.add_argument("--chebi-obo", type=Path, required=True)
    ap.add_argument("--drug-master", type=Path, required=True)
    ap.add_argument("--chembl-db", type=Path, required=True)
    ap.add_argument("--registry-3d", type=Path, default=None)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--hierarchy-report", type=Path, required=True)
    ap.add_argument("--coverage-report", type=Path, required=True)
    ap.add_argument("--max-rows", type=int, default=None)
    args = ap.parse_args()

    created_at = utc_now()
    fetch_date = utc_today()

    for p in [args.xref, args.chebi_obo, args.drug_master, args.chembl_db]:
        if not p.exists():
            raise SystemExit(f"[ERROR] missing required input: {p}")

    xref_rows, xref_stats = load_xref(args.xref, args.max_rows)
    target_inchikeys = {r.inchikey for r in xref_rows}

    # load inputs
    dbid_to_groups, drugbank_report = load_drugbank_groups(args.drug_master)
    zinc_by_ik, zinc_report = load_zinc_from_registry(args.registry_3d)
    atc_by_chembl, atc_report = load_atc_map_from_chembl(args.chembl_db)
    ik_to_chebi, chebi_name, chebi_parents, chebi_parse_report = parse_chebi_obo(args.chebi_obo, target_inchikeys)

    memo_path: Dict[str, List[str]] = {}

    header = [
        "inchikey",
        "chebi_id",
        "chebi_name",
        "chebi_ancestor_ids",
        "chebi_paths",
        "atc_code",
        "atc_level1_code",
        "atc_level1_description",
        "atc_level2_code",
        "atc_level2_description",
        "drugbank_group",
        "zinc_id",
        "zinc_label",
        "semantic_tags",
        "chebi_source",
        "atc_source",
        "drugbank_source",
        "zinc_source",
        "fetch_date",
        "source_version",
    ]

    out_rows: List[Dict[str, str]] = []
    counts = Counter()
    hierarchy_samples: List[Dict[str, Any]] = []

    for r in xref_rows:
        ik = r.inchikey
        counts["rows"] += 1

        chebi_ids = sorted(ik_to_chebi.get(ik, set()))
        chebi_names: List[str] = []
        chebi_ancestors: Set[str] = set()
        chebi_paths: List[str] = []

        for cid in chebi_ids:
            chebi_names.append(chebi_name.get(cid, ""))
            path = shortest_path_to_root(cid, chebi_parents, memo_path, set())
            if len(path) > 1:
                chebi_ancestors.update(path[1:])
            chebi_paths.append(">".join(path))

            if len(hierarchy_samples) < 30:
                hierarchy_samples.append(
                    {
                        "inchikey": ik,
                        "chebi_id": cid,
                        "path": path,
                        "path_length": len(path),
                        "leaf_name": chebi_name.get(cid, ""),
                        "parent_name": chebi_name.get(path[1], "") if len(path) > 1 else "",
                    }
                )

        if chebi_ids:
            counts["rows_with_chebi"] += 1
            if chebi_ancestors:
                counts["rows_with_chebi_expandable"] += 1

        # ATC from all linked ChEMBL IDs
        atc_codes: Set[str] = set()
        atc_l1: Set[str] = set()
        atc_l1_desc: Set[str] = set()
        atc_l2: Set[str] = set()
        atc_l2_desc: Set[str] = set()

        for c in r.chembl_ids:
            m = atc_by_chembl.get(c, {})
            atc_codes.update(m.get("atc_code", set()))
            atc_l1.update(m.get("atc_level1_code", set()))
            atc_l1_desc.update(m.get("atc_level1_desc", set()))
            atc_l2.update(m.get("atc_level2_code", set()))
            atc_l2_desc.update(m.get("atc_level2_desc", set()))

        if atc_codes:
            counts["rows_with_atc"] += 1

        # DrugBank groups
        groups: Set[str] = set()
        for dbid in r.drugbank_ids:
            groups.update(dbid_to_groups.get(dbid, set()))

        if groups:
            counts["rows_with_drugbank_group"] += 1

        # ZINC labels (optional input)
        zinc_id = normalize(zinc_by_ik.get(ik, {}).get("zinc_id", ""))
        zinc_labels: Set[str] = set()
        if zinc_id:
            zinc_labels.add("zinc_registered")
        if normalize(zinc_by_ik.get(ik, {}).get("zinc_3d_available", "0")) == "1":
            zinc_labels.add("zinc_3d_available")

        if zinc_labels:
            counts["rows_with_zinc_label"] += 1

        # semantic tag aggregation
        semantic_tags: Set[str] = set()
        semantic_tags.update({f"chebi:{x}" for x in chebi_ids})
        semantic_tags.update({f"atc:{x}" for x in atc_codes})
        semantic_tags.update({f"group:{x}" for x in groups})
        semantic_tags.update({f"zinc:{x}" for x in zinc_labels})
        if not semantic_tags:
            semantic_tags.add("unclassified")

        source_version_tokens = [
            f"ChEBI:{chebi_parse_report.get('data_version', 'unknown')}",
            f"ChEMBL:{atc_report.get('chembl_release', 'CHEMBL_UNKNOWN')}",
            drugbank_report.get("version", "DrugBank:unknown"),
            f"ZINC:{'molecule_3d_registry_v1' if zinc_report.get('available') else 'unavailable'}",
        ]

        out_rows.append(
            {
                "inchikey": ik,
                "chebi_id": join_values(chebi_ids) or "NA",
                "chebi_name": join_values(chebi_names) or "NA",
                "chebi_ancestor_ids": join_values(chebi_ancestors) or "NA",
                "chebi_paths": "|".join(sorted(set(chebi_paths))) if chebi_paths else "NA",
                "atc_code": join_values(atc_codes) or "NA",
                "atc_level1_code": join_values(atc_l1) or "NA",
                "atc_level1_description": join_values(atc_l1_desc) or "NA",
                "atc_level2_code": join_values(atc_l2) or "NA",
                "atc_level2_description": join_values(atc_l2_desc) or "NA",
                "drugbank_group": join_values(groups) or "NA",
                "zinc_id": zinc_id or "NA",
                "zinc_label": join_values(zinc_labels) or "NA",
                "semantic_tags": join_values(semantic_tags),
                "chebi_source": "ChEBI_obo",
                "atc_source": "ChEMBL_atc_classification",
                "drugbank_source": "DrugBank_groups",
                "zinc_source": "molecule_3d_registry_v1" if zinc_report.get("available") else "NA",
                "fetch_date": fetch_date,
                "source_version": join_values(source_version_tokens),
            }
        )

    rows_written = write_tsv(args.out, out_rows, header)

    hierarchy_report = {
        "name": "molecule_semantic_tags_v1.chebi_hierarchy",
        "created_at": created_at,
        "chebi_parse": chebi_parse_report,
        "rows_with_chebi": int(counts["rows_with_chebi"]),
        "rows_with_chebi_expandable": int(counts["rows_with_chebi_expandable"]),
        "expandable_rate_among_chebi_rows": (
            counts["rows_with_chebi_expandable"] / counts["rows_with_chebi"]
            if counts["rows_with_chebi"]
            else 0.0
        ),
        "samples": hierarchy_samples,
    }

    coverage_report = {
        "name": "molecule_semantic_tags_v1.coverage",
        "created_at": created_at,
        "rows": rows_written,
        "metrics": {
            "rows_with_chebi": int(counts["rows_with_chebi"]),
            "rows_with_atc": int(counts["rows_with_atc"]),
            "rows_with_drugbank_group": int(counts["rows_with_drugbank_group"]),
            "rows_with_zinc_label": int(counts["rows_with_zinc_label"]),
            "rows_with_chebi_expandable": int(counts["rows_with_chebi_expandable"]),
            "chebi_rate": (counts["rows_with_chebi"] / rows_written) if rows_written else 0.0,
            "atc_rate": (counts["rows_with_atc"] / rows_written) if rows_written else 0.0,
            "drugbank_group_rate": (counts["rows_with_drugbank_group"] / rows_written) if rows_written else 0.0,
            "zinc_label_rate": (counts["rows_with_zinc_label"] / rows_written) if rows_written else 0.0,
        },
        "non_empty_rates": {
            col: {
                "rate": (sum(1 for row in out_rows if normalize(row.get(col, "")) != "") / rows_written)
                if rows_written
                else 0.0
            }
            for col in [
                "chebi_source",
                "atc_source",
                "drugbank_source",
                "zinc_source",
                "source_version",
            ]
        },
        "sources": {
            "chebi": f"ChEBI_obo:{chebi_parse_report.get('data_version', 'unknown')}",
            "atc": f"ChEMBL:{atc_report.get('chembl_release', 'CHEMBL_UNKNOWN')}",
            "drugbank": drugbank_report.get("version", "DrugBank:unknown"),
            "zinc": "molecule_3d_registry_v1" if zinc_report.get("available") else "NA",
        },
    }

    build_report = {
        "name": "molecule_semantic_tags_v1.build",
        "created_at": created_at,
        "sample_mode": args.max_rows is not None,
        "max_rows": args.max_rows,
        "inputs": {
            "xref": str(args.xref),
            "chebi_obo": str(args.chebi_obo),
            "drug_master": str(args.drug_master),
            "chembl_db": str(args.chembl_db),
            "registry_3d": str(args.registry_3d) if args.registry_3d else None,
            "registry_3d_exists": bool(args.registry_3d and args.registry_3d.exists()),
        },
        "output": str(args.out),
        "metrics": {
            "rows_written": rows_written,
            **xref_stats,
            **coverage_report["metrics"],
        },
        "chebi_parse": chebi_parse_report,
        "atc": atc_report,
        "drugbank": drugbank_report,
        "zinc": zinc_report,
        "coverage_report": str(args.coverage_report),
        "hierarchy_report": str(args.hierarchy_report),
    }

    args.hierarchy_report.parent.mkdir(parents=True, exist_ok=True)
    args.hierarchy_report.write_text(json.dumps(hierarchy_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    args.coverage_report.parent.mkdir(parents=True, exist_ok=True)
    args.coverage_report.write_text(json.dumps(coverage_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(build_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"[OK] build -> {args.out} (rows={rows_written})")
    print(
        "[OK] mapped rates: chebi={:.4f} atc={:.4f} drugbank_group={:.4f} zinc_label={:.4f}".format(
            coverage_report["metrics"]["chebi_rate"],
            coverage_report["metrics"]["atc_rate"],
            coverage_report["metrics"]["drugbank_group_rate"],
            coverage_report["metrics"]["zinc_label_rate"],
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
