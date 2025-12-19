#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
from pathlib import Path
from datetime import date

IDMAP = Path("data/raw/rna/rnacentral/id_mapping.tsv")  # 也支持 .gz
MRNA_MASTER = Path("data/output/rna_master_mrna_v1.fixed.tsv")
OUT_XREF = Path("data/output/rna_xref_mrna_enst_urs_v1.tsv")
OUT_REP = Path("data/output/mrna_enst_urs_coverage.txt")

FETCH_DATE = date.today().strftime("%Y-%m-%d")

if not IDMAP.exists():
    raise SystemExit(f"[ERROR] missing {IDMAP}")
if not MRNA_MASTER.exists():
    raise SystemExit(f"[ERROR] missing {MRNA_MASTER}")

def open_maybe_gz(p: Path):
    if str(p).endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p, "r")

# 1) 读取 mRNA master 里的 ENST_9606 集合
enst_set = set()
with open(MRNA_MASTER, "r") as f:
    header = f.readline().rstrip("\n").split("\t")
    rna_id_idx = header.index("rna_id")
    rna_type_idx = header.index("rna_type")

    for line in f:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if parts[rna_type_idx].lower() != "mrna":
            continue
        rid = parts[rna_id_idx]
        if rid.startswith("ENST") and rid.endswith("_9606"):
            enst_set.add(rid)

print(f"[OK] ENST in mRNA master: {len(enst_set):,}")

# 2) 扫描 id_mapping.tsv（无表头，固定列位）
# 观察到的格式：URS | database | external_id | taxid | ...（后面可忽略）
pairs = {}      # ENST_9606 -> URS_9606（若多对一，保留第一个）
multi_map = 0
scanned = 0
kept_tax = 0
kept_db = 0
kept_enst = 0
kept_in_master = 0

with open_maybe_gz(IDMAP) as f:
    for line in f:
        if not line.strip():
            continue
        scanned += 1
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue

        urs_raw = parts[0]
        db = parts[1].upper()
        ext = parts[2]
        tax = parts[3]

        if tax != "9606":
            continue
        kept_tax += 1

        if db not in ("ENSEMBL", "ENSEMBL_GENCODE", "GENCODE"):
            continue
        kept_db += 1

        if not ext.startswith("ENST"):
            continue
        kept_enst += 1

        ext_nover = ext.split(".")[0]
        enst_id = f"{ext_nover}_9606"
        if enst_id not in enst_set:
            continue
        kept_in_master += 1

        urs = urs_raw if urs_raw.endswith("_9606") else f"{urs_raw}_9606"
        if enst_id in pairs and pairs[enst_id] != urs:
            multi_map += 1
        else:
            pairs[enst_id] = urs

mapped = len(pairs)
total = len(enst_set)
pct = (mapped / total * 100.0) if total else 0.0

# 3) 输出 xref
OUT_XREF.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_XREF, "w") as g:
    g.write("\t".join([
        "rna_id", "rna_type", "taxon_id",
        "xref_db", "xref_id", "xref_level",
        "source", "fetch_date", "source_version"
    ]) + "\n")
    for enst_id, urs in sorted(pairs.items()):
        g.write("\t".join([
            enst_id, "mrna", "9606",
            "RNAcentral", urs, "urs_id",
            "RNAcentral", FETCH_DATE, "RNAcentral:25"
        ]) + "\n")

# 4) coverage 报告
with open(OUT_REP, "w") as g:
    g.write("mRNA ENST -> URS coverage report\n")
    g.write("="*60 + "\n")
    g.write(f"mRNA ENST total: {total}\n")
    g.write(f"Mapped to URS:  {mapped}\n")
    g.write(f"Coverage:      {pct:.2f}%\n")
    g.write(f"ENST mapped to multiple URS (observed): {multi_map}\n\n")
    g.write(f"Scan stats:\n")
    g.write(f"  lines_scanned: {scanned}\n")
    g.write(f"  pass_taxid_9606: {kept_tax}\n")
    g.write(f"  pass_db_ensembl: {kept_db}\n")
    g.write(f"  pass_ext_ENST: {kept_enst}\n")
    g.write(f"  pass_in_mRNA_master: {kept_in_master}\n\n")
    g.write(f"id_mapping: {IDMAP}\n")
    g.write(f"mRNA master: {MRNA_MASTER}\n")
    g.write(f"xref out: {OUT_XREF}\n")

print(f"[OK] xref written: {OUT_XREF} (rows={mapped:,})")
print(f"[OK] report: {OUT_REP} (coverage={pct:.2f}%)")
print(f"[WARN] multi-map events observed: {multi_map}")