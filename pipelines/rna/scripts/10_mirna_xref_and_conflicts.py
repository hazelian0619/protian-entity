#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
from pathlib import Path
from datetime import date
from collections import defaultdict

INP = Path("data/processed/rna/mirna_with_gene.tsv")
OUT_XREF = Path("data/output/rna_xref_mirna_v1.tsv")

OUT_MB2URS = Path("data/output/mirna_conflict_mirbase_to_urs.tsv")
OUT_URS2MB = Path("data/output/mirna_conflict_urs_to_mirbase.tsv")
OUT_MIMAT2URS = Path("data/output/mirna_conflict_mimat_to_urs.tsv")
OUT_URS2MIMAT = Path("data/output/mirna_conflict_urs_to_mimat.tsv")

FETCH_DATE = date.today().strftime("%Y-%m-%d")
SOURCE = "RNAcentral;miRBase"
SOURCE_VERSION = "RNAcentral:25;miRBase:22.1"
TAXON_ID = "9606"

if not INP.exists():
    raise SystemExit(f"[ERROR] missing {INP}")

def split_line(line: str):
    # 兼容 tab/空格混合
    return line.rstrip("\n").split("\t") if "\t" in line else re.split(r"\s+", line.strip())

with open(INP, "r") as f:
    header = split_line(f.readline())
    # 兼容你这个文件可能只有前6列
    # 预期列名：urs mirbase_id mimat_id sequence sequence_len gene_symbol ...
    idx = {name: header.index(name) for name in header}
    need = ["urs", "mirbase_id", "mimat_id"]
    for k in need:
        if k not in idx:
            raise SystemExit(f"[ERROR] column not found: {k}, header={header}")

    rows = []
    for line in f:
        if not line.strip():
            continue
        parts = split_line(line)
        if len(parts) < len(header):
            continue
        urs = parts[idx["urs"]]
        mirbase_id = parts[idx["mirbase_id"]]
        mimat_id = parts[idx["mimat_id"]]
        rows.append((urs, mirbase_id, mimat_id))

# 1) 输出 xref（一个 URS 对应两条 xref：hsa-miR-xxx 和 MIMATxxxx）
OUT_XREF.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_XREF, "w") as g:
    g.write("\t".join([
        "rna_id","rna_type","taxon_id",
        "xref_db","xref_id","xref_level",
        "source","fetch_date","source_version"
    ]) + "\n")

    for urs, mirbase_id, mimat_id in rows:
        g.write("\t".join([urs, "mirna", TAXON_ID, "miRBase", mirbase_id, "mirbase_id",
                           SOURCE, FETCH_DATE, SOURCE_VERSION]) + "\n")
        g.write("\t".join([urs, "mirna", TAXON_ID, "miRBase", mimat_id, "mimat_id",
                           SOURCE, FETCH_DATE, SOURCE_VERSION]) + "\n")

# 2) 冲突审计：统计多对多
mirbase2urs = defaultdict(set)
urs2mirbase = defaultdict(set)
mimat2urs = defaultdict(set)
urs2mimat = defaultdict(set)

for urs, mirbase_id, mimat_id in rows:
    if mirbase_id:
        mirbase2urs[mirbase_id].add(urs)
        urs2mirbase[urs].add(mirbase_id)
    if mimat_id:
        mimat2urs[mimat_id].add(urs)
        urs2mimat[urs].add(mimat_id)

def write_conflict(outp: Path, left2right: dict, left_name: str, right_name: str):
    with open(outp, "w") as g:
        g.write("\t".join([left_name, f"{right_name}_count", f"{right_name}_list"]) + "\n")
        for left, rights in sorted(left2right.items(), key=lambda x: (len(x[1]), x[0]), reverse=True):
            if len(rights) > 1:
                g.write("\t".join([left, str(len(rights)), ";".join(sorted(rights))]) + "\n")

write_conflict(OUT_MB2URS, mirbase2urs, "mirbase_id", "urs")
write_conflict(OUT_URS2MB, urs2mirbase, "urs", "mirbase_id")
write_conflict(OUT_MIMAT2URS, mimat2urs, "mimat_id", "urs")
write_conflict(OUT_URS2MIMAT, urs2mimat, "urs", "mimat_id")

print(f"[OK] xref: {OUT_XREF}")
print(f"[OK] conflicts: {OUT_MB2URS}, {OUT_URS2MB}, {OUT_MIMAT2URS}, {OUT_URS2MIMAT}")