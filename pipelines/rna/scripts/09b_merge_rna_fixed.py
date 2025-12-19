#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import sys

MIRNA_FILE = Path("data/output/rna_master_mirna_v1.fixed.tsv")
MRNA_FILE  = Path("data/output/rna_master_mrna_v1.fixed.tsv")
OUT_FILE   = Path("data/output/rna_master_v1.tsv")

for p in [MIRNA_FILE, MRNA_FILE]:
    if not p.exists():
        print(f"[ERROR] missing: {p}")
        sys.exit(1)

def read_tsv(p: Path):
    with open(p, 'r') as f:
        header = f.readline().rstrip('\n').split('\t')
        rows = [line.rstrip('\n').split('\t') for line in f if line.strip()]
    return header, rows

h1, r1 = read_tsv(MIRNA_FILE)
h2, r2 = read_tsv(MRNA_FILE)

print(f"[OK] miRNA: {len(r1)}")
print(f"[OK] mRNA : {len(r2)}")

if h1 != h2:
    print("[ERROR] header mismatch")
    print(h1)
    print(h2)
    sys.exit(1)

all_rows = r1 + r2
rna_ids = [x[0] for x in all_rows]
if len(set(rna_ids)) != len(rna_ids):
    print("[ERROR] duplicated rna_id found")
    sys.exit(1)

with open(OUT_FILE, 'w') as g:
    g.write('\t'.join(h1) + '\n')
    for row in sorted(all_rows, key=lambda x: x[0]):
        g.write('\t'.join(row) + '\n')

print(f"[OK] merged -> {OUT_FILE}")
print(f"[OK] total rows: {len(all_rows)}")