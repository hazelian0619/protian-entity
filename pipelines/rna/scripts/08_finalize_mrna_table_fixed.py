#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
from pathlib import Path
from datetime import date

def _find_project_root(start: Path) -> Path:
    # Prefer explicit override for reproducible runs
    env_root = os.environ.get("KG_ROOT")
    if env_root:
        return Path(env_root).expanduser().resolve()

    # Otherwise, walk upwards until we find the repo root
    for parent in [start] + list(start.parents):
        if (parent / ".git").exists() and (parent / "data").exists():
            return parent

    # Fallback: keep previous behavior as last resort
    return start.parent.parent


project_root = _find_project_root(Path(__file__).resolve())

# 切换到项目根目录（确保相对路径有效）
os.chdir(project_root)
print("B5: 重新生成 mRNA 表（修复 TSV 格式）")
print("=" * 70)

INPUT_FILE = Path("data/processed/rna/mrna_with_sequence.tsv")
OUTPUT_FILE = Path("data/output/rna_master_mrna_v1.tsv")

if not INPUT_FILE.exists():
    print(f"[ERROR] 找不到: {INPUT_FILE}")
    sys.exit(1)

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

# 读取
records = []
with open(INPUT_FILE, 'r') as f:
    header = f.readline().strip().split('\t')
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 8:
            records.append({
                'rna_id': parts[0],
                'ensembl_transcript_id': parts[1],
                'ensembl_gene_id': parts[2],
                'gene_symbol': parts[3],
                'hgnc_id': parts[4],
                'ncbi_gene_id': parts[5],
                'sequence': parts[6],
                'sequence_len': parts[7]
            })

print(f"[OK] 读取 {len(records)} 条记录")

# 标准化
FETCH_DATE = date.today().strftime('%Y-%m-%d')
standardized = []

for rec in records:
    standardized.append({
        'rna_id': rec['rna_id'],
        'rna_type': 'mrna',
        'rna_name': rec['ensembl_transcript_id'],
        'sequence': rec['sequence'],
        'sequence_len': rec['sequence_len'],
        'taxon_id': '9606',
        'symbol': rec['gene_symbol'],
        'source': 'RNAcentral;Ensembl',
        'fetch_date': FETCH_DATE,
        'source_version': 'RNAcentral:25;Ensembl:113',
        'hgnc_id': rec['hgnc_id'],
        'ensembl_gene_id': rec['ensembl_gene_id'],
        'ncbi_gene_id': rec['ncbi_gene_id'],
        'mirbase_id': '',
        'ensembl_transcript_id': rec['ensembl_transcript_id'],
        'rfam_id': '',
        'secondary_structure': '',
        'pdb_ids': ''
    })

# 写入
FIELDS = ['rna_id','rna_type','rna_name','sequence','sequence_len','taxon_id',
          'symbol','source','fetch_date','source_version','hgnc_id',
          'ensembl_gene_id','ncbi_gene_id','mirbase_id','ensembl_transcript_id',
          'rfam_id','secondary_structure','pdb_ids']

with open(OUTPUT_FILE, 'w') as f:
    f.write('\t'.join(FIELDS) + '\n')
    for rec in sorted(standardized, key=lambda x: x['rna_id']):
        values = [str(rec[k]) if rec[k] else '' for k in FIELDS]
        f.write('\t'.join(values) + '\n')

print(f"[OK] 已写入 {len(standardized)} 条")
print(f"✅ 输出: {OUTPUT_FILE}")