#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import re
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
print("A5: 重新生成 miRNA 表（修复 TSV 格式）")
print("=" * 70)

INPUT_FILE = Path("data/processed/rna/mirna_with_gene.tsv")
OUTPUT_FILE = Path("data/output/rna_master_mirna_v1.tsv")

if not INPUT_FILE.exists():
    print(f"[ERROR] 找不到: {INPUT_FILE}")
    sys.exit(1)

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

# 读取（用正则处理空白）
records = []
with open(INPUT_FILE, 'r') as f:
    header = f.readline()
    for line in f:
        parts = re.split(r'\s+', line.strip())
        if len(parts) >= 6:
            records.append({
                'urs': parts[0],
                'mirbase_id': parts[1],
                'mimat_id': parts[2],
                'sequence': parts[3],
                'sequence_len': parts[4],
                'gene_symbol': parts[5],
                'hgnc_id': parts[6] if len(parts) > 6 else '',
                'ensembl_gene_id': parts[7] if len(parts) > 7 else '',
                'ncbi_gene_id': parts[8] if len(parts) > 8 else ''
            })

print(f"[OK] 读取 {len(records)} 条记录")

# 标准化
FETCH_DATE = date.today().strftime('%Y-%m-%d')
standardized = []

for rec in records:
    standardized.append({
        'rna_id': rec['urs'],
        'rna_type': 'mirna',
        'rna_name': rec['mirbase_id'],
        'sequence': rec['sequence'],
        'sequence_len': rec['sequence_len'],
        'taxon_id': '9606',
        'symbol': rec['gene_symbol'],
        'source': 'RNAcentral;miRBase',
        'fetch_date': FETCH_DATE,
        'source_version': 'RNAcentral:25;miRBase:22.1',
        'hgnc_id': rec['hgnc_id'],
        'ensembl_gene_id': rec['ensembl_gene_id'],
        'ncbi_gene_id': rec['ncbi_gene_id'],
        'mirbase_id': rec['mimat_id'],
        'ensembl_transcript_id': '',
        'rfam_id': '',
        'secondary_structure': '',
        'pdb_ids': ''
    })

# 写入（确保用 tab 分隔）
FIELDS = ['rna_id','rna_type','rna_name','sequence','sequence_len','taxon_id',
          'symbol','source','fetch_date','source_version','hgnc_id',
          'ensembl_gene_id','ncbi_gene_id','mirbase_id','ensembl_transcript_id',
          'rfam_id','secondary_structure','pdb_ids']

with open(OUTPUT_FILE, 'w') as f:
    f.write('\t'.join(FIELDS) + '\n')
    for rec in sorted(standardized, key=lambda x: x['rna_id']):
        # 用 list comprehension 确保顺序
        values = [str(rec[k]) if rec[k] else '' for k in FIELDS]
        f.write('\t'.join(values) + '\n')

print(f"[OK] 已写入 {len(standardized)} 条")
print(f"✅ 输出: {OUTPUT_FILE}")