#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
终极版：从最早的正确数据重新构建 RNA Master v1
"""

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
print("=" * 70)
print("终极版：重建 RNA Master v1")
print("=" * 70)

SEED_GENES = Path("data/processed/rna/seed_genes.tsv")
MIRNA_URS = Path("data/processed/rna/mirna_urs_mapping.tsv")
MRNA_SEQ = Path("data/processed/rna/mrna_with_sequence.tsv")
OUTPUT = Path("data/output/rna_master_v1.tsv")

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

FETCH_DATE = date.today().strftime('%Y-%m-%d')

# ============================================================
# 第1步：读取 seed_genes（用于补充 gene 信息）
# ============================================================
print("[STEP 1] 读取 seed_genes...")

seed_genes = {}
with open(SEED_GENES, 'r') as f:
    f.readline()  # 跳过表头
    for line in f:
        parts = re.split(r'\s+', line.strip())
        if len(parts) >= 4:
            ensembl_gene_id = parts[2].split('.')[0]
            seed_genes[ensembl_gene_id] = {
                'symbol': parts[0],
                'hgnc_id': parts[1],
                'ncbi_gene_id': parts[3]
            }

print(f"[OK] seed_genes: {len(seed_genes)}\n")

# ============================================================
# 第2步：处理 miRNA
# ============================================================
print("[STEP 2] 处理 miRNA...")

mirna_records = []
with open(MIRNA_URS, 'r') as f:
    f.readline()  # 跳过表头
    for line in f:
        parts = re.split(r'\s+', line.strip())
        if len(parts) >= 5:
            mirna_records.append({
                'rna_id': parts[0],
                'rna_type': 'mirna',
                'rna_name': parts[1],
                'sequence': parts[3],
                'sequence_len': parts[4],
                'taxon_id': '9606',
                'symbol': parts[5] if len(parts) > 5 else '',
                'source': 'RNAcentral;miRBase',
                'fetch_date': FETCH_DATE,
                'source_version': 'RNAcentral:25;miRBase:22.1',
                'hgnc_id': '',
                'ensembl_gene_id': '',
                'ncbi_gene_id': '',
                'mirbase_id': parts[2],
                'ensembl_transcript_id': '',
                'rfam_id': '',
                'secondary_structure': '',
                'pdb_ids': ''
            })

print(f"[OK] miRNA: {len(mirna_records)}\n")

# ============================================================
# 第3步：处理 mRNA
# ============================================================
print("[STEP 3] 处理 mRNA...")

mrna_records = []
with open(MRNA_SEQ, 'r') as f:
    header = f.readline().strip()
    
    for line in f:
        # mrna_with_sequence.tsv 格式：
        # rna_id, ensembl_transcript_id, ensembl_gene_id, gene_symbol, hgnc_id, ncbi_gene_id, sequence, sequence_len
        parts = line.strip().split('\t')
        
        if len(parts) >= 8:
            # 标准 tab 分隔
            rna_id = parts[0]
            transcript_id = parts[1]
            gene_id = parts[2]
            symbol = parts[3]
            hgnc = parts[4]
            ncbi = parts[5]
            sequence = parts[6]
            seq_len = parts[7]
        else:
            # 可能有粘连，用正则+启发式修复
            parts = re.split(r'\s+', line.strip(), maxsplit=7)
            if len(parts) < 8:
                continue
            
            rna_id = parts[0]
            transcript_id = parts[1]
            gene_id = parts[2]
            symbol = parts[3]
            hgnc = parts[4]
            
            # ncbi_gene_id 和 sequence 可能粘在一起
            rest = parts[5] + parts[6] if len(parts) > 6 else parts[5]
            
            # 尝试分离：找到第一个字母（序列开始）
            seq_start = 0
            for i, c in enumerate(rest):
                if c in 'ACGUN':
                    seq_start = i
                    break
            
            if seq_start > 0:
                ncbi = rest[:seq_start]
                sequence = rest[seq_start:]
            else:
                ncbi = ''
                sequence = rest
            
            seq_len = str(len(sequence))
        
        mrna_records.append({
            'rna_id': rna_id,
            'rna_type': 'mRNA',
            'rna_name': transcript_id,
            'sequence': sequence,
            'sequence_len': seq_len,
            'taxon_id': '9606',
            'symbol': symbol,
            'source': 'RNAcentral;Ensembl',
            'fetch_date': FETCH_DATE,
            'source_version': 'RNAcentral:25;Ensembl:113',
            'hgnc_id': hgnc,
            'ensembl_gene_id': gene_id,
            'ncbi_gene_id': ncbi,
            'mirbase_id': '',
            'ensembl_transcript_id': transcript_id,
            'rfam_id': '',
            'secondary_structure': '',
            'pdb_ids': ''
        })

print(f"[OK] mRNA: {len(mrna_records)}\n")

# ============================================================
# 第4步：合并 + 写入
# ============================================================
print("[STEP 4] 合并并写入...")

all_records = mirna_records + mrna_records

FIELDS = ['rna_id','rna_type','rna_name','sequence','sequence_len','taxon_id',
          'symbol','source','fetch_date','source_version','hgnc_id',
          'ensembl_gene_id','ncbi_gene_id','mirbase_id','ensembl_transcript_id',
          'rfam_id','secondary_structure','pdb_ids']

with open(OUTPUT, 'w') as f:
    f.write('\t'.join(FIELDS) + '\n')
    for rec in sorted(all_records, key=lambda x: x['rna_id']):
        values = [str(rec[k]) if rec[k] else '' for k in FIELDS]
        f.write('\t'.join(values) + '\n')

print(f"[OK] 已写入 {len(all_records):,} 条\n")

# ============================================================
# 第5步：QA
# ============================================================
print("=" * 70)
print("QA")
print("=" * 70)

rna_types = {}
for rec in all_records:
    rt = rec['rna_type']
    rna_types[rt] = rna_types.get(rt, 0) + 1

print("\n[统计] RNA 类型:")
for rt in sorted(rna_types.keys()):
    count = rna_types[rt]
    pct = count / len(all_records) * 100
    print(f"  {rt:10s}: {count:7,d} ({pct:5.1f}%)")

lengths = [int(rec['sequence_len']) for rec in all_records if rec['sequence_len'].isdigit()]
print(f"\n[统计] 序列长度:")
print(f"  最短: {min(lengths):,} nt")
print(f"  最长: {max(lengths):,} nt")
print(f"  平均: {sum(lengths)//len(lengths):,} nt")

print("\n" + "=" * 70)
print("🎉🎉🎉 RNA Master v1 完成！")
print("=" * 70)
print(f"✅ 输出: {OUTPUT}")
print(f"✅ 总记录: {len(all_records):,}")
print(f"✅   miRNA: {len(mirna_records):,}")
print(f"✅   mRNA: {len(mrna_records):,}")
print("\n🚀 下一步: 小分子实体")