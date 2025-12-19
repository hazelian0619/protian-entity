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
print("A5: 生成标准化 miRNA 子表")
print("=" * 70)

INPUT_FILE = Path("data/processed/rna/mirna_with_gene.tsv")
OUTPUT_FILE = Path("data/output/rna_master_mirna_v1.tsv")

if not INPUT_FILE.exists():
    print(f"[ERROR] 找不到: {INPUT_FILE}")
    sys.exit(1)

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

# 先诊断文件
print("\n[诊断] 检查输入文件...")
with open(INPUT_FILE, 'r') as f:
    first_line = f.readline()
    print(f"  表头: {first_line.strip()}")
    print(f"  列数: {len(first_line.strip().split(chr(9)))}")
    
    second_line = f.readline()
    if second_line:
        print(f"  第1行数据: {second_line.strip()[:100]}...")
        print(f"  第1行列数: {len(second_line.strip().split(chr(9)))}")

# 重新读取
records = []
with open(INPUT_FILE, 'r') as f:
    header = f.readline()
    line_num = 0
    for line in f:
        line_num += 1
        parts = line.strip().split('\t')
        
        # 调试：如果前10行解析失败，显示详情
        if line_num <= 10 and len(parts) < 9:
            print(f"  [WARN] 第{line_num}行列数不足: {len(parts)} 列")
            print(f"         内容: {line.strip()[:80]}")
        
        if len(parts) >= 9:
            records.append({
                'urs': parts[0],
                'mirbase_id': parts[1],
                'mimat_id': parts[2],
                'sequence': parts[3],
                'sequence_len': parts[4],
                'gene_symbol': parts[5] if len(parts) > 5 else '',
                'hgnc_id': parts[6] if len(parts) > 6 else '',
                'ensembl_gene_id': parts[7] if len(parts) > 7 else '',
                'ncbi_gene_id': parts[8] if len(parts) > 8 else ''
            })

print(f"\n[OK] 成功读取 {len(records)} 条记录")

if len(records) == 0:
    print("[ERROR] 没有读取到任何记录，请检查文件格式")
    sys.exit(1)

# 继续标准化...
FETCH_DATE = date.today().strftime('%Y-%m-%d')
SOURCE = "RNAcentral;miRBase"
SOURCE_VERSION = "RNAcentral:25;miRBase:22.1"

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
        'source': SOURCE,
        'fetch_date': FETCH_DATE,
        'source_version': SOURCE_VERSION,
        'hgnc_id': rec['hgnc_id'],
        'ensembl_gene_id': rec['ensembl_gene_id'],
        'ncbi_gene_id': rec['ncbi_gene_id'],
        'mirbase_id': rec['mimat_id'],
        'ensembl_transcript_id': '',
        'rfam_id': '',
        'secondary_structure': '',
        'pdb_ids': ''
    })

FIELDS = ['rna_id','rna_type','rna_name','sequence','sequence_len','taxon_id',
          'symbol','source','fetch_date','source_version','hgnc_id',
          'ensembl_gene_id','ncbi_gene_id','mirbase_id','ensembl_transcript_id',
          'rfam_id','secondary_structure','pdb_ids']

with open(OUTPUT_FILE, 'w') as f:
    f.write('\t'.join(FIELDS) + '\n')
    for rec in sorted(standardized, key=lambda x: x['rna_id']):
        f.write('\t'.join(str(rec[k]) for k in FIELDS) + '\n')

print(f"[OK] 已写入 {len(standardized)} 条记录")
print(f"[OK] 输出: {OUTPUT_FILE}")
print("\n✅ A5 完成！")