#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
======================================================================
B3+B4 合并版：直接从 Ensembl cDNA FASTA 提取 mRNA（绕过 RNAcentral）
======================================================================

新策略（因为 RNAcentral 映射不完整）：
1. 直接解析 Ensembl cDNA FASTA
2. 从 header 里提取：transcript_id (ENST), gene_id (ENSG), gene_symbol
3. 筛选：gene_id 在 seed_genes 里的
4. 输出：mrna_with_sequence.tsv（包含序列）

Header 格式示例：
>ENST00000622028.1 cdna chromosome:GRCh38:21:10649400:10649835:-1 gene:ENSG00000277282.1 gene_biotype:IG_V_gene transcript_biotype:IG_V_gene gene_symbol:IGHV1OR21-1 ...
"""

import os
import sys
import re
from pathlib import Path

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
print("B3+B4: 从 Ensembl cDNA 直接提取 mRNA")
print("=" * 70)

SEED_FILE = Path("data/processed/rna/seed_genes.tsv")
FASTA_FILE = Path("data/raw/rna/ensembl/Homo_sapiens.GRCh38.cdna.all.fa")
OUTPUT_FILE = Path("data/processed/rna/mrna_with_sequence.tsv")

# 检查文件
for p in [SEED_FILE, FASTA_FILE]:
    if not p.exists():
        print(f"[ERROR] 找不到文件: {p}")
        sys.exit(1)
    print(f"[OK] 输入: {p}")
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
print(f"[OK] 输出: {OUTPUT_FILE}\n")

# ============================================================
# 第1步：读取 seed_genes，建立白名单
# ============================================================
print("[STEP 1] 读取种子基因集...")

seed_genes = {}  # {ensembl_gene_id: {symbol, hgnc_id, ncbi_gene_id}}

with open(SEED_FILE, 'r') as f:
    header = f.readline().strip().split('\t')
    idx_symbol = header.index('gene_symbol')
    idx_hgnc   = header.index('hgnc_id')
    idx_ens    = header.index('ensembl_gene_id')
    idx_ncbi   = header.index('ncbi_gene_id')
    
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 4:
            continue
        ens = parts[idx_ens].split('.')[0]  # 去版本号
        if not ens:
            continue
        seed_genes[ens] = {
            'gene_symbol': parts[idx_symbol],
            'hgnc_id': parts[idx_hgnc],
            'ncbi_gene_id': parts[idx_ncbi]
        }

print(f"[OK] 种子基因数: {len(seed_genes)}\n")

# ============================================================
# 第2步：解析 Ensembl cDNA FASTA
# ============================================================
print("[STEP 2] 解析 Ensembl cDNA FASTA...")
print("[INFO] 这一步会读取完整 FASTA，约需 1-2 分钟...\n")

records = []
current_transcript_id = None
current_gene_id = None
current_gene_symbol = None
current_seq = []
line_count = 0
kept_count = 0

with open(FASTA_FILE, 'r') as f:
    for line in f:
        line_count += 1
        if line_count % 100000 == 0:
            print(f"  已处理 {line_count // 1000}k 行... 保留 {kept_count} 条转录本")
        
        line = line.strip()
        
        if line.startswith('>'):
            # 保存上一个序列（如果基因在白名单）
            if current_transcript_id and current_gene_id and current_gene_id in seed_genes:
                seq = ''.join(current_seq).upper().replace('T', 'U')
                info = seed_genes[current_gene_id]
                records.append({
                    'ensembl_transcript_id': current_transcript_id,
                    'ensembl_gene_id': current_gene_id,
                    'gene_symbol': info['gene_symbol'],
                    'hgnc_id': info['hgnc_id'],
                    'ncbi_gene_id': info['ncbi_gene_id'],
                    'sequence': seq,
                    'sequence_len': len(seq)
                })
                kept_count += 1
            
            # 解析新 header
            # 格式：>ENST00000622028.1 cdna ... gene:ENSG00000277282.1 ... gene_symbol:IGHV1OR21-1 ...
            
            # 提取 transcript ID（第一个空格前）
            parts = line[1:].split()
            current_transcript_id = parts[0] if parts else None
            
            # 提取 gene ID（gene:ENSG...）
            gene_match = re.search(r'gene:(ENSG\d+)\.?\d*', line)
            current_gene_id = gene_match.group(1) if gene_match else None
            
            # 提取 gene_symbol（gene_symbol:XXX）
            symbol_match = re.search(r'gene_symbol:(\S+)', line)
            current_gene_symbol = symbol_match.group(1) if symbol_match else None
            
            # 重置序列缓存
            current_seq = []
        else:
            # 序列行
            current_seq.append(line)
    
    # 处理最后一个序列
    if current_transcript_id and current_gene_id and current_gene_id in seed_genes:
        seq = ''.join(current_seq).upper().replace('T', 'U')
        info = seed_genes[current_gene_id]
        records.append({
            'ensembl_transcript_id': current_transcript_id,
            'ensembl_gene_id': current_gene_id,
            'gene_symbol': info['gene_symbol'],
            'hgnc_id': info['hgnc_id'],
            'ncbi_gene_id': info['ncbi_gene_id'],
            'sequence': seq,
            'sequence_len': len(seq)
        })
        kept_count += 1

print(f"\n[OK] FASTA 解析完成")
print(f"[OK] 总行数: {line_count}")
print(f"[OK] 提取到 {kept_count} 条转录本（来自 seed_genes）\n")

if kept_count == 0:
    print("[ERROR] 没有提取到任何转录本，请检查：")
    print("  1. seed_genes.tsv 的 ensembl_gene_id 是否正确")
    print("  2. Ensembl cDNA FASTA 文件是否完整")
    sys.exit(1)

# ============================================================
# 第3步：生成 rna_id（URS 或临时 ID）
# ============================================================
print("[STEP 3] 生成 rna_id...")

# 简化版：直接用 ENST_9606 作为 rna_id（后续可以与 RNAcentral 对齐）
for rec in records:
    transcript_base = rec['ensembl_transcript_id'].split('.')[0]  # 去版本号
    rec['rna_id'] = f"{transcript_base}_9606"

print(f"[OK] 已生成 rna_id\n")

# ============================================================
# 第4步：写入输出
# ============================================================
print("[STEP 4] 写入输出...")

with open(OUTPUT_FILE, 'w') as f:
    f.write("rna_id\tensembl_transcript_id\tensembl_gene_id\tgene_symbol\thgnc_id\tncbi_gene_id\tsequence\tsequence_len\n")
    for rec in sorted(records, key=lambda x: x['rna_id']):
        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            rec['rna_id'],
            rec['ensembl_transcript_id'],
            rec['ensembl_gene_id'],
            rec['gene_symbol'],
            rec['hgnc_id'],
            rec['ncbi_gene_id'],
            rec['sequence'],
            rec['sequence_len']
        ))

print(f"[OK] 已写入 {len(records)} 条记录")

# ============================================================
# 第5步：QA 检查
# ============================================================
print("\n" + "=" * 70)
print("QA 检查")
print("=" * 70)

unique_genes = len(set(r['ensembl_gene_id'] for r in records))
lengths = [r['sequence_len'] for r in records]

print(f"[统计] 转录本总数: {len(records)}")
print(f"[统计] 唯一基因数: {unique_genes}")
print(f"[统计] 平均每个基因: {len(records) / unique_genes:.1f} 个转录本")
print(f"[统计] 序列长度: {min(lengths)}-{max(lengths)} nt (平均 {sum(lengths)/len(lengths):.0f})")

print(f"\n[样例] 前5条记录:")
for rec in records[:5]:
    print(f"  {rec['rna_id']:30s} {rec['gene_symbol']:10s} {rec['sequence'][:30]}... (len={rec['sequence_len']})")

print("\n" + "=" * 70)
print("✅ B3+B4 完成！")
print("=" * 70)
print(f"输出: {OUTPUT_FILE}")
print(f"记录数: {len(records)}")
print("\n下一步: B5 标准化为最终 mRNA 子表")