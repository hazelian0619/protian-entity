#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
B1: 提取种子基因集（从 protein_master_v6）
"""

import os
import sys
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
print("B1: 提取种子基因集")
print("=" * 70)

# 可能的文件位置（按优先级）
possible_paths = [
    Path("data/raw/protein_master_v6_clean.tsv"),                                      # 方案A/B 创建后的位置
    Path("/Users/pluviophile/graph/1025/data/processed/protein_master_v6_clean.tsv"),  # 旧项目直接引用
    Path("data/processed/protein_master_v6_clean.tsv"),
    Path("protein_master_v6_clean.tsv")
]

INPUT_FILE = None
for p in possible_paths:
    if p.exists():
        INPUT_FILE = p
        print(f"[OK] 找到文件: {p}")
        break

if INPUT_FILE is None:
    print("[ERROR] 找不到 protein_master_v6_clean.tsv")
    print("\n请执行以下命令之一：")
    print("\n【方案A：软链接（推荐）】")
    print("  mkdir -p data/raw")
    print("  ln -s /Users/pluviophile/graph/1025/data/processed/protein_master_v6_clean.tsv \\")
    print("        data/raw/protein_master_v6_clean.tsv")
    print("\n【方案B：复制文件】")
    print("  mkdir -p data/raw")
    print("  cp /Users/pluviophile/graph/1025/data/processed/protein_master_v6_clean.tsv \\")
    print("     data/raw/protein_master_v6_clean.tsv")
    sys.exit(1)

OUTPUT_FILE = Path("data/processed/rna/seed_genes.tsv")
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

print(f"[OK] 输出: {OUTPUT_FILE}\n")

# ============================================================
# 读取蛋白表
# ============================================================
print("[STEP 1] 读取蛋白表...")

genes = {}

with open(INPUT_FILE, 'r') as f:
    header_line = f.readline().strip()
    headers = header_line.split('\t')
    
    print(f"[INFO] 列数: {len(headers)}")
    print(f"[INFO] 前10列: {', '.join(headers[:10])}\n")
    
    # 智能识别列名（可能叫 gene_symbol 或 symbol）
    idx_symbol = None
    for i, h in enumerate(headers):
        if 'symbol' in h.lower() and 'gene' in h.lower():
            idx_symbol = i
            break
    if idx_symbol is None:
        for i, h in enumerate(headers):
            if h.lower() == 'symbol':
                idx_symbol = i
                break
    
    idx_hgnc = next((i for i, h in enumerate(headers) if 'hgnc' in h.lower()), None)
    idx_ensembl = next((i for i, h in enumerate(headers) if 'ensembl_gene' in h.lower()), None)
    idx_ncbi = next((i for i, h in enumerate(headers) if 'ncbi_gene' in h.lower() or 'entrez' in h.lower()), None)
    
    print(f"[INFO] 识别到的列索引:")
    print(f"  symbol: {idx_symbol} ({headers[idx_symbol] if idx_symbol is not None else 'N/A'})")
    print(f"  hgnc: {idx_hgnc} ({headers[idx_hgnc] if idx_hgnc is not None else 'N/A'})")
    print(f"  ensembl: {idx_ensembl} ({headers[idx_ensembl] if idx_ensembl is not None else 'N/A'})")
    print(f"  ncbi: {idx_ncbi} ({headers[idx_ncbi] if idx_ncbi is not None else 'N/A'})\n")
    
    if idx_ensembl is None:
        print("[ERROR] 必须有 ensembl_gene_id 列")
        sys.exit(1)
    
    for line_num, line in enumerate(f, 1):
        parts = line.strip().split('\t')
        if len(parts) <= max(filter(None, [idx_symbol, idx_hgnc, idx_ensembl, idx_ncbi])):
            continue
        
        ensembl = parts[idx_ensembl] if idx_ensembl is not None else ''
        
        # 跳过空值和占位符
        if not ensembl or ensembl in ['', 'NA', 'N/A', 'None', '-']:
            continue
        
        if ensembl not in genes:
            genes[ensembl] = {
                'symbol': parts[idx_symbol] if idx_symbol is not None else '',
                'hgnc_id': parts[idx_hgnc] if idx_hgnc is not None else '',
                'ensembl_gene_id': ensembl,
                'ncbi_gene_id': parts[idx_ncbi] if idx_ncbi is not None else ''
            }

print(f"[OK] 提取到 {len(genes)} 个唯一基因\n")

if len(genes) < 10000:
    print(f"[WARN] 基因数量偏少（{len(genes)}），预期应该 >15000")
    print("       请检查 protein_master_v6 文件是否正确")

# ============================================================
# 写入输出
# ============================================================
print("[STEP 2] 写入种子基因表...")

with open(OUTPUT_FILE, 'w') as f:
    f.write("gene_symbol\thgnc_id\tensembl_gene_id\tncbi_gene_id\n")
    for gene_info in sorted(genes.values(), key=lambda x: x['ensembl_gene_id']):
        f.write("{}\t{}\t{}\t{}\n".format(
            gene_info['symbol'],
            gene_info['hgnc_id'],
            gene_info['ensembl_gene_id'],
            gene_info['ncbi_gene_id']
        ))

print(f"[OK] 已写入 {len(genes)} 个基因")

# ============================================================
# QA
# ============================================================
ensembl_coverage = 100.0  # 我们用 ensembl 去重，所以一定是100%
symbol_coverage = sum(1 for g in genes.values() if g['symbol']) / len(genes) * 100
hgnc_coverage = sum(1 for g in genes.values() if g['hgnc_id']) / len(genes) * 100

print("\n" + "=" * 70)
print("QA 检查")
print("=" * 70)
print(f"✅ ensembl_gene_id 非空率: {ensembl_coverage:.1f}%")
print(f"{'✅' if symbol_coverage >= 90 else '⚠️ '} gene_symbol 非空率: {symbol_coverage:.1f}%")
print(f"{'✅' if hgnc_coverage >= 70 else '⚠️ '} hgnc_id 非空率: {hgnc_coverage:.1f}%")

print("\n" + "=" * 70)
print("✅ B1 完成！")
print("=" * 70)
print(f"输出: {OUTPUT_FILE}")
print(f"种子基因数: {len(genes)}")
print("\n下一步: B2 下载 Ensembl cDNA 文件")