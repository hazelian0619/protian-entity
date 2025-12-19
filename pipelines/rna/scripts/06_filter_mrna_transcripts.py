
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
======================================================================
B3: 从 RNAcentral id_mapping 筛选 seed genes 对应的转录本 URS（修复版）
======================================================================

真实数据格式（每行 6 列，制表符分隔）：
  列1: URS（RNAcentral 主键）
  列2: database（ENSEMBL 或 ENSEMBL_GENCODE）
  列3: external_id（只有 ENST，转录本 ID）
  列4: taxon_id（9606=人类）
  列5: rna_type（mRNA/lncRNA 等，可能不准）
  列6: gene_id（ENSG，带版本号，这才是关键）

筛选逻辑：
  1. taxon_id = 9606（只要人类）
  2. database 是 ENSEMBL 或 ENSEMBL_GENCODE
  3. 第6列的 ENSG（去掉版本号）在 seed_genes.tsv 里

输出：
  mrna_urs_list.tsv
  列：urs, ensembl_transcript_id, ensembl_gene_id, gene_symbol, hgnc_id, ncbi_gene_id
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
print("B3: 筛选 mRNA 转录本 URS（修复版）")
print("=" * 70)

SEED_FILE = Path("data/processed/rna/seed_genes.tsv")
XREF_FILE = Path("data/raw/rna/rnacentral/id_mapping.tsv")
OUTPUT_FILE = Path("data/processed/rna/mrna_urs_list.tsv")

# 检查文件
for p in [SEED_FILE, XREF_FILE]:
    if not p.exists():
        print(f"[ERROR] 找不到文件: {p}")
        sys.exit(1)
    print(f"[OK] 输入: {p}")
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
print(f"[OK] 输出: {OUTPUT_FILE}\n")

# ============================================================
# 第1步：读取种子基因集，建立白名单
# ============================================================
print("[STEP 1] 读取种子基因集...")

seed_genes = {}  # {ensembl_gene_id_无版本号: {symbol, hgnc_id, ncbi_gene_id}}

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
        
        ens = parts[idx_ens]
        if not ens:
            continue
        
        # 去掉版本号（ENSG00000141510.8 → ENSG00000141510）
        ens_base = ens.split('.')[0]
        
        seed_genes[ens_base] = {
            'gene_symbol': parts[idx_symbol],
            'hgnc_id': parts[idx_hgnc],
            'ncbi_gene_id': parts[idx_ncbi],
            'ensembl_gene_id_full': ens  # 保留原始（带版本号）
        }

print(f"[OK] 种子基因数: {len(seed_genes)}")
print(f"[样例] 前3个基因: {list(seed_genes.keys())[:3]}\n")

# ============================================================
# 第2步：扫描 id_mapping.tsv，筛选符合条件的转录本
# ============================================================
print("[STEP 2] 扫描 RNAcentral id_mapping.tsv...")
print("[INFO] 这一步会处理 1.8 亿行，需要几分钟...\n")

# 真实格式（制表符分隔，6列）：
# URS, database, external_id, taxon_id, rna_type, gene_id

records = []
line_count = 0
kept_count = 0
warning_no_gene = 0  # 统计：有多少行缺第6列

with open(XREF_FILE, 'r') as f:
    for line in f:
        line_count += 1
        
        if line_count % 1000000 == 0:
            print(f"  已处理 {line_count // 1000000} 百万行... 保留 {kept_count} 条")
        
        parts = line.strip().split('\t')
        
        # 必须有至少 5 列
        if len(parts) < 5:
            continue
        
        urs       = parts[0]
        database  = parts[1]
        external  = parts[2]
        taxon_id  = parts[3]
        rna_type  = parts[4] if len(parts) > 4 else ''
        gene_id   = parts[5] if len(parts) > 5 else ''  # 第6列：基因 ID
        
        # ==================== 筛选条件 ====================
        
        # 条件1：只要人类
        if taxon_id != '9606':
            continue
        
        # 条件2：只要 Ensembl 相关的库
        if database not in ('ENSEMBL', 'ENSEMBL_GENCODE'):
            continue
        
        # 条件3：external_id 必须是 ENST 开头（转录本）
        if not external.startswith('ENST'):
            continue
        
        # 条件4：第6列必须有 gene_id（ENSG）
        if not gene_id or not gene_id.startswith('ENSG'):
            warning_no_gene += 1
            continue
        
        # 条件5：gene_id（去版本号）必须在 seed_genes 里
        gene_id_base = gene_id.split('.')[0]  # ENSG00000141510.8 → ENSG00000141510
        
        if gene_id_base not in seed_genes:
            continue
        
        # ==================== 通过筛选，记录 ====================
        
        info = seed_genes[gene_id_base]
        full_urs = urs + "_9606"
        
        records.append({
            'urs': full_urs,
            'ensembl_transcript_id': external,  # ENST...
            'ensembl_gene_id': gene_id_base,    # ENSG...（无版本号）
            'gene_symbol': info['gene_symbol'],
            'hgnc_id': info['hgnc_id'],
            'ncbi_gene_id': info['ncbi_gene_id']
        })
        kept_count += 1

print(f"\n[OK] 扫描完成，总行数: {line_count}")
print(f"[OK] 筛选出 {kept_count} 条转录本映射\n")

if warning_no_gene > 0:
    print(f"[INFO] 跳过了 {warning_no_gene} 条（缺第6列 gene_id）")

if kept_count == 0:
    print("[ERROR] 仍然没有筛选出任何转录本，请检查：")
    print("  1. seed_genes.tsv 的 ensembl_gene_id 格式")
    print("  2. id_mapping.tsv 第6列是否真的有 ENSG")
    sys.exit(1)

# ============================================================
# 第3步：去重（同一个 URS 可能有多个 database 来源）
# ============================================================
print("[STEP 3] 去重...")

# 按 urs 去重（保留第一个）
unique = {}
for rec in records:
    urs = rec['urs']
    if urs not in unique:
        unique[urs] = rec

print(f"[OK] 去重前: {len(records)} 条")
print(f"[OK] 去重后: {len(unique)} 条\n")

records = list(unique.values())

# ============================================================
# 第4步：写入输出
# ============================================================
print("[STEP 4] 写入输出...")

with open(OUTPUT_FILE, 'w') as f:
    f.write("urs\tensembl_transcript_id\tensembl_gene_id\tgene_symbol\thgnc_id\tncbi_gene_id\n")
    for rec in sorted(records, key=lambda x: x['urs']):
        f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
            rec['urs'],
            rec['ensembl_transcript_id'],
            rec['ensembl_gene_id'],
            rec['gene_symbol'],
            rec['hgnc_id'],
            rec['ncbi_gene_id']
        ))

print(f"[OK] 已写入 {len(records)} 条记录到: {OUTPUT_FILE}")

# ============================================================
# 第5步：QA 检查
# ============================================================
print("\n" + "=" * 70)
print("QA 检查")
print("=" * 70)

# 统计：有多少个唯一的基因
unique_genes = len(set(r['ensembl_gene_id'] for r in records))
print(f"[统计] 唯一转录本数: {len(records)}")
print(f"[统计] 唯一基因数: {unique_genes}")
print(f"[统计] 平均每个基因: {len(records) / unique_genes:.1f} 个转录本")

# 样例
print(f"\n[样例] 前5条记录:")
for rec in records[:5]:
    print(f"  {rec['urs']:30s} {rec['ensembl_transcript_id']:20s} {rec['gene_symbol']:10s}")

print("\n" + "=" * 70)
print("✅ B3 完成！")
print("=" * 70)
print(f"输出: {OUTPUT_FILE}")
print(f"记录数: {len(records)}")
print(f"\n下一步: B4 从 Ensembl cDNA 补序列")