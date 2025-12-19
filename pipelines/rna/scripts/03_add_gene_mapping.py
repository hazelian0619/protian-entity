#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
========================================
步骤 A4: 补充 gene 映射
========================================
输入文件1: data/processed/rna/mirna_urs_mapping.tsv (A3的输出)
输入文件2: data/raw/rna/rnacentral/id_mapping.tsv (RNAcentral xref)
输出文件: data/processed/rna/mirna_with_gene.tsv

功能：
1. 读取带 URS 的 miRNA 列表
2. 从 RNAcentral xref 中找到 URS -> HGNC/Ensembl/NCBI 的映射
3. 补充 gene_symbol, hgnc_id, ensembl_gene_id, ncbi_gene_id
4. 输出完整的 miRNA 表（包含 gene 信息）
"""

import os
import sys
from pathlib import Path
from collections import defaultdict

# ============================================================
# 第1步：自动定位项目根目录
# ============================================================
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
print("A4: 补充 gene 映射")
print("=" * 70)
print(f"[INFO] 项目根目录: {project_root}\n")

# ============================================================
# 第2步：定义输入输出路径
# ============================================================
INPUT_MIRNA = Path("data/processed/rna/mirna_urs_mapping.tsv")
INPUT_XREF = Path("data/raw/rna/rnacentral/id_mapping.tsv")
OUTPUT_TSV = Path("data/processed/rna/mirna_with_gene.tsv")

for fpath in [INPUT_MIRNA, INPUT_XREF]:
    if not fpath.exists():
        print(f"[ERROR] 找不到文件: {fpath}")
        sys.exit(1)
    print(f"[OK] 输入文件存在: {fpath}")

OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
print(f"[OK] 输出路径: {OUTPUT_TSV}\n")

# ============================================================
# 第3步：读取 A3 的 miRNA 列表
# ============================================================
print("[STEP 1] 读取 miRNA 列表...")

mirna_dict = {}  # {urs: {mirbase_id, mimat_id, sequence, ...}}

with open(INPUT_MIRNA, 'r') as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 5:
            continue
        urs, mirbase_id, mimat_id, sequence, sequence_len = parts
        
        mirna_dict[urs] = {
            'mirbase_id': mirbase_id,
            'mimat_id': mimat_id,
            'sequence': sequence,
            'sequence_len': int(sequence_len),
            # 预设 gene 字段为空
            'gene_symbol': '',
            'hgnc_id': '',
            'ensembl_gene_id': '',
            'ncbi_gene_id': ''
        }

print(f"[OK] 读取到 {len(mirna_dict)} 个 miRNA\n")

# ============================================================
# 第4步：从 RNAcentral xref 提取 gene 映射
# ============================================================
print("[STEP 2] 从 RNAcentral xref 提取 gene 映射...")
print("[INFO] 这一步需要 1-2 分钟（扫描 11GB 文件）...\n")

# 我们需要找这些数据库的映射：
# - HGNC: HGNC:xxxxx (基因官方 ID)
# - Ensembl: ENSG... (Ensembl 基因 ID)
# - NCBI Gene: 数字 (NCBI 基因 ID)

urs_to_hgnc = {}       # {urs: hgnc_id}
urs_to_ensembl = {}    # {urs: ensembl_gene_id}
urs_to_ncbi = {}       # {urs: ncbi_gene_id}
urs_to_symbol = {}     # {urs: gene_symbol}

line_count = 0
relevant_count = 0

with open(INPUT_XREF, 'r') as f:
    for line in f:
        line_count += 1
        
        if line_count % 1000000 == 0:
            print(f"  已处理 {line_count // 1000000} 百万行...")
        
        parts = line.strip().split('\t')
        if len(parts) < 5:
            continue
        
        urs = parts[0]
        database = parts[1]
        external_id = parts[2]
        taxon_id = parts[3]
        
        # 只处理人类 + 我们已有的 URS
        full_urs = urs if urs.endswith('_9606') else urs + '_9606'
        if taxon_id != '9606' or full_urs not in mirna_dict:
            continue
        
        relevant_count += 1
        
        # 根据数据库类型提取不同的 ID
        if database == 'HGNC':
            # HGNC ID 格式: HGNC:31586
            if external_id.startswith('HGNC:'):
                urs_to_hgnc[full_urs] = external_id
                # 尝试从 external_id 的其他部分提取 symbol
                # 注意：id_mapping.tsv 可能没有 symbol，需要后续补充
        
        elif database == 'ENA' and 'gene' in external_id.lower():
            # 某些情况下，ENA 条目会包含 gene symbol
            # 格式可能是: gene:MIR21 或类似
            pass  # 暂时跳过，因为格式不统一
        
        elif database == 'ENSEMBL':
            # Ensembl 基因 ID 格式: ENSG00000...
            if external_id.startswith('ENSG'):
                urs_to_ensembl[full_urs] = external_id
        
        elif database == 'NCBI_GENE' or database == 'NCBI':
            # NCBI Gene ID 是纯数字
            if external_id.isdigit():
                urs_to_ncbi[full_urs] = external_id

print(f"[OK] 处理完成，总计 {line_count} 行")
print(f"[OK] 找到与我们 miRNA 相关的映射 {relevant_count} 条")
print(f"[统计] HGNC 映射: {len(urs_to_hgnc)}")
print(f"[统计] Ensembl 映射: {len(urs_to_ensembl)}")
print(f"[统计] NCBI 映射: {len(urs_to_ncbi)}\n")

# ============================================================
# 第5步：合并 gene 信息到 miRNA 表
# ============================================================
print("[STEP 3] 合并 gene 信息...")

for urs, info in mirna_dict.items():
    if urs in urs_to_hgnc:
        info['hgnc_id'] = urs_to_hgnc[urs]
    if urs in urs_to_ensembl:
        info['ensembl_gene_id'] = urs_to_ensembl[urs]
    if urs in urs_to_ncbi:
        info['ncbi_gene_id'] = urs_to_ncbi[urs]
    
    # 尝试从 mirbase_id 推断 gene_symbol
    # 例如: hsa-miR-21-5p -> MIR21
    mirbase_id = info['mirbase_id']
    if mirbase_id.startswith('hsa-'):
        # 提取 miR-xxx 部分
        parts = mirbase_id.replace('hsa-', '').split('-')
        if len(parts) >= 2 and parts[0] in ['miR', 'let']:
            if parts[0] == 'let':
                # hsa-let-7a-5p -> MIRLET7A
                info['gene_symbol'] = 'MIRLET' + parts[1].upper()
            else:
                # hsa-miR-21-5p -> MIR21
                info['gene_symbol'] = 'MIR' + parts[1].upper()

# 统计 gene 映射覆盖率
has_hgnc = sum(1 for v in mirna_dict.values() if v['hgnc_id'])
has_ensembl = sum(1 for v in mirna_dict.values() if v['ensembl_gene_id'])
has_ncbi = sum(1 for v in mirna_dict.values() if v['ncbi_gene_id'])
has_symbol = sum(1 for v in mirna_dict.values() if v['gene_symbol'])
has_any_gene = sum(1 for v in mirna_dict.values() if v['hgnc_id'] or v['ensembl_gene_id'] or v['ncbi_gene_id'])

print(f"[统计] Gene 映射覆盖率:")
print(f"  hgnc_id 非空: {has_hgnc} ({has_hgnc/len(mirna_dict)*100:.1f}%)")
print(f"  ensembl_gene_id 非空: {has_ensembl} ({has_ensembl/len(mirna_dict)*100:.1f}%)")
print(f"  ncbi_gene_id 非空: {has_ncbi} ({has_ncbi/len(mirna_dict)*100:.1f}%)")
print(f"  gene_symbol 非空: {has_symbol} ({has_symbol/len(mirna_dict)*100:.1f}%)")
print(f"  至少有一个 gene ID: {has_any_gene} ({has_any_gene/len(mirna_dict)*100:.1f}%)\n")

# ============================================================
# 第6步：写入输出文件
# ============================================================
print("[STEP 4] 写入输出文件...")

with open(OUTPUT_TSV, 'w') as f:
    f.write("urs\tmirbase_id\tmimat_id\tsequence\tsequence_len\t")
    f.write("gene_symbol\thgnc_id\tensembl_gene_id\tncbi_gene_id\n")
    
    for urs in sorted(mirna_dict.keys()):
        info = mirna_dict[urs]
        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            urs,
            info['mirbase_id'],
            info['mimat_id'],
            info['sequence'],
            info['sequence_len'],
            info['gene_symbol'],
            info['hgnc_id'],
            info['ensembl_gene_id'],
            info['ncbi_gene_id']
        ))

print(f"[OK] 已写入 {len(mirna_dict)} 条记录到: {OUTPUT_TSV}\n")

# ============================================================
# 第7步：质量检查
# ============================================================
print("=" * 70)
print("质量检查报告")
print("=" * 70)

# 显示样例（有 gene 的 + 没 gene 的各几个）
with_gene = [v for v in mirna_dict.values() if v['hgnc_id'] or v['ensembl_gene_id']]
without_gene = [v for v in mirna_dict.values() if not (v['hgnc_id'] or v['ensembl_gene_id'])]

print("\n[样例] 有 gene 映射的 miRNA (前3个):")
for info in with_gene[:3]:
    print(f"  {info['mirbase_id']:20s} -> {info['gene_symbol']:10s} {info['hgnc_id']:15s}")

print("\n[样例] 无 gene 映射的 miRNA (前3个):")
for info in without_gene[:3]:
    print(f"  {info['mirbase_id']:20s} -> (无 gene 映射)")

print(f"\n[检查] miRNA 的 gene 映射情况:")
print(f"  有 gene 信息: {len(with_gene)} ({len(with_gene)/len(mirna_dict)*100:.1f}%)")
print(f"  无 gene 信息: {len(without_gene)} ({len(without_gene)/len(mirna_dict)*100:.1f}%)")

if len(with_gene) / len(mirna_dict) > 0.3:
    print(f"  ✅ gene 覆盖率合理（miRNA 的 gene 映射本身就稀疏）")
else:
    print(f"  ⚠️  gene 覆盖率较低，但对 miRNA 来说可接受")

# ============================================================
# 第8步：最终验收
# ============================================================
print("\n" + "=" * 70)
print("✅✅ A4 步骤完成！")
print("=" * 70)
print(f"输出文件: {OUTPUT_TSV}")
print(f"记录数量: {len(mirna_dict)}")
print(f"Gene 覆盖率: {has_any_gene/len(mirna_dict)*100:.1f}%")
print(f"\n下一步: 运行 A5 脚本（生成标准化 miRNA 子表）")