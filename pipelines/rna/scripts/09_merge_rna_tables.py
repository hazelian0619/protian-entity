#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
======================================================================
M1-M3: 合并 miRNA + mRNA 子表，生成最终 RNA 主表
======================================================================

任务：
1. 读取 rna_master_mirna_v1.tsv 和 rna_master_mrna_v1.tsv
2. 检查列对齐（必须完全一致）
3. 纵向拼接（concat）
4. 全局去重检查（rna_id 唯一性）
5. 全局 QA
6. 输出 rna_master_v1.tsv
"""

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
print("=" * 70)
print("M1-M3: 合并 RNA 子表")
print("=" * 70)

MIRNA_FILE = Path("data/output/rna_master_mirna_v1.tsv")
MRNA_FILE = Path("data/output/rna_master_mrna_v1.tsv")
OUTPUT_FILE = Path("data/output/rna_master_v1.tsv")
QA_REPORT = Path("data/output/rna_master_v1_qa.txt")

# 检查文件
for p in [MIRNA_FILE, MRNA_FILE]:
    if not p.exists():
        print(f"[ERROR] 找不到文件: {p}")
        sys.exit(1)
    print(f"[OK] 输入: {p}")
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
print(f"[OK] 输出: {OUTPUT_FILE}\n")

# ============================================================
# M1: 读取两个子表
# ============================================================
print("[STEP M1] 读取子表...")

mirna_records = []
with open(MIRNA_FILE, 'r') as f:
    mirna_header = f.readline().strip().split('\t')
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == len(mirna_header):
            mirna_records.append(parts)

mrna_records = []
with open(MRNA_FILE, 'r') as f:
    mrna_header = f.readline().strip().split('\t')
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == len(mrna_header):
            mrna_records.append(parts)

print(f"[OK] miRNA 记录: {len(mirna_records)}")
print(f"[OK] mRNA 记录: {len(mrna_records)}")
print(f"[OK] 合计: {len(mirna_records) + len(mrna_records)}\n")

# ============================================================
# M2: 检查列对齐
# ============================================================
print("[STEP M2] 检查列对齐...")

if mirna_header != mrna_header:
    print("[ERROR] 两个子表的列不一致！")
    print(f"miRNA 列: {mirna_header}")
    print(f"mRNA 列: {mrna_header}")
    sys.exit(1)

print(f"[OK] 列完全一致（{len(mirna_header)} 列）")
print(f"[OK] 列名: {', '.join(mirna_header[:8])}...\n")

# ============================================================
# M3: 合并记录
# ============================================================
print("[STEP M3] 合并记录...")

all_records = mirna_records + mrna_records
print(f"[OK] 合并后记录数: {len(all_records)}\n")

# ============================================================
# M4: 全局去重检查
# ============================================================
print("[STEP M4] 全局去重检查...")

# rna_id 是第1列（索引0）
rna_ids = [rec[0] for rec in all_records]
unique_rna_ids = set(rna_ids)

if len(rna_ids) != len(unique_rna_ids):
    duplicates = len(rna_ids) - len(unique_rna_ids)
    print(f"[ERROR] 发现 {duplicates} 个重复 rna_id！")
    
    # 找出重复的 ID
    from collections import Counter
    dup_ids = [id for id, count in Counter(rna_ids).items() if count > 1]
    print(f"[ERROR] 重复的 rna_id 示例: {dup_ids[:5]}")
    sys.exit(1)

print(f"[OK] rna_id 唯一性检查通过（{len(unique_rna_ids)} 个唯一 ID）\n")

# ============================================================
# M5: 写入最终文件
# ============================================================
print("[STEP M5] 写入最终文件...")

with open(OUTPUT_FILE, 'w') as f:
    # 写表头
    f.write('\t'.join(mirna_header) + '\n')
    
    # 写记录（按 rna_id 排序）
    for rec in sorted(all_records, key=lambda x: x[0]):
        f.write('\t'.join(rec) + '\n')

print(f"[OK] 已写入 {len(all_records)} 条记录\n")

# ============================================================
# M6: 全局 QA
# ============================================================
print("=" * 70)
print("QA 检查")
print("=" * 70)

# 统计 rna_type 分布
rna_type_idx = mirna_header.index('rna_type')
rna_type_counts = {}
for rec in all_records:
    rna_type = rec[rna_type_idx]
    rna_type_counts[rna_type] = rna_type_counts.get(rna_type, 0) + 1

print("\n[统计] RNA 类型分布:")
for rna_type, count in sorted(rna_type_counts.items()):
    pct = count / len(all_records) * 100
    print(f"  {rna_type:10s}: {count:7d} ({pct:5.1f}%)")

# 统计序列长度
seq_idx = mirna_header.index('sequence')
seq_len_idx = mirna_header.index('sequence_len')
lengths = [int(rec[seq_len_idx]) for rec in all_records if rec[seq_len_idx].isdigit()]

print(f"\n[统计] 序列长度分布:")
print(f"  最短: {min(lengths)} nt")
print(f"  最长: {max(lengths)} nt")
print(f"  平均: {sum(lengths) / len(lengths):.0f} nt")
print(f"  中位数: {sorted(lengths)[len(lengths)//2]} nt")

# 统计 gene 覆盖
gene_idx = mirna_header.index('ensembl_gene_id')
unique_genes = set(rec[gene_idx] for rec in all_records if rec[gene_idx])
print(f"\n[统计] 基因覆盖:")
print(f"  唯一基因数: {len(unique_genes)}")

# QA 检查
qa_results = []

def qa_check(name, passed, expected="100%"):
    status = "✅" if passed else "❌"
    qa_results.append(f"[{status}] {name}: {expected}")
    print(f"{status} {name:50s} ({expected})")

print(f"\n[检查] 核心指标:")
qa_check("rna_id 唯一性", len(rna_ids) == len(unique_rna_ids))
qa_check("rna_type 非空", all(rec[rna_type_idx] for rec in all_records))
qa_check("sequence 非空", all(rec[seq_idx] for rec in all_records))
qa_check("记录总数符合预期", 
         290000 <= len(all_records) <= 300000,
         "290k-300k")

# 保存 QA 报告
with open(QA_REPORT, 'w') as f:
    f.write("RNA Master v1 QA 报告\n")
    f.write("=" * 70 + "\n\n")
    f.write(f"生成时间: {date.today()}\n")
    f.write(f"输出文件: {OUTPUT_FILE}\n\n")
    f.write(f"记录总数: {len(all_records)}\n")
    f.write(f"miRNA: {len(mirna_records)}\n")
    f.write(f"mRNA: {len(mrna_records)}\n\n")
    f.write("RNA 类型分布:\n")
    for rna_type, count in sorted(rna_type_counts.items()):
        f.write(f"  {rna_type}: {count}\n")
    f.write(f"\n唯一基因数: {len(unique_genes)}\n")
    f.write(f"序列长度: {min(lengths)}-{max(lengths)} nt\n\n")
    f.write("QA 检查结果:\n")
    for result in qa_results:
        f.write(result + "\n")

print(f"\n[OK] QA 报告: {QA_REPORT}")

# ============================================================
# 完成
# ============================================================
print("\n" + "=" * 70)
print("🎉🎉🎉 RNA Master v1 完成！")
print("=" * 70)
print(f"✅ 输出: {OUTPUT_FILE}")
print(f"✅ 记录数: {len(all_records):,}")
print(f"✅ miRNA: {len(mirna_records):,}")
print(f"✅ mRNA: {len(mrna_records):,}")
print(f"✅ 唯一基因: {len(unique_genes):,}")
print("\n下一步: 小分子实体（Molecule Master v1）")