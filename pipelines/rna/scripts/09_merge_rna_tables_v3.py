#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M1-M3: 合并 RNA 子表（简化版，纯 tab 分隔）
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

for p in [MIRNA_FILE, MRNA_FILE]:
    if not p.exists():
        print(f"[ERROR] 找不到: {p}")
        sys.exit(1)

print(f"[OK] 输入: {MIRNA_FILE}")
print(f"[OK] 输入: {MRNA_FILE}")
print(f"[OK] 输出: {OUTPUT_FILE}\n")

# ============================================================
# M1: 读取（纯 tab 分隔）
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

print(f"[OK] miRNA 记录: {len(mirna_records):,}")
print(f"[OK] mRNA 记录: {len(mrna_records):,}")
print(f"[OK] 合计: {len(mirna_records) + len(mrna_records):,}\n")

if len(mirna_records) == 0 or len(mrna_records) == 0:
    print("[ERROR] 读取失败，记录数为 0")
    print(f"miRNA 表头列数: {len(mirna_header)}")
    print(f"mRNA 表头列数: {len(mrna_header)}")
    print("请检查：")
    print("  1. 文件是否用纯 tab 分隔")
    print("  2. 每行列数是否一致")
    sys.exit(1)

# ============================================================
# M2: 检查列对齐
# ============================================================
print("[STEP M2] 检查列对齐...")

if mirna_header != mrna_header:
    print("[ERROR] 两表列名不一致")
    print(f"miRNA: {mirna_header}")
    print(f"mRNA: {mrna_header}")
    sys.exit(1)

print(f"[OK] 列一致（{len(mirna_header)} 列）\n")

# ============================================================
# M3: 合并
# ============================================================
print("[STEP M3] 合并...")

all_records = mirna_records + mrna_records
print(f"[OK] 合并后: {len(all_records):,} 条\n")

# ============================================================
# M4: 去重检查
# ============================================================
print("[STEP M4] 去重检查...")

rna_ids = [rec[0] for rec in all_records]
unique_ids = set(rna_ids)

if len(rna_ids) != len(unique_ids):
    dup = len(rna_ids) - len(unique_ids)
    print(f"[ERROR] 发现 {dup} 个重复")
    sys.exit(1)

print(f"[OK] rna_id 唯一性通过（{len(unique_ids):,} 个）\n")

# ============================================================
# M5: 写入
# ============================================================
print("[STEP M5] 写入...")

with open(OUTPUT_FILE, 'w') as f:
    f.write('\t'.join(mirna_header) + '\n')
    for rec in sorted(all_records, key=lambda x: x[0]):
        f.write('\t'.join(rec) + '\n')

print(f"[OK] 已写入 {len(all_records):,} 条\n")

# ============================================================
# M6: QA
# ============================================================
print("=" * 70)
print("QA 检查")
print("=" * 70)

rna_type_idx = mirna_header.index('rna_type')
seq_len_idx = mirna_header.index('sequence_len')

rna_types = {}
for rec in all_records:
    rt = rec[rna_type_idx]
    rna_types[rt] = rna_types.get(rt, 0) + 1

print("\n[统计] RNA 类型:")
for rt in sorted(rna_types.keys()):
    count = rna_types[rt]
    pct = count / len(all_records) * 100
    print(f"  {rt:10s}: {count:7,d} ({pct:5.1f}%)")

lengths = [int(rec[seq_len_idx]) for rec in all_records if rec[seq_len_idx].isdigit()]
if lengths:
    print(f"\n[统计] 序列长度:")
    print(f"  最短: {min(lengths):,} nt")
    print(f"  最长: {max(lengths):,} nt")
    print(f"  平均: {sum(lengths)//len(lengths):,} nt")

print("\n" + "=" * 70)
print("🎉🎉🎉 RNA Master v1 完成！")
print("=" * 70)
print(f"✅ 输出: {OUTPUT_FILE}")
print(f"✅ 总记录: {len(all_records):,}")
print(f"✅   miRNA: {len(mirna_records):,}")
print(f"✅   mRNA: {len(mrna_records):,}")
print(f"✅ 唯一 ID: {len(unique_ids):,}")
print("\n🚀 下一步: 小分子实体（Molecule Master v1）")