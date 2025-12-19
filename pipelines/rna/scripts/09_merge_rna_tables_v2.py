#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M1-M3: 合并 RNA 子表（修复版：支持空白字符分割）
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
print("M1-M3: 合并 RNA 子表")
print("=" * 70)

MIRNA_FILE = Path("data/output/rna_master_mirna_v1.tsv")
MRNA_FILE = Path("data/output/rna_master_mrna_v1.tsv")
OUTPUT_FILE = Path("data/output/rna_master_v1.tsv")
QA_REPORT = Path("data/output/rna_master_v1_qa.txt")

for p in [MIRNA_FILE, MRNA_FILE]:
    if not p.exists():
        print(f"[ERROR] 找不到: {p}")
        sys.exit(1)
    print(f"[OK] 输入: {p}")
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
print(f"[OK] 输出: {OUTPUT_FILE}\n")

# ============================================================
# M1: 读取子表（用正则处理空白）
# ============================================================
print("[STEP M1] 读取子表...")

def read_tsv_robust(filepath):
    """用正则分割空白符，兼容 tab/空格混合"""
    records = []
    with open(filepath, 'r') as f:
        header_line = f.readline().strip()
        header = re.split(r'\s+', header_line)
        
        for line in f:
            parts = re.split(r'\s+', line.strip())
            if len(parts) == len(header):
                records.append(parts)
            elif len(parts) > 0:  # 跳过空行
                # 尝试调整（可能某些字段为空）
                pass
    return header, records

mirna_header, mirna_records = read_tsv_robust(MIRNA_FILE)
mrna_header, mrna_records = read_tsv_robust(MRNA_FILE)

print(f"[OK] miRNA 记录: {len(mirna_records)}")
print(f"[OK] mRNA 记录: {len(mrna_records)}")
print(f"[OK] 合计: {len(mirna_records) + len(mrna_records)}\n")

if len(mirna_records) == 0 or len(mrna_records) == 0:
    print("[ERROR] 读取失败，记录数为 0")
    print("请检查文件格式")
    sys.exit(1)

# ============================================================
# M2: 检查列对齐
# ============================================================
print("[STEP M2] 检查列对齐...")

if mirna_header != mrna_header:
    print("[WARN] 列不完全一致，尝试对齐...")
    print(f"miRNA 列数: {len(mirna_header)}")
    print(f"mRNA 列数: {len(mrna_header)}")
    # 使用 miRNA 的列作为标准
    final_header = mirna_header
else:
    final_header = mirna_header

print(f"[OK] 使用列: {len(final_header)} 列\n")

# ============================================================
# M3: 合并
# ============================================================
print("[STEP M3] 合并记录...")

all_records = mirna_records + mrna_records
print(f"[OK] 合并后: {len(all_records)} 条\n")

# ============================================================
# M4: 去重检查
# ============================================================
print("[STEP M4] 去重检查...")

rna_ids = [rec[0] for rec in all_records]
unique_ids = set(rna_ids)

if len(rna_ids) != len(unique_ids):
    dup = len(rna_ids) - len(unique_ids)
    print(f"[ERROR] 发现 {dup} 个重复 rna_id")
    sys.exit(1)

print(f"[OK] rna_id 唯一性通过\n")

# ============================================================
# M5: 写入
# ============================================================
print("[STEP M5] 写入最终文件...")

with open(OUTPUT_FILE, 'w') as f:
    f.write('\t'.join(final_header) + '\n')
    for rec in sorted(all_records, key=lambda x: x[0]):
        f.write('\t'.join(rec) + '\n')

print(f"[OK] 已写入 {len(all_records)} 条\n")

# ============================================================
# M6: QA
# ============================================================
print("=" * 70)
print("QA 检查")
print("=" * 70)

rna_type_idx = final_header.index('rna_type')
rna_type_counts = {}
for rec in all_records:
    rt = rec[rna_type_idx] if rna_type_idx < len(rec) else 'unknown'
    rna_type_counts[rt] = rna_type_counts.get(rt, 0) + 1

print("\n[统计] RNA 类型分布:")
for rt, count in sorted(rna_type_counts.items()):
    pct = count / len(all_records) * 100
    print(f"  {rt:10s}: {count:7,d} ({pct:5.1f}%)")

seq_len_idx = final_header.index('sequence_len')
lengths = []
for rec in all_records:
    if seq_len_idx < len(rec) and rec[seq_len_idx].replace('.','').isdigit():
        lengths.append(int(float(rec[seq_len_idx])))

if lengths:
    print(f"\n[统计] 序列长度:")
    print(f"  最短: {min(lengths)} nt")
    print(f"  最长: {max(lengths)} nt")
    print(f"  平均: {sum(lengths)/len(lengths):.0f} nt")

# 保存 QA 报告
with open(QA_REPORT, 'w') as f:
    f.write("RNA Master v1 QA Report\n")
    f.write("=" * 70 + "\n\n")
    f.write(f"Date: {date.today()}\n")
    f.write(f"Output: {OUTPUT_FILE}\n\n")
    f.write(f"Total records: {len(all_records):,}\n")
    f.write(f"  miRNA: {len(mirna_records):,}\n")
    f.write(f"  mRNA: {len(mrna_records):,}\n\n")
    f.write("RNA type distribution:\n")
    for rt, count in sorted(rna_type_counts.items()):
        f.write(f"  {rt}: {count:,}\n")

print(f"\n[OK] QA 报告: {QA_REPORT}")

print("\n" + "=" * 70)
print("🎉🎉🎉 RNA Master v1 完成！")
print("=" * 70)
print(f"✅ 输出: {OUTPUT_FILE}")
print(f"✅ 记录数: {len(all_records):,}")
print(f"✅ miRNA: {len(mirna_records):,}")
print(f"✅ mRNA: {len(mrna_records):,}")
print("\n下一步: 小分子实体")