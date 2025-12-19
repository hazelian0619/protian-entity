#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
========================================
步骤 A2: 从 miRBase 提取人类 miRNA
========================================
输入文件: data/raw/rna/mirbase/mature.fa
输出文件: data/processed/rna/hsa_mirna_list.tsv

功能：
1. 读取 miRBase 的 mature.fa 文件（包含所有物种的 miRNA）
2. 只保留人类（hsa- 开头）的 miRNA
3. 提取序列并转换成 RNA 格式（T->U）
4. 去重（如果同一个 ID 有多个序列，保留最长的）
5. 输出标准 TSV 格式
"""

import os
import sys
from pathlib import Path

# ============================================================
# 第1步：自动定位项目根目录（无论在哪里运行都能找到）
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
print(f"[INFO] 项目根目录: {project_root}")
print(f"[INFO] 当前工作目录: {os.getcwd()}\n")

# ============================================================
# 第2步：定义输入输出路径
# ============================================================
INPUT_FASTA = Path("data/raw/rna/mirbase/mature.fa")
OUTPUT_TSV = Path("data/processed/rna/hsa_mirna_list.tsv")

# 检查输入文件是否存在
if not INPUT_FASTA.exists():
    print(f"[ERROR] 找不到输入文件: {INPUT_FASTA}")
    print(f"[ERROR] 当前目录: {os.getcwd()}")
    print(f"[ERROR] 请确认文件是否在正确位置")
    sys.exit(1)

# 创建输出目录（如果不存在）
OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
print(f"[OK] 输入文件存在: {INPUT_FASTA}")
print(f"[OK] 输出路径: {OUTPUT_TSV}\n")

# ============================================================
# 第3步：解析 FASTA 文件
# ============================================================
print("[STEP 1] 开始解析 FASTA 文件...")

records = []              # 存储所有人类 miRNA 记录
current_id = None         # 当前读取的 miRNA 名称（如 hsa-miR-21）
current_mimat = None      # 当前读取的 MIMAT ID（如 MIMAT0000076）
current_seq = []          # 当前读取的序列片段

with open(INPUT_FASTA, 'r') as f:
    for line_num, line in enumerate(f, 1):
        line = line.strip()
        
        if line.startswith('>'):
            # 遇到新的序列头（以 > 开头）
            
            # 先保存上一个序列（如果是人类的）
            if current_id and current_id.startswith('hsa-'):
                seq = ''.join(current_seq).upper().replace('T', 'U')  # DNA->RNA
                records.append({
                    'mirbase_id': current_id,      # 如 hsa-miR-21-5p
                    'mimat_id': current_mimat,     # 如 MIMAT0000076
                    'sequence': seq,               # RNA 序列（AUGC）
                    'sequence_len': len(seq)       # 序列长度
                })
            
            # 解析新的序列头
            # 格式示例: >hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p
            parts = line[1:].split()  # 去掉 > 号，按空格拆分
            current_id = parts[0] if len(parts) > 0 else None
            current_mimat = parts[1] if len(parts) > 1 else None
            current_seq = []  # 重置序列缓存
            
        else:
            # 普通行，属于序列内容
            current_seq.append(line)
    
    # 处理文件最后一个序列
    if current_id and current_id.startswith('hsa-'):
        seq = ''.join(current_seq).upper().replace('T', 'U')
        records.append({
            'mirbase_id': current_id,
            'mimat_id': current_mimat,
            'sequence': seq,
            'sequence_len': len(seq)
        })

print(f"[OK] 初步提取到 {len(records)} 条人类 miRNA 记录\n")

# ============================================================
# 第4步：去重（按 mirbase_id，保留序列最长的）
# ============================================================
print("[STEP 2] 去重处理...")

unique = {}
duplicate_count = 0

for rec in records:
    mid = rec['mirbase_id']
    if mid not in unique:
        # 第一次遇到这个 ID
        unique[mid] = rec
    else:
        # 发现重复，比较序列长度
        duplicate_count += 1
        if len(rec['sequence']) > len(unique[mid]['sequence']):
            # 新序列更长，替换
            print(f"  [WARN] 发现重复 ID: {mid}，保留更长的序列")
            unique[mid] = rec

records = list(unique.values())
print(f"[OK] 去重后剩余 {len(records)} 条记录")
if duplicate_count > 0:
    print(f"[INFO] 发现并处理了 {duplicate_count} 个重复 ID\n")
else:
    print(f"[OK] 没有重复记录\n")

# ============================================================
# 第5步：写入 TSV 文件
# ============================================================
print("[STEP 3] 写入输出文件...")

with open(OUTPUT_TSV, 'w') as f:
    # 写表头
    f.write("mirbase_id\tmimat_id\tsequence\tsequence_len\n")
    
    # 写数据（按 ID 排序）
    for rec in sorted(records, key=lambda x: x['mirbase_id']):
        f.write("{}\t{}\t{}\t{}\n".format(
            rec['mirbase_id'],
            rec['mimat_id'],
            rec['sequence'],
            rec['sequence_len']
        ))

print(f"[OK] 已写入 {len(records)} 条记录到: {OUTPUT_TSV}\n")

# ============================================================
# 第6步：质量检查
# ============================================================
print("=" * 60)
print("质量检查报告")
print("=" * 60)

# 检查1: 显示前5个样例
print("\n[样例] 前5条记录:")
for rec in list(records)[:5]:
    print("  {:25s} {:15s} {}... (长度={})".format(
        rec['mirbase_id'],
        rec['mimat_id'],
        rec['sequence'][:30],
        rec['sequence_len']
    ))

# 检查2: 序列字符合法性
print("\n[检查] 序列字符合法性:")
illegal = set()
for rec in records:
    chars = set(rec['sequence']) - {'A', 'C', 'G', 'U', 'N'}
    illegal.update(chars)

if illegal:
    print(f"  ❌ 发现非法字符: {illegal}")
else:
    print(f"  ✅ 通过（所有序列只含 A/C/G/U/N）")

# 检查3: 数量合理性
print("\n[检查] 记录数量:")
print(f"  总计: {len(records)} 条")
if len(records) < 2000:
    print(f"  ❌ 数量过少（预期 2000-3500）")
    sys.exit(1)
elif len(records) > 3500:
    print(f"  ❌ 数量过多（预期 2000-3500）")
    sys.exit(1)
else:
    print(f"  ✅ 数量合理（预期范围: 2000-3500）")

# 检查4: 序列长度分布
lengths = [rec['sequence_len'] for rec in records]
print("\n[统计] 序列长度分布:")
print(f"  最短: {min(lengths)} nt")
print(f"  最长: {max(lengths)} nt")
print(f"  平均: {sum(lengths)/len(lengths):.1f} nt")

# ============================================================
# 第7步：最终验收
# ============================================================
print("\n" + "=" * 60)
print("✅✅ A2 步骤完成！")
print("=" * 60)
print(f"输出文件: {OUTPUT_TSV}")
print(f"记录数量: {len(records)}")
print(f"\n下一步: 运行 A3 脚本（映射到 RNAcentral URS）")
