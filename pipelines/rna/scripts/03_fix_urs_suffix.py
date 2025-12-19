#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
快速修复：给 URS 加上 _9606 后缀
"""
import os
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
INPUT = Path("data/processed/rna/mirna_urs_mapping.tsv")
OUTPUT = Path("data/processed/rna/mirna_urs_mapping_fixed.tsv")

print("[修复] 给 URS 加上 _9606 后缀...")

with open(INPUT, 'r') as f_in, open(OUTPUT, 'w') as f_out:
    header = f_in.readline()
    f_out.write(header)  # 表头不变
    
    count = 0
    for line in f_in:
        parts = line.strip().split('\t')
        if len(parts) >= 5:
            urs = parts[0]
            # 如果 URS 不以 _9606 结尾，加上后缀
            if not urs.endswith('_9606'):
                parts[0] = urs + '_9606'
            f_out.write('\t'.join(parts) + '\n')
            count += 1

print(f"[OK] 已修复 {count} 条记录")
print(f"[OK] 输出: {OUTPUT}")

# 替换原文件
OUTPUT.replace(INPUT)
print(f"[OK] 已覆盖原文件: {INPUT}")

# 验证
with open(INPUT, 'r') as f:
    f.readline()  # 跳过表头
    first_line = f.readline()
    urs = first_line.split('\t')[0]
    if urs.endswith('_9606'):
        print(f"[验证] ✅ URS 格式正确: {urs}")
    else:
        print(f"[验证] ❌ URS 格式错误: {urs}")