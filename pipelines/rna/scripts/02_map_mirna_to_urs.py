#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
========================================
步骤 A3: 将 miRNA 映射到 RNAcentral URS
========================================
输入文件1: data/processed/rna/hsa_mirna_list.tsv (A2的输出)
输入文件2: data/raw/rna/rnacentral/id_mapping.tsv (RNAcentral xref)
输出文件: data/processed/rna/mirna_urs_mapping.tsv

功能：
1. 读取 A2 提取的人类 miRNA 列表（含 MIMAT ID）
2. 从 RNAcentral xref 中找到 MIMAT -> URS 的映射
3. 处理一对多关系（一个 MIMAT 可能对应多个 URS，选择最合适的）
4. 输出带 URS 的 miRNA 表
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
print("A3: 映射 miRNA 到 RNAcentral URS")
print("=" * 70)
print(f"[INFO] 项目根目录: {project_root}")
print(f"[INFO] 当前工作目录: {os.getcwd()}\n")

# ============================================================
# 第2步：定义输入输出路径
# ============================================================
INPUT_MIRNA = Path("data/processed/rna/hsa_mirna_list.tsv")
INPUT_XREF = Path("data/raw/rna/rnacentral/id_mapping.tsv")
OUTPUT_TSV = Path("data/processed/rna/mirna_urs_mapping.tsv")

# 检查文件
for fpath in [INPUT_MIRNA, INPUT_XREF]:
    if not fpath.exists():
        print(f"[ERROR] 找不到文件: {fpath}")
        sys.exit(1)
    print(f"[OK] 输入文件存在: {fpath}")

OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
print(f"[OK] 输出路径: {OUTPUT_TSV}\n")

# ============================================================
# 第3步：读取 A2 的 miRNA 列表
# ============================================================
print("[STEP 1] 读取 miRNA 列表...")

mirna_dict = {}  # {mimat_id: {'mirbase_id': ..., 'sequence': ..., 'sequence_len': ...}}

with open(INPUT_MIRNA, 'r') as f:
    header = f.readline()  # 跳过表头
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 4:
            continue
        mirbase_id, mimat_id, sequence, sequence_len = parts
        mirna_dict[mimat_id] = {
            'mirbase_id': mirbase_id,
            'sequence': sequence,
            'sequence_len': int(sequence_len)
        }

print(f"[OK] 读取到 {len(mirna_dict)} 个 miRNA（按 MIMAT ID 索引）")
print(f"[样例] {list(mirna_dict.keys())[:3]}\n")

# ============================================================
# 第4步：从 RNAcentral xref 提取 MIMAT -> URS 映射
# ============================================================
print("[STEP 2] 从 RNAcentral xref 提取 MIRBASE 映射...")
print("[INFO] 这一步可能需要 1-2 分钟（文件 11GB）...\n")

# id_mapping.tsv 格式（6列）：
# URS0000000001   ENA     GU786683.1:1..200:rRNA  77133   rRNA    ""
# 列：URS, database, external_id, taxon_id, rna_type, ...

mimat_to_urs = defaultdict(list)  # {mimat_id: [(urs, rna_type), ...]}
line_count = 0
mirbase_count = 0

with open(INPUT_XREF, 'r') as f:
    for line in f:
        line_count += 1
        
        # 每处理 100 万行显示进度
        if line_count % 1000000 == 0:
            print(f"  已处理 {line_count // 1000000} 百万行...")
        
        parts = line.strip().split('\t')
        if len(parts) < 5:
            continue
        
        urs_raw = parts[0]
        urs = urs_raw if urs_raw.endswith('_9606') else urs_raw + '_9606'
        database = parts[1].upper()
        external_id = parts[2]
        taxon_id = parts[3]
        rna_type = parts[4]
        
        # 只保留：数据库=MIRBASE 且 物种=9606（人类）
        if database == 'MIRBASE' and taxon_id == '9606':
            mirbase_count += 1
            # external_id 可能是 MIMAT0000062 或带其他信息
            # 提取 MIMAT 部分（前缀是 MIMAT）
            mimat_id = external_id.split(':')[0] if ':' in external_id else external_id
            
            if mimat_id in mirna_dict:
                mimat_to_urs[mimat_id].append((urs, rna_type))

print(f"[OK] 处理完成，总计 {line_count} 行")
print(f"[OK] 找到 {mirbase_count} 个人类 MIRBASE 映射")
print(f"[OK] 其中 {len(mimat_to_urs)} 个 MIMAT 能匹配上我们的 miRNA 列表\n")

# ============================================================
# 第5步：处理一对多映射（选择策略）
# ============================================================
print("[STEP 3] 处理一对多映射...")

# 统计映射情况
one_to_one = 0    # 一个 MIMAT 对应一个 URS
one_to_many = 0   # 一个 MIMAT 对应多个 URS
no_mapping = 0    # 没有找到 URS

final_mapping = {}  # {mimat_id: urs}

for mimat_id in mirna_dict:
    if mimat_id not in mimat_to_urs:
        # 没有映射
        no_mapping += 1
        continue
    
    urs_list = mimat_to_urs[mimat_id]
    
    if len(urs_list) == 1:
        # 一对一映射
        one_to_one += 1
        final_mapping[mimat_id] = urs_list[0][0]
    else:
        # 一对多映射：选择第一个（通常是最常用的）
        # 也可以选择 URS 最长的，或者标记为 canonical 的
        one_to_many += 1
        final_mapping[mimat_id] = urs_list[0][0]
        
        # 如果有多个，记录警告
        if len(urs_list) > 2:
            print(f"  [WARN] {mimat_id} 有 {len(urs_list)} 个 URS，选择第一个: {urs_list[0][0]}")

print(f"[统计] 一对一映射: {one_to_one}")
print(f"[统计] 一对多映射: {one_to_many}")
print(f"[统计] 无映射: {no_mapping}")
print(f"[OK] 最终映射数量: {len(final_mapping)}\n")

# ============================================================
# 第6步：生成输出文件
# ============================================================
print("[STEP 4] 生成输出文件...")

output_records = []

for mimat_id, urs in final_mapping.items():
    mirna_info = mirna_dict[mimat_id]
    output_records.append({
        'urs': urs,
        'mirbase_id': mirna_info['mirbase_id'],
        'mimat_id': mimat_id,
        'sequence': mirna_info['sequence'],
        'sequence_len': mirna_info['sequence_len']
    })

# 按 URS 排序
output_records.sort(key=lambda x: x['urs'])

with open(OUTPUT_TSV, 'w') as f:
    f.write("urs\tmirbase_id\tmimat_id\tsequence\tsequence_len\n")
    for rec in output_records:
        f.write("{}\t{}\t{}\t{}\t{}\n".format(
            rec['urs'],
            rec['mirbase_id'],
            rec['mimat_id'],
            rec['sequence'],
            rec['sequence_len']
        ))

print(f"[OK] 已写入 {len(output_records)} 条记录到: {OUTPUT_TSV}\n")

# ============================================================
# 第7步：质量检查
# ============================================================
print("=" * 70)
print("质量检查报告")
print("=" * 70)

# 检查1: 显示前5个样例
print("\n[样例] 前5条记录:")
for rec in output_records[:5]:
    print("  URS: {:<25s} miRNA: {:<20s} 序列: {}... (len={})".format(
        rec['urs'],
        rec['mirbase_id'],
        rec['sequence'][:30],
        rec['sequence_len']
    ))

# 检查2: URS 格式检查
print("\n[检查] URS 格式:")
invalid_urs = [rec['urs'] for rec in output_records if not rec['urs'].endswith('_9606')]
if invalid_urs:
    print(f"  ❌ 发现 {len(invalid_urs)} 个 URS 不以 _9606 结尾")
    print(f"  示例: {invalid_urs[:3]}")
else:
    print(f"  ✅ 所有 URS 都以 _9606 结尾")

# 检查3: 映射覆盖率
coverage = len(final_mapping) / len(mirna_dict) * 100
print(f"\n[统计] 映射覆盖率:")
print(f"  输入 miRNA 总数: {len(mirna_dict)}")
print(f"  成功映射数量: {len(final_mapping)}")
print(f"  覆盖率: {coverage:.1f}%")

if coverage < 85:
    print(f"  ⚠️  覆盖率偏低（建议 >90%），可能需要检查数据源版本")
elif coverage < 95:
    print(f"  ✅ 覆盖率良好（可接受）")
else:
    print(f"  ✅ 覆盖率优秀")

# 检查4: 未映射的 miRNA（记录到日志）
if no_mapping > 0:
    unmapped_file = Path("data/processed/rna/mirna_unmapped.txt")
    with open(unmapped_file, 'w') as f:
        f.write("# 未找到 URS 映射的 miRNA（共 {} 个）\n".format(no_mapping))
        for mimat_id, info in mirna_dict.items():
            if mimat_id not in final_mapping:
                f.write("{}\t{}\n".format(mimat_id, info['mirbase_id']))
    print(f"\n[INFO] 未映射的 miRNA 列表已保存到: {unmapped_file}")

# ============================================================
# 第8步：最终验收
# ============================================================
print("\n" + "=" * 70)
print("✅✅ A3 步骤完成！")
print("=" * 70)
print(f"输出文件: {OUTPUT_TSV}")
print(f"记录数量: {len(output_records)}")
print(f"映射覆盖率: {coverage:.1f}%")
print(f"\n下一步: 运行 A4 脚本（补充 gene 映射）")