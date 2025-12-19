#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
B5: 标准化 mRNA 子表（按 DATA_DICTIONARY）

输入: data/processed/rna/mrna_with_sequence.tsv
输出: data/output/rna_master_mrna_v1.tsv
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
print("B5: 标准化 mRNA 子表")
print("=" * 70)

INPUT_FILE = Path("data/processed/rna/mrna_with_sequence.tsv")
OUTPUT_FILE = Path("data/output/rna_master_mrna_v1.tsv")
QA_REPORT = Path("data/output/mrna_qa_report.txt")

if not INPUT_FILE.exists():
    print(f"[ERROR] 找不到: {INPUT_FILE}")
    sys.exit(1)

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
print(f"[OK] 输入: {INPUT_FILE}")
print(f"[OK] 输出: {OUTPUT_FILE}\n")

# 固定值
FETCH_DATE = date.today().strftime('%Y-%m-%d')
SOURCE = "RNAcentral;Ensembl"
SOURCE_VERSION = "RNAcentral:25;Ensembl:113"
TAXON_ID = 9606
RNA_TYPE = "mrna"

# ============================================================
# 读取数据
# ============================================================
print("[STEP 1] 读取数据...")

records = []
with open(INPUT_FILE, 'r') as f:
    header = f.readline().strip().split('\t')
    for line_num, line in enumerate(f, 1):
        parts = line.strip().split('\t')
        if len(parts) < 8:
            continue
        
        # 修复 ncbi_gene_id（从浮点数转整数）
        ncbi_gene_id = parts[5]
        if ncbi_gene_id and ncbi_gene_id != '':
            try:
                ncbi_gene_id = str(int(float(ncbi_gene_id)))
            except:
                ncbi_gene_id = ''
        
        records.append({
            'rna_id': parts[0],
            'ensembl_transcript_id': parts[1],
            'ensembl_gene_id': parts[2],
            'gene_symbol': parts[3],
            'hgnc_id': parts[4],
            'ncbi_gene_id': ncbi_gene_id,
            'sequence': parts[6],
            'sequence_len': int(parts[7])
        })

print(f"[OK] 读取到 {len(records)} 条记录\n")

# ============================================================
# 标准化字段
# ============================================================
print("[STEP 2] 标准化字段...")

standardized = []
for rec in records:
    standardized.append({
        # Tier 1
        'rna_id': rec['rna_id'],
        'rna_type': RNA_TYPE,
        'rna_name': rec['ensembl_transcript_id'],  # 用 ENST 作为 name
        'sequence': rec['sequence'],
        'sequence_len': rec['sequence_len'],
        'taxon_id': TAXON_ID,
        'symbol': rec['gene_symbol'],
        'source': SOURCE,
        'fetch_date': FETCH_DATE,
        'source_version': SOURCE_VERSION,
        'hgnc_id': rec['hgnc_id'],
        'ensembl_gene_id': rec['ensembl_gene_id'],
        'ncbi_gene_id': rec['ncbi_gene_id'],
        # Tier 2
        'mirbase_id': '',  # mRNA 没有 mirbase_id
        'ensembl_transcript_id': rec['ensembl_transcript_id'],
        'rfam_id': '',
        'secondary_structure': '',
        'pdb_ids': ''
    })

print(f"[OK] 标准化完成\n")

# ============================================================
# 写入输出
# ============================================================
print("[STEP 3] 写入输出...")

FIELDS = [
    'rna_id', 'rna_type', 'rna_name', 'sequence', 'sequence_len', 'taxon_id',
    'symbol', 'source', 'fetch_date', 'source_version',
    'hgnc_id', 'ensembl_gene_id', 'ncbi_gene_id',
    'mirbase_id', 'ensembl_transcript_id', 'rfam_id',
    'secondary_structure', 'pdb_ids'
]

with open(OUTPUT_FILE, 'w') as f:
    f.write('\t'.join(FIELDS) + '\n')
    for rec in sorted(standardized, key=lambda x: x['rna_id']):
        f.write('\t'.join(str(rec[k]) for k in FIELDS) + '\n')

print(f"[OK] 已写入 {len(standardized)} 条记录\n")

# ============================================================
# QA 检查
# ============================================================
print("=" * 70)
print("QA 检查")
print("=" * 70)

qa_results = []
all_pass = True

def qa_check(name, condition, expected="100%"):
    global all_pass
    passed = condition
    status = "✅" if passed else "❌"
    qa_results.append(f"[{status}] {name}: {expected}")
    print(f"{status} {name:50s} ({expected})")
    if not passed:
        all_pass = False
    return passed

print("\n[检查组1] 主键")
qa_check("rna_id 非空率", all(r['rna_id'] for r in standardized))
qa_check("rna_id 唯一性", len(set(r['rna_id'] for r in standardized)) == len(standardized))
qa_check("rna_id 格式", all(r['rna_id'].endswith('_9606') for r in standardized))

print("\n[检查组2] 序列")
qa_check("sequence 非空率", all(r['sequence'] for r in standardized))
qa_check("sequence 字符合法", all(set(r['sequence']) <= {'A','C','G','U','N'} for r in standardized))
qa_check("sequence_len 准确", all(len(r['sequence']) == r['sequence_len'] for r in standardized))

print("\n[检查组3] 版本")
qa_check("source_version 非空", all(r['source_version'] for r in standardized))
qa_check("fetch_date 非空", all(r['fetch_date'] for r in standardized))

print("\n[检查组4] gene 映射")
symbol_pct = sum(1 for r in standardized if r['symbol']) / len(standardized) * 100
ensembl_pct = sum(1 for r in standardized if r['ensembl_gene_id']) / len(standardized) * 100
qa_check(f"symbol 非空率", symbol_pct >= 99, f">99% (实际{symbol_pct:.1f}%)")
qa_check(f"ensembl_gene_id 非空率", ensembl_pct >= 99, f">99% (实际{ensembl_pct:.1f}%)")

print("\n[检查组5] 数据量")
unique_genes = len(set(r['ensembl_gene_id'] for r in standardized if r['ensembl_gene_id']))
qa_check("记录数量", len(standardized) >= 100000, f">100k (实际{len(standardized)})")
qa_check("唯一基因数", unique_genes >= 18000, f">18k (实际{unique_genes})")

# 统计
lengths = [r['sequence_len'] for r in standardized]
print(f"\n[统计]")
print(f"  转录本总数: {len(standardized)}")
print(f"  唯一基因数: {unique_genes}")
print(f"  平均每基因: {len(standardized) / unique_genes:.1f} 个转录本")
print(f"  序列长度: {min(lengths)}-{max(lengths)} nt (平均 {sum(lengths)/len(lengths):.0f})")

# 保存 QA 报告
with open(QA_REPORT, 'w') as f:
    f.write("mRNA 子表 QA 报告\n")
    f.write("=" * 70 + "\n")
    f.write(f"生成时间: {FETCH_DATE}\n")
    f.write(f"记录数量: {len(standardized)}\n\n")
    for r in qa_results:
        f.write(r + "\n")

print(f"\n[OK] QA 报告: {QA_REPORT}")

print("\n" + "=" * 70)
if all_pass:
    print("🎉🎉 B5 完成！mRNA 管线全部完成！")
else:
    print("⚠️  B5 完成但有警告")
print("=" * 70)
print(f"✅ 输出: {OUTPUT_FILE}")
print(f"✅ 记录: {len(standardized)} 条")
print(f"✅ 字段: {len(FIELDS)} 个")
print("\n下一步: M1 合并 miRNA + mRNA")