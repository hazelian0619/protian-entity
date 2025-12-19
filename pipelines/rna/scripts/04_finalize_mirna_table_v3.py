#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
print("A5: 生成标准化 miRNA 子表")
print("=" * 70)

INPUT_FILE = Path("data/processed/rna/mirna_with_gene.tsv")
OUTPUT_FILE = Path("data/output/rna_master_mirna_v1.tsv")
QA_REPORT = Path("data/output/mirna_qa_report.txt")

if not INPUT_FILE.exists():
    print(f"[ERROR] 找不到: {INPUT_FILE}")
    sys.exit(1)

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

print(f"[OK] 输入: {INPUT_FILE}")
print(f"[OK] 输出: {OUTPUT_FILE}\n")

# ============================================================
# 读取数据（用空白字符分割，不限tab）
# ============================================================
print("[STEP 1] 读取数据...")

records = []
with open(INPUT_FILE, 'r') as f:
    header = f.readline()  # 跳过表头
    
    for line_num, line in enumerate(f, 1):
        # 用正则分割：任意空白字符（tab、空格、多个空格）
        parts = re.split(r'\s+', line.strip())
        
        # 预期：urs, mirbase_id, mimat_id, sequence, sequence_len, gene_symbol, hgnc_id, ensembl_gene_id, ncbi_gene_id
        if len(parts) < 6:
            if line_num <= 5:
                print(f"  [WARN] 第{line_num}行字段不足（{len(parts)}列）: {line.strip()[:60]}...")
            continue
        
        records.append({
            'urs': parts[0],
            'mirbase_id': parts[1],
            'mimat_id': parts[2],
            'sequence': parts[3],
            'sequence_len': parts[4],
            'gene_symbol': parts[5] if len(parts) > 5 else '',
            'hgnc_id': parts[6] if len(parts) > 6 else '',
            'ensembl_gene_id': parts[7] if len(parts) > 7 else '',
            'ncbi_gene_id': parts[8] if len(parts) > 8 else ''
        })

print(f"[OK] 成功读取 {len(records)} 条记录\n")

if len(records) == 0:
    print("[ERROR] 没有读取到任何记录")
    sys.exit(1)

# ============================================================
# 标准化字段
# ============================================================
print("[STEP 2] 标准化字段...")

FETCH_DATE = date.today().strftime('%Y-%m-%d')
SOURCE = "RNAcentral;miRBase"
SOURCE_VERSION = "RNAcentral:25;miRBase:22.1"

standardized = []
for rec in records:
    standardized.append({
        # Tier 1
        'rna_id': rec['urs'],
        'rna_type': 'mirna',
        'rna_name': rec['mirbase_id'],
        'sequence': rec['sequence'],
        'sequence_len': rec['sequence_len'],
        'taxon_id': '9606',
        'symbol': rec['gene_symbol'],
        'source': SOURCE,
        'fetch_date': FETCH_DATE,
        'source_version': SOURCE_VERSION,
        'hgnc_id': rec['hgnc_id'],
        'ensembl_gene_id': rec['ensembl_gene_id'],
        'ncbi_gene_id': rec['ncbi_gene_id'],
        # Tier 2
        'mirbase_id': rec['mimat_id'],
        'ensembl_transcript_id': '',
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
    'symbol', 'source', 'fetch_date', 'source_version', 'hgnc_id',
    'ensembl_gene_id', 'ncbi_gene_id', 'mirbase_id', 'ensembl_transcript_id',
    'rfam_id', 'secondary_structure', 'pdb_ids'
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
print("QA 质量检查")
print("=" * 70)

qa_results = []
all_pass = True

def qa_check(name, condition, expected="100%", critical=True):
    global all_pass
    passed = condition
    status = "✅ PASS" if passed else ("❌ FAIL" if critical else "⚠️  WARN")
    qa_results.append(f"[{status}] {name}: {expected}")
    print(f"{status:12s} {name:40s} ({expected})")
    if not passed and critical:
        all_pass = False
    return passed

print("\n[检查组1] 主键")
qa_check("rna_id 非空率", all(r['rna_id'] for r in standardized))
qa_check("rna_id 唯一性", len(set(r['rna_id'] for r in standardized)) == len(standardized))
qa_check("rna_id 格式（_9606结尾）", all(r['rna_id'].endswith('_9606') for r in standardized))

print("\n[检查组2] 序列")
qa_check("sequence 非空率", all(r['sequence'] for r in standardized))
qa_check("sequence 字符合法", all(set(r['sequence']) <= {'A','C','G','U','N'} for r in standardized))
qa_check("sequence_len 准确", all(len(r['sequence']) == int(r['sequence_len']) for r in standardized))

print("\n[检查组3] 版本")
qa_check("source_version 非空", all(r['source_version'] for r in standardized))
qa_check("fetch_date 非空", all(r['fetch_date'] for r in standardized))

print("\n[检查组4] gene映射（非强制）")
symbol_pct = sum(1 for r in standardized if r['symbol']) / len(standardized) * 100
qa_check(f"symbol非空率", symbol_pct >= 50, f">50% (实际{symbol_pct:.1f}%)", critical=False)

print("\n[检查组5] 数据量")
qa_check("记录数量", 2000 <= len(standardized) <= 3500, "2000-3500")

# 统计
lengths = [int(r['sequence_len']) for r in standardized]
print(f"\n[统计]")
print(f"  记录总数: {len(standardized)}")
print(f"  序列长度: {min(lengths)}-{max(lengths)} nt (平均 {sum(lengths)/len(lengths):.1f})")
print(f"  symbol覆盖: {symbol_pct:.1f}%")

# 保存QA报告
with open(QA_REPORT, 'w') as f:
    f.write("miRNA子表 QA报告\n")
    f.write("=" * 70 + "\n")
    f.write(f"生成时间: {FETCH_DATE}\n")
    f.write(f"记录数量: {len(standardized)}\n\n")
    for r in qa_results:
        f.write(r + "\n")

print(f"\n[OK] QA报告: {QA_REPORT}")

# ============================================================
# 最终验收
# ============================================================
print("\n" + "=" * 70)
if all_pass:
    print("🎉🎉 A5完成！工程师A的miRNA管线全部完成！")
else:
    print("⚠️  A5完成但有警告")
print("=" * 70)
print(f"✅ 输出: {OUTPUT_FILE}")
print(f"✅ 记录: {len(standardized)} 条")
print(f"✅ 字段: {len(FIELDS)} 个")
print("\n下一步: 工程师B完成mRNA管线 或 直接合并")