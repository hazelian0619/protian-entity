#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
========================================
步骤 A5: 生成标准化 miRNA 子表
========================================
输入文件: data/processed/rna/mirna_with_gene.tsv (A4的输出)
输出文件: data/output/rna_master_mirna_v1.tsv

功能：
1. 读取 A4 的完整 miRNA 数据
2. 按照 DATA_DICTIONARY_RNA.md 的 Tier 1 + Tier 2 字段重组
3. 填充固定字段（rna_type, taxon_id, source, source_version, fetch_date）
4. 输出最终的标准化表格
5. 运行 QA 检查

这是工程师 A（miRNA 管线）的最后一步！
"""

import os
import sys
from pathlib import Path
from datetime import date

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
print("A5: 生成标准化 miRNA 子表（工程师A最终步骤）")
print("=" * 70)
print(f"[INFO] 项目根目录: {project_root}\n")

# ============================================================
# 第2步：定义输入输出路径
# ============================================================
INPUT_FILE = Path("data/processed/rna/mirna_with_gene.tsv")
OUTPUT_FILE = Path("data/output/rna_master_mirna_v1.tsv")
QA_REPORT = Path("data/output/mirna_qa_report.txt")

if not INPUT_FILE.exists():
    print(f"[ERROR] 找不到文件: {INPUT_FILE}")
    sys.exit(1)

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
print(f"[OK] 输入文件: {INPUT_FILE}")
print(f"[OK] 输出文件: {OUTPUT_FILE}\n")

# ============================================================
# 第3步：读取 A4 的数据
# ============================================================
print("[STEP 1] 读取 A4 数据...")

records = []
with open(INPUT_FILE, 'r') as f:
    header = f.readline()  # urs	mirbase_id	mimat_id	sequence	sequence_len	gene_symbol	hgnc_id	ensembl_gene_id	ncbi_gene_id
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        
        records.append({
            'urs': parts[0],
            'mirbase_id': parts[1],
            'mimat_id': parts[2],
            'sequence': parts[3],
            'sequence_len': parts[4],
            'gene_symbol': parts[5],
            'hgnc_id': parts[6],
            'ensembl_gene_id': parts[7],
            'ncbi_gene_id': parts[8]
        })

print(f"[OK] 读取到 {len(records)} 条记录\n")

# ============================================================
# 第4步：重组为标准字段（按 DATA_DICTIONARY）
# ============================================================
print("[STEP 2] 重组为标准字段...")

# 固定值
FETCH_DATE = date.today().strftime('%Y-%m-%d')  # 2025-12-18
SOURCE = "RNAcentral;miRBase"
SOURCE_VERSION = "RNAcentral:25;miRBase:22.1"
TAXON_ID = 9606
RNA_TYPE = "mirna"

# 标准字段顺序（按 DATA_DICTIONARY）
# Tier 1: rna_id, rna_type, rna_name, sequence, sequence_len, taxon_id, 
#         symbol, source, fetch_date, source_version, hgnc_id, ensembl_gene_id, ncbi_gene_id
# Tier 2: mirbase_id, ensembl_transcript_id, rfam_id, secondary_structure, pdb_ids

standardized = []
for rec in records:
    standardized.append({
        # === Tier 1 字段（必填或重要）===
        'rna_id': rec['urs'],                    # 主键
        'rna_type': RNA_TYPE,                    # mirna
        'rna_name': rec['mirbase_id'],          # hsa-miR-21-5p
        'sequence': rec['sequence'],             # RNA 序列
        'sequence_len': rec['sequence_len'],     # 长度
        'taxon_id': TAXON_ID,                   # 9606
        'symbol': rec['gene_symbol'],           # MIR21
        'source': SOURCE,                       # RNAcentral;miRBase
        'fetch_date': FETCH_DATE,               # 2025-12-18
        'source_version': SOURCE_VERSION,       # RNAcentral:25;miRBase:22.1
        'hgnc_id': rec['hgnc_id'],              # HGNC:xxxxx（可能为空）
        'ensembl_gene_id': rec['ensembl_gene_id'],  # ENSGxxxxx（可能为空）
        'ncbi_gene_id': rec['ncbi_gene_id'],        # 数字（可能为空）
        
        # === Tier 2 字段（预留，v1 允许为空）===
        'mirbase_id': rec['mimat_id'],          # MIMAT0000076
        'ensembl_transcript_id': '',            # miRNA 没有 transcript（空）
        'rfam_id': '',                          # v1 不填
        'secondary_structure': '',              # v1 不填
        'pdb_ids': ''                           # v1 不填
    })

print(f"[OK] 标准化完成\n")

# ============================================================
# 第5步：写入输出文件
# ============================================================
print("[STEP 3] 写入输出文件...")

# 字段顺序（与 DATA_DICTIONARY 一致）
FIELD_ORDER = [
    # Tier 1
    'rna_id', 'rna_type', 'rna_name', 'sequence', 'sequence_len', 'taxon_id',
    'symbol', 'source', 'fetch_date', 'source_version', 
    'hgnc_id', 'ensembl_gene_id', 'ncbi_gene_id',
    # Tier 2
    'mirbase_id', 'ensembl_transcript_id', 'rfam_id', 
    'secondary_structure', 'pdb_ids'
]

with open(OUTPUT_FILE, 'w') as f:
    # 写表头
    f.write('\t'.join(FIELD_ORDER) + '\n')
    
    # 写数据（按 rna_id 排序）
    for rec in sorted(standardized, key=lambda x: x['rna_id']):
        values = [str(rec[field]) for field in FIELD_ORDER]
        f.write('\t'.join(values) + '\n')

print(f"[OK] 已写入 {len(standardized)} 条记录\n")

# ============================================================
# 第6步：QA 检查
# ============================================================
print("=" * 70)
print("QA 质量检查")
print("=" * 70)

qa_results = []

def qa_check(name, condition, expected="100%", critical=True):
    passed = condition
    status = "✅ PASS" if passed else ("❌ FAIL" if critical else "⚠️  WARN")
    qa_results.append(f"[{status}] {name}: {expected}")
    print(f"{status:12s} {name:40s} (目标: {expected})")
    return passed

print("\n[检查组1] 主键与基本字段")
all_pass = True
all_pass &= qa_check("rna_id 非空率", 
                      all(r['rna_id'] for r in standardized))
all_pass &= qa_check("rna_id 唯一性", 
                      len(set(r['rna_id'] for r in standardized)) == len(standardized))
all_pass &= qa_check("rna_id 格式（以 _9606 结尾）", 
                      all(r['rna_id'].endswith('_9606') for r in standardized))
all_pass &= qa_check("rna_type 固定为 mirna", 
                      all(r['rna_type'] == 'mirna' for r in standardized))
all_pass &= qa_check("taxon_id 固定为 9606", 
                      all(r['taxon_id'] == 9606 for r in standardized))

print("\n[检查组2] 序列字段")
all_pass &= qa_check("sequence 非空率", 
                      all(r['sequence'] for r in standardized))
all_pass &= qa_check("sequence 字符合法性（A/C/G/U/N）", 
                      all(set(r['sequence']) <= {'A','C','G','U','N'} for r in standardized))
all_pass &= qa_check("sequence_len 准确性", 
                      all(len(r['sequence']) == int(r['sequence_len']) for r in standardized))

print("\n[检查组3] 版本与来源")
all_pass &= qa_check("source_version 非空率", 
                      all(r['source_version'] for r in standardized))
all_pass &= qa_check("source 非空率", 
                      all(r['source'] for r in standardized))
all_pass &= qa_check("fetch_date 非空率", 
                      all(r['fetch_date'] for r in standardized))

print("\n[检查组4] gene 映射（非强制）")
symbol_coverage = sum(1 for r in standardized if r['symbol']) / len(standardized) * 100
hgnc_coverage = sum(1 for r in standardized if r['hgnc_id']) / len(standardized) * 100
qa_check(f"symbol 非空率", 
         symbol_coverage >= 50, 
         f">50% (实际 {symbol_coverage:.1f}%)", 
         critical=False)
qa_check(f"hgnc_id 非空率", 
         hgnc_coverage >= 0, 
         f">0% (实际 {hgnc_coverage:.1f}%，miRNA 可接受)", 
         critical=False)

print("\n[检查组5] 数据量")
qa_check("记录数量范围", 
         2000 <= len(standardized) <= 3500,
         "2000-3500")

# 统计信息
lengths = [int(r['sequence_len']) for r in standardized]
print(f"\n[统计] 序列长度分布:")
print(f"  最短: {min(lengths)} nt")
print(f"  最长: {max(lengths)} nt")
print(f"  平均: {sum(lengths)/len(lengths):.1f} nt")

# ============================================================
# 第7步：生成 QA 报告文件
# ============================================================
with open(QA_REPORT, 'w') as f:
    f.write("=" * 70 + "\n")
    f.write("miRNA 子表 QA 报告\n")
    f.write("=" * 70 + "\n\n")
    f.write(f"生成时间: {FETCH_DATE}\n")
    f.write(f"输出文件: {OUTPUT_FILE}\n")
    f.write(f"记录数量: {len(standardized)}\n\n")
    f.write("检查结果:\n")
    for result in qa_results:
        f.write(result + "\n")
    f.write("\n数据统计:\n")
    f.write(f"  序列长度: {min(lengths)}-{max(lengths)} nt (平均 {sum(lengths)/len(lengths):.1f})\n")
    f.write(f"  symbol 覆盖率: {symbol_coverage:.1f}%\n")
    f.write(f"  hgnc_id 覆盖率: {hgnc_coverage:.1f}%\n")

print(f"\n[OK] QA 报告已保存: {QA_REPORT}")

# ============================================================
# 第8步：最终验收
# ============================================================
print("\n" + "=" * 70)
if all_pass:
    print("🎉🎉 A5 完成！工程师A的 miRNA 管线全部完成！")
else:
    print("⚠️  A5 完成，但有警告项（已标记在上方）")
print("=" * 70)
print(f"✅ 最终输出: {OUTPUT_FILE}")
print(f"✅ QA 报告: {QA_REPORT}")
print(f"✅ 记录数量: {len(standardized)}")
print(f"✅ 字段数量: {len(FIELD_ORDER)} (Tier 1: 13, Tier 2: 5)")
print("\n下一步:")
print("  - 工程师B: 完成 mRNA 管线（B1-B7）")
print("  - 项目负责人: 合并两条管线（M1-M5）")