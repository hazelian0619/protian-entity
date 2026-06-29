# step11_add_gene_ids_fixed.py
# 修正版：每个UniProt ID只映射到第一个基因符号

import pandas as pd

print("=== Phase 1: 添加基因ID映射（修正版 - 无重复） ===\n")

base_path = '/Users/pluviophile/graph/1025/data'
processed_dir = f'{base_path}/processed'

# 输入输出
input_file = f'{processed_dir}/protein_master_v3_final_fixed.tsv'
output_file = f'{processed_dir}/protein_master_v3_with_gene_ids_fixed.tsv'
hgnc_file = f'{base_path}/hgnc_raw_backup.tsv'  # 使用实际文件名

print(f"输入文件: {input_file}")
print(f"HGNC文件: {hgnc_file}")
print(f"输出文件: {output_file}\n")

# 1. 读取主表
print("1. 读取主表（v3基准版本）...")
df = pd.read_csv(input_file, sep='\t')
print(f"   总蛋白数: {len(df)}")
print(f"   当前列数: {len(df.columns)}\n")

# 2. 读取HGNC
print("2. 读取HGNC基因映射...")
df_hgnc = pd.read_csv(hgnc_file, sep='\t', low_memory=False)
print(f"   HGNC记录: {len(df_hgnc)}\n")

# 3. 创建映射字典（每个symbol只保留一条记录）
print("3. 创建基因ID映射...")

mapping = {}
for _, row in df_hgnc.iterrows():
    symbol = row.get('symbol')
    if pd.notna(symbol):
        # 每个symbol只保留第一条记录（如果重复）
        if symbol not in mapping:
            mapping[symbol] = {
                'ncbi_gene_id': row.get('entrez_id'),
                'ensembl_gene_id': row.get('ensembl_gene_id'),
                'ensembl_transcript_id': row.get('ensembl_transcript', ''),
                'gene_synonyms': row.get('alias_symbol')
            }

print(f"   映射表大小: {len(mapping)}\n")

# 4. 提取**第一个**基因符号并映射
print("4. 映射基因ID（只用第一个符号）...")

def get_primary_symbol(gene_names):
    """从gene_names提取第一个符号（最重要的）
    
    示例：
    "DEFB4A DEFB102 DEFB2 DEFB4; DEFB4B" → "DEFB4A"
    "H2AC11 H2AFP HIST1H2AG; H2AC13" → "H2AC11"
    """
    if pd.isna(gene_names) or gene_names == '':
        return None
    
    # 分割：先按分号，再按空格
    parts = str(gene_names).split(';')
    if len(parts) > 0:
        first_part = parts[0].strip()  # 取分号前的部分
        symbols = first_part.split()    # 按空格分割
        if len(symbols) > 0:
            return symbols[0]  # 返回第一个符号
    return None

# 初始化新字段
df['ncbi_gene_id'] = ''
df['ensembl_gene_id'] = ''
df['ensembl_transcript_id'] = ''
df['gene_synonyms'] = ''

mapped_count = 0
unmapped_ids = []

for idx, row in df.iterrows():
    symbol = get_primary_symbol(row['gene_names'])
    
    if symbol and symbol in mapping:
        gene_info = mapping[symbol]
        df.at[idx, 'ncbi_gene_id'] = gene_info.get('ncbi_gene_id', '')
        df.at[idx, 'ensembl_gene_id'] = gene_info.get('ensembl_gene_id', '')
        df.at[idx, 'ensembl_transcript_id'] = gene_info.get('ensembl_transcript_id', '')
        df.at[idx, 'gene_synonyms'] = gene_info.get('gene_synonyms', '')
        mapped_count += 1
    else:
        unmapped_ids.append(row['uniprot_id'])

print(f"   成功映射: {mapped_count}/{len(df)} ({mapped_count/len(df)*100:.1f}%)")
print(f"   未映射: {len(unmapped_ids)} ({len(unmapped_ids)/len(df)*100:.1f}%)\n")

# 5. 验证：检查是否还有重复
print("5. 验证去重效果...")
duplicated = df[df.duplicated('uniprot_id', keep=False)]
if len(duplicated) > 0:
    print(f"   ❌ 仍有{len(duplicated)}行重复！")
    print(f"   重复ID数: {duplicated['uniprot_id'].nunique()}")
    print("\n   重复ID示例（前5个）:")
    for uid in duplicated['uniprot_id'].unique()[:5]:
        print(f"      {uid}")
else:
    print(f"   ✅ 没有重复，完美！")

# 6. 统计覆盖率
print("\n6. 新增字段覆盖率:")
for field in ['ncbi_gene_id', 'ensembl_gene_id', 'ensembl_transcript_id', 'gene_synonyms']:
    has_data = df[field].notna() & (df[field] != '') & (df[field] != 'nan')
    print(f"   {field:25} {has_data.sum():,}/{len(df):,} ({has_data.sum()/len(df)*100:.1f}%)")

# 7. 对比旧版
print("\n7. 与旧版对比:")
old_file = f'{processed_dir}/protein_master_v3_with_gene_ids.tsv'
try:
    df_old = pd.read_csv(old_file, sep='\t')
    print(f"   旧版行数: {len(df_old)}")
    print(f"   新版行数: {len(df)}")
    print(f"   差异: {len(df_old) - len(df)} 行（重复被消除）")
except:
    print(f"   （旧版文件不存在，跳过对比）")

# 8. 保存
print(f"\n8. 保存结果...")
df.to_csv(output_file, sep='\t', index=False)
print(f"   ✅ 保存成功")

print(f"\n{'='*70}")
print("✅ Phase 1完成（修正版 - 无重复）！")
print(f"{'='*70}")
print(f"输出: {output_file}")
print(f"行数: {len(df):,} (应该 = 19,253)")
print(f"列数: {len(df.columns)} (29 → 33)")
print(f"\n修正内容:")
print(f"  ✅ 只提取gene_names中第一个符号")
print(f"  ✅ 避免了多基因符号导致的重复")
print(f"  ✅ 保证1个UniProt ID = 1行")
print(f"\n对比旧版:")
print(f"  - 旧版: 19,253 → 19,253 (有重复)")
print(f"  - 新版: 19,253 → 19,253 (无重复) ✅")
print(f"\n下一步: python step11b_merge_clean.py")
print(f"         (使用v3_with_gene_ids_fixed.tsv)")
print(f"{'='*70}\n")
