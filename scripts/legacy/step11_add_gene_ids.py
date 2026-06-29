# step11_add_gene_ids.py
# 位置：/Users/pluviophile/graph/1025/scripts/step11_add_gene_ids.py

import pandas as pd

print("=== Phase 1: 添加基因ID映射 ===\n")

# 路径配置
base_path = '/Users/pluviophile/graph/1025/data'
processed_dir = f'{base_path}/processed'

# 输入：干净的v3
input_file = f'{processed_dir}/protein_master_v3_final_fixed.tsv'
# 输出：带基因ID的v3
output_file = f'{processed_dir}/protein_master_v3_with_gene_ids.tsv'

# HGNC数据
hgnc_file = f'{processed_dir}/hgnc_core.tsv'

# 1. 读取主表
print("1. 读取主表（v3基准版本）...")
df = pd.read_csv(input_file, sep='\t')
print(f"   总蛋白数: {len(df)}")
print(f"   当前列数: {len(df.columns)}\n")

# 2. 读取HGNC
print("2. 读取HGNC基因映射...")
df_hgnc = pd.read_csv(hgnc_file, sep='\t')
print(f"   HGNC记录: {len(df_hgnc)}\n")

# 3. 创建映射字典
print("3. 创建基因ID映射...")

# 从symbol映射到IDs
mapping = {}
for _, row in df_hgnc.iterrows():
    symbol = row.get('symbol')
    if pd.notna(symbol):
        mapping[symbol] = {
            'ncbi_gene_id': row.get('entrez_id'),
            'ensembl_gene_id': row.get('ensembl_gene_id'),
            'ensembl_transcript_id': row.get('ensembl_transcript'),  # 如果有
            'gene_synonyms': row.get('alias_symbol')
        }

print(f"   映射表大小: {len(mapping)}\n")

# 4. 提取symbol并映射
print("4. 映射基因ID...")

def get_primary_symbol(gene_names):
    """从gene_names提取primary symbol（第一个）"""
    if pd.isna(gene_names) or gene_names == '':
        return None
    return str(gene_names).split()[0]

df['primary_symbol'] = df['gene_names'].apply(get_primary_symbol)

# 映射
mapped_count = 0
for idx, row in df.iterrows():
    symbol = row['primary_symbol']
    if symbol and symbol in mapping:
        gene_info = mapping[symbol]
        df.at[idx, 'ncbi_gene_id'] = gene_info.get('ncbi_gene_id')
        df.at[idx, 'ensembl_gene_id'] = gene_info.get('ensembl_gene_id')
        df.at[idx, 'ensembl_transcript_id'] = gene_info.get('ensembl_transcript_id')
        df.at[idx, 'gene_synonyms'] = gene_info.get('gene_synonyms')
        mapped_count += 1

# 删除临时列
df = df.drop(columns=['primary_symbol'])

print(f"   成功映射: {mapped_count}/{len(df)} ({mapped_count/len(df)*100:.1f}%)\n")

# 5. 统计
print("5. 新增字段覆盖率:")
for field in ['ncbi_gene_id', 'ensembl_gene_id', 'ensembl_transcript_id', 'gene_synonyms']:
    has_data = df[field].notna().sum()
    print(f"   {field:25} {has_data:,}/{len(df):,} ({has_data/len(df)*100:.1f}%)")

# 6. 保存
print(f"\n6. 保存结果...")
df.to_csv(output_file, sep='\t', index=False)

print(f"\n{'='*60}")
print("✅ Phase 1完成！")
print(f"{'='*60}")
print(f"输入: {input_file}")
print(f"输出: {output_file}")
print(f"列数: {len(df.columns)} (29 → 33)")
print(f"\n新增字段:")
print(f"  - ncbi_gene_id")
print(f"  - ensembl_gene_id")
print(f"  - ensembl_transcript_id")
print(f"  - gene_synonyms")
print(f"\n下一步: 等待step10完成后，运行step11b_merge_clean.py")
