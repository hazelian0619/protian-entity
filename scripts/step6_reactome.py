# step6_reactome.py - Reactome通路成员提取
import pandas as pd
from datetime import datetime

print("=== 步骤6: Reactome通路成员 ===")

# 路径
base_path = '/Users/pluviophile/graph/1025/data'
raw_dir = f'{base_path}/raw'
processed_dir = f'{base_path}/processed'
reactome_file = f'{raw_dir}/UniProt2Reactome_All_Levels.txt'
master_file = f'{processed_dir}/protein_master_v2.tsv'
output_file = f'{processed_dir}/pathway_members.tsv'

# 1. 加载主表（获取我们的蛋白ID列表）
print("1. 加载protein_master_v2.tsv...")
df_master = pd.read_csv(master_file, sep='\t')
our_proteins = set(df_master['uniprot_id'])
print(f"   主表: {len(our_proteins)}个蛋白")

# 2. 加载Reactome原始文件（6列）
print("2. 加载Reactome文件...")
df_reactome = pd.read_csv(
    reactome_file, 
    sep='\t',
    names=['uniprot_id', 'pathway_id', 'url', 'pathway_name', 'evidence_code', 'species'],
    dtype=str
)
print(f"   原始行数: {len(df_reactome):,} (全物种)")

# 3. 过滤人类通路
print("3. 过滤人类通路...")
df_human = df_reactome[df_reactome['species'] == 'Homo sapiens'].copy()
print(f"   人类通路: {len(df_human):,}行")

# 4. 只保留主表中的蛋白
print("4. 保留主表蛋白...")
df_filtered = df_human[df_human['uniprot_id'].isin(our_proteins)].copy()
print(f"   匹配主表: {len(df_filtered):,}行")

# 5. 去重（同一蛋白-通路对可能有多条证据）
print("5. 去重处理...")
df_unique = df_filtered.drop_duplicates(subset=['uniprot_id', 'pathway_id'])
print(f"   去重后: {len(df_unique):,}行")

# 6. 清理并添加元数据
print("6. 数据清理...")
df_output = df_unique[['uniprot_id', 'pathway_id', 'pathway_name', 'species']].copy()
df_output['source'] = 'Reactome'
df_output['fetch_date'] = datetime.now().strftime('%Y-%m-%d')

# 7. 输出
df_output.to_csv(output_file, sep='\t', index=False)
print(f"\n✅ 输出: {output_file}")
print(f"   总行数: {len(df_output):,}")
print(f"   唯一蛋白数: {df_output['uniprot_id'].nunique()}")
print(f"   唯一通路数: {df_output['pathway_id'].nunique()}")

# 8. 示例
print(f"\n   示例前3行:")
sample = df_output.head(3)[['uniprot_id', 'pathway_id', 'pathway_name']]
for _, row in sample.iterrows():
    name_short = row['pathway_name'][:40] + '...' if len(row['pathway_name']) > 40 else row['pathway_name']
    print(f"   {row['uniprot_id']} → {row['pathway_id']} ({name_short})")

# 9. 统计
avg_pathways = len(df_output) / df_output['uniprot_id'].nunique()
print(f"\n   平均每个蛋白参与 {avg_pathways:.1f} 条通路")

# 10. Top通路
print(f"\n   参与蛋白最多的前5条通路:")
top_pathways = df_output['pathway_name'].value_counts().head(5)
for pathway, count in top_pathways.items():
    pathway_short = pathway[:45] + '...' if len(pathway) > 45 else pathway
    print(f"   - {pathway_short} ({count}个蛋白)")

print("\n=== 步骤6完成！Reactome通路ready ===")
print(f"完成度: 85% → 93% (+8%)")
