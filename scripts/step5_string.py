# step5_string.py - STRING蛋白交互网络提取（v12.0）
import pandas as pd
from datetime import datetime
import os

print("=== 步骤5: STRING蛋白交互网络 (v12.0) ===")

# 路径
base_path = '/Users/pluviophile/graph/1025/data'
raw_dir = f'{base_path}/raw'
processed_dir = f'{base_path}/processed'
string_file = f'{raw_dir}/9606.protein.links.v12.0.txt'
master_file = f'{processed_dir}/protein_master_v2.tsv'
output_file = f'{processed_dir}/protein_edges.tsv'

# 1. 加载主表（建立映射）
print("1. 加载protein_master_v2.tsv...")
df_master = pd.read_csv(master_file, sep='\t')
print(f"   主表: {len(df_master)}个蛋白")

# 2. 建立ENSEMBL→UniProt映射
print("2. 建立映射表（ENSEMBL→UniProt）...")
mapping = {}
for _, row in df_master.iterrows():
    uniprot_id = row['uniprot_id']
    string_ids = str(row.get('string_ids', ''))
    if pd.notna(string_ids) and string_ids != '' and string_ids != 'nan':
        # string_ids格式: "9606.ENSP00000269305; 9606.ENSP00000269306"
        for sid in string_ids.split(';'):
            sid = sid.strip()
            if sid.startswith('9606.'):
                mapping[sid] = uniprot_id
            elif sid:  # 如果没有9606前缀，加上
                mapping[f'9606.{sid}'] = uniprot_id
print(f"   映射: {len(mapping)}个ENSEMBL ID")

# 3. 加载STRING原始文件
print("3. 加载STRING文件（~11M行，需要2-3分钟）...")
df_string = pd.read_csv(string_file, sep=' ')
print(f"   原始行数: {len(df_string):,}")

# 4. 过滤高可信交互（score≥400）
print("4. 过滤高可信交互（score≥400）...")
df_filtered = df_string[df_string['combined_score'] >= 400].copy()
print(f"   高可信: {len(df_filtered):,}行")

# 5. 映射到UniProt ID
print("5. 映射到UniProt ID（保留在主表内的蛋白）...")
edges = []
mapped_count = 0
for idx, row in df_filtered.iterrows():
    p1 = row['protein1']
    p2 = row['protein2']
    score = row['combined_score']
    
    # 两个蛋白都要在映射表里
    if p1 in mapping and p2 in mapping:
        uniprot1 = mapping[p1]
        uniprot2 = mapping[p2]
        edges.append({
            'source_uniprot': uniprot1,
            'target_uniprot': uniprot2,
            'combined_score': score,
            'source': 'STRING_v12.0',
            'fetch_date': datetime.now().strftime('%Y-%m-%d')
        })
        mapped_count += 1
    
    # 每100万行打印进度
    if (idx + 1) % 1000000 == 0:
        print(f"   处理进度: {idx+1:,}/{len(df_filtered):,} ({mapped_count:,}条已映射)")

df_edges = pd.DataFrame(edges)
print(f"   映射成功: {len(df_edges):,}条边")

# 6. 去重（无向图，A-B和B-A算同一条边）
print("6. 去重处理...")
df_edges['pair'] = df_edges.apply(
    lambda x: tuple(sorted([x['source_uniprot'], x['target_uniprot']])), 
    axis=1
)
df_edges = df_edges.drop_duplicates(subset='pair').drop(columns='pair')
print(f"   去重后: {len(df_edges):,}条边")

# 7. 输出
df_edges.to_csv(output_file, sep='\t', index=False)
print(f"\n✅ 输出: {output_file}")
print(f"   总边数: {len(df_edges):,}")
print(f"   平均score: {df_edges['combined_score'].mean():.1f}")
print(f"\n   示例前3行:")
print(df_edges.head(3).to_string(index=False))

# 8. 统计
unique_proteins = set(df_edges['source_uniprot']).union(set(df_edges['target_uniprot']))
print(f"\n   涉及蛋白数: {len(unique_proteins)}/{len(df_master)} ({len(unique_proteins)/len(df_master)*100:.1f}%)")

print("\n=== 步骤5完成！STRING交互网络ready ===")
print(f"完成度: 70% → 85% (+15%)")
