# step11b_merge_clean_fixed.py
# 修正版：处理重复ID问题

import pandas as pd

print("=== Step 11b: 合并清理（修正版 - 处理重复） ===\n")

base_path = '/Users/pluviophile/graph/1025/data/processed'

v5_file = f'{base_path}/protein_master_v5_final.tsv'
v3_gene_file = f'{base_path}/protein_master_v3_with_gene_ids.tsv'
output_file = f'{base_path}/protein_master_v6_clean.tsv'

# 1. 读取v5
print("1. 读取step10输出（v5）...")
df_v5 = pd.read_csv(v5_file, sep='\t')
print(f"   总行数: {len(df_v5)}")

# 去重v5（保留第一个）
v5_dup = df_v5[df_v5.duplicated('uniprot_id', keep=False)]
if len(v5_dup) > 0:
    print(f"   ⚠️ 发现{len(v5_dup)}行重复")
    df_v5 = df_v5.drop_duplicates(subset='uniprot_id', keep='first')
    print(f"   ✅ 去重后: {len(df_v5)}行")

# 2. 读取v3_with_gene_ids
print("\n2. 读取Phase 1输出（v3+基因ID）...")
df_v3 = pd.read_csv(v3_gene_file, sep='\t')
print(f"   总行数: {len(df_v3)}")

# 去重v3（保留第一个）
v3_dup = df_v3[df_v3.duplicated('uniprot_id', keep=False)]
if len(v3_dup) > 0:
    print(f"   ⚠️ 发现{len(v3_dup)}行重复")
    print(f"   重复的ID数: {v3_dup['uniprot_id'].nunique()}")
    
    # 显示重复示例
    print("\n   重复示例（前5个）:")
    for uid in v3_dup['uniprot_id'].unique()[:5]:
        dup_rows = df_v3[df_v3['uniprot_id'] == uid]
        print(f"      {uid}: {len(dup_rows)}行")
        # 检查差异
        if len(dup_rows) == 2:
            row1 = dup_rows.iloc[0]
            row2 = dup_rows.iloc[1]
            diff_cols = []
            for col in dup_rows.columns:
                if str(row1[col]) != str(row2[col]):
                    diff_cols.append(col)
            if diff_cols:
                print(f"         差异字段: {diff_cols[:3]}")
    
    df_v3 = df_v3.drop_duplicates(subset='uniprot_id', keep='first')
    print(f"\n   ✅ 去重后: {len(df_v3)}行")

# 3. 提取step10更新的字段
print("\n3. 提取step10更新的4个字段...")
function_fields = ['function', 'go_biological_process', 
                   'go_molecular_function', 'go_cellular_component']

df_updates = df_v5[['uniprot_id'] + function_fields].copy()

# 4. 删除v3中的旧字段
print("\n4. 删除v3中的旧功能字段...")
for field in function_fields:
    if field in df_v3.columns:
        df_v3 = df_v3.drop(columns=[field])
        print(f"   删除: {field}")

# 5. 合并（现在应该1:1）
print("\n5. 合并更新...")
df_final = df_v3.merge(df_updates, on='uniprot_id', how='left')

print(f"   合并后行数: {len(df_final)}")
print(f"   预期行数: {len(df_v3)}")

if len(df_final) != len(df_v3):
    print(f"   ❌ 仍有问题！行数: {len(df_final)} vs {len(df_v3)}")
else:
    print(f"   ✅ 行数正确！")

# 6. 统计
print("\n6. 最终字段覆盖率:")

gene_id_fields = ['ncbi_gene_id', 'ensembl_gene_id', 'gene_synonyms']
print("\n   基因ID字段:")
for field in gene_id_fields:
    if field in df_final.columns:
        has_data = df_final[field].notna() & (df_final[field] != '')
        print(f"      {field:25} {has_data.sum():,}/{len(df_final):,} ({has_data.sum()/len(df_final)*100:.1f}%)")

print("\n   功能注释字段:")
for field in function_fields:
    if field in df_final.columns:
        has_data = df_final[field].notna() & (df_final[field] != '')
        print(f"      {field:30} {has_data.sum():,}/{len(df_final):,} ({has_data.sum()/len(df_final)*100:.1f}%)")

# 7. 保存
print(f"\n7. 保存...")
df_final.to_csv(output_file, sep='\t', index=False)

print(f"\n{'='*70}")
print("✅ 完成！")
print(f"{'='*70}")
print(f"输出: {output_file}")
print(f"行数: {len(df_final):,} (应该 = 19,135)")
print(f"列数: {len(df_final.columns)}")
print(f"\n注意: 去重后行数从19,253 → 19,135")
print(f"原因: Phase 1引入了66个重复ID（184行）")
print(f"{'='*70}\n")
