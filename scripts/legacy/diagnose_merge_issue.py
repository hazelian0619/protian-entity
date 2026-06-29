# diagnose_merge_issue.py
import pandas as pd

print("=== 诊断merge行数异常 ===\n")

base_path = '/Users/pluviophile/graph/1025/data/processed'

# 1. 读取v5和v3
print("1. 读取文件...")
df_v5 = pd.read_csv(f'{base_path}/protein_master_v5_final.tsv', sep='\t')
df_v3 = pd.read_csv(f'{base_path}/protein_master_v3_with_gene_ids.tsv', sep='\t')

print(f"   v5行数: {len(df_v5)}")
print(f"   v3行数: {len(df_v3)}\n")

# 2. 检查v5是否有重复ID
print("2. 检查v5重复ID...")
v5_duplicates = df_v5[df_v5.duplicated('uniprot_id', keep=False)]
if len(v5_duplicates) > 0:
    print(f"   ❌ v5有{len(v5_duplicates)}行重复ID！")
    print(f"   重复的ID数: {v5_duplicates['uniprot_id'].nunique()}")
    print("\n   示例重复:")
    for uid in v5_duplicates['uniprot_id'].unique()[:5]:
        dup_rows = df_v5[df_v5['uniprot_id'] == uid]
        print(f"      {uid}: {len(dup_rows)}行")
else:
    print(f"   ✅ v5没有重复ID")

# 3. 检查v3是否有重复ID
print("\n3. 检查v3重复ID...")
v3_duplicates = df_v3[df_v3.duplicated('uniprot_id', keep=False)]
if len(v3_duplicates) > 0:
    print(f"   ❌ v3有{len(v3_duplicates)}行重复ID！")
else:
    print(f"   ✅ v3没有重复ID")

# 4. ID差异分析
print("\n4. ID差异分析...")
v5_ids = set(df_v5['uniprot_id'])
v3_ids = set(df_v3['uniprot_id'])

print(f"   v5唯一ID: {len(v5_ids)}")
print(f"   v3唯一ID: {len(v3_ids)}")
print(f"   共同ID: {len(v5_ids & v3_ids)}")
print(f"   v5独有ID: {len(v5_ids - v3_ids)}")
print(f"   v3独有ID: {len(v3_ids - v5_ids)}")

# 5. 模拟merge
print("\n5. 模拟merge结果...")
df_test = df_v3.merge(df_v5[['uniprot_id', 'function']], 
                      on='uniprot_id', how='left', suffixes=('_old', ''))
print(f"   merge后行数: {len(df_test)}")
print(f"   预期行数: {len(df_v3)}")
if len(df_test) > len(df_v3):
    print(f"   ❌ 行数膨胀: +{len(df_test) - len(df_v3)}")
