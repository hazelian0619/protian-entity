# diagnose_v3_original.py
import pandas as pd

print("=== 诊断v3原始数据 ===\n")

base_path = '/Users/pluviophile/graph/1025/data/processed'
v3_file = f'{base_path}/protein_master_v3_final_fixed.tsv'

# 读取v3
df = pd.read_csv(v3_file, sep='\t')
print(f"v3总行数: {len(df)}\n")

# 检查v3是否有重复uniprot_id
duplicated = df[df.duplicated('uniprot_id', keep=False)]

if len(duplicated) > 0:
    print(f"❌ v3原始数据就有{len(duplicated)}行重复！")
    print(f"   重复ID数: {duplicated['uniprot_id'].nunique()}\n")
    
    print("重复ID示例（前5个）:")
    for uid in duplicated['uniprot_id'].unique()[:5]:
        dup_rows = df[df['uniprot_id'] == uid]
        print(f"  {uid}: {len(dup_rows)}行")
        print(f"    gene_names: {dup_rows.iloc[0]['gene_names']}")
else:
    print("✅ v3原始数据没有重复")
