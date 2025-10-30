# fix_alphafold_urls.py - 修正URL版本号
import pandas as pd

print("=== 修正AlphaFold URL版本号 ===\n")

base_path = '/Users/pluviophile/graph/1025/data/processed'
input_file = f'{base_path}/protein_master_v3_final.tsv'
output_file = f'{base_path}/protein_master_v3_final_fixed.tsv'

# 读取
df = pd.read_csv(input_file, sep='\t')

# 修正URL（v4 → v6）
df['alphafold_pdb_url'] = df['alphafold_pdb_url'].str.replace('v4.pdb', 'v6.pdb')

# 保存
df.to_csv(output_file, sep='\t', index=False)

print(f"✅ 修正完成！")
print(f"输入: {input_file}")
print(f"输出: {output_file}")

# 验证
v4_count = df['alphafold_pdb_url'].str.contains('v4', na=False).sum()
v6_count = df['alphafold_pdb_url'].str.contains('v6', na=False).sum()

print(f"\n验证:")
print(f"  v4 URLs: {v4_count}")
print(f"  v6 URLs: {v6_count}")

if v4_count == 0:
    print("\n  ✅ 所有URL已更新为v6")
