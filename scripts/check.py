# check.py - 诊断疾病信息
import pandas as pd  # ⬅️ 添加这一行！

df = pd.read_csv('/Users/pluviophile/graph/1025/data/processed/protein_master_v3_final_fixed.tsv', sep='\t')

# 有疾病信息的
has_disease = df[df['diseases'].notna() & (df['diseases'] != '')]

print(f"有疾病信息: {len(has_disease)}/{len(df)} ({len(has_disease)/len(df)*100:.1f}%)")
print("\n示例10个:")
print(has_disease[['uniprot_id', 'gene_names', 'diseases']].head(10).to_string())

print("\n疾病信息长度分布:")
print(has_disease['diseases'].str.len().describe())

print("\n最长的3个疾病信息:")
longest = has_disease.nlargest(3, has_disease['diseases'].str.len())
for _, row in longest.iterrows():
    print(f"\n{row['uniprot_id']} ({row['gene_names']}):")
    print(f"  长度: {len(row['diseases'])} 字符")
    print(f"  内容: {row['diseases'][:200]}...")
