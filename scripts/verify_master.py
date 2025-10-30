# verify_master.py - 验证protein_master.tsv质量
import pandas as pd
import os

master_file = '/Users/pluviophile/graph/1025/data/processed/protein_master.tsv'

print("=== 验证protein_master.tsv ===")
# 加载
df = pd.read_csv(master_file, sep='\t')

# 指标1：行数
print(f"✅ 总行数: {len(df)} (应该~19000)")

# 指标2：空值率
seq_empty = df['sequence'].isna().sum()
func_empty = df['function'].isna().sum()
dis_empty = df['diseases'].isna().sum()
print(f"✅ 空序列: {seq_empty} ({seq_empty/len(df)*100:.2f}%, 应该<1%)")
print(f"✅ 空功能: {func_empty} ({func_empty/len(df)*100:.2f}%, 应该<20%)")
print(f"✅ 空疾病: {dis_empty} ({dis_empty/len(df)*100:.2f}%, 应该~50%)")

# 指标3：关键蛋白
key_proteins = ['P04637', 'P05231', 'P01308', 'P35222', 'Q9Y6K9']
found = df[df['uniprot_id'].isin(key_proteins)]
print(f"✅ 关键蛋白: {len(found)}/5 (TP53/IL6/INS/CTNNB1/NFKB2)")

# 示例
print("\n示例: TP53 (P04637)")
tp53 = df[df['uniprot_id'] == 'P04637']
if not tp53.empty:
    print(f"  Symbol: {tp53['symbol'].values[0]}")
    print(f"  Seq长度: {tp53['sequence_len'].values[0]}")
    print(f"  有Function: {'是' if not pd.isna(tp53['function'].values[0]) else '否'}")
    print(f"  有Diseases: {'是' if not pd.isna(tp53['diseases'].values[0]) else '否'}")
else:
    print("  ❌ 未找到TP53")

print("\n=== 验证完成 ===")
