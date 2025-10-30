# diagnose_hgnc_transcript.py
import pandas as pd

print("=== HGNC Transcript字段诊断 ===\n")

hgnc_file = '/Users/pluviophile/graph/1025/data/processed/hgnc_core.tsv'

# 1. 读取HGNC
print("1. 读取HGNC数据...")
df_hgnc = pd.read_csv(hgnc_file, sep='\t', low_memory=False)
print(f"   总记录: {len(df_hgnc)}")
print(f"   总字段: {len(df_hgnc.columns)}\n")

# 2. 查找包含transcript的字段
print("2. 查找包含'transcript'的字段...")
transcript_cols = [col for col in df_hgnc.columns if 'transcript' in col.lower()]

if transcript_cols:
    print(f"   ✅ 找到 {len(transcript_cols)} 个相关字段:")
    for col in transcript_cols:
        non_null = df_hgnc[col].notna().sum()
        print(f"      - {col:40} 非空: {non_null:,}/{len(df_hgnc):,} ({non_null/len(df_hgnc)*100:.1f}%)")
    
    print("\n3. 示例数据（前10行）:")
    print(df_hgnc[['symbol', 'hgnc_id'] + transcript_cols].head(10).to_string())
else:
    print("   ❌ 没有找到包含'transcript'的字段")
    print("   说明：HGNC可能不提供transcript ID")

# 3. 查看所有Ensembl相关字段
print("\n4. 所有Ensembl相关字段:")
ensembl_cols = [col for col in df_hgnc.columns if 'ensembl' in col.lower()]
for col in ensembl_cols:
    non_null = df_hgnc[col].notna().sum()
    print(f"   {col:40} 非空: {non_null:,}/{len(df_hgnc):,} ({non_null/len(df_hgnc)*100:.1f}%)")

# 4. 示例：看看某个基因的完整信息
print("\n5. 示例：TP53的HGNC信息")
tp53 = df_hgnc[df_hgnc['symbol'] == 'TP53']
if not tp53.empty:
    ensembl_fields = {col: tp53.iloc[0][col] for col in ensembl_cols}
    print("   Ensembl相关字段:")
    for k, v in ensembl_fields.items():
        print(f"      {k}: {v}")
else:
    print("   未找到TP53")
