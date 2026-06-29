# quick_check_v6.py
import pandas as pd

df = pd.read_csv('/Users/pluviophile/graph/1025/data/processed/protein_master_v6_clean.tsv', sep='\t')

print(f"总行数: {len(df)}")
print(f"总列数: {len(df.columns)}")
print(f"\ngene_synonyms统计:")
has_syn = df['gene_synonyms'].notna() & (df['gene_synonyms'] != '') & (df['gene_synonyms'] != 'nan')
print(f"  有数据: {has_syn.sum()}/{len(df)} ({has_syn.sum()/len(df)*100:.1f}%)")

print(f"\n示例（有synonyms的）:")
with_syn = df[has_syn].head(5)
print(with_syn[['uniprot_id', 'symbol', 'gene_synonyms']].to_string())
