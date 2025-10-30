# find_missing.py - 找出缺失的39个UniProt ID
import pandas as pd
import json

# 加载
seed_file = '/Users/pluviophile/graph/1025/data/processed/uniprot_seed.json'
master_file = '/Users/pluviophile/graph/1025/data/processed/protein_master.tsv'

with open(seed_file, 'r') as f:
    seed_list = json.load(f)
df_seed = pd.DataFrame(seed_list)

df_master = pd.read_csv(master_file, sep='\t')

# 找缺失
seed_ids = set(df_seed['uniprot_id'])
master_ids = set(df_master['uniprot_id'])
missing_ids = sorted(list(seed_ids - master_ids))

print(f"=== 缺失ID分析 ===")
print(f"种子总数: {len(seed_ids)}")
print(f"主表行数: {len(master_ids)}")
print(f"缺失数量: {len(missing_ids)}")
print(f"\n缺失ID列表:")
for i, uid in enumerate(missing_ids, 1):
    # 找对应基因名
    symbol = df_seed[df_seed['uniprot_id'] == uid]['symbol'].values[0]
    print(f"{i}. {uid} ({symbol})")

# 保存到文件
with open('/Users/pluviophile/graph/1025/data/processed/missing_ids.txt', 'w') as f:
    f.write('\n'.join(missing_ids))
print(f"\n已保存到: missing_ids.txt")
