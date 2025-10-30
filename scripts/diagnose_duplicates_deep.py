# diagnose_duplicates_deep.py
# 深度分析重复ID的原因和解决方案

import pandas as pd

print("=== 重复ID深度诊断 ===\n")

base_path = '/Users/pluviophile/graph/1025/data/processed'

# 1. 读取v3_with_gene_ids
print("1. 读取v3_with_gene_ids...")
df = pd.read_csv(f'{base_path}/protein_master_v3_with_gene_ids.tsv', sep='\t')
print(f"   总行数: {len(df)}\n")

# 2. 找出重复的ID
print("2. 识别重复ID...")
duplicated = df[df.duplicated('uniprot_id', keep=False)]
print(f"   重复行数: {len(duplicated)}")
print(f"   重复ID数: {duplicated['uniprot_id'].nunique()}\n")

# 3. 详细分析每个重复ID
print("3. 详细分析重复ID（前10个）:\n")
print("="*100)

for uid in duplicated['uniprot_id'].unique()[:10]:
    dup_rows = df[df['uniprot_id'] == uid]
    print(f"\n【{uid}】 - {len(dup_rows)}行重复")
    print(f"  gene_names: {dup_rows.iloc[0]['gene_names']}")
    
    # 对比差异
    print("\n  差异对比:")
    print(f"  {'字段':<30} | {'第1行':<30} | {'第2行':<30}")
    print(f"  {'-'*30}-+-{'-'*30}-+-{'-'*30}")
    
    row1 = dup_rows.iloc[0]
    row2 = dup_rows.iloc[1] if len(dup_rows) > 1 else row1
    
    # 检查关键字段
    key_fields = ['hgnc_id', 'symbol', 'ncbi_gene_id', 'ensembl_gene_id', 'gene_synonyms']
    has_diff = False
    
    for field in key_fields:
        val1 = str(row1.get(field, ''))[:30]
        val2 = str(row2.get(field, ''))[:30]
        if val1 != val2:
            print(f"  {field:<30} | {val1:<30} | {val2:<30} ⚠️")
            has_diff = True
        else:
            print(f"  {field:<30} | {val1:<30} | (相同)")
    
    if not has_diff:
        print("  ⚠️ 所有检查字段都相同，是完全重复！")
    
    print("\n" + "="*100)

# 4. 统计分析
print("\n4. 重复类型统计:\n")

# 检查是否有完全重复（所有字段都一样）
complete_dup_count = 0
partial_dup_count = 0

for uid in duplicated['uniprot_id'].unique():
    dup_rows = df[df['uniprot_id'] == uid]
    
    # 检查是否所有字段都相同
    if len(dup_rows) == 2:
        row1_values = dup_rows.iloc[0].to_dict()
        row2_values = dup_rows.iloc[1].to_dict()
        
        diff_fields = [k for k in row1_values.keys() 
                      if str(row1_values[k]) != str(row2_values[k])]
        
        if len(diff_fields) == 0:
            complete_dup_count += 1
        else:
            partial_dup_count += 1

print(f"  完全重复（所有字段都相同）: {complete_dup_count}个ID")
print(f"  部分重复（有字段差异）:     {partial_dup_count}个ID")

# 5. 分析HGNC映射问题
print("\n5. HGNC映射分析:\n")

# 读取HGNC数据
hgnc_file = f'/Users/pluviophile/graph/1025/data/hgnc_core.tsv'
df_hgnc = pd.read_csv(hgnc_file, sep='\t', low_memory=False)

print(f"  HGNC总记录: {len(df_hgnc)}")

# 检查HGNC中是否有重复的symbol
hgnc_dup_symbols = df_hgnc[df_hgnc.duplicated('symbol', keep=False)]
if len(hgnc_dup_symbols) > 0:
    print(f"  ⚠️ HGNC中有{len(hgnc_dup_symbols)}行symbol重复")
    print(f"  重复的symbol数: {hgnc_dup_symbols['symbol'].nunique()}")
    
    print("\n  HGNC重复symbol示例（前5个）:")
    for sym in hgnc_dup_symbols['symbol'].unique()[:5]:
        sym_rows = df_hgnc[df_hgnc['symbol'] == sym]
        print(f"    {sym}: {len(sym_rows)}条记录")
        if 'status' in df_hgnc.columns:
            statuses = sym_rows['status'].value_counts().to_dict()
            print(f"      状态: {statuses}")
else:
    print("  ✅ HGNC中没有symbol重复")

# 6. 建议的解决方案
print("\n" + "="*100)
print("6. 解决方案建议:\n")

if complete_dup_count > 0:
    print(f"  ✅ 方案A: 删除{complete_dup_count}个完全重复的ID（简单去重）")

if partial_dup_count > 0:
    print(f"  ⚠️ 方案B: 处理{partial_dup_count}个部分重复的ID:")
    print("     - 如果HGNC有status字段，优先保留'Approved'")
    print("     - 如果没有status，保留ncbi_gene_id非空的")
    print("     - 最后保留第一条")

print("\n  📋 方案C: 修正Phase 1映射逻辑，避免重复产生")

print("\n" + "="*100)
print("\n运行完成！请查看上述分析，决定使用哪个方案。")

