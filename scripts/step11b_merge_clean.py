# step11b_merge_clean.py
# 位置：/Users/pluviophile/graph/1025/scripts/step11b_merge_clean.py
# 
# 功能：合并step10的功能更新 + Phase 1的基因ID，清理v4引入的空字段

import pandas as pd

print("=== Step 11b: 合并清理 - step10更新 + Phase 1基因ID ===\n")

# 路径配置
base_path = '/Users/pluviophile/graph/1025/data/processed'

# 输入文件
v5_file = f'{base_path}/protein_master_v5_final.tsv'        # step10输出
v3_gene_file = f'{base_path}/protein_master_v3_with_gene_ids.tsv'  # Phase 1输出

# 输出文件
output_file = f'{base_path}/protein_master_v6_clean.tsv'

print("=" * 70)
print("📋 文件路径")
print("=" * 70)
print(f"v5 (step10输出):      {v5_file}")
print(f"v3+基因ID (Phase 1):  {v3_gene_file}")
print(f"输出 (v6 clean):      {output_file}\n")

# 1. 读取v5（step10的更新）
print("1. 读取step10输出（v5）...")
try:
    df_v5 = pd.read_csv(v5_file, sep='\t')
    print(f"   ✅ 总蛋白数: {len(df_v5)}")
    print(f"   ✅ 列数: {len(df_v5.columns)}")
except Exception as e:
    print(f"   ❌ 读取失败: {e}")
    exit(1)

# 2. 读取v3_with_gene_ids（干净版本）
print("\n2. 读取Phase 1输出（v3+基因ID）...")
try:
    df_clean = pd.read_csv(v3_gene_file, sep='\t')
    print(f"   ✅ 总蛋白数: {len(df_clean)}")
    print(f"   ✅ 列数: {len(df_clean.columns)}")
except Exception as e:
    print(f"   ❌ 读取失败: {e}")
    exit(1)

# 3. 验证数据一致性
print("\n3. 验证数据一致性...")
if len(df_v5) != len(df_clean):
    print(f"   ⚠️ 警告：行数不一致 (v5: {len(df_v5)}, v3: {len(df_clean)})")
else:
    print(f"   ✅ 行数一致: {len(df_v5)}")

# 检查uniprot_id是否匹配
v5_ids = set(df_v5['uniprot_id'])
v3_ids = set(df_clean['uniprot_id'])
common_ids = v5_ids & v3_ids
print(f"   ✅ 共同ID数: {len(common_ids)}/{len(df_clean)}")

# 4. 提取step10更新的字段
print("\n4. 提取step10更新的4个字段...")
function_fields = ['function', 'go_biological_process', 
                   'go_molecular_function', 'go_cellular_component']

# 检查字段是否存在
missing_fields = [f for f in function_fields if f not in df_v5.columns]
if missing_fields:
    print(f"   ⚠️ v5中缺少字段: {missing_fields}")
    print("   使用v5中存在的字段...")
    function_fields = [f for f in function_fields if f in df_v5.columns]

print(f"   提取字段: {function_fields}")
df_updates = df_v5[['uniprot_id'] + function_fields].copy()

# 5. 统计更新覆盖率
print("\n5. step10更新统计:")
for field in function_fields:
    if field in df_updates.columns:
        has_data = df_updates[field].notna() & (df_updates[field] != '')
        print(f"   {field:30} {has_data.sum():,}/{len(df_updates):,} ({has_data.sum()/len(df_updates)*100:.1f}%)")

# 6. 更新到干净版本
print("\n6. 合并更新到干净版本...")

# 先删除v3_with_gene_ids中的这4个字段（如果存在）
for field in function_fields:
    if field in df_clean.columns:
        df_clean = df_clean.drop(columns=[field])
        print(f"   删除旧字段: {field}")

# 合并更新
print("   执行合并...")
df_final = df_clean.merge(df_updates, on='uniprot_id', how='left')

print(f"   ✅ 合并完成")

# 7. 验证最终结果
print("\n7. 验证最终数据...")
print(f"   行数: {len(df_final)} (应该 = {len(df_clean)})")
print(f"   列数: {len(df_final.columns)}")

if len(df_final) != len(df_clean):
    print(f"   ⚠️ 警告：合并后行数变化！")
else:
    print(f"   ✅ 行数保持不变")

# 8. 统计最终字段覆盖率
print("\n8. 最终字段覆盖率:")

# 基因ID字段
gene_id_fields = ['ncbi_gene_id', 'ensembl_gene_id', 'gene_synonyms']
print("\n   基因ID字段 (Phase 1):")
for field in gene_id_fields:
    if field in df_final.columns:
        has_data = df_final[field].notna() & (df_final[field] != '')
        print(f"      {field:25} {has_data.sum():,}/{len(df_final):,} ({has_data.sum()/len(df_final)*100:.1f}%)")

# 功能注释字段
print("\n   功能注释字段 (step10更新):")
for field in function_fields:
    if field in df_final.columns:
        has_data = df_final[field].notna() & (df_final[field] != '')
        print(f"      {field:30} {has_data.sum():,}/{len(df_final):,} ({has_data.sum()/len(df_final)*100:.1f}%)")

# 9. 列出所有字段
print(f"\n9. v6_clean字段列表 ({len(df_final.columns)}列):")
for i, col in enumerate(df_final.columns, 1):
    print(f"   {i:2}. {col}")

# 10. 保存
print(f"\n10. 保存干净的v6...")
try:
    df_final.to_csv(output_file, sep='\t', index=False)
    print(f"    ✅ 保存成功")
except Exception as e:
    print(f"    ❌ 保存失败: {e}")
    exit(1)

# 11. 最终总结
print(f"\n{'='*70}")
print("✅ 合并清理完成！")
print(f"{'='*70}")
print(f"输出文件: {output_file}")
print(f"行数: {len(df_final):,}")
print(f"列数: {len(df_final.columns)}")
print(f"\n数据来源:")
print(f"  ✅ v3原始29列 (基础信息)")
print(f"  ✅ Phase 1的4列 (基因ID映射)")
print(f"  ✅ step10的4列 (功能+GO注释更新)")
print(f"\n预期列数: 33列 (29 + 4基因ID)")
print(f"实际列数: {len(df_final.columns)}列")

if len(df_final.columns) == 33:
    print("\n🎉 列数正确！数据结构完整！")
elif len(df_final.columns) > 33:
    print(f"\n⚠️ 列数多于预期，可能有冗余字段")
    extra_cols = len(df_final.columns) - 33
    print(f"   多出{extra_cols}列，请检查")
else:
    print(f"\n⚠️ 列数少于预期，可能有字段缺失")

# 12. 示例数据
print(f"\n示例数据（前2行）:")
print("="*70)
example_cols = ['uniprot_id', 'gene_names', 'ncbi_gene_id', 'ensembl_gene_id', 
                'function', 'go_biological_process']
available_cols = [c for c in example_cols if c in df_final.columns]
print(df_final[available_cols].head(2).to_string())

print(f"\n{'='*70}")
print("下一步: python step12_add_uniprot_metadata.py")
print(f"{'='*70}\n")
