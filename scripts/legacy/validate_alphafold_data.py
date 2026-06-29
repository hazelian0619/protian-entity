# validate_alphafold_data.py - 验证数据完整性
import pandas as pd

print("=== AlphaFold数据完整性验证 ===\n")

base_path = '/Users/pluviophile/graph/1025/data/processed'

# 1. 读取数据
df_master = pd.read_csv(f'{base_path}/protein_master_v3_final.tsv', sep='\t')
df_quality = pd.read_csv(f'{base_path}/alphafold_quality.tsv', sep='\t')

print(f"主表: {len(df_master)} 个蛋白")
print(f"附表: {len(df_quality)} 个蛋白\n")

# 2. 检查mean_plddt匹配
print("=== 检查1: mean_plddt值匹配 ===")

# 测试几个示例
test_ids = ['A0A087X1C5', 'A0A0B4J2F0', 'A0A0K2S4Q6', 'A0A024R1R8']

for uid in test_ids:
    # 从主表获取
    master_row = df_master[df_master['uniprot_id'] == uid]
    
    if len(master_row) == 0:
        print(f"  ⚠️ {uid}: 不在主表中")
        continue
    
    master_plddt = master_row.iloc[0]['alphafold_mean_plddt']
    
    # 从附表获取
    quality_row = df_quality[df_quality['uniprot_id'] == uid]
    
    if len(quality_row) == 0:
        print(f"  ❌ {uid}: 主表有(pLDDT={master_plddt})，附表缺失")
    else:
        quality_plddt = quality_row.iloc[0]['mean_plddt']
        
        if abs(master_plddt - quality_plddt) < 0.01:
            print(f"  ✅ {uid}: 匹配 (主表={master_plddt}, 附表={quality_plddt})")
        else:
            print(f"  ❌ {uid}: 不匹配！(主表={master_plddt}, 附表={quality_plddt})")

# 3. 检查URL版本号
print(f"\n=== 检查2: URL版本号 ===")
sample_urls = df_master['alphafold_pdb_url'].dropna().head(3)
for url in sample_urls:
    if 'v4' in url:
        print(f"  ⚠️ URL包含v4: {url}")
    elif 'v6' in url:
        print(f"  ✅ URL包含v6: {url}")
    else:
        print(f"  ❓ URL版本不明: {url}")

# 4. 检查覆盖率
print(f"\n=== 检查3: 覆盖率统计 ===")

has_af = df_master['has_alphafold'].sum()
has_plddt = df_master['alphafold_mean_plddt'].notna().sum()

print(f"  has_alphafold=True: {has_af}/{len(df_master)} ({has_af/len(df_master)*100:.1f}%)")
print(f"  有mean_plddt值: {has_plddt}/{len(df_master)} ({has_plddt/len(df_master)*100:.1f}%)")

# 5. 检查缺失的蛋白
print(f"\n=== 检查4: 缺失数据分析 ===")

# 主表中标记有AlphaFold但没有pLDDT的
has_af_no_plddt = df_master[(df_master['has_alphafold'] == True) & 
                             (df_master['alphafold_mean_plddt'].isna())]

print(f"  标记有AlphaFold但无pLDDT: {len(has_af_no_plddt)}")
if len(has_af_no_plddt) > 0:
    print(f"  示例ID: {has_af_no_plddt['uniprot_id'].head(5).tolist()}")

# 6. 质量分布对比
print(f"\n=== 检查5: 质量分布 ===")

print(f"  附表质量分布:")
quality_dist = df_quality['quality_grade'].value_counts()
for grade, count in quality_dist.items():
    print(f"    {grade}: {count} ({count/len(df_quality)*100:.1f}%)")

print(f"\n  主表pLDDT统计:")
print(f"    平均值: {df_master['alphafold_mean_plddt'].mean():.2f}")
print(f"    中位数: {df_master['alphafold_mean_plddt'].median():.2f}")
print(f"    最小值: {df_master['alphafold_mean_plddt'].min():.2f}")
print(f"    最大值: {df_master['alphafold_mean_plddt'].max():.2f}")

# 7. 关键问题检查
print(f"\n=== 检查6: 关键问题 ===")

# URL版本问题
v4_urls = df_master['alphafold_pdb_url'].str.contains('v4', na=False).sum()
v6_urls = df_master['alphafold_pdb_url'].str.contains('v6', na=False).sum()

if v4_urls > 0:
    print(f"  ⚠️ 发现{v4_urls}个v4 URL（应该是v6）")
else:
    print(f"  ✅ 所有URL都是v6版本")

# 8. 总结
print(f"\n=== 总结 ===")

issues = []

if v4_urls > 0:
    issues.append("URL版本号错误（v4应改为v6）")

if len(has_af_no_plddt) > 0:
    issues.append(f"{len(has_af_no_plddt)}个蛋白标记有AlphaFold但无pLDDT值")

if has_plddt < has_af:
    issues.append("pLDDT覆盖率低于has_alphafold覆盖率")

if len(issues) == 0:
    print("  ✅ 数据完整性验证通过！")
else:
    print("  发现以下问题:")
    for i, issue in enumerate(issues, 1):
        print(f"    {i}. {issue}")
