# step8b_update_master_table.py - 更新主表（添加AlphaFold字段）
import pandas as pd
import os

print("=== 步骤8b: 更新主表AlphaFold信息 ===")

base_path = '/Users/pluviophile/graph/1025/data'
processed_dir = f'{base_path}/processed'
master_file = f'{processed_dir}/protein_master_v3.tsv'
quality_file = f'{processed_dir}/alphafold_quality.tsv'
output_file = f'{processed_dir}/protein_master_v3_final.tsv'

# 1. 读取主表
print("1. 读取主表...")
df_master = pd.read_csv(master_file, sep='\t')
print(f"   主表: {len(df_master)} 个蛋白")

# 2. 读取附表
if not os.path.exists(quality_file):
    print(f"\n❌ 未找到附表: {quality_file}")
    print("请先运行: python step8a_extract_alphafold_quality.py")
    exit(1)

print("2. 读取附表...")
df_quality = pd.read_csv(quality_file, sep='\t')
print(f"   附表: {len(df_quality)} 个蛋白")

# 3. 清理主表（删除旧的AlphaFold字段）
print("3. 清理旧字段...")
columns_to_drop = ['alphafold_version', 'alphafold_note']
existing_columns = [col for col in columns_to_drop if col in df_master.columns]
if existing_columns:
    df_master = df_master.drop(columns=existing_columns)
    print(f"   删除: {', '.join(existing_columns)}")

# 重命名
if 'alphafold_url' in df_master.columns:
    df_master = df_master.rename(columns={'alphafold_url': 'alphafold_entry_url'})
    print("   重命名: alphafold_url → alphafold_entry_url")

# 4. 添加新字段（如果不存在）
print("4. 添加新字段...")
new_fields = {
    'has_alphafold': None,
    'alphafold_id': None,
    'alphafold_entry_url': None,
    'alphafold_pdb_url': None,
    'alphafold_mean_plddt': None
}

for field, default in new_fields.items():
    if field not in df_master.columns:
        df_master[field] = default

# 5. 生成AlphaFold信息
print("5. 生成AlphaFold信息...")

# 创建quality字典（快速查找）
quality_dict = df_quality.set_index('uniprot_id')['mean_plddt'].to_dict()

exists_count = 0
not_found_count = 0

for idx, row in df_master.iterrows():
    uniprot_id = row['uniprot_id']
    
    # 提取主ID（去掉异构体后缀）
    main_id = uniprot_id.split('-')[0] if '-' in str(uniprot_id) else uniprot_id
    
    # 检查是否在AlphaFold数据库中
    if main_id in quality_dict:
        # 有质量数据 = 在AlphaFold数据库中
        df_master.at[idx, 'has_alphafold'] = True
        df_master.at[idx, 'alphafold_id'] = f"AF-{main_id}-F1"
        df_master.at[idx, 'alphafold_entry_url'] = f"https://alphafold.ebi.ac.uk/entry/{main_id}"
        df_master.at[idx, 'alphafold_pdb_url'] = f"https://alphafold.ebi.ac.uk/files/AF-{main_id}-F1-model_v4.pdb"
        df_master.at[idx, 'alphafold_mean_plddt'] = quality_dict[main_id]
        exists_count += 1
    else:
        # 不在数据库中
        df_master.at[idx, 'has_alphafold'] = False
        not_found_count += 1

# 6. 保存
df_master.to_csv(output_file, sep='\t', index=False)

print(f"\n✅ 完成！")
print(f"\n输出文件:")
print(f"  主表: {output_file}")
print(f"  附表: {quality_file} (已存在)")

print(f"\n统计信息:")
print(f"  有AlphaFold: {exists_count}/{len(df_master)} ({exists_count/len(df_master)*100:.1f}%)")
print(f"  无AlphaFold: {not_found_count}/{len(df_master)} ({not_found_count/len(df_master)*100:.1f}%)")

print(f"\n新增字段（主表）:")
print(f"  - has_alphafold (是否有AlphaFold预测)")
print(f"  - alphafold_id (模型ID)")
print(f"  - alphafold_entry_url (网页链接)")
print(f"  - alphafold_pdb_url (PDB文件下载链接)")
print(f"  - alphafold_mean_plddt (平均置信度)")

print(f"\n附表字段 ({quality_file}):")
print(f"  - uniprot_id")
print(f"  - mean_plddt (平均置信度)")
print(f"  - median_plddt (中位数)")
print(f"  - min_plddt (最小值)")
print(f"  - max_plddt (最大值)")
print(f"  - very_high_conf_pct (>90分比例)")
print(f"  - high_conf_pct (70-90分比例)")
print(f"  - low_conf_pct (50-70分比例)")
print(f"  - very_low_conf_pct (<50分比例)")
print(f"  - quality_grade (质量等级)")
print(f"  - residue_count (残基数)")
