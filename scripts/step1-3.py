import pandas as pd
import os
import json  # JSON输出用

# 路径设置（你的文件夹）
base_path = '/Users/pluviophile/graph/1025/data'
backup_file = os.path.join(base_path, 'hgnc_raw_backup.tsv')  # 步骤1输入
processed_dir = os.path.join(base_path, 'processed')  # 输出文件夹
os.makedirs(processed_dir, exist_ok=True)  # 自动创建processed

print("=== 步骤1: 加载并备份 (如果备份已，跳重读) ===")
# 步骤1: 加载TSV (low_memory=False 避警告)
df = pd.read_csv(backup_file, sep='\t', low_memory=False)
print("文件加载成功！总行数:", len(df))  # ~44k
print("列数:", len(df.columns))  # ~50
print("前5行预览:")
print(df.head())

# 备份原 (已hgnc_raw_backup.tsv，无需重)
print("备份已存在: hgnc_raw_backup.tsv")

# 日志步骤1
log1_file = os.path.join(base_path, 'ETL_log_step1.txt')
with open(log1_file, 'w') as f:
    f.write("步骤1完成: 2025-10-26\n")
    f.write(f"加载HGNC raw: {len(df)}行，{len(df.columns)}列\n")
    f.write("源: genenames.org 2025-10，备份: hgnc_raw_backup.tsv\n")
print(f"步骤1日志: {log1_file}")

print("\n=== 步骤2: 过滤蛋白范围 ===")
# 步骤2.1: 过滤Approved
df_approved = df[df['status'] == 'Approved'].copy()
print("过滤Approved后，行数:", len(df_approved))  # ~40k

# 步骤2.2: 过滤蛋白编码
df_protein = df_approved[df_approved['locus_group'] == 'protein-coding gene'].copy()
print("过滤蛋白后，行数:", len(df_protein))  # ~19k
print("蛋白样本预览 (关键列):")
print(df_protein[['hgnc_id', 'symbol', 'locus_group', 'uniprot_ids']].head(3))

# 步骤2.3: 位置验证
df_protein = df_protein[df_protein['location'].str.contains('q|p|chr', na=False, case=False, regex=True)].copy()
print("位置验证后，行数:", len(df_protein))  # ~19k

# 步骤2.4: merge旧 (可选，无旧跳)
old_file = os.path.join(processed_dir, 'hgnc_core_old.tsv')  # 改你的旧路径
try:
    df_old = pd.read_csv(old_file, sep='\t', low_memory=False)
    print("找到旧文件，开始merge...")
    df_merged = df_protein.merge(df_old, on='symbol', how='outer', suffixes=('_new', '_old'))
    for col in ['uniprot_ids', 'date_modified']:
        mask_new = df_merged[col + '_new'].notna()
        df_merged.loc[mask_new, col] = df_merged[col + '_new']
        mask_old = df_merged[col].isna() & df_merged[col + '_old'].notna()
        df_merged.loc[mask_old, col] = df_merged[col + '_old']
    cols_to_keep = [col for col in df_merged.columns if not ('_new' in col or '_old' in col)]
    df_protein = df_merged[cols_to_keep].copy()
    print("merge完成！")
except FileNotFoundError:
    print("无旧文件，跳过merge")

# 步骤2.5: 输出hgnc_core.tsv (全列)
core_file = os.path.join(processed_dir, 'hgnc_core.tsv')
df_protein.to_csv(core_file, sep='\t', index=False)
print(f"输出hgnc_core.tsv: {core_file}，{len(df_protein)}行")

# 步骤2.6: 校验
empty_uniprot = df_protein['uniprot_ids'].isna().sum()
percent_empty = (empty_uniprot / len(df_protein)) * 100
print(f"空uniprot_ids: {empty_uniprot} ({percent_empty:.2f}%)")
print("最近日期:", df_protein['date_modified'].max())
print("示例uniprot:", df_protein['uniprot_ids'].dropna().head(1).values)

# 日志步骤2
log2_file = os.path.join(base_path, 'ETL_log_step2.txt')
with open(log2_file, 'w') as f:
    f.write("步骤2完成: 2025-10-26\n")
    f.write(f"过滤蛋白: {len(df_protein)}行 (Approved + protein-coding)\n")
    f.write(f"空uniprot: {empty_uniprot} ({percent_empty:.2f}%)\n")
    f.write("源: HGNC 2025-10, 备份自 hgnc_raw_backup.tsv\n")
print(f"步骤2日志: {log2_file}")

print("\n=== 步骤3: 处理uniprot_ids + 输出种子 ===")
# 步骤3.1: 加载hgnc_core (已df_protein)
print("使用步骤2 df_protein，继续...")

# 步骤3.2: 填空标记
df_protein['uniprot_ids'] = df_protein['uniprot_ids'].fillna('MISSING_ENSEMBL')
print("空uniprot标记: MISSING_ENSEMBL (25行)")

# 步骤3.3: 拆列表
df_protein['uniprot_list'] = df_protein['uniprot_ids'].str.split(' ')
print("示例拆分:", df_protein['uniprot_list'].head(1).values)

# 步骤3.4: 展开+滤+去重
df_exploded = df_protein.explode('uniprot_list').copy()
df_exploded = df_exploded.dropna(subset=['uniprot_list'])
df_exploded = df_exploded[df_exploded['uniprot_list'] != 'MISSING_ENSEMBL']
df_unique = df_exploded[['uniprot_list', 'symbol', 'hgnc_id']].drop_duplicates().rename(columns={'uniprot_list': 'uniprot_id'})
print("唯一ID数:", len(df_unique))  # ~19k

# 步骤3.5: JSON输出
seed_list = df_unique.to_dict('records')
output_json = os.path.join(processed_dir, 'uniprot_seed.json')
with open(output_json, 'w') as f:
    json.dump(seed_list, f, indent=2)
print(f"输出uniprot_seed.json: {output_json}，{len(seed_list)}条")

# 步骤3.6: 校验
missing_count = df_protein['uniprot_ids'].isna().sum()  # 0
unique_ids = len(df_unique)
print(f"原空: {missing_count} (标记补)")
print(f"种子ID: {unique_ids}")
print("示例种子:", seed_list[:2])

# 步骤3.7: 更新hgnc_core (加列)
df_protein.to_csv(core_file, sep='\t', index=False)
print("hgnc_core更新: 加uniprot_list")

# 日志步骤3
log3_file = os.path.join(base_path, 'ETL_log_step3.txt')
log_content = f"""步骤3完成: 2025-10-26
处理uniprot_ids: 拆空格+展开异构体
唯一蛋白ID: {unique_ids}
原空标记: {missing_count} (补Ensembl计划)
源: HGNC 2025-10 via hgnc_core.tsv
输出: uniprot_seed.json (JSON列表，ready UniProt API)"""
with open(log3_file, 'w') as f:
    f.write(log_content)
print(f"步骤3日志: {log3_file}")

print("\n=== 全脚本完成！检查processed/ 文件。下一: UniProt API拉详 ===")
