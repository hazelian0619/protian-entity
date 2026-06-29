# step9_diseases_disgenet.py
# 位置：/Users/pluviophile/graph/1025/scripts/step9_diseases_disgenet.py

import pandas as pd
import requests
import gzip
from io import BytesIO

print("=== 步骤9: 补充疾病信息 (DisGeNET) ===\n")

# 路径配置
base_path = '/Users/pluviophile/graph/1025/data'
processed_dir = f'{base_path}/processed'
master_file = f'{processed_dir}/protein_master_v3_final_fixed.tsv'
output_file = f'{processed_dir}/protein_master_v4_diseases.tsv'

# 1. 读取主表
print("1. 读取主表...")
df = pd.read_csv(master_file, sep='\t')
print(f"   总蛋白数: {len(df)}")

# 统计当前疾病信息
has_disease = ~(df['diseases'].isna() | (df['diseases'] == ''))
print(f"   已有疾病信息: {has_disease.sum()}/{len(df)} ({has_disease.sum()/len(df)*100:.1f}%)\n")

# 2. 下载DisGeNET数据
print("2. 下载DisGeNET基因-疾病关联...")
print("   URL: https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz")
print("   注意：文件约50MB，下载需要2-5分钟\n")

disgenet_url = "https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz"

try:
    response = requests.get(disgenet_url, timeout=300)
    
    if response.status_code == 200:
        print("   ✅ 下载成功，解析中...")
        
        # 解压并读取
        with gzip.open(BytesIO(response.content), 'rt') as f:
            df_disgenet = pd.read_csv(f, sep='\t')
        
        print(f"   DisGeNET总记录: {len(df_disgenet)}")
        
        # 筛选人类基因（NCBI Taxonomy ID = 9606）
        if 'NID' in df_disgenet.columns:
            df_disgenet = df_disgenet[df_disgenet['NID'] == 9606]
            print(f"   人类基因关联: {len(df_disgenet)}")
        
    else:
        print(f"   ❌ 下载失败: HTTP {response.status_code}")
        print("   跳过DisGeNET更新")
        df_disgenet = None
        
except Exception as e:
    print(f"   ❌ 错误: {str(e)}")
    print("   跳过DisGeNET更新")
    df_disgenet = None

# 3. 匹配并更新
if df_disgenet is not None:
    print("\n3. 匹配基因-疾病关联...")
    
    # 创建基因-疾病映射（按基因符号）
    gene_disease_map = {}
    
    for _, row in df_disgenet.iterrows():
        gene = row.get('geneName', '')
        disease = row.get('diseaseName', '')
        score = row.get('score', 0)
        
        # 只保留高置信度关联
        if gene and disease and score > 0.3:
            if gene not in gene_disease_map:
                gene_disease_map[gene] = set()
            gene_disease_map[gene].add(disease)
    
    print(f"   有疾病关联的基因: {len(gene_disease_map)}")
    
    # 更新主表
    updated_count = 0
    
    for idx, row in df.iterrows():
        # 如果已有疾病信息，跳过
        if pd.notna(row['diseases']) and str(row['diseases']).strip() != '':
            continue
        
        # 从gene_names提取基因符号
        gene_names_str = str(row.get('gene_names', ''))
        if not gene_names_str or gene_names_str == 'nan':
            continue
        
        gene_names = gene_names_str.split()
        
        # 尝试匹配
        diseases_found = set()
        for gene in gene_names:
            if gene in gene_disease_map:
                diseases_found.update(gene_disease_map[gene])
        
        if diseases_found:
            # 去重并合并（最多10个疾病）
            diseases_list = list(diseases_found)[:10]
            df.at[idx, 'diseases'] = '; '.join(diseases_list)
            updated_count += 1
    
    print(f"   ✅ 补充了 {updated_count} 个蛋白的疾病信息\n")

# 4. 保存
print("4. 保存结果...")
df.to_csv(output_file, sep='\t', index=False)

# 5. 最终统计
has_disease_new = ~(df['diseases'].isna() | (df['diseases'] == ''))

print(f"\n{'='*60}")
print("✅ 完成！")
print(f"{'='*60}")
print(f"输出文件: {output_file}")
print(f"\n疾病信息覆盖率:")
print(f"  更新前: {has_disease.sum():,}/{len(df):,} ({has_disease.sum()/len(df)*100:.1f}%)")
print(f"  更新后: {has_disease_new.sum():,}/{len(df):,} ({has_disease_new.sum()/len(df)*100:.1f}%)")
print(f"  增加: +{has_disease_new.sum() - has_disease.sum():,} 个蛋白")
print(f"\n下一步: python step10_function_go_complete.py")
