# step9_diseases_disgenet_api.py
# 位置：/Users/pluviophile/graph/1025/scripts/step9_diseases_disgenet_api.py

import pandas as pd
import requests
import time
import json
from tqdm import tqdm

print("=== 步骤9: 补充疾病信息 (DisGeNET API) ===\n")

# ⚠️ 替换为你的API Key（从DisGeNET控制台获取）
DISGENET_API_KEY = "YOUR_API_KEY_HERE"

# 路径配置
base_path = '/Users/pluviophile/graph/1025/data'
processed_dir = f'{base_path}/processed'
master_file = f'{processed_dir}/protein_master_v4_diseases.tsv'  # step9输出的
output_file = f'{processed_dir}/protein_master_v4_diseases_final.tsv'
cache_file = f'{base_path}/disgenet_cache.json'

# 1. 读取主表
print("1. 读取主表...")
df = pd.read_csv(master_file, sep='\t')
print(f"   总蛋白数: {len(df)}")

has_disease = ~(df['diseases'].isna() | (df['diseases'] == ''))
print(f"   已有疾病信息: {has_disease.sum()}/{len(df)} ({has_disease.sum()/len(df)*100:.1f}%)\n")

# 2. 提取需要查询的基因
print("2. 识别需要补充的基因...")
genes_to_query = set()

for idx, row in df.iterrows():
    # 跳过已有疾病信息的
    if pd.notna(row['diseases']) and str(row['diseases']).strip() != '':
        continue
    
    gene_names = str(row.get('gene_names', '')).split()
    genes_to_query.update(gene_names)

# 移除空值和'nan'
genes_to_query = {g for g in genes_to_query if g and g != 'nan'}
print(f"   需要查询的基因: {len(genes_to_query)}\n")

# 3. 加载缓存（如果有）
gene_disease_cache = {}
if os.path.exists(cache_file):
    with open(cache_file, 'r') as f:
        gene_disease_cache = json.load(f)
    print(f"3. 加载缓存: {len(gene_disease_cache)} 个基因\n")
else:
    print("3. 无缓存文件\n")

# 4. 查询DisGeNET API
print("4. 查询DisGeNET API...")
print(f"   API URL: https://api.disgenet.com/api/v1/gda/gene")
print(f"   API Key: {DISGENET_API_KEY[:10]}...\n")

def query_disgenet_gene(gene_symbol, api_key):
    """查询DisGeNET API获取基因-疾病关联"""
    
    # 如果已在缓存，返回缓存
    if gene_symbol in gene_disease_cache:
        return gene_disease_cache[gene_symbol]
    
    url = f"https://api.disgenet.com/api/v1/gda/gene"
    
    headers = {
        "Authorization": f"Bearer {api_key}",
        "accept": "application/json"
    }
    
    params = {
        "gene": gene_symbol,
        "source": "CURATED",
        "min_score": 0.3,
        "limit": 20  # 每个基因最多20个疾病
    }
    
    try:
        response = requests.get(url, headers=headers, params=params, timeout=15)
        
        if response.status_code == 200:
            data = response.json()
            
            diseases = []
            if isinstance(data, list):
                for item in data:
                    disease_name = item.get('disease_name')
                    if disease_name:
                        diseases.append(disease_name)
            
            # 缓存结果
            gene_disease_cache[gene_symbol] = diseases
            return diseases
        
        elif response.status_code == 401:
            print(f"\n   ❌ API认证失败！")
            print(f"   请检查API Key: {api_key[:20]}...")
            return None
        
        elif response.status_code == 404:
            # 基因不存在，缓存空结果
            gene_disease_cache[gene_symbol] = []
            return []
        
        else:
            print(f"\n   ⚠️ HTTP {response.status_code}: {gene_symbol}")
            return []
        
    except Exception as e:
        print(f"\n   ❌ 错误 ({gene_symbol}): {str(e)}")
        return []

# 查询所有基因
queried = 0
failed = 0
auth_failed = False

for gene in tqdm(list(genes_to_query), desc="查询进度"):
    result = query_disgenet_gene(gene, DISGENET_API_KEY)
    
    if result is None:  # 认证失败
        auth_failed = True
        break
    
    if result:
        queried += 1
    
    # 避免API限流：每个请求间隔0.5秒
    time.sleep(0.5)
    
    # 每查询100个基因保存一次缓存
    if (queried + failed) % 100 == 0:
        with open(cache_file, 'w') as f:
            json.dump(gene_disease_cache, f)

# 保存最终缓存
with open(cache_file, 'w') as f:
    json.dump(gene_disease_cache, f)

if auth_failed:
    print("\n❌ API认证失败，请检查API Key")
    print("修改脚本中的 DISGENET_API_KEY 变量")
    exit(1)

print(f"\n   查询结果: {queried} 个基因有疾病关联")

# 5. 更新主表
print("\n5. 更新主表...")
updated_count = 0

for idx, row in df.iterrows():
    # 跳过已有疾病信息的
    if pd.notna(row['diseases']) and str(row['diseases']).strip() != '':
        continue
    
    gene_names = str(row.get('gene_names', '')).split()
    
    diseases_found = set()
    for gene in gene_names:
        if gene in gene_disease_cache:
            diseases_found.update(gene_disease_cache[gene])
    
    if diseases_found:
        # 最多保留10个疾病
        diseases_list = list(diseases_found)[:10]
        df.at[idx, 'diseases'] = '; '.join(diseases_list)
        updated_count += 1

print(f"   ✅ 补充了 {updated_count} 个蛋白的疾病信息")

# 6. 保存
df.to_csv(output_file, sep='\t', index=False)

has_disease_new = ~(df['diseases'].isna() | (df['diseases'] == ''))

print(f"\n{'='*60}")
print("✅ 完成！")
print(f"{'='*60}")
print(f"输出文件: {output_file}")
print(f"缓存文件: {cache_file}")
print(f"\n疾病信息覆盖率:")
print(f"  更新前: {has_disease.sum():,}/{len(df):,} ({has_disease.sum()/len(df)*100:.1f}%)")
print(f"  更新后: {has_disease_new.sum():,}/{len(df):,} ({has_disease_new.sum()/len(df)*100:.1f}%)")
print(f"  增加: +{has_disease_new.sum() - has_disease.sum():,} 个蛋白")

print(f"\n注意:")
print(f"  - 缓存已保存，下次运行会跳过已查询的基因")
print(f"  - 如果中断，重新运行会从断点继续")
