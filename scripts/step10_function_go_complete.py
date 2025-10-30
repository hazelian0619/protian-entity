# step10_function_go_complete.py
# 位置：/Users/pluviophile/graph/1025/scripts/step10_function_go_complete.py

import pandas as pd
import requests
import time
from tqdm import tqdm

print("=== 步骤10: 补充功能和GO注释 ===\n")

# 路径配置
base_path = '/Users/pluviophile/graph/1025/data'
processed_dir = f'{base_path}/processed'
master_file = f'{processed_dir}/protein_master_v4_diseases.tsv'
output_file = f'{processed_dir}/protein_master_v5_final.tsv'

# 1. 读取主表
print("1. 读取主表...")
df = pd.read_csv(master_file, sep='\t')
print(f"   总蛋白数: {len(df)}\n")

# 2. 识别缺失的蛋白
print("2. 识别需要补充的蛋白...")

needs_update = (
    (df['function'].isna() | (df['function'] == '')) |
    (df['go_biological_process'].isna() | (df['go_biological_process'] == '')) |
    (df['go_molecular_function'].isna() | (df['go_molecular_function'] == '')) |
    (df['go_cellular_component'].isna() | (df['go_cellular_component'] == ''))
)

to_update = df[needs_update]['uniprot_id'].tolist()
print(f"   需要更新: {len(to_update)} 个蛋白\n")

if len(to_update) == 0:
    print("✅ 所有蛋白的功能和GO注释已完整！")
    df.to_csv(output_file, sep='\t', index=False)
    exit(0)

# 3. 从UniProt API获取数据
print("3. 从UniProt REST API获取数据...")
print(f"   预计时间: {len(to_update) * 0.3 / 60:.0f}-{len(to_update) * 0.5 / 60:.0f} 分钟\n")

def fetch_uniprot_data(uniprot_id):
    """从UniProt REST API获取蛋白信息"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    
    try:
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            
            result = {
                'function': None,
                'go_bp': [],
                'go_mf': [],
                'go_cc': []
            }
            
            # 提取function
            if 'comments' in data:
                for comment in data['comments']:
                    if comment.get('commentType') == 'FUNCTION':
                        texts = comment.get('texts', [])
                        if texts:
                            result['function'] = texts[0].get('value', '')
                            break
            
            # 提取GO术语
            if 'uniProtKBCrossReferences' in data:
                for ref in data['uniProtKBCrossReferences']:
                    if ref.get('database') == 'GO':
                        go_id = ref.get('id')
                        properties = ref.get('properties', [])
                        
                        go_term = None
                        go_aspect = None
                        
                        for prop in properties:
                            if prop.get('key') == 'GoTerm':
                                go_term = prop.get('value')
                            elif prop.get('key') == 'GoEvidenceType':
                                # 可以用来判断aspect，但通常GoTerm已包含P:/F:/C:前缀
                                pass
                        
                        if go_term and go_id:
                            # 提取aspect和term
                            if 'P:' in go_term:
                                term = go_term.replace('P:', '').strip()
                                result['go_bp'].append(f"{term} [{go_id}]")
                            elif 'F:' in go_term:
                                term = go_term.replace('F:', '').strip()
                                result['go_mf'].append(f"{term} [{go_id}]")
                            elif 'C:' in go_term:
                                term = go_term.replace('C:', '').strip()
                                result['go_cc'].append(f"{term} [{go_id}]")
            
            return result
        
        return None
        
    except Exception as e:
        return None

# 4. 批量更新
updated = 0
failed = 0
checkpoint_interval = 100

for i, uniprot_id in enumerate(tqdm(to_update, desc="更新进度"), 1):
    data = fetch_uniprot_data(uniprot_id)
    
    if data:
        idx = df[df['uniprot_id'] == uniprot_id].index[0]
        
        # 更新function（如果当前为空）
        if data['function']:
            if pd.isna(df.at[idx, 'function']) or str(df.at[idx, 'function']).strip() == '':
                df.at[idx, 'function'] = data['function']
        
        # 更新GO术语（如果当前为空）
        if data['go_bp']:
            if pd.isna(df.at[idx, 'go_biological_process']) or str(df.at[idx, 'go_biological_process']).strip() == '':
                df.at[idx, 'go_biological_process'] = '; '.join(data['go_bp'])
        
        if data['go_mf']:
            if pd.isna(df.at[idx, 'go_molecular_function']) or str(df.at[idx, 'go_molecular_function']).strip() == '':
                df.at[idx, 'go_molecular_function'] = '; '.join(data['go_mf'])
        
        if data['go_cc']:
            if pd.isna(df.at[idx, 'go_cellular_component']) or str(df.at[idx, 'go_cellular_component']).strip() == '':
                df.at[idx, 'go_cellular_component'] = '; '.join(data['go_cc'])
        
        updated += 1
    else:
        failed += 1
    
    # 避免API限流
    time.sleep(0.2)
    
    # 每100个保存checkpoint
    if i % checkpoint_interval == 0:
        df.to_csv(f'{processed_dir}/checkpoint_function_{i}.tsv', sep='\t', index=False)

# 5. 保存最终结果
df.to_csv(output_file, sep='\t', index=False)

# 6. 最终统计
print(f"\n{'='*60}")
print("✅ 完成！")
print(f"{'='*60}")
print(f"输出文件: {output_file}")
print(f"\n更新统计:")
print(f"  成功: {updated:,}")
print(f"  失败: {failed:,}")

print(f"\n最终覆盖率:")
for field in ['function', 'go_biological_process', 'go_molecular_function', 'go_cellular_component']:
    has_data = ~(df[field].isna() | (df[field] == ''))
    print(f"  {field:30} {has_data.sum():,}/{len(df):,} ({has_data.sum()/len(df)*100:.1f}%)")

print(f"\n{'='*60}")
print("🎉 所有数据补充完成！")
print(f"{'='*60}")
