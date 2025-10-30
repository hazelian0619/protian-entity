# supplement_missing.py - 补充缺失的UniProt ID（银牌方案）
import pandas as pd
import requests
import time
import json
from datetime import datetime

print("=== 银牌方案：补充缺失ID ===")

# 路径
base_path = '/Users/pluviophile/graph/1025/data'
processed_dir = f'{base_path}/processed'
master_file = f'{processed_dir}/protein_master.tsv'
seed_file = f'{processed_dir}/uniprot_seed.json'
missing_file = f'{processed_dir}/missing_ids.txt'

# 加载
df_master = pd.read_csv(master_file, sep='\t')
with open(seed_file, 'r') as f:
    df_seed = pd.DataFrame(json.load(f))
with open(missing_file, 'r') as f:
    missing_ids = [line.strip() for line in f if line.strip()]

print(f"当前主表: {len(df_master)}行")
print(f"待补充: {len(missing_ids)} ID")

# 单独查询UniProt API
supplement_entries = []
failed_ids = []

for i, uid in enumerate(missing_ids, 1):
    url = f'https://rest.uniprot.org/uniprotkb/{uid}.json'
    print(f"[{i}/{len(missing_ids)}] 查询 {uid}...", end=' ')
    
    try:
        response = requests.get(url, headers={'Accept': 'application/json'}, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            
            # 提取字段（同步步骤4的结构）
            entry = {
                'uniprot_id': data.get('primaryAccession', uid),
                'entry_name': data.get('uniProtkbId', ''),
                'protein_name': data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                'gene_names': ' '.join([g.get('geneName', {}).get('value', '') for g in data.get('genes', [])]),
                'sequence': data.get('sequence', {}).get('value', ''),
                'sequence_len': len(data.get('sequence', {}).get('value', '')),
                'mass': data.get('sequence', {}).get('molWeight', 0),
                'function': '; '.join([c.get('texts', [{}])[0].get('value', '') for c in data.get('comments', []) if c.get('commentType') == 'FUNCTION']),
                'subcellular_location': '; '.join([c.get('texts', [{}])[0].get('value', '') for c in data.get('comments', []) if c.get('commentType') == 'SUBCELLULAR LOCATION']),
                'ptms': '; '.join([c.get('texts', [{}])[0].get('value', '') for c in data.get('comments', []) if c.get('commentType') == 'PTM']),
                'diseases': '; '.join([c.get('disease', {}).get('diseaseId', '') for c in data.get('comments', []) if c.get('commentType') == 'DISEASE']),
                'isoforms': str(data.get('comments', [{}])[0].get('isoforms', [])) if data.get('comments') else '',
                'go_biological_process': '; '.join([g.get('id', '') for g in data.get('uniProtKBCrossReferences', []) if g.get('database') == 'GO' and 'P:' in g.get('properties', {}).get('GoTerm', '')]),
                'go_molecular_function': '; '.join([g.get('id', '') for g in data.get('uniProtKBCrossReferences', []) if g.get('database') == 'GO' and 'F:' in g.get('properties', {}).get('GoTerm', '')]),
                'go_cellular_component': '; '.join([g.get('id', '') for g in data.get('uniProtKBCrossReferences', []) if g.get('database') == 'GO' and 'C:' in g.get('properties', {}).get('GoTerm', '')]),
                'domains': '',
                'pdb_ids': '; '.join([x.get('id', '') for x in data.get('uniProtKBCrossReferences', []) if x.get('database') == 'PDB']),
                'string_ids': '; '.join([x.get('id', '') for x in data.get('uniProtKBCrossReferences', []) if x.get('database') == 'STRING']),
                'keywords': '; '.join([k.get('name', '') for k in data.get('keywords', [])]),
                'date_modified': data.get('entryAudit', {}).get('lastAnnotationUpdateDate', ''),
                'source': 'UniProtKB/Swiss-Prot',
                'fetch_date': datetime.now().strftime('%Y-%m-%d')
            }
            
            # 添加HGNC信息（从seed merge）
            seed_row = df_seed[df_seed['uniprot_id'] == uid]
            if not seed_row.empty:
                entry['hgnc_id'] = seed_row.iloc[0]['hgnc_id']
                entry['symbol'] = seed_row.iloc[0]['symbol']
            
            supplement_entries.append(entry)
            print(f"✅ 成功")
        
        elif response.status_code == 404:
            failed_ids.append((uid, '404 Not Found'))
            print(f"❌ 404（已废弃或不存在）")
        
        else:
            failed_ids.append((uid, f'{response.status_code}'))
            print(f"❌ 错误{response.status_code}")
    
    except Exception as e:
        failed_ids.append((uid, str(e)[:50]))
        print(f"❌ 异常: {e}")
    
    time.sleep(0.5)  # 防限

print(f"\n=== 补充完成 ===")
print(f"成功: {len(supplement_entries)}")
print(f"失败: {len(failed_ids)}")

# 合并到主表
if supplement_entries:
    df_supplement = pd.DataFrame(supplement_entries)
    df_master_v2 = pd.concat([df_master, df_supplement], ignore_index=True)
    
    # 保存新版
    output_file = f'{processed_dir}/protein_master_v2.tsv'
    df_master_v2.to_csv(output_file, sep='\t', index=False)
    print(f"\n新版主表: {output_file}")
    print(f"总行数: {len(df_master_v2)} (原{len(df_master)} + 补{len(supplement_entries)})")
    
    # 覆盖率
    coverage = len(df_master_v2) / len(df_seed) * 100
    print(f"覆盖率: {coverage:.2f}%")
else:
    print("\n⚠️ 无成功补充，未生成新文件")

# 记录失败
if failed_ids:
    with open(f'{processed_dir}/failed_ids.txt', 'w') as f:
        f.write("UniProt ID\t原因\n")
        for uid, reason in failed_ids:
            f.write(f"{uid}\t{reason}\n")
    print(f"失败记录: failed_ids.txt")

print("\n=== 银牌方案完成 ===")
