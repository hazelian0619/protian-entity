# supplement_missing_fixed.py - 修复复合ID问题
import pandas as pd
import requests
import time
import json
from datetime import datetime

print("=== 银牌方案修复版：处理复合ID ===")

base_path = '/Users/pluviophile/graph/1025/data'
processed_dir = f'{base_path}/processed'
master_file = f'{processed_dir}/protein_master.tsv'
seed_file = f'{processed_dir}/uniprot_seed.json'
missing_file = f'{processed_dir}/missing_ids.txt'

df_master = pd.read_csv(master_file, sep='\t')
with open(seed_file, 'r') as f:
    df_seed = pd.DataFrame(json.load(f))
with open(missing_file, 'r') as f:
    missing_ids = [line.strip() for line in f if line.strip()]

print(f"当前主表: {len(df_master)}行")
print(f"待补充: {len(missing_ids)} ID（包含复合ID）")

supplement_entries = []
failed_ids = []

for i, uid_raw in enumerate(missing_ids, 1):
    # 处理复合ID：取第一个
    uid = uid_raw.split('|')[0]  # A6ZKI3|O15255 → A6ZKI3
    
    print(f"[{i}/{len(missing_ids)}] {uid_raw} → 查询{uid}...", end=' ')
    
    url = f'https://rest.uniprot.org/uniprotkb/{uid}.json'
    
    try:
        response = requests.get(url, headers={'Accept': 'application/json'}, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            
            # 提取基本字段
            sequence = data.get('sequence', {}).get('value', '')
            protein_desc = data.get('proteinDescription', {})
            rec_name = protein_desc.get('recommendedName', {})
            protein_name = rec_name.get('fullName', {}).get('value', '') if rec_name else ''
            
            # 构建entry（保留原始复合ID作为uniprot_id）
            entry = {
                'uniprot_id': uid_raw,  # 保留原始"A6ZKI3|O15255"
                'entry_name': data.get('uniProtkbId', ''),
                'protein_name': protein_name,
                'gene_names': ' '.join([g.get('geneName', {}).get('value', '') for g in data.get('genes', [])]),
                'sequence': sequence,
                'sequence_len': len(sequence),
                'mass': data.get('sequence', {}).get('molWeight', 0),
                'function': '',
                'subcellular_location': '',
                'ptms': '',
                'diseases': '',
                'isoforms': '',
                'go_biological_process': '',
                'go_molecular_function': '',
                'go_cellular_component': '',
                'domains': '',
                'pdb_ids': '',
                'string_ids': '',
                'keywords': '; '.join([k.get('name', '') for k in data.get('keywords', [])]),
                'date_modified': data.get('entryAudit', {}).get('lastAnnotationUpdateDate', ''),
                'source': 'UniProtKB (补充)',
                'fetch_date': datetime.now().strftime('%Y-%m-%d')
            }
            
            # 添加HGNC信息
            seed_row = df_seed[df_seed['uniprot_id'] == uid_raw]
            if not seed_row.empty:
                entry['hgnc_id'] = seed_row.iloc[0]['hgnc_id']
                entry['symbol'] = seed_row.iloc[0]['symbol']
            
            supplement_entries.append(entry)
            print(f"✅ {len(sequence)}aa")
        
        elif response.status_code == 404:
            failed_ids.append((uid_raw, f'{uid} 404'))
            print(f"❌ 404")
        else:
            failed_ids.append((uid_raw, f'{uid} {response.status_code}'))
            print(f"❌ {response.status_code}")
    
    except Exception as e:
        failed_ids.append((uid_raw, str(e)[:50]))
        print(f"❌ 异常")
    
    time.sleep(0.5)

print(f"\n=== 补充完成 ===")
print(f"✅ 成功: {len(supplement_entries)}/39")
print(f"❌ 失败: {len(failed_ids)}/39")

if supplement_entries:
    df_supplement = pd.DataFrame(supplement_entries)
    df_master_v2 = pd.concat([df_master, df_supplement], ignore_index=True)
    
    output_file = f'{processed_dir}/protein_master_v2.tsv'
    df_master_v2.to_csv(output_file, sep='\t', index=False)
    
    coverage = len(df_master_v2) / 19253 * 100
    print(f"\n📊 新版主表: {output_file}")
    print(f"   总行数: {len(df_master_v2)} (原{len(df_master)} + 补{len(supplement_entries)})")
    print(f"   覆盖率: {coverage:.2f}%")
else:
    print("\n⚠️ 仍无成功，保持原版protein_master.tsv")

if failed_ids:
    with open(f'{processed_dir}/failed_ids_v2.txt', 'w') as f:
        for uid, reason in failed_ids:
            f.write(f"{uid}\t{reason}\n")

print("\n=== 修复版完成 ===")
