# step4_bulk_final.py - 处理下载的完整UniProt TSV
import pandas as pd
import os
import json
from datetime import datetime

# 路径
base_path = '/Users/pluviophile/graph/1025/data'
raw_dir = os.path.join(base_path, 'raw')
processed_dir = os.path.join(base_path, 'processed')
uniprot_tsv = os.path.join(raw_dir, 'uniprot_swissprot_full.tsv')  # 你下载的文件
seed_file = os.path.join(processed_dir, 'uniprot_seed.json')
master_file = os.path.join(processed_dir, 'protein_master.tsv')

print("=== 步骤4批量: 加载UniProt TSV ===")
# 加载bulk TSV
df_uniprot = pd.read_csv(uniprot_tsv, sep='\t', low_memory=False)
print(f"加载UniProt: {len(df_uniprot)}行, {len(df_uniprot.columns)}列")

# 提取Schema字段（UniProt列名→项目字段）
df_master = pd.DataFrame({
    'uniprot_id': df_uniprot['Entry'],
    'entry_name': df_uniprot['Entry Name'],
    'protein_name': df_uniprot['Protein names'],
    'gene_names': df_uniprot['Gene Names'],
    'sequence': df_uniprot['Sequence'],
    'sequence_len': df_uniprot['Length'],
    'mass': df_uniprot['Mass'],
    'function': df_uniprot['Function [CC]'],
    'subcellular_location': df_uniprot['Subcellular location [CC]'],
    'ptms': df_uniprot['Post-translational modification'],
    'diseases': df_uniprot['Involvement in disease'],
    'isoforms': df_uniprot['Alternative products (isoforms)'],
    'go_biological_process': df_uniprot['Gene Ontology (biological process)'],
    'go_molecular_function': df_uniprot['Gene Ontology (molecular function)'],
    'go_cellular_component': df_uniprot['Gene Ontology (cellular component)'],
    'domains': df_uniprot['Domain [FT]'],
    'pdb_ids': df_uniprot['PDB'],
    'string_ids': df_uniprot['STRING'],
    'keywords': df_uniprot['Keywords'],
    'date_modified': df_uniprot['Date of last modification'],
    'source': 'UniProtKB/Swiss-Prot',
    'fetch_date': datetime.now().strftime('%Y-%m-%d')
})

# merge种子（HGNC info）
with open(seed_file, 'r') as f:
    df_seed = pd.DataFrame(json.load(f))
df_master = df_master.merge(df_seed[['uniprot_id', 'hgnc_id', 'symbol']], 
                             on='uniprot_id', how='inner')  # inner=只保留种子里的

print(f"merge种子: {len(df_master)}行（match HGNC）")

# 输出
df_master.to_csv(master_file, sep='\t', index=False)
print(f"输出: {master_file}，{len(df_master)}行")

# 校验
seq_empty = df_master['sequence'].isna().sum()
func_empty = df_master['function'].isna().sum()
print(f"空序列: {seq_empty} (<1%)")
print(f"空功能: {func_empty} (~10-20%正常)")
print("示例 (P04217):")
print(df_master[df_master['uniprot_id'] == 'P04217'][['uniprot_id', 'protein_name', 'sequence_len', 'function']].to_string())

# 日志
log_file = os.path.join(base_path, 'ETL_log_step4_bulk.txt')
with open(log_file, 'w') as f:
    f.write(f"步骤4批量: 2025-10-26 {datetime.now().strftime('%H:%M')}\n")
    f.write(f"源: UniProt Swiss-Prot bulk TSV (59MB)\n")
    f.write(f"输出: {len(df_master)}行（match {len(df_seed)} HGNC种子）\n")
    f.write(f"空序列: {seq_empty}, 空功能: {func_empty}\n")
print(f"日志: {log_file}")
print("=== 完成！主表protein_master.tsv ready ===")
