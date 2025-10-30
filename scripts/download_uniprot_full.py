# download_uniprot_full.py - 下载完整UniProt TSV（含Sequence等所有字段）
import requests
import gzip
import shutil

# UniProt REST API - Swiss-Prot人类蛋白，指定所有重要字段
fields = [
    'accession', 'id', 'reviewed', 'protein_name', 'gene_names', 'organism_name',
    'length', 'sequence', 'mass', 'cc_function', 'cc_subcellular_location',
    'cc_ptm', 'cc_disease', 'cc_alternative_products', 'go_p', 'go_f', 'go_c',
    'ft_domain', 'xref_pdb', 'xref_string', 'date_modified', 'keyword'
]
fields_param = ','.join(fields)
url = f"https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv&query=(proteome:UP000005640)+AND+(reviewed:true)&fields={fields_param}"

output_gz = "/Users/pluviophile/graph/1025/data/raw/uniprot_swissprot_full.tsv.gz"
output = output_gz.replace('.gz', '')

print("下载中...（~500MB，5-10分钟）")
response = requests.get(url, stream=True)
with open(output_gz, 'wb') as f:
    for chunk in response.iter_content(chunk_size=8192):
        f.write(chunk)
print(f"下载完成: {output_gz}")

print("解压中...")
with gzip.open(output_gz, 'rb') as f_in:
    with open(output, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
print(f"完成！文件: {output}（应该~500MB，含Sequence列）")

# 验证
import pandas as pd
df = pd.read_csv(output, sep='\t', nrows=5)
print("文件列数:", len(df.columns))
print("列名预览:", df.columns.tolist()[:10])
if 'Sequence' in df.columns:
    print("✅ 成功！包含Sequence列")
else:
    print("❌ 警告：缺少Sequence列")
