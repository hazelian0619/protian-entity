# step7_ptm_fixed.py - PhosphoSitePlus PTM位点提取（修正版）
import pandas as pd
from datetime import datetime
import re

print("=== 步骤7: PhosphoSitePlus PTM位点（修正版） ===")

# 路径
base_path = '/Users/pluviophile/graph/1025/data'
raw_dir = f'{base_path}/raw'
processed_dir = f'{base_path}/processed'
ptm_file = f'{raw_dir}/Phosphorylation_site_dataset'
master_file = f'{processed_dir}/protein_master_v2.tsv'
output_file = f'{processed_dir}/ptm_sites.tsv'

# 1. 加载主表
print("1. 加载protein_master_v2.tsv...")
df_master = pd.read_csv(master_file, sep='\t')
our_proteins = set(df_master['uniprot_id'])
print(f"   主表: {len(our_proteins)}个蛋白")

# 2. 加载PTM文件（跳过版权说明，找到表头行）
print("2. 加载PhosphoSitePlus文件...")
try:
    # 读取前几行找表头
    with open(ptm_file, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    # 找到包含GENE的行作为表头
    header_line = 0
    for i, line in enumerate(lines):
        if line.startswith('GENE'):
            header_line = i
            break
    
    print(f"   表头在第 {header_line + 1} 行")
    
    # 从表头行开始读取
    df_ptm = pd.read_csv(ptm_file, sep='\t', skiprows=header_line)
    print(f"   原始行数: {len(df_ptm):,} (全物种)")
    print(f"   列名: {df_ptm.columns.tolist()}")
    
except Exception as e:
    print(f"   ❌ 文件读取失败: {e}")
    exit(1)

# 3. 过滤人类数据（注意：是'human'不是'Homo sapiens'）
print("3. 过滤人类PTM...")
df_human = df_ptm[df_ptm['ORGANISM'] == 'human'].copy()
print(f"   人类PTM: {len(df_human):,}行")

# 4. 提取UniProt ID（去掉异构体后缀）
print("4. 提取UniProt ID...")
df_human['uniprot_id'] = df_human['ACC_ID'].str.split('-').str[0]

# 5. 只保留主表中的蛋白
print("5. 保留主表蛋白...")
df_filtered = df_human[df_human['uniprot_id'].isin(our_proteins)].copy()
print(f"   匹配主表: {len(df_filtered):,}行")

# 6. 解析MOD_RSD列（如S119-p → S, 119, Phosphorylation）
print("6. 解析修饰位点...")

def parse_mod_rsd(mod_rsd):
    """
    解析MOD_RSD格式
    S119-p → (S, 119, Phosphorylation)
    T2-p → (T, 2, Phosphorylation)
    K382-ac → (K, 382, Acetylation)
    """
    if pd.isna(mod_rsd) or mod_rsd == '':
        return None, None, None
    
    # 匹配模式：字母+数字+修饰类型
    match = re.match(r'([A-Z])(\d+)-(\w+)', str(mod_rsd))
    if match:
        residue = match.group(1)
        position = match.group(2)
        mod_type_code = match.group(3)
        
        # 修饰类型映射
        mod_type_map = {
            'p': 'Phosphorylation',
            'ub': 'Ubiquitination',
            'ac': 'Acetylation',
            'me': 'Methylation',
            'me1': 'Mono-methylation',
            'me2': 'Di-methylation',
            'me3': 'Tri-methylation',
            'sm': 'Sumoylation',
            'gl': 'O-GlcNAcylation'
        }
        ptm_type = mod_type_map.get(mod_type_code, mod_type_code)
        
        return residue, position, ptm_type
    return None, None, None

# 应用解析函数
df_filtered[['residue', 'position', 'ptm_type']] = df_filtered['MOD_RSD'].apply(
    lambda x: pd.Series(parse_mod_rsd(x))
)

# 去掉解析失败的行
df_filtered = df_filtered.dropna(subset=['residue', 'position', 'ptm_type'])
print(f"   解析成功: {len(df_filtered):,}行")

# 7. 构建输出DataFrame
print("7. 数据清理...")
df_output = pd.DataFrame({
    'uniprot_id': df_filtered['uniprot_id'],
    'gene_name': df_filtered['GENE'],
    'ptm_type': df_filtered['ptm_type'],
    'residue': df_filtered['residue'],
    'position': df_filtered['position'],
    'site_sequence': df_filtered['SITE_+/-7_AA'],  # ±7个氨基酸序列
    'domain': df_filtered['DOMAIN'],
    'lit_references': df_filtered['LT_LIT'].fillna(0).astype(int),  # 文献数
    'ms_references': df_filtered['MS_LIT'].fillna(0).astype(int),   # 质谱文献数
    'source': 'PhosphoSitePlus',
    'fetch_date': datetime.now().strftime('%Y-%m-%d')
})

# 8. 去重
print("8. 去重处理...")
df_output = df_output.drop_duplicates(subset=['uniprot_id', 'residue', 'position'])
print(f"   去重后: {len(df_output):,}行")

# 9. 输出
df_output.to_csv(output_file, sep='\t', index=False)
print(f"\n✅ 输出: {output_file}")
print(f"   总行数: {len(df_output):,}")
print(f"   唯一蛋白数: {df_output['uniprot_id'].nunique()}")
print(f"   唯一基因数: {df_output['gene_name'].nunique()}")

# 10. 示例
print(f"\n   示例前5行:")
sample = df_output.head(5)[['uniprot_id', 'gene_name', 'ptm_type', 'residue', 'position']]
for _, row in sample.iterrows():
    print(f"   {row['uniprot_id']} ({row['gene_name']}) → {row['ptm_type']} {row['residue']}{row['position']}")

# 11. 统计
print(f"\n   修饰类型分布:")
ptm_type_counts = df_output['ptm_type'].value_counts()
for ptm_type, count in ptm_type_counts.head(10).items():
    print(f"   - {ptm_type}: {count:,}个位点")

avg_sites = len(df_output) / df_output['uniprot_id'].nunique()
print(f"\n   平均每个蛋白 {avg_sites:.1f} 个修饰位点")

# Top修饰蛋白
print(f"\n   修饰位点最多的前5个蛋白:")
top_proteins = df_output.groupby(['gene_name', 'uniprot_id']).size().sort_values(ascending=False).head(5)
for (gene, uniprot), count in top_proteins.items():
    print(f"   - {gene} ({uniprot}): {count}个位点")

print("\n=== 步骤7完成！PTM位点ready ===")
print(f"完成度: 93% → 95%+ (+2-5%)")
