# step8a_extract_alphafold_quality_v6.py - 适配v6版本
import pandas as pd
import tarfile
import gzip
import os
from tqdm import tqdm

print("=== 步骤8a: 提取AlphaFold质量数据 (v6版本) ===")

base_path = '/Users/pluviophile/graph/1025/data'
processed_dir = f'{base_path}/processed'
tar_file = f'{base_path}/UP000005640_9606_HUMAN_v6.tar'
output_file = f'{processed_dir}/alphafold_quality.tsv'

# 1. 检查文件
if not os.path.exists(tar_file):
    print(f"❌ 未找到文件: {tar_file}")
    exit(1)

print(f"✅ 找到压缩包: {tar_file}")
print(f"   文件大小: {os.path.getsize(tar_file) / 1024 / 1024 / 1024:.2f} GB")

# 2. 解析tar包中的.pdb.gz文件
print("\n开始提取质量数据...")
print("注意：v6版本文件是.gz压缩，需要解压后提取pLDDT")
print("预计时间：15-25分钟\n")

quality_data = []
error_count = 0

try:
    with tarfile.open(tar_file, 'r') as tar:
        # 获取所有.pdb.gz文件
        members = [m for m in tar.getmembers() 
                  if m.name.endswith('.pdb.gz') and m.isfile()]
        
        print(f"找到 {len(members)} 个PDB.GZ文件\n")
        
        for member in tqdm(members, desc="解析进度"):
            try:
                # 提取UniProt ID
                # 文件名: AF-A0A024R1R8-F1-model_v6.pdb.gz
                filename = member.name.split('/')[-1]
                parts = filename.split('-')
                
                if len(parts) < 2:
                    continue
                    
                uniprot_id = parts[1]
                
                # 提取并解压文件
                f = tar.extractfile(member)
                if f is None:
                    continue
                
                # 解压gzip
                pdb_content = gzip.decompress(f.read()).decode('utf-8')
                
                # 从B-factor列提取pLDDT分数
                plddt_scores = []
                
                for line in pdb_content.split('\n'):
                    if line.startswith('ATOM'):
                        try:
                            # PDB格式：B-factor在第61-66列
                            plddt = float(line[60:66].strip())
                            plddt_scores.append(plddt)
                        except:
                            continue
                
                if not plddt_scores:
                    continue
                
                # 计算统计数据
                mean_plddt = sum(plddt_scores) / len(plddt_scores)
                median_plddt = sorted(plddt_scores)[len(plddt_scores) // 2]
                min_plddt = min(plddt_scores)
                max_plddt = max(plddt_scores)
                
                # 按置信度分类
                very_high = sum(1 for s in plddt_scores if s >= 90)
                high = sum(1 for s in plddt_scores if 70 <= s < 90)
                low = sum(1 for s in plddt_scores if 50 <= s < 70)
                very_low = sum(1 for s in plddt_scores if s < 50)
                total = len(plddt_scores)
                
                # 质量分级
                if mean_plddt >= 90:
                    quality_grade = 'Excellent'
                elif mean_plddt >= 70:
                    quality_grade = 'Good'
                elif mean_plddt >= 50:
                    quality_grade = 'Medium'
                else:
                    quality_grade = 'Poor'
                
                # 保存数据
                quality_data.append({
                    'uniprot_id': uniprot_id,
                    'mean_plddt': round(mean_plddt, 2),
                    'median_plddt': round(median_plddt, 2),
                    'min_plddt': round(min_plddt, 2),
                    'max_plddt': round(max_plddt, 2),
                    'very_high_conf_pct': round(very_high / total * 100, 1),
                    'high_conf_pct': round(high / total * 100, 1),
                    'low_conf_pct': round(low / total * 100, 1),
                    'very_low_conf_pct': round(very_low / total * 100, 1),
                    'quality_grade': quality_grade,
                    'residue_count': total,
                    'alphafold_version': 'v6'
                })
                
            except Exception as e:
                error_count += 1
                if error_count <= 10:
                    print(f"\n   警告: {member.name} - {str(e)}")
                continue
    
    # 3. 保存结果
    if not quality_data:
        print("\n❌ 未提取到任何数据！")
        exit(1)
    
    df_quality = pd.DataFrame(quality_data)
    df_quality.to_csv(output_file, sep='\t', index=False)
    
    print(f"\n✅ 提取完成！")
    print(f"\n输出信息:")
    print(f"   文件: {output_file}")
    print(f"   总行数: {len(df_quality)}")
    print(f"   错误数: {error_count}")
    
    # 4. 统计信息
    print(f"\n质量分布:")
    quality_counts = df_quality['quality_grade'].value_counts()
    for grade in ['Excellent', 'Good', 'Medium', 'Poor']:
        count = quality_counts.get(grade, 0)
        if count > 0:
            print(f"   {grade}: {count} ({count/len(df_quality)*100:.1f}%)")
    
    print(f"\n平均pLDDT统计:")
    print(f"   平均值: {df_quality['mean_plddt'].mean():.2f}")
    print(f"   中位数: {df_quality['mean_plddt'].median():.2f}")
    print(f"   最小值: {df_quality['mean_plddt'].min():.2f}")
    print(f"   最大值: {df_quality['mean_plddt'].max():.2f}")
    
    print(f"\n高置信度残基统计:")
    print(f"   >90分平均比例: {df_quality['very_high_conf_pct'].mean():.1f}%")
    print(f"   70-90分平均比例: {df_quality['high_conf_pct'].mean():.1f}%")
    
    # 5. 保存checkpoint（每5000个）
    print(f"\n提示：可以删除原始tar包节省空间")
    print(f"   提取的数据已保存到: {output_file}")
    
except Exception as e:
    print(f"\n❌ 处理失败: {str(e)}")
    import traceback
    traceback.print_exc()
    exit(1)

print(f"\n下一步:")
print(f"   运行: python step8b_update_master_table.py")
