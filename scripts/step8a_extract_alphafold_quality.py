# step8a_extract_alphafold_quality.py - 修改版（支持v4-v6）
import pandas as pd
import tarfile
import json
import os
import glob
from tqdm import tqdm

print("=== 步骤8a: 提取AlphaFold质量数据 ===")

base_path = '/Users/pluviophile/graph/1025/data'
processed_dir = f'{base_path}/processed'
output_file = f'{processed_dir}/alphafold_quality.tsv'

# 1. 自动查找tar文件（支持v4-v6）
print("1. 查找AlphaFold压缩包...")
tar_pattern = f'{base_path}/UP000005640_9606_HUMAN_v*.tar'
tar_files = glob.glob(tar_pattern)

if not tar_files:
    print(f"❌ 未找到文件: {tar_pattern}")
    print("\n请确认文件位置：")
    print(f"  当前目录: {base_path}")
    print(f"  期待文件: UP000005640_9606_HUMAN_v4.tar 或 v5.tar 或 v6.tar")
    exit(1)

# 选择最新版本
tar_file = sorted(tar_files)[-1]
print(f"✅ 找到压缩包: {tar_file}")
print(f"   文件大小: {os.path.getsize(tar_file) / 1024 / 1024 / 1024:.2f} GB")

# 提取版本号
version = tar_file.split('_v')[-1].replace('.tar', '')
print(f"   AlphaFold版本: v{version}")

# 2. 解析tar包中的confidence JSON文件
print("\n2. 开始提取质量数据...")
print("   注意：这个过程需要10-20分钟，请耐心等待")

quality_data = []
error_count = 0

try:
    with tarfile.open(tar_file, 'r') as tar:
        # 获取所有confidence JSON文件
        members = [m for m in tar.getmembers() 
                  if 'confidence' in m.name and m.name.endswith('.json')]
        
        print(f"   找到 {len(members)} 个置信度文件")
        
        for member in tqdm(members, desc="   解析进度"):
            try:
                # 提取UniProt ID
                # 文件名格式: AF-P04637-F1-confidence_v6.json
                filename = member.name.split('/')[-1]
                parts = filename.split('-')
                
                if len(parts) >= 2:
                    uniprot_id = parts[1]
                else:
                    continue
                
                # 读取JSON
                f = tar.extractfile(member)
                if f is None:
                    continue
                    
                data = json.load(f)
                
                # 提取pLDDT数组
                plddt_scores = data.get('plddt', [])
                
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
                    'alphafold_version': f'v{version}'
                })
                
            except Exception as e:
                error_count += 1
                if error_count <= 10:
                    print(f"      错误: {member.name} - {str(e)}")
                continue
    
    # 3. 转换为DataFrame并保存
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
    print(f"   AlphaFold版本: v{version}")
    
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
    
except Exception as e:
    print(f"\n❌ 处理失败: {str(e)}")
    import traceback
    traceback.print_exc()
    exit(1)

print(f"\n下一步:")
print(f"   运行: python step8b_update_master_table.py")
