import pandas as pd
import numpy as np
import argparse # 导入argparse模块
import sys

# 定义分类级别列名
taxonomy_columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

# 读取文件
def read_files(taxonomy_file, otutab_file):
    """读取分类信息文件和OTU特征表。"""
    # 读取分类信息文件
    # 假设文件是以制表符分隔
    taxonomy_df = pd.read_csv(taxonomy_file, sep='\t')
    
    # 读取OTU特征表
    # 假设文件是以制表符分隔
    otutab_df = pd.read_csv(otutab_file, sep='\t')
    
    return taxonomy_df, otutab_df

# 合并文件
def merge_tables(taxonomy_df, otutab_df):
    """
    合并分类信息和OTU丰度表
    """
    # 1. 自动检测并统一 ID 列名
    # 定义可能的候选列名
    id_variants = ['OTU_ID', 'OTUID', '#OTU ID', '#OTU_ID', 'ID']
    
    # 统一 taxonomy_df 的 ID 列
    for col in id_variants:
        if col in taxonomy_df.columns:
            taxonomy_df = taxonomy_df.rename(columns={col: '#OTU_ID'})
            break
            
    # 统一 otutab_df 的 ID 列
    for col in id_variants:
        if col in otutab_df.columns:
            otutab_df = otutab_df.rename(columns={col: '#OTU_ID'})
            break

    # 2. 检查 ID 列是否存在
    if '#OTU_ID' not in taxonomy_df.columns or '#OTU_ID' not in otutab_df.columns:
        print(f"Error: 找不到 ID 列。Taxonomy 列名: {list(taxonomy_df.columns)}, OTU Table 列名: {list(otutab_df.columns)}")
        sys.exit(1)
    
    # 3. 执行合并
    merged_df = pd.merge(taxonomy_df, otutab_df, on='#OTU_ID', how='right')
    
    # 4. 重新排列列顺序
    columns_order = ['#OTU_ID'] + taxonomy_columns + [col for col in merged_df.columns 
                                                    if col not in ['#OTU_ID'] + taxonomy_columns]
    merged_df = merged_df[columns_order]
    
    return merged_df

def main():
    """主函数，处理命令行参数并执行文件合并操作。"""
    # 1. 设置命令行参数解析器
    parser = argparse.ArgumentParser(
        description="将分类信息文件和OTU特征表合并，并输出结果。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # 添加参数
    parser.add_argument(
        '-t', '--taxonomy', 
        type=str, 
        required=True, 
        help="分类信息文件的路径 (e.g., algae_taxonomy_ASV.txt)"
    )
    parser.add_argument(
        '-o', '--otutab', 
        type=str, 
        required=True, 
        help="OTU特征表的路径 (e.g., algae_ASV_otutab_rare.txt)"
    )
    parser.add_argument(
        '-out', '--output', 
        type=str, 
        required=True, 
        help="合并后结果文件的输出路径 (e.g., merged_algae_ASV_taxonomy.txt)"
    )
    
    # 解析参数
    args = parser.parse_args()
    
    taxonomy_file = args.taxonomy
    otutab_file = args.otutab
    output_file = args.output
    
    print(f"输入分类文件: {taxonomy_file}")
    print(f"输入OTU表文件: {otutab_file}")
    print(f"输出文件: {output_file}")
    print("-" * 30)

    try:
        # 读取文件
        print("正在读取文件...")
        taxonomy_df, otutab_df = read_files(taxonomy_file, otutab_file)
        
        # 合并表格
        print("正在合并表格...")
        merged_df = merge_tables(taxonomy_df, otutab_df)
        
        # 保存结果
        print("正在保存结果...")
        # 保存为制表符分隔文件，不包含索引
        merged_df.to_csv(output_file, sep='\t', index=False)
        
        # 输出一些基本统计信息
        print("\n🎉 处理完成!")
        print(f"总行数 (ASV/OTU 数量): {len(merged_df)}")
        print(f"总列数: {len(merged_df.columns)}")
        # 减去 #OTU ID (1) 和 taxonomy_columns (7) 
        print(f"样本数: {len(merged_df.columns) - len(taxonomy_columns) - 1}")
        
        # 检查是否有缺失的分类信息
        missing_taxonomy = merged_df[taxonomy_columns].isna().sum()
        print("\n每个分类级别的缺失值数量:")
        print(missing_taxonomy)
        
    except FileNotFoundError:
        print("\n❌ 错误: 找不到指定的输入文件。请检查文件路径是否正确。")
    except Exception as e:
        print(f"\n❌ 处理过程中出现错误: {str(e)}")

if __name__ == "__main__":
    main()