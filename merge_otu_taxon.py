#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
OTU到Taxon转换脚本
将OTU表转换为基于分类学级别的taxon丰度表
"""

import pandas as pd
import numpy as np
import argparse
import sys
import os
import logging
from datetime import datetime

# 设置日志
def setup_logger(log_file):
    logger = logging.getLogger('otu_to_taxon')
    logger.setLevel(logging.INFO)
    
    # 创建文件处理器
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    
    # 创建控制台处理器
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    # 创建格式化程序
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    # 将处理器添加到记录器
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger

def convert_otu_to_taxon(input_file, output_file, log_file):
    """
    将OTU表转换为taxon表，按照分类学级别合并
    
    参数:
    -----------
    input_file : str
        包含物种注释的OTU表输入文件路径
    output_file : str
        taxon表保存路径
    log_file : str
        日志文件保存路径
    """
    # 设置日志
    logger = setup_logger(log_file)
    logger.info(f"开始处理OTU表: {input_file}")
    
    # 读取OTU表
    try:
        otu_table = pd.read_csv(input_file, sep='\t')
        logger.info(f"成功读取OTU表，包含{otu_table.shape[0]}个OTU和{otu_table.shape[1]}列")
    except Exception as e:
        logger.error(f"读取输入文件错误: {e}")
        sys.exit(1)
    
    # 检查表格是否具有预期结构
    if len(otu_table.columns) <= 8:
        logger.error(f"错误: 输入表格列数少于9列")
        sys.exit(1)
    
    # 获取分类学列名和样本列名
    taxonomy_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    sample_cols = list(otu_table.columns[8:])
    
    logger.info(f"分类学级别: {', '.join(taxonomy_levels)}")
    logger.info(f"样本数量: {len(sample_cols)}")
    
    # 创建空的taxon表和收集合并记录的字典
    taxon_table = []
    taxon_records = {}
    
    # 从最精确的分类级别开始向上迭代
    for level_idx, level in enumerate(reversed(taxonomy_levels)):
        level_idx = 6 - level_idx  # 转换为实际列索引(从0开始)
        level_col = level_idx + 1  # OTU表中的实际列索引(从1开始)
        
        logger.info(f"处理分类学级别: {level} (列索引: {level_col})")
        
        # 筛选尚未分配到更精确级别的OTU
        if level_idx == 6:  # Species级别
            unassigned_otus = otu_table.copy()
        else:
            assigned_otu_ids = [record for records in taxon_records.values() for record in records]
            unassigned_otus = otu_table[~otu_table['#OTU_ID'].isin(assigned_otu_ids)]
        
        logger.info(f"在{level}级别有{len(unassigned_otus)}个未分配的OTU")
        
        # 分组并收集非'Unassigned'的分类单元
        for taxon_name, group in unassigned_otus.groupby(level):
            if taxon_name != 'Unassigned' and not pd.isna(taxon_name):
                # 收集这个taxon包含的OTU
                otu_ids = group['#OTU_ID'].tolist()
                taxon_records[f"{level}:{taxon_name}"] = otu_ids
                
                # 汇总样本计数
                abundance = group[sample_cols].sum()
                
                # 将此taxon添加到结果表中
                taxon_entry = {
                    'Classification': level,
                    'Taxon': taxon_name
                }
                for col in sample_cols:
                    taxon_entry[col] = abundance[col]
                
                taxon_table.append(taxon_entry)
                
                logger.info(f"在{level}级别合并了{len(otu_ids)}个OTU到taxon '{taxon_name}'")
                logger.info(f"  包含的OTU: {', '.join(otu_ids)}")
    
    # 将结果转换为DataFrame并保存
    result = pd.DataFrame(taxon_table)
    
    # 按Classification级别排序（Kingdom, Phylum, Class...）
    level_order = {level: i for i, level in enumerate(taxonomy_levels)}
    result['level_order'] = result['Classification'].map(level_order)
    result = result.sort_values('level_order').drop('level_order', axis=1)
    
    # 保存结果
    result.to_csv(output_file, sep='\t', index=False)
    logger.info(f"已将taxon表保存到{output_file}")
    
    # 记录合并统计信息
    logger.info(f"合并统计信息:")
    for level in taxonomy_levels:
        level_count = result[result['Classification'] == level].shape[0]
        logger.info(f"  {level}级别: {level_count}个taxon")
    
    total_otus = sum(len(otus) for otus in taxon_records.values())
    logger.info(f"总共合并了{total_otus}个OTU为{len(taxon_table)}个taxon")
    
    # 记录未分配的OTU
    all_assigned_otus = [otu for otus in taxon_records.values() for otu in otus]
    unassigned_otus = set(otu_table['#OTU_ID']) - set(all_assigned_otus)
    if unassigned_otus:
        logger.warning(f"有{len(unassigned_otus)}个OTU未被分配到任何taxon")
        logger.warning(f"未分配的OTU: {', '.join(list(unassigned_otus))}")
    
    logger.info("处理完成")
    return result, taxon_records

def main():
    parser = argparse.ArgumentParser(description='将OTU表转换为基于分类学级别的taxon丰度表')
    parser.add_argument('-i', '--input', required=True, help='输入OTU表文件（包含物种分类信息）')
    parser.add_argument('-o', '--output', required=True, help='输出taxon丰度表文件')
    parser.add_argument('-l', '--log', help='日志文件路径（默认: 与输出文件同名加.log后缀）')
    
    args = parser.parse_args()
    
    # 如果未指定日志文件，使用默认名称
    if not args.log:
        args.log = os.path.splitext(args.output)[0] + '.log'
    
    convert_otu_to_taxon(args.input, args.output, args.log)

if __name__ == "__main__":
    main()