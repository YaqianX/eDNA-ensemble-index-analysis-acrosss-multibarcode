#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OTU / Taxon 特征表相对丰度标准化脚本  (v2.0)
=============================================

适用输入:
  - 纯数值 OTU 表 (第一列为 #OTU_ID，后续列全为样品计数)
  - merge_otu_taxon.py 输出的 taxon 丰度表
    (第一列 Classification 为索引，第二列 Taxon 为字符串，其余为样品计数)
  - 任何含混合列（分类学注释 + 数值计数）的特征表

处理逻辑:
  1. 读取文件，以第一列为行索引
  2. 自动检测并隔离非数值列（如 Taxon、Kingdom 等）
  3. 仅对数值列进行列求和归一化 (0-1 或 百分比)
  4. 输出时将非数值列原样置于数值列之前

用法示例:
  python otutab_normalize.py -i taxon_abundance.txt -o taxon_normalized.txt
  python otutab_normalize.py -i otutab.txt -o normalized.txt -m relative
  python otutab_normalize.py -i taxon_abundance.txt -m percentage -f '%.6f'
"""

import pandas as pd
import numpy as np
import sys
import os
import argparse


def normalize_otu_table(input_file, output_file=None, separator='\t',
                        output_format='%.8f', method='relative'):
    """
    对 OTU / Taxon 特征表进行相对丰度标准化。

    自动处理含非数值列（分类学注释）的混合格式表格：
    非数值列在标准化过程中被隔离保留，最终原样写入输出文件。

    Parameters
    ----------
    input_file    : str   输入文件路径
    output_file   : str   输出文件路径（None 时自动生成）
    separator     : str   字段分隔符，默认制表符
    output_format : str   数值输出格式，默认 '%.8f'
    method        : str   'relative'（0-1）或 'percentage'（0-100）

    Returns
    -------
    pandas.DataFrame  标准化后的完整数据框（含非数值列），失败时返回 None
    """

    try:
        print(f"正在读取文件: {input_file}")
        print(f"使用分隔符: {'制表符' if separator == chr(9) else repr(separator)}")

        raw_table = pd.read_csv(input_file, sep=separator, index_col=0)
        print(f"原始数据维度: {raw_table.shape[0]} 行 × {raw_table.shape[1]} 列")
        print(f"所有列名: {list(raw_table.columns)}")

        # ── 分离数值列与非数值列 ──────────────────────────────────────────
        numeric_cols = [c for c in raw_table.columns
                        if pd.api.types.is_numeric_dtype(raw_table[c])]
        meta_cols    = [c for c in raw_table.columns
                        if not pd.api.types.is_numeric_dtype(raw_table[c])]

        if meta_cols:
            print(f"\n检测到非数值（注释）列，将跳过标准化并原样保留: {meta_cols}")

        if not numeric_cols:
            print("错误: 未检测到任何数值列，无法进行标准化。")
            return None

        otu_numeric = raw_table[numeric_cols].astype(float)
        otu_meta    = raw_table[meta_cols] if meta_cols else None

        print(f"参与标准化的样品列 ({len(numeric_cols)} 个): {numeric_cols}")

        if (otu_numeric < 0).any().any():
            print("警告: 数值列中包含负值，请确认输入数据是否正确。")

        # ── 列求和 ────────────────────────────────────────────────────────
        sample_sums = otu_numeric.sum(axis=0)
        print(f"\n各样品的总计数:")
        for sample, total in sample_sums.items():
            print(f"  {sample}: {total:,.2f}")

        zero_samples = sample_sums[sample_sums == 0]
        if len(zero_samples) > 0:
            print(f"\n警告: 以下样品总计数为 0，相对丰度将全部填充为 0: "
                  f"{list(zero_samples.index)}")

        # ── 标准化 ────────────────────────────────────────────────────────
        print(f"\n正在进行 {method} 标准化...")
        normalized_numeric = otu_numeric.div(sample_sums, axis=1)

        if method == 'percentage':
            normalized_numeric = normalized_numeric * 100
            print("已转换为百分比形式 (×100)")

        normalized_numeric = normalized_numeric.fillna(0)

        # ── 验证 ──────────────────────────────────────────────────────────
        expected = 100.0 if method == 'percentage' else 1.0
        print(f"\n标准化后各样品总和（期望值 ≈ {expected}）:")
        for sample, s in normalized_numeric.sum(axis=0).items():
            tol = 1e-3 if method == 'percentage' else 1e-6
            flag = "" if abs(s - expected) < tol or sample_sums[sample] == 0 else " ← 偏差超出容差，请检查"
            print(f"  {sample}: {s:.8f}{flag}")

        # ── 统计信息 ──────────────────────────────────────────────────────
        print(f"\n数据统计:")
        print(f"  标准化前数值范围: {otu_numeric.min().min():.4f} – "
              f"{otu_numeric.max().max():.4f}")
        print(f"  标准化后数值范围: {normalized_numeric.min().min():.8f} – "
              f"{normalized_numeric.max().max():.8f}")

        # ── 拼合非数值列 ──────────────────────────────────────────────────
        if otu_meta is not None:
            output_table = pd.concat([otu_meta, normalized_numeric], axis=1)
            print(f"\n非数值列 {meta_cols} 已原样置于输出表格首列。")
        else:
            output_table = normalized_numeric

        # ── 确定输出路径 ──────────────────────────────────────────────────
        if output_file is None:
            base = os.path.splitext(input_file)[0]
            suffix = "_percentage" if method == 'percentage' else "_normalized"
            output_file = f"{base}{suffix}.txt"

        print(f"\n正在保存标准化表到: {output_file}")
        output_table.to_csv(output_file, sep=separator, float_format=output_format)
        print("标准化完成！")

        return output_table

    except FileNotFoundError:
        print(f"错误: 找不到输入文件 '{input_file}'")
        return None
    except pd.errors.EmptyDataError:
        print(f"错误: 输入文件 '{input_file}' 为空")
        return None
    except Exception as e:
        print(f"处理过程中出现错误: {str(e)}")
        import traceback
        traceback.print_exc()
        return None


def display_top_taxa(output_table, top_n=10):
    """显示数值列均值最高的前 top_n 个分类群。"""
    if output_table is None:
        return
    numeric_cols = [c for c in output_table.columns
                    if pd.api.types.is_numeric_dtype(output_table[c])]
    if not numeric_cols:
        return
    mean_abund = output_table[numeric_cols].mean(axis=1).sort_values(ascending=False)
    print(f"\n数值列均值最高的前 {top_n} 个分类群:")
    print("-" * 50)
    for i, (idx, val) in enumerate(mean_abund.head(top_n).items(), 1):
        print(f"  {i:2d}. {idx}: {val:.8f}")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='OTU / Taxon 特征表相对丰度标准化工具 (v2.0)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
适用场景:
  - 纯 OTU 计数表（传统格式）
  - merge_otu_taxon.py 输出的含 Classification/Taxon 列的混合格式表

使用示例:
  python otutab_normalize.py -i taxon_abundance.txt -o taxon_normalized.txt
  python otutab_normalize.py -i otutab.txt -m relative -f '%.8f'
  python otutab_normalize.py -i merged_table.txt -m percentage -q
        """
    )
    parser.add_argument('-i', '--input', required=True,
                        help='输入特征表文件路径')
    parser.add_argument('-o', '--output',
                        help='输出文件路径（默认: 输入文件名_normalized.txt）')
    parser.add_argument('-s', '--separator', default='\t',
                        choices=['\t', ',', ';', ' '],
                        help='字段分隔符（默认: 制表符）')
    parser.add_argument('-m', '--method', default='relative',
                        choices=['relative', 'percentage'],
                        help='标准化方法: relative(0-1) 或 percentage(0-100)（默认: relative）')
    parser.add_argument('-f', '--format', default='%.8f',
                        help='数值输出格式（默认: %%.8f）')
    parser.add_argument('-t', '--top', type=int, default=10,
                        help='显示均值最高的前 N 个分类群（默认: 10）')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='静默模式，仅显示错误信息')
    parser.add_argument('--version', action='version', version='OTU标准化工具 v2.0')
    return parser.parse_args()


def main():
    args = parse_arguments()

    if not args.quiet:
        print("=" * 60)
        print("OTU / Taxon 特征表相对丰度标准化工具  v2.0")
        print("=" * 60)
        print(f"输入文件  : {args.input}")
        print(f"输出文件  : {args.output if args.output else '自动生成'}")
        print(f"分隔符    : {'制表符' if args.separator == chr(9) else repr(args.separator)}")
        print(f"标准化方法: {args.method}")
        print(f"数值格式  : {args.format}")
        print("-" * 60)

    if not os.path.exists(args.input):
        print(f"错误: 输入文件 '{args.input}' 不存在")
        return 1

    output_table = normalize_otu_table(
        input_file=args.input,
        output_file=args.output,
        separator=args.separator,
        output_format=args.format,
        method=args.method
    )

    if output_table is not None and not args.quiet:
        display_top_taxa(output_table, args.top)
        return 0
    return 0 if output_table is not None else 1


if __name__ == "__main__":
    sys.exit(main())
