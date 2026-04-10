#!/usr/bin/env python3
"""
eDNA Ensemble Index Pipeline
=============================
实现从归一化相对丰度表到多条形码合并 eDNA ensemble index 的完整流程。

参考方法：
  Kelly et al. (2019), Djurhuus et al. (2020)

流程：
  Step 1a : 按样品前缀分组 → 各基质组内均值相对丰度
  Step 1b : ÷ 该分类群在该条形码所有基质中的最大值 → eDNA index
  Step 2  : 跨条形码合并 → 共有分类群取非零值均值 → ensemble index

用法示例
--------
# 单条形码（仅输出 Step1 结果）
python edna_ensemble_pipeline.py \\
    --inputs 18S:data/18s_normalize.txt \\
    --matrix-map "B:Biofilm,M:Membrane,S:Sediment" \\
    --outdir results/

# 多条形码合并（完整流程）
python edna_ensemble_pipeline.py \\
    --inputs 18S:data/18s_normalize.txt COI:data/co1_normalize.txt \\
    --matrix-map "B:Biofilm,M:Membrane,S:Sediment" \\
    --exclude-matrix "COI:Biofilm" \\
    --outdir results/

# 从已有的 Step1 文件直接进行 Step2 合并
python edna_ensemble_pipeline.py \\
    --step1-inputs 18S:data/18s_step1.txt COI:data/co1_step1.txt \\
    --exclude-matrix "COI:Biofilm" \\
    --outdir results/

# 自定义输出文件名（三种方式可组合使用）
python edna_ensemble_pipeline.py \\
    --inputs 18S:data/18s_normalize.txt COI:data/co1_normalize.txt \\
    --matrix-map "B:Biofilm,M:Membrane,S:Sediment" \\
    --exclude-matrix "COI:Biofilm" \\
    --prefix benthic_2024 \\
    --output-merged benthic_ensemble_final.txt \\
    --out-step1a "18S:benthic_18s_mean.txt,COI:benthic_co1_mean.txt" \\
    --out-step1b "18S:benthic_18s_index.txt,COI:benthic_co1_index.txt" \\
    --outdir results/
"""

import argparse
import sys
import os
import textwrap
from pathlib import Path

import pandas as pd
import numpy as np


# ═══════════════════════════════════════════════════════════════════════════════
# 工具函数
# ═══════════════════════════════════════════════════════════════════════════════

def log(msg: str, level: str = "INFO"):
    prefix = {"INFO": "[ INFO ]", "WARN": "[ WARN ]", "STEP": "[ STEP ]",
               "OK":   "[  OK  ]", "ERR":  "[ ERR  ]"}.get(level, "[      ]")
    print(f"{prefix} {msg}", flush=True)


def read_table(path: str, index_col: str = "Taxon") -> pd.DataFrame:
    """读取 TSV/CSV 文件，自动判断分隔符，以 index_col 为行索引。"""
    path = str(path)
    sep = "\t" if path.endswith(".txt") or path.endswith(".tsv") else ","
    df = pd.read_csv(path, sep=sep, index_col=0)

    # 若第一列不叫 Taxon 但存在 Taxon 列，则重设
    if index_col in df.columns:
        df = df.set_index(index_col)

    # 删除非数值列（如 Classification）
    non_numeric = [c for c in df.columns if not pd.api.types.is_numeric_dtype(df[c])]
    if non_numeric:
        log(f"  忽略非数值列: {non_numeric}", "WARN")
        df = df.drop(columns=non_numeric)

    return df.astype(float)


def parse_matrix_map(matrix_map_str: str) -> dict:
    """
    解析样品前缀 → 基质名称的映射字符串。
    格式: "B:Biofilm,M:Membrane,S:Sediment"
    支持多字符前缀，如 "Bio:Biofilm,Mem:Membrane,Sed:Sediment"
    """
    mapping = {}
    for item in matrix_map_str.split(","):
        item = item.strip()
        if ":" not in item:
            raise ValueError(f"格式错误，期望 prefix:name，得到: '{item}'")
        prefix, name = item.split(":", 1)
        mapping[prefix.strip()] = name.strip()
    return mapping


def infer_matrix(col: str, prefix_map: dict) -> str | None:
    """根据列名前缀推断基质类型，前缀越长优先级越高。"""
    matched = None
    best_len = 0
    for prefix, name in prefix_map.items():
        if col.startswith(prefix) and len(prefix) > best_len:
            matched = name
            best_len = len(prefix)
    return matched


def parse_exclude(exclude_list: list) -> dict:
    """
    解析需排除的 (条形码, 基质) 组合。
    格式: ["COI:Biofilm", "COI:Membrane"]
    返回: {"COI": {"Biofilm", "Membrane"}}
    """
    result = {}
    for item in (exclude_list or []):
        if ":" not in item:
            raise ValueError(f"--exclude-matrix 格式错误，期望 BARCODE:Matrix，得到: '{item}'")
        bc, mat = item.split(":", 1)
        result.setdefault(bc.strip(), set()).add(mat.strip())
    return result


def parse_name_map(name_str: str) -> dict:
    """
    解析 BARCODE:filename 格式的自定义文件名字符串。
    格式: "18S:my_18s_mean.txt,COI:my_co1_mean.txt"
    返回: {"18S": "my_18s_mean.txt", "COI": "my_co1_mean.txt"}
    """
    result = {}
    if not name_str:
        return result
    for item in name_str.split(","):
        item = item.strip()
        if ":" not in item:
            raise ValueError(f"文件名映射格式错误，期望 BARCODE:filename，得到: '{item}'")
        bc, fname = item.split(":", 1)
        result[bc.strip()] = fname.strip()
    return result


def resolve_filename(outdir: Path, bc: str, step: str,
                     prefix: str, name_map: dict, default: str) -> Path:
    """
    确定某个输出文件的最终路径，优先级：
      1. name_map 中指定的文件名（--out-step1a / --out-step1b）
      2. prefix + default 自动拼接
      3. 纯 default（无 prefix）
    """
    if bc in name_map:
        return outdir / name_map[bc]
    if prefix:
        return outdir / f"{prefix}_{default}"
    return outdir / default




def step1a_group_mean(df: pd.DataFrame, prefix_map: dict,
                      matrix_order: list) -> pd.DataFrame:
    """
    输入: 归一化相对丰度表（列 = 样品，行 = 分类群）
    输出: 各基质均值相对丰度（列 = 基质名称，行 = 分类群）
    """
    groups = {mat: [] for mat in matrix_order}
    unmatched = []

    for col in df.columns:
        mat = infer_matrix(col, prefix_map)
        if mat is None:
            unmatched.append(col)
        elif mat in groups:
            groups[mat].append(col)

    if unmatched:
        log(f"  以下样品列无法匹配任何前缀，已跳过: {unmatched}", "WARN")

    result = {}
    for mat in matrix_order:
        cols = groups[mat]
        if not cols:
            log(f"  基质 '{mat}' 未找到任何匹配样品列，填充 0。", "WARN")
            result[mat] = pd.Series(0.0, index=df.index)
        else:
            result[mat] = df[cols].mean(axis=1)
            log(f"  {mat}: {len(cols)} 个样品 → 均值完成")

    return pd.DataFrame(result)


# ═══════════════════════════════════════════════════════════════════════════════
# Step 1b：计算 eDNA index（组内归一化）
# ═══════════════════════════════════════════════════════════════════════════════

def step1b_edna_index(mean_df: pd.DataFrame,
                      exclude_matrices: set = None) -> pd.DataFrame:
    """
    输入: 均值相对丰度表（列 = 基质，行 = 分类群）
    输出: eDNA index 表（各值除以该分类群在有效基质中的最大值）

    exclude_matrices: 该条形码中不参与归一化计算的基质集合（如 COI 的 Biofilm）
    """
    exclude_matrices = exclude_matrices or set()

    # 仅用有效基质列计算最大值
    valid_cols = [c for c in mean_df.columns if c not in exclude_matrices]
    if not valid_cols:
        raise ValueError("所有基质列均被排除，无法计算 eDNA index。")

    row_max = mean_df[valid_cols].max(axis=1)

    # 最大值为 0 的分类群（全基质未检出）保持为 0，避免除以零
    index_df = mean_df.copy()
    nonzero_mask = row_max > 0
    index_df.loc[nonzero_mask] = (
        mean_df.loc[nonzero_mask].div(row_max[nonzero_mask], axis=0)
    )

    zero_taxa = mean_df.index[~nonzero_mask].tolist()
    if zero_taxa:
        log(f"  {len(zero_taxa)} 个分类群在所有有效基质中均为 0，保留为 0。", "WARN")

    # 排除的基质列强制归零（该条形码不适用该基质）
    for mat in exclude_matrices:
        if mat in index_df.columns:
            index_df[mat] = 0.0

    return index_df


# ═══════════════════════════════════════════════════════════════════════════════
# Step 2：跨条形码合并 → ensemble index
# ═══════════════════════════════════════════════════════════════════════════════

def step2_merge(barcode_dfs: dict, matrix_order: list,
                exclude_map: dict = None) -> pd.DataFrame:
    """
    输入:
        barcode_dfs : {barcode_name: eDNA_index_df}
        matrix_order: 基质顺序列表
        exclude_map : {barcode: {matrix, ...}} 该条形码不覆盖的基质
    输出:
        合并后的 ensemble index 前体表，附带 Status 列
    """
    exclude_map = exclude_map or {}
    all_taxa = sorted(set().union(*[set(df.index) for df in barcode_dfs.values()]))
    barcode_names = list(barcode_dfs.keys())

    records = []

    for taxon in all_taxa:
        row = {}
        detected_in = [bc for bc in barcode_names if taxon in barcode_dfs[bc].index]

        for mat in matrix_order:
            # 收集各条形码在该基质的非零值
            # 条件：条形码检测到该分类群 且 该基质不在该条形码的排除列表中
            valid_vals = []
            for bc in detected_in:
                if mat in (exclude_map.get(bc, set())):
                    continue  # 该条形码不覆盖该基质，跳过
                val = barcode_dfs[bc].loc[taxon, mat] if mat in barcode_dfs[bc].columns else 0.0
                if val > 0:
                    valid_vals.append(val)

            row[mat] = np.mean(valid_vals) if valid_vals else 0.0

        # Status 标注
        if len(detected_in) > 1:
            sources = "+".join(sorted(detected_in))
            row["Status"] = f"Merged_Average_(Sources:{sources})"
        elif len(detected_in) == 1:
            row["Status"] = f"Unique_to_{detected_in[0]}"
        else:
            row["Status"] = "Not_detected"

        row["Taxon"] = taxon
        records.append(row)

    result = pd.DataFrame(records).set_index("Taxon")
    result = result[["Status"] + matrix_order]
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# 主流程
# ═══════════════════════════════════════════════════════════════════════════════

def run_full_pipeline(args):
    """从原始归一化相对丰度表开始的完整流程。"""
    prefix_map = parse_matrix_map(args.matrix_map)
    matrix_order = list(dict.fromkeys(prefix_map.values()))
    exclude_map = parse_exclude(args.exclude_matrix)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    name_map_1a = parse_name_map(args.out_step1a)
    name_map_1b = parse_name_map(args.out_step1b)
    prefix = args.prefix or ""

    # 解析输入 "BARCODE:path" 格式
    inputs = {}
    for item in args.inputs:
        if ":" not in item:
            log(f"--inputs 格式错误，期望 BARCODE:path，得到: '{item}'", "ERR")
            sys.exit(1)
        bc, path = item.split(":", 1)
        inputs[bc.strip()] = path.strip()

    barcode_step1 = {}

    for bc, path in inputs.items():
        log(f"处理条形码: {bc}  ←  {path}", "STEP")

        # Step 1a
        log("  Step 1a: 读取归一化相对丰度表...", "INFO")
        raw_df = read_table(path)
        log(f"  读取完成: {len(raw_df)} 个分类群 × {len(raw_df.columns)} 个样品")

        mean_df = step1a_group_mean(raw_df, prefix_map, matrix_order)
        log(f"  分组均值完成: 基质列 = {list(mean_df.columns)}")

        out1a = resolve_filename(outdir, bc, "1a", prefix, name_map_1a,
                                 f"{bc}_step1a_matrix_mean.txt")
        mean_df.to_csv(out1a, sep="\t", float_format="%.8e")
        log(f"  Step 1a 输出 → {out1a}", "OK")

        # Step 1b
        log("  Step 1b: 计算 eDNA index...", "INFO")
        excl_mats = exclude_map.get(bc, set())
        if excl_mats:
            log(f"  该条形码排除基质（不参与归一化）: {excl_mats}", "WARN")

        index_df = step1b_edna_index(mean_df, exclude_matrices=excl_mats)

        out1b = resolve_filename(outdir, bc, "1b", prefix, name_map_1b,
                                 f"{bc}_step1b_eDNA_index.txt")
        index_df.to_csv(out1b, sep="\t", float_format="%.8e")
        log(f"  Step 1b 输出 → {out1b}", "OK")

        barcode_step1[bc] = index_df

    # Step 2（多条形码时合并）
    if len(barcode_step1) > 1:
        log("Step 2: 跨条形码合并 → ensemble index", "STEP")
        merged = step2_merge(barcode_step1, matrix_order, exclude_map)
        _save_and_report(merged, outdir, matrix_order,
                         output_merged=args.output_merged, prefix=prefix)
    else:
        log("仅有单条形码，跳过 Step 2 合并。", "WARN")


def run_from_step1(args):
    """从已有 Step1 文件直接进行 Step2 合并。"""
    exclude_map = parse_exclude(args.exclude_matrix)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    barcode_step1 = {}
    matrix_order = None

    for item in args.step1_inputs:
        if ":" not in item:
            log(f"--step1-inputs 格式错误，期望 BARCODE:path，得到: '{item}'", "ERR")
            sys.exit(1)
        bc, path = item.split(":", 1)
        bc, path = bc.strip(), path.strip()

        log(f"读取 Step1 文件: {bc}  ←  {path}", "STEP")
        df = read_table(path)
        barcode_step1[bc] = df.fillna(0).astype(float)

        # 以第一个文件的列顺序为准
        if matrix_order is None:
            matrix_order = list(df.columns)
        log(f"  {len(df)} 个分类群 | 基质列: {list(df.columns)}", "INFO")

    log("Step 2: 跨条形码合并 → ensemble index", "STEP")
    merged = step2_merge(barcode_step1, matrix_order, exclude_map)
    _save_and_report(merged, outdir, matrix_order,
                     output_merged=args.output_merged, prefix=args.prefix or "")


def _save_and_report(merged: pd.DataFrame, outdir: Path, matrix_order: list,
                     output_merged: str = None, prefix: str = ""):
    """保存合并结果并打印统计摘要。"""
    if output_merged:
        out2 = outdir / output_merged
    elif prefix:
        out2 = outdir / f"{prefix}_ensemble_index_merged.txt"
    else:
        out2 = outdir / "ensemble_index_merged.txt"

    merged.to_csv(out2, sep="\t", float_format="%.8e")
    log(f"ensemble index 输出 → {out2}", "OK")

    # 统计摘要
    log("─" * 50, "INFO")
    log(f"合并后总分类群数: {len(merged)}", "INFO")
    for status, grp in merged.groupby("Status"):
        log(f"  {status}: {len(grp)} 个", "INFO")

    # 共有分类群详情
    shared = merged[merged["Status"].str.startswith("Merged_Average")]
    if not shared.empty:
        log("─" * 50, "INFO")
        log("共有分类群合并详情:", "INFO")
        for taxon, row in shared.iterrows():
            vals = "  ".join(f"{m}={row[m]:.4e}" for m in matrix_order)
            log(f"  {taxon}: {vals}", "INFO")


# ═══════════════════════════════════════════════════════════════════════════════
# 命令行入口
# ═══════════════════════════════════════════════════════════════════════════════

def build_parser():
    parser = argparse.ArgumentParser(
        prog="edna_ensemble_pipeline.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            eDNA Ensemble Index Pipeline
            ────────────────────────────
            从归一化相对丰度表（或已有 Step1 文件）计算多条形码
            eDNA ensemble index，并输出合并结果。

            参考方法: Kelly et al. (2019), Djurhuus et al. (2020)
        """),
        epilog=textwrap.dedent("""\
            示例:
              # 完整流程（从归一化表开始）
              python edna_ensemble_pipeline.py \\
                  --inputs 18S:18s_normalize.txt COI:co1_normalize.txt \\
                  --matrix-map "B:Biofilm,M:Membrane,S:Sediment" \\
                  --exclude-matrix COI:Biofilm \\
                  --outdir results/

              # 添加统一前缀
              python edna_ensemble_pipeline.py \\
                  --inputs 18S:18s_normalize.txt COI:co1_normalize.txt \\
                  --matrix-map "B:Biofilm,M:Membrane,S:Sediment" \\
                  --exclude-matrix COI:Biofilm \\
                  --prefix benthic_2024 \\
                  --outdir results/

              # 完全自定义所有文件名
              python edna_ensemble_pipeline.py \\
                  --inputs 18S:18s_normalize.txt COI:co1_normalize.txt \\
                  --matrix-map "B:Biofilm,M:Membrane,S:Sediment" \\
                  --exclude-matrix COI:Biofilm \\
                  --out-step1a "18S:18s_mean.txt,COI:co1_mean.txt" \\
                  --out-step1b "18S:18s_index.txt,COI:co1_index.txt" \\
                  --output-merged final_ensemble.txt \\
                  --outdir results/

              # 仅 Step2 合并（从已有 Step1 文件开始）
              python edna_ensemble_pipeline.py \\
                  --step1-inputs 18S:18s_step1.txt COI:co1_step1.txt \\
                  --exclude-matrix COI:Biofilm \\
                  --output-merged final_ensemble.txt \\
                  --outdir results/
        """)
    )

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--inputs", nargs="+", metavar="BARCODE:PATH",
        help=(
            "完整流程入口：归一化相对丰度表，格式为 BARCODE:文件路径。\n"
            "可传入多个条形码，例如：\n"
            "  --inputs 18S:18s.txt COI:co1.txt"
        )
    )
    mode.add_argument(
        "--step1-inputs", nargs="+", metavar="BARCODE:PATH",
        help=(
            "Step2 入口：已计算好的 eDNA index Step1 文件，格式同 --inputs。\n"
            "列应为基质名称（Biofilm / Membrane / Sediment）。"
        )
    )

    parser.add_argument(
        "--matrix-map", metavar="STR",
        default="B:Biofilm,M:Membrane,S:Sediment",
        help=(
            "样品列前缀 → 基质名称映射（仅 --inputs 模式需要）。\n"
            "格式: 前缀1:基质名1,前缀2:基质名2,...\n"
            "默认: 'B:Biofilm,M:Membrane,S:Sediment'"
        )
    )
    parser.add_argument(
        "--exclude-matrix", nargs="*", metavar="BARCODE:Matrix", default=[],
        help=(
            "指定某条形码不覆盖的基质（不参与该基质的归一化和均值计算）。\n"
            "格式: BARCODE:基质名，可多次指定，例如：\n"
            "  --exclude-matrix COI:Biofilm"
        )
    )
    parser.add_argument(
        "--outdir", "-o", metavar="DIR", default="edna_output",
        help="输出目录（不存在时自动创建）。默认: edna_output/"
    )

    # ── 输出文件命名选项 ──────────────────────────────────────────────────────
    naming = parser.add_argument_group(
        "输出文件命名（三种方式可单独或组合使用，优先级：指定文件名 > prefix > 默认名）"
    )
    naming.add_argument(
        "--prefix", metavar="STR", default="",
        help=(
            "为所有输出文件添加统一前缀。\n"
            "例如 --prefix benthic_2024 会生成:\n"
            "  benthic_2024_18S_step1a_matrix_mean.txt\n"
            "  benthic_2024_ensemble_index_merged.txt 等"
        )
    )
    naming.add_argument(
        "--output-merged", metavar="FILENAME", default="",
        help=(
            "指定最终合并文件（Step 2）的完整文件名（不含目录）。\n"
            "优先级高于 --prefix。\n"
            "例如: --output-merged benthic_final_ensemble.txt"
        )
    )
    naming.add_argument(
        "--out-step1a", metavar="BARCODE:FILENAME[,...]", default="",
        help=(
            "指定各条形码 Step 1a 输出文件名，格式: BARCODE:文件名。\n"
            "多个条形码用逗号分隔，优先级高于 --prefix。\n"
            "例如: --out-step1a \"18S:18s_mean.txt,COI:co1_mean.txt\""
        )
    )
    naming.add_argument(
        "--out-step1b", metavar="BARCODE:FILENAME[,...]", default="",
        help=(
            "指定各条形码 Step 1b 输出文件名，格式同 --out-step1a。\n"
            "例如: --out-step1b \"18S:18s_index.txt,COI:co1_index.txt\""
        )
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    log("eDNA Ensemble Index Pipeline 启动", "INFO")
    log(f"输出目录: {args.outdir}", "INFO")
    log("─" * 50, "INFO")

    try:
        if args.inputs:
            run_full_pipeline(args)
        else:
            run_from_step1(args)
    except Exception as e:
        log(f"运行出错: {e}", "ERR")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    log("─" * 50, "INFO")
    log("全部完成。", "OK")


if __name__ == "__main__":
    main()