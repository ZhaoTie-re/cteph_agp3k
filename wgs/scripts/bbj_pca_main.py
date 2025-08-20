"""
模块名称：BBJ PCA 主脚本（bbj_pca_main.py）

概述：
    本脚本用于调用 `prepare_bbj_pca_inputs` 和 `plot_pca_pairwise_pdf` 两个主要函数，自动完成BBJ数据的PCA输入文件准备及主成分分析结果的成对散点图PDF绘制。

适用场景：
    适用于需要对BBJ基因型数据进行主成分分析（PCA），并可视化主要成分成对关系的科研与分析流程。

输入参数说明（对应 argparse 参数）：
    --bbj_bed_prefix   ：BBJ BED文件的前缀（必需）
    --high_ld          ：高LD区域文件路径（必需）
    --output_prefix    ：输出文件前缀（可选，默认"bbj.maf"）
    --threads          ：线程数（可选，默认16）

输出结果：
    生成PCA主成分成对散点图的PDF文件（如 "bbj_pca_pairwise_plots.pdf"）。
"""

# 导入所需模块
import argparse
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt

# 动态添加脚本路径，便于加载工具模块
script_path = os.path.abspath('/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/scripts')
if script_path not in sys.path:
    sys.path.append(script_path)

import importlib
import bbj_projection_tools

#
# 强制重新加载模块（适用于开发调试阶段，确保修改能即时生效）
importlib.reload(bbj_projection_tools)

from bbj_projection_tools import (
    prepare_bbj_pca_inputs,
    plot_pca_pairwise_pdf,
)


def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="Run BBJ PCA workflow and plot results.")
    parser.add_argument("--bbj_bed_prefix", type=str, required=True, help="Prefix of BBJ BED files")  # BBJ BED文件前缀
    parser.add_argument("--high_ld", type=str, required=True, help="Path to high LD regions file")    # 高LD区域文件路径
    parser.add_argument("--output_prefix", type=str, default="bbj.maf", help="Output prefix (default: bbj.maf)")  # 输出文件前缀
    parser.add_argument("--threads", type=int, default=16, help="Number of threads (default: 16)")    # 线程数
    args = parser.parse_args()

    # 调用PCA输入准备函数
    no_ld_prefix, prune_in, prune_out, eigenvec_path, eigenval_path = prepare_bbj_pca_inputs(
        bbj_bed_prefix=args.bbj_bed_prefix,
        output_prefix=args.output_prefix,
        high_ld=args.high_ld,
        threads=args.threads,
    )

    # 绘制PCA成对主成分散点图并保存为PDF
    pdf_path = plot_pca_pairwise_pdf(eigenvec_path, eigenval_path, "bbj_pca_pairwise_plots.pdf")
    print("PDF saved to:", pdf_path)


if __name__ == "__main__":
    main()
