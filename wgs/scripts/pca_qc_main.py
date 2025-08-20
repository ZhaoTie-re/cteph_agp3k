#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PCA 质控主脚本
---------------------------------
功能：
1) 基于提供的 PLINK 二进制基因型（--bfile 前缀）进行变体筛选与 LD 修剪，并计算 PCA。
2) 使用计算得到的 eigenvec/eigenval 生成成对（PC1 vs PC2, PC3 vs PC4, ...）散点图 PDF，按病例/对照分组上色。

使用示例：
python pca_qc_main.py \
  --bed_prefix /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/13.run_variant_qc/cteph_agp3k.sqc.vqc \
  --high_ld /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/high-LD-regions-hg38-GRCh38.txt \
  --output_prefix cteph_agp3k \
  --threads 32 \
  --maf_threshold 0.05 \
  --case_prefix PHOM \
  --case_name CTEPH \
  --control_name AGP3K \
  --plink2_path /home/b/b37974/plink2

注意：
- 该脚本依赖同目录下的 `pca_qc_tools.py`，其中应当实现：
  - run_prune_and_pca(bed_prefix, maf_threshold, high_ld, threads, output_prefix, plink2_path) -> 返回
    (no_high_ld_prefix, prune_in, prune_out, eigenvec_file, eigenval_file, eigenvec_allele_file)
  - plot_pca_pairwise_pdf(eigenvec_file, eigenval_file, case_prefix, case_name, control_name, output_pdf) -> 返回 pdf 路径
- 在无显示环境的服务器上作图，使用 matplotlib 的 'Agg' 后端。

作者: ZHAO TIE
"""

from __future__ import annotations

import argparse
import os
import sys
import logging

# 必须在导入 pyplot 之前设置后端
import matplotlib
matplotlib.use("Agg")

# 将脚本目录加入 sys.path 以便导入本地工具模块
SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.append(SCRIPT_DIR)

# 本地工具
import importlib
import pca_qc_tools  # type: ignore
# 开发调试时可保留一次 reload，生产环境中不会有副作用
importlib.reload(pca_qc_tools)
from pca_qc_tools import run_prune_and_pca, plot_pca_pairwise_pdf  # type: ignore


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="运行 PLINK2 进行 LD 修剪与 PCA，并输出成对 PC 散点图 PDF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bed_prefix", required=True,
                        help="输入 PLINK 二进制文件前缀（不带扩展名）")
    parser.add_argument("--high_ld", required=True,
                        help="高 LD 区域 BED/TXT（兼容 plink --exclude range）")
    parser.add_argument("--output_prefix", required=True,
                        help="输出文件名前缀")
    parser.add_argument("--threads", type=int, default=32,
                        help="线程数（传递给 PLINK2）")
    parser.add_argument("--maf_threshold", type=float, default=0.05,
                        help="MAF 过滤阈值")
    parser.add_argument("--case_prefix", required=True,
                        help="病例样本 IID 的前缀，用于自动分组标记")
    parser.add_argument("--case_name", default="CASE",
                        help="病例组显示名称")
    parser.add_argument("--control_name", default="CTRL",
                        help="对照组显示名称")
    parser.add_argument("--plink2_path", default="/home/b/b37974/plink2",
                        help="plink2 可执行文件路径")
    parser.add_argument("--output_pdf", default=None,
                        help="可选：输出 PDF 路径；默认使用 <output_prefix>.pca_pairwise.pdf")
    parser.add_argument("--verbose", action="store_true", help="输出调试信息")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='[%(levelname)s] %(message)s')

    bed_prefix: str = args.bed_prefix
    high_ld: str = args.high_ld
    output_prefix: str = args.output_prefix
    threads: int = args.threads
    maf_threshold: float = args.maf_threshold
    case_prefix: str = args.case_prefix
    case_name: str = args.case_name
    control_name: str = args.control_name
    plink2_path: str = args.plink2_path
    output_pdf: str = args.output_pdf or f"{output_prefix}.pca_pairwise.pdf"

    # Step 1: 运行 LD 修剪与 PCA
    try:
        (
            no_high_ld_prefix,
            prune_in,
            prune_out,
            eigenvec_file,
            eigenval_file,
            _eigenvec_allele_file,
        ) = run_prune_and_pca(
            bed_prefix=bed_prefix,
            maf_threshold=maf_threshold,
            high_ld=high_ld,
            threads=threads,
            output_prefix=output_prefix,
            plink2_path=plink2_path,
        )
    except Exception:
        logging.exception("运行 LD 修剪与 PCA 失败")
        sys.exit(1)

    # Existence checks for eigenvec/eigenval
    if not os.path.exists(str(eigenvec_file)) or not os.path.exists(str(eigenval_file)):
        logging.error("找不到 eigenvec/eigenval 文件，请检查上游步骤或参数")
        sys.exit(1)

    # Step 2: 生成成对 PC 的散点图 PDF
    try:
        pdf_path: str = plot_pca_pairwise_pdf(
            eigenvec_file=str(eigenvec_file),
            eigenval_file=str(eigenval_file),
            case_prefix=case_prefix,
            case_name=case_name,
            control_name=control_name,
            output_pdf=output_pdf,
        )
    except Exception:
        logging.exception("绘制 PCA 散点图失败")
        sys.exit(1)

    logging.info("[PCA QC 完成]")
    logging.info(f"输入: --bed_prefix {bed_prefix}")
    logging.info(f"高 LD 区域: {high_ld}")
    logging.info(f"MAF 阈值: {maf_threshold}")
    logging.info(f"线程数: {threads}")
    logging.info(f"中间输出前缀(去高 LD): {no_high_ld_prefix}")
    logging.info(f"eigenvec: {eigenvec_file}")
    logging.info(f"eigenval: {eigenval_file}")
    logging.info(f"PDF 输出: {pdf_path}")


if __name__ == "__main__":
    main()