#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
平台敏感性位点过滤的按染色体（per-chromosome）CLI。

流程步骤：
  1) build_aligned_geno_mats（PLINK vs VCF）→ pre_path, post_path
  2) harmonize_gt_matrices → out_pre, out_post
  3) reorder_pre_to_post → out_pre_order, out_post_order
  4) build_sample_group_info(info_xls, sample_list)
  5) compute_coverage_transition_counts（按染色体命名输出）
  6) summarize_coverage_transition_significance（按染色体命名输出）

日志：终端打印简要进度；详细日志由 geno_miss_bias_tools.* 内部函数各自写入。
作者：ZHAO TIE
"""

import os
import sys
import argparse
import importlib
import time
import re

# 确保工具模块可被导入
SCRIPT_DIR = os.path.abspath('/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/scripts')
if SCRIPT_DIR not in sys.path:
    sys.path.append(SCRIPT_DIR)

import geno_miss_bias_tools  # noqa: E402
importlib.reload(geno_miss_bias_tools)  # dev convenience

from geno_miss_bias_tools import (  # noqa: E402
    build_aligned_geno_mats,
    harmonize_gt_matrices,
    reorder_pre_to_post,
    build_sample_group_info,
    compute_coverage_transition_counts,
    summarize_coverage_transition_significance,
)


def _ts():
    return time.strftime('%Y-%m-%d %H:%M:%S')


def info(msg: str):
    print(f"[{_ts()}][INFO] {msg}")


def warn(msg: str):
    print(f"[{_ts()}][WARN] {msg}")


def err(msg: str):
    print(f"[{_ts()}][ERROR] {msg}", file=sys.stderr)


def _normalize_chr(chr_arg: str) -> str:
    """接受 '22' 或 'chr22'，返回纯数字 '22'。格式不符则抛出 SystemExit(2)。"""
    if chr_arg is None:
        err('必须提供 --chr，形如 22 或 chr22')
        sys.exit(2)
    m = re.fullmatch(r'(?:chr)?(\d+)', chr_arg.strip(), flags=re.IGNORECASE)
    if not m:
        err('参数 --chr 只能为 22 或 chr22 这样的格式（仅数字染色体）。')
        sys.exit(2)
    return m.group(1)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description='平台敏感性位点过滤（按染色体）。',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # 必需的输入输出
    p.add_argument('--plink-prefix', required=True,
                   help='PLINK bed/bim/fam 前缀（QC 后）。')
    p.add_argument('--vcf', required=True,
                   help='VCF.GZ 路径（QC 前）。需要 .tbi 索引。')
    p.add_argument('--info-xls', required=True,
                   help='样本信息 Excel（至少包含 ID 与 Target DP (JHRPv4)）。')

    # 染色体及计算相关
    p.add_argument('--chr', required=True,
                   help='染色体标签，仅接受形如 22 或 chr22。输出文件名中将使用该染色体号。')
    p.add_argument('--n-chunks', type=int, default=12,
                   help='将变体切分为 N 块以并行导出/查询。')
    p.add_argument('--chunk-size', type=int, default=None,
                   help='替代 --n-chunks：按固定变体数切块（可选）。')
    p.add_argument('--max-parallel', type=int, default=8,
                   help='同时处理的最大并行块数。')
    p.add_argument('--plink-threads', type=int, default=8,
                   help='plink2 导出线程数。')
    p.add_argument('--bcftools-threads', type=int, default=8,
                   help='bcftools view/query 线程数。')
    p.add_argument('--keep-temp', action='store_true',
                   help='保留中间分块文件（默认不保留）。')

    # 统计参数
    p.add_argument('--coverage-labels', default='15x,30x',
                   help='覆盖度分组标签，格式 "A,B"；需与 df_info 归一化结果一致。')
    p.add_argument('--n-resamples', type=int, default=9999,
                   help='2×3 Fisher 蒙特卡罗重采样次数。')
    p.add_argument('--rng', type=int, default=42,
                   help='Monte Carlo 随机种子。')

    # 其他
    p.add_argument('--work-dir', default=os.getcwd(),
                   help='工作目录（默认：当前目录）。')

    return p.parse_args()


def main():
    args = parse_args()
    chrom_norm = _normalize_chr(args.chr)
    info(f"已规范化染色体标签为: {chrom_norm}")

    os.makedirs(args.work_dir, exist_ok=True)
    os.chdir(args.work_dir)

    cov_labels = tuple([s.strip() for s in args.coverage_labels.split(',')])
    if len(cov_labels) != 2:
        err('参数 --coverage-labels 必须包含两个以逗号分隔的标签，例如 15x,30x')
        sys.exit(2)

    info('步骤1：构建对齐的基因型矩阵（build_aligned_geno_mats）...')
    pre_path, post_path = build_aligned_geno_mats(
        plink_prefix=args.plink_prefix,
        vcf_path=args.vcf,
        chrom=chrom_norm,
        n_chunks=args.n_chunks,
        chunk_size=args.chunk_size,
        max_parallel=args.max_parallel,
        plink_threads=args.plink_threads,
        bcftools_threads=args.bcftools_threads,
        return_paths_only=True,
        keep_temp=args.keep_temp,
    )
    info(f"pre 路径={pre_path}")
    info(f"post 路径={post_path}")

    info('步骤2：统一矩阵格式（harmonize_gt_matrices）...')
    out_pre, out_post = harmonize_gt_matrices(pre_mt=pre_path, post_mt=post_path)
    info(f"统一后 pre={out_pre}")
    info(f"统一后 post={out_post}")

    info('步骤3：按照 post 顺序重排 pre（reorder_pre_to_post）...')
    out_pre_order, out_post_order = reorder_pre_to_post(pre_mt=out_pre, post_mt=out_post)
    info(f"重排后 pre={out_pre_order}")
    info(f"重排后 post={out_post_order}")

    # 从重排后的 pre 矩阵表头提取样本列表
    info('步骤4：构建样本分组信息（build_sample_group_info）...')
    with open(out_pre_order, 'r') as f:
        header = f.readline().rstrip('\n').split('\t')
    sample_list = sorted(list(set(header[1:])))
    df_info = build_sample_group_info(info_xls=args.info_xls, sample_list=sample_list)
    info(f"df_info 行数={len(df_info)}")

    info('步骤5：计算覆盖度分组的转移计数（compute_coverage_transition_counts）...')
    out_tsv = compute_coverage_transition_counts(
        out_pre_order=out_pre_order,
        out_post_order=out_post_order,
        df_info=df_info,
        coverage_labels=cov_labels,
        chrom=chrom_norm,
    )
    info(f"覆盖度转移统计={out_tsv}")

    info('步骤6：三层显著性汇总（summarize_coverage_transition_significance）...')
    summary_path = summarize_coverage_transition_significance(
        transitions_tsv=out_tsv,
        coverage_labels=cov_labels,
        n_resamples=args.n_resamples,
        rng=args.rng,
        chrom=chrom_norm,
    )
    info(f"统计汇总={summary_path}")

    info('完成。')


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        err('用户中断')
        sys.exit(130)
    except SystemExit:
        raise
    except Exception as e:
        err(f"致命错误: {e}")
        sys.exit(1)
