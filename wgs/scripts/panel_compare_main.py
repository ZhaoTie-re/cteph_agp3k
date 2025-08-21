#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本名称：panel_compare_main.py
=================================

【概述】
本脚本串联完成三步流程：
1) 运行 PLINK2 变异质控汇总（`run_plink2_variant_qc`），生成 `*.variant_qc_summary.tsv`；
2) 基于 ToMMo 面板（VCF）为上述结果追加列：`IN_TOMMO, TOMMO_AAF, TOMMO_FILTER`（`run_plink2_variant_qc_with_tommo`）；
3) 生成面板比较的多页 PDF 报告（`plot_tommo_panel_compare_pdf`）。

【适用场景】
- WGS/WES/芯片的质控后，需与 ToMMo（或其他外部 VCF 面板）对齐，快速产出统计表与可视化报告。

【输入与输出】
- 输入：
  - `--bed_prefix`：PLINK 二进制基因型前缀（.bed/.bim/.fam 同前缀）。
  - `--tommo_vcf_path`：ToMMo 面板的 bgzip 压缩 VCF 路径（需有 .tbi 索引）。
- 输出：
  - `{output_prefix}.variant_qc_summary.tsv`
  - `{output_prefix}.variant_qc_summary.variant_qc_with_tommo.tsv`
  - `{output_prefix}.variant_qc_with_tommo.panel_compare.pdf`

【性能参数】
- `--threads`：传给 `run_plink2_variant_qc` 与 `bcftools view --threads` 的线程数（建议与分配 CPU 一致）。
- `--chunk_size`：流式合并 ToMMo 注释时的分块大小（行数；内存紧张可调小）。
- `--max_workers`：按染色体并发的上限（I/O 友好建议 2~4；SSD 可适度提高）。

【作者】
- 作者（Author）：ZHAO TIE
- 维护（Maintainer）：ZHAO TIE

【使用示例】
- 典型：
  ```bash
  python panel_compare_main.py \
    --bed_prefix \
      /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/16.rm_maf0_vmiss1_repeat/cteph_agp3k.sqc.vqc.bbj_sample_keep.rm_maf0_vmiss1 \
    --output_prefix cteph_agp3k \
    --threads 32 \
    --tommo_vcf_path \
      /LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz \
    --chunk_size 500000 \
    --max_workers 4
  ```

【注意事项】
- ToMMo VCF 的染色体命名需与 `VARIANT_ID` 一致（例如 `chr1` vs `1`）。
- 若 ToMMo 的 AF 字段并非 `INFO/AF`，请在 `panel_compare_tools.py` 中调整格式串（已留有注释）。
- 无显示环境建议使用非交互后端（脚本已设置）。
"""

import argparse
import os
import sys
import matplotlib
matplotlib.use('Agg')  # 无显示环境下渲染 PDF/PNG
import matplotlib.pyplot as plt  # noqa: F401  (部分绘图函数内部需要)

# 将当前脚本目录加入 sys.path，确保本地模块可导入
HERE = os.path.abspath(os.path.dirname(__file__))
if HERE not in sys.path:
    sys.path.append(HERE)

# 导入你自己的工具模块
import importlib
import variant_qc_calculator
import panel_compare_tools

# 可在开发阶段强制 reload（正式环境可关闭）
# importlib.reload(variant_qc_calculator)
# importlib.reload(panel_compare_tools)

from variant_qc_calculator import run_plink2_variant_qc
from panel_compare_tools import (
    run_plink2_variant_qc_with_tommo,
    plot_tommo_panel_compare_pdf,
)


def parse_args() -> argparse.Namespace:
    """解析命令行参数。"""
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "运行 PLINK2 变异质控 → 追加 ToMMo 注释 → 生成 PDF 报告 的一体化脚本。\n"
            "作者: ZHAO TIE"
        ),
    )
    p.add_argument(
        "--bed_prefix", required=True,
        help="PLINK 二进制基因型前缀（.bed/.bim/.fam 同前缀）",
    )
    p.add_argument(
        "--output_prefix", required=True,
        help="输出文件前缀（会用于三步产物的文件名）",
    )
    p.add_argument(
        "--threads", type=int, default=16,
        help="线程数：用于 PLINK2 与 bcftools view --threads",
    )
    p.add_argument(
        "--tommo_vcf_path", required=True,
        help="ToMMo 面板 VCF（.vcf.gz），需配套 .tbi 索引",
    )
    p.add_argument(
        "--chunk_size", type=int, default=500_000,
        help="合并 ToMMo 注释时分块大小（行数）",
    )
    p.add_argument(
        "--max_workers", type=int, default=4,
        help="按染色体并行的最大并发数",
    )
    p.add_argument(
        "--skip_pdf", action="store_true",
        help="仅生成带有 ToMMo 注释的 TSV，不绘制 PDF",
    )
    p.add_argument(
        "--regions_chunk_lines", type=int, default=50_000,
        help="每个染色体再切分的chunk行数，用于进度条控制",
    )
    return p.parse_args()


def _progress(msg: str):
    print(f"[panel_compare_main] {msg}", file=sys.stderr, flush=True)


def main():
    args = parse_args()

    bed_prefix     = args.bed_prefix
    output_prefix  = args.output_prefix
    threads        = int(args.threads)
    tommo_vcf_path = args.tommo_vcf_path
    chunk_size     = int(args.chunk_size)
    max_workers    = int(args.max_workers)
    regions_chunk_lines = int(args.regions_chunk_lines)

    _progress("Step1: 运行 PLINK2 变异质控汇总 …")
    variant_qc_summary = run_plink2_variant_qc(
        bed_prefix=bed_prefix,
        output_prefix=output_prefix,
        threads=threads,
    )
    _progress(f"产出: {variant_qc_summary}")

    _progress("Step2: ToMMo 批量注释并合并三列（IN_TOMMO/TOMMO_AAF/TOMMO_FILTER） …")
    out_tsv = run_plink2_variant_qc_with_tommo(
        variant_qc_summary=variant_qc_summary,
        tommo_vcf_path=tommo_vcf_path,
        threads=threads,         # 传给 bcftools view --threads
        chunk_size=chunk_size,   # 流式合并，内存友好
        max_workers=max_workers, # 染色体并行数（I/O 友好 2~4）
        regions_chunk_lines=regions_chunk_lines,
    )
    _progress(f"产出: {out_tsv}")

    if args.skip_pdf:
        _progress("Step3: 跳过 PDF 生成（--skip_pdf 打开）")
        print(out_tsv)
        return

    _progress("Step3: 生成面板比较 PDF 报告 …")
    pdf_path = plot_tommo_panel_compare_pdf(
        variant_qc_with_tommo=out_tsv,
    )
    _progress(f"产出: {pdf_path}")

    # 终端输出，便于流水线收集
    print(out_tsv)
    print("Saved:", pdf_path)


if __name__ == "__main__":
    main()