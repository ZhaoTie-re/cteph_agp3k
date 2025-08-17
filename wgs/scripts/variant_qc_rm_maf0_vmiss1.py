"""
脚本功能：运行 Plink2 对输入 Plink 数据进行变异层面 QC 分析，并排除 MAF=0 或 VMISS=1 的变异。
步骤包括：
  1. 运行 Plink2 计算变异的频率、缺失率、HWE等
  2. 筛选出低质量变异（MAF=0 或 VMISS=1）
  3. 使用 Plink2 排除这些变异，输出新的 bed 文件

作者: ZHAO TIE
"""

import argparse
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="执行变异层面QC并删除MAF=0或缺失率=1的变异")
parser.add_argument('--script_path', type=str, required=True, help="variant_qc_calculator.py 所在路径")
parser.add_argument('--bed_prefix', type=str, required=True, help="输入 Plink 数据前缀")
parser.add_argument('--output_prefix', type=str, required=True, help="输出文件前缀")
parser.add_argument('--threads', type=int, default=32, help="并行线程数")
args = parser.parse_args()

script_path = os.path.abspath(args.script_path)
if script_path not in sys.path:
    sys.path.append(script_path)

bed_prefix = args.bed_prefix
output_prefix = args.output_prefix
threads = args.threads


import importlib
import variant_qc_calculator

# 强制重新加载模块（适用于开发调试阶段）
importlib.reload(variant_qc_calculator)

from variant_qc_calculator import (
    init_globals_for_chunk,
    process_chunk,
    run_plink2_variant_qc,
    extract_maf0_or_vmiss1_variants_streaming,
    run_plink2_exclude_variants,
)

if __name__ == "__main__":
    variant_qc_summary = run_plink2_variant_qc(
        bed_prefix=bed_prefix,
        output_prefix=output_prefix,
        threads=threads,
    )

    maf0_or_vmiss1_flags_tsv, maf0_or_vmiss1_ids_tsv = extract_maf0_or_vmiss1_variants_streaming(
        input_file=variant_qc_summary, chunksize=1000000
    )

    run_plink2_exclude_variants(
        bed_prefix=bed_prefix,
        output_prefix=f"{output_prefix}.rm_maf0_vmiss1",
        exclude_variants_file=maf0_or_vmiss1_ids_tsv,
        threads=threads
    )
