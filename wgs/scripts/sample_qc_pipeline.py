"""
本脚本用于根据样本QC标记文件筛选通过质量控制的样本，并使用PLINK2工具从原始PLINK文件中提取这些样本的子集。
参数说明：
    --sample_qc_flags: 样本QC标记CSV文件路径，包含样本的各种QC通过与否标记。
    --bed_prefix: 输入的PLINK文件前缀（不含.bed/.bim/.fam后缀）。
    --out_prefix: 输出文件的前缀。
    --mode: 筛选模式，选择使用 PASS_SMISS 还是 PASS_MEAN_DP 作为过滤条件，默认"smiss"。
    --include_pass_pi_hat: 是否包含 PASS_PI_HAT 过滤条件，默认为False。

作者: ZHAO TIE
"""

import os
import subprocess
import tempfile
import argparse
import pandas as pd

# 解析命令行参数
parser = argparse.ArgumentParser(description="Run PLINK2 extraction for QC-passed samples")
parser.add_argument("--sample_qc_flags", required=True, help="Path to sample QC flags CSV")
parser.add_argument("--bed_prefix", required=True, help="Prefix to input PLINK BED/BIM/FAM files")
parser.add_argument("--out_prefix", required=True, help="Prefix for output files")
parser.add_argument("--mode", choices=["smiss", "meandp"], default="smiss", help="QC mode to filter samples")
parser.add_argument("--include_pass_pi_hat", action="store_true", help="Whether to include PASS_PI_HAT filter")
args = parser.parse_args()

def filter_pass_samples(sample_qc_flags_path, mode="meandp", include_pass_pi_hat=True):
    """
    根据指定模式筛选通过QC的样本。

    参数:
        sample_qc_flags_path (str): 样本QC标记文件（CSV）的路径。
        mode (str): 模式选择，"smiss" 表示使用 PASS_SMISS，"meandp" 表示使用 PASS_MEAN_DP。
        include_pass_pi_hat (bool): 是否要求同时满足 PASS_PI_HAT 条件，默认为 True。

    返回:
        pd.DataFrame: 包含通过筛选条件的 #FID 和 IID 的 DataFrame。
    """
    # 读取QC标记数据
    df = pd.read_csv(sample_qc_flags_path)

    # 根据mode选择基础筛选条件
    if mode == "smiss":
        cond = (df["PASS_SMISS"] == True) & (df["PASS_HET_F"] == True)
    elif mode == "meandp":
        cond = (df["PASS_MEAN_DP"] == True) & (df["PASS_HET_F"] == True)
    else:
        raise ValueError("mode 必须为 'smiss' 或 'meandp'")

    # 若需要额外加上 PI_HAT 条件
    if include_pass_pi_hat:
        cond = cond & (df["PASS_PI_HAT"] == True)

    return df.loc[cond, ["#FID", "IID"]].reset_index(drop=True)

def run_plink2_extract_samples(sample_qc_flags_path, bed_prefix, out_prefix,
                               mode="smiss", include_pass_pi_hat=True, plink2_path="plink2"):
    """
    从 PLINK 文件中提取通过QC的样本子集。

    参数:
        sample_qc_flags_path (str): 样本QC标记CSV路径
        bed_prefix (str): 原始 PLINK 文件前缀（不含 .bed/.bim/.fam）
        out_prefix (str): 输出文件前缀
        mode (str): 筛选模式，"smiss" 或 "meandp"
        include_pass_pi_hat (bool): 是否筛选 PASS_PI_HAT 样本
        plink2_path (str): plink2 执行路径
    """
    # 获取符合条件的样本
    pass_df = filter_pass_samples(sample_qc_flags_path, mode, include_pass_pi_hat)

    # 临时写入样本列表
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp:
        pass_df.to_csv(tmp.name, sep="\t", header=False, index=False)
        keep_file = tmp.name

    # 构造 PLINK2 命令
    cmd = [
        plink2_path,
        "--bfile", bed_prefix,
        "--keep", keep_file,
        "--make-bed",
        "--out", out_prefix,
        "--threads", "16",  # 根据需要调整线程数
    ]

    # 执行PLINK2命令
    try:
        subprocess.run(cmd, check=True)
        print(f"提取完成，输出文件前缀为：{out_prefix}")
    finally:
        os.remove(keep_file)  # 清理临时文件

# 主程序入口，执行样本提取流程
run_plink2_extract_samples(
    sample_qc_flags_path=args.sample_qc_flags,
    bed_prefix=args.bed_prefix,
    out_prefix=args.out_prefix,
    mode=args.mode,
    include_pass_pi_hat=args.include_pass_pi_hat,
    plink2_path='/home/b/b37974/plink2'
)
