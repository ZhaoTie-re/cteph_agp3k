# %%
"""
sample_qc_main.py

本脚本用于运行样本质量控制（QC）流程。
通过输入相关文件和参数，计算样本的多项QC指标，生成统计结果和可视化图表，
并输出样本QC标记以供后续分析使用。

参数说明：
--info_file           : 输入的样本信息xlsx文件路径
--bed_prefix          : 输入的BED文件前缀
--high_ld             : 高连锁不平衡区域文件路径
--case_prefix         : 病例样本前缀
--out_prefix          : 输出文件前缀
--script_path         : 脚本目录路径
--pi_threshold        : PI_HAT阈值，默认0.2
--dp_robust_z_threshold: DP的稳健Z值阈值，默认-3.0
--het_threshold       : HET F阈值，默认"5sd"
--smiss_threshold     : 样本缺失率阈值，默认0.10
--threads             : 使用线程数，默认8

作者: ZHAO TIE
"""

import argparse
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt

# 解析命令行参数
parser = argparse.ArgumentParser(description="Run sample QC pipeline")
parser.add_argument("--info_file", required=True, help="Path to info xlsx file")
parser.add_argument("--bed_prefix", required=True, help="Input BED file prefix")
parser.add_argument("--high_ld", required=True, help="High LD regions file")
parser.add_argument("--case_prefix", required=True, help="Prefix for case samples")
parser.add_argument("--out_prefix", required=True, help="Prefix for output files")
parser.add_argument("--script_path", required=True, help="Path to script directory")
parser.add_argument("--pi_threshold", type=float, default=0.2, help="PI_HAT threshold")
parser.add_argument("--dp_robust_z_threshold", type=float, default=-3.0, help="DP robust Z threshold")
parser.add_argument("--het_threshold", type=str, default="5sd", help="HET F threshold")
parser.add_argument("--smiss_threshold", type=float, default=0.10, help="Sample missing rate threshold")
parser.add_argument("--threads", type=int, default=8, help="Number of threads to use (e.g., 4 or 8)")

args = parser.parse_args()

# 设置matplotlib默认样式，避免与其他样式冲突
plt.style.use('default')

# 添加脚本路径，确保能够找到自定义模块
script_path = os.path.abspath(args.script_path)
if script_path not in sys.path:
    sys.path.append(script_path)

# %%
# 提取参数变量，方便后续调用
threads = args.threads                              # 使用线程数
info_file = args.info_file                          # 样本信息文件路径
bed_prefix = args.bed_prefix                        # 输入BED文件前缀
high_ld = args.high_ld                              # 高LD区域文件
case_prefix = args.case_prefix                      # 病例样本前缀
out_prefix = args.out_prefix                        # 输出文件前缀
pi_threshold = args.pi_threshold                    # PI_HAT阈值
dp_robust_z_threshold = args.dp_robust_z_threshold  # DP稳健Z阈值
het_threshold = args.het_threshold                  # HET F阈值
smiss_threshold = args.smiss_threshold              # 样本缺失率阈值

# %%
if __name__ == "__main__":
    # 导入自定义模块中的函数
    from sample_qc_calculator import (
        update_fam_and_pheno,
        run_sample_qc_summary,
        compute_robust_z_by_group
    )
    from sample_qc_flags import (
        plot_pi_hat_distribution,
        plot_kinship_network_and_prune,
        plot_hetf_vs_meandp_table,
        plot_hetf_vs_smiss_table,
        generate_sample_qc_flags
    )

    # =====================
    # 数据准备阶段
    # =====================
    # 更新FAM和PHENO文件，准备QC所需数据
    update_fam_and_pheno(
        info_file=info_file,
        bed_prefix=bed_prefix,
        case_prefix=case_prefix,
        out_prefix=out_prefix,
    )

    # 运行样本QC汇总计算，包括高LD区域处理等
    run_sample_qc_summary(
        info_file=info_file,
        bed_prefix=out_prefix,
        high_ld_file=high_ld,
        threads=threads
    )

    # 读取样本QC汇总结果
    sample_qc_summary = pd.read_csv(f"{out_prefix}.sample_qc_summary.csv")

    # =====================
    # QC计算阶段
    # =====================
    # 计算MEAN_DP的稳健Z值，左尾检验
    sample_qc_summary_update = compute_robust_z_by_group(
        df=sample_qc_summary,
        group_col="TARGET_DP",
        value_col="MEAN_DP",
        z_col_name="ROBUST_Z_DP",
        p_col_name="P_ROBUST_Z_DP",
        fdr_col_name="FDR_ROBUST_Z_DP",
        tail="left"
    )

    # 保存更新后的样本QC汇总结果
    sample_qc_summary_update.to_csv(f"{out_prefix}.sample_qc_summary.csv", index=False)

    # 读取PI_HAT文件，仅保留超过阈值的样本对
    pihat_df = pd.read_csv(
        f"{out_prefix}.pi_hat.csv",
        usecols=["FID1", "IID1", "FID2", "IID2", "PI_HAT"],
        dtype={
            "FID1": "category",
            "IID1": "category",
            "FID2": "category",
            "IID2": "category",
            "PI_HAT": "float32"
        }
    )
    pihat_df = pihat_df[pihat_df["PI_HAT"] > pi_threshold].reset_index(drop=True)

    # =====================
    # 图表生成阶段
    # =====================
    # 绘制PI_HAT分布图
    plot_pi_hat_distribution(
        pihat_df=pihat_df,
        pi_threshold=pi_threshold,
        case_prefix=case_prefix,
        out_prefix=out_prefix
    )

    # 绘制亲缘关系网络图并进行样本剔除，返回剔除样本列表
    pi_outlier = plot_kinship_network_and_prune(
        pihat_df=pihat_df,
        sample_qc_summary=sample_qc_summary_update,
        pi_threshold=pi_threshold,
        case_prefix=case_prefix,
        out_prefix=out_prefix
    )

    # 绘制HET F与平均DP的关系表格
    plot_hetf_vs_meandp_table(
        sample_qc_summary=sample_qc_summary_update,
        dp_robust_z_threshold=dp_robust_z_threshold,
        het_threshold=het_threshold,
        case_prefix=case_prefix,
        out_prefix=out_prefix
    )

    # 绘制HET F与样本缺失率的关系表格
    plot_hetf_vs_smiss_table(
        sample_qc_summary=sample_qc_summary_update,
        smiss_threshold=smiss_threshold,
        het_threshold=het_threshold,
        case_prefix=case_prefix,
        out_prefix=out_prefix
    )

    # =====================
    # 标记输出阶段
    # =====================
    # 生成样本QC标记文件，包含各类QC异常标记
    generate_sample_qc_flags(
        sample_qc_summary=sample_qc_summary_update,
        smiss_threshold=smiss_threshold,
        dp_robust_z_threshold=dp_robust_z_threshold,
        het_threshold=het_threshold,
        pi_outlier=pi_outlier,
        out_prefix=out_prefix
    )
