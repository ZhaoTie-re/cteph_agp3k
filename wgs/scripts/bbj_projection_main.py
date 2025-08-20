# -*- coding: utf-8 -*-
"""
模块名称：BBJ 投影主脚本（bbj_projection_main.py）

【概述】
本脚本串联两个工具函数：
1) run_bbj_projection：对“我的数据”按 BBJ 的 prune 集合对齐，与 BBJ 合并后，依据 BBJ PCA 权重进行投影，产出 *.sscore。
2) plot_projection_pairwise_pdf：基于 *.sscore 中的 PC*_AVG 列，按分组（BBJ / 对照 / 病例）绘制 PC 成对散点图 PDF。

【适用场景】
- 已完成 QC 的 cohort，需将其投影到 BBJ 的 PCA 空间进行群体结构比对与可视化。

【参数（args）】
- --my_bed_prefix：我的 PLINK 前缀（QC 后）。
- --bbj_bed_prefix：BBJ 的 PLINK 前缀（已完成 bbj_prepare 与 bbj_pca）。
- --bbj_prune_in：BBJ 的 prune.in 文件路径。
- --my_prefix_out：我的数据在本流程中的输出基名（影响输出文件名前缀）。
- --bbj_prefix_out：BBJ 数据在本流程中的输出基名。
- --bbj_pca_acount：BBJ PCA acount 文件路径（plink2 --read-freq 使用）。
- --bbj_pca_eigenvec_allele：BBJ PCA eigenvec.allele 文件路径（plink2 --score 使用）。
- --threads：plink2 的线程数。
- --case_prefix：病例样本 IID 前缀（用于识别病例）。
- --case_name：病例在图例中的名称。
- --control_name：对照在图例中的名称。

【输出】
- 投影得分：{my_prefix_out}.{bbj_prefix_out}.projection.sscore
- 投影散点图 PDF：默认文件名 `bbj_projection_pairwise_plots.pdf`（可在函数中传入 output_pdf 自定义）

【注意】
- 本脚本不会修改内部工具的命令格式（cmd），仅对流程进行封装与调度。

【作者】
- ZHAO TIE
"""

import argparse
import sys
import os
from pathlib import Path

# 可视化相关导入（仅用于类型与可执行调用）
import matplotlib.pyplot as plt  # noqa: F401  # 防止某些环境下的懒加载/后端警告

# 将工具模块路径加入 sys.path（根据你的目录结构调整）
SCRIPT_DIR = Path(__file__).resolve().parent
TOOLS_DIR = SCRIPT_DIR  # 工具与主脚本在同一目录；若变动，请修改此行
if str(TOOLS_DIR) not in sys.path:
    sys.path.append(str(TOOLS_DIR))

# 导入自定义工具
import importlib
import bbj_projection_tools
# 开发调试时重新加载，确保修改即时生效；生产可移除
importlib.reload(bbj_projection_tools)

from bbj_projection_tools import run_bbj_projection, plot_projection_pairwise_pdf


def parse_args() -> argparse.Namespace:
    """参数解析：提供默认值，并允许命令行覆盖。"""
    parser = argparse.ArgumentParser(
        description="BBJ 投影与可视化：将 cohort 投影至 BBJ PCA 空间并绘制 PC 成对散点图"
    )
    # —— 投影阶段 ——
    parser.add_argument(
        "--my_bed_prefix", type=str, required=False,
        default="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/14.run_pca/cteph_agp3k.no_high_ld",
        help="我的 PLINK 前缀（QC 后）"
    )
    parser.add_argument(
        "--bbj_bed_prefix", type=str, required=False,
        default="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/bbj_projection/02.bbj_pca/bbj.maf.no_high_ld",
        help="BBJ 的 PLINK 前缀（已完成 bbj_prepare & bbj_pca）"
    )
    parser.add_argument(
        "--bbj_prune_in", type=str, required=False,
        default="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/bbj_projection/02.bbj_pca/bbj.maf.no_high_ld.prune.in",
        help="BBJ 的 prune.in 文件路径"
    )
    parser.add_argument(
        "--my_prefix_out", type=str, required=False, default="cteph_agp3k",
        help="我的数据在流程中的输出基名"
    )
    parser.add_argument(
        "--bbj_prefix_out", type=str, required=False, default="bbj",
        help="BBJ 数据在流程中的输出基名"
    )
    parser.add_argument(
        "--bbj_pca_acount", type=str, required=False,
        default="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/bbj_projection/02.bbj_pca/bbj.maf.no_high_ld.prune.pca.acount",
        help="plink2 --read-freq 使用的 BBJ PCA acount 文件路径"
    )
    parser.add_argument(
        "--bbj_pca_eigenvec_allele", type=str, required=False,
        default="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/bbj_projection/02.bbj_pca/bbj.maf.no_high_ld.prune.pca.eigenvec.allele",
        help="plink2 --score 使用的 BBJ PCA eigenvec.allele 文件路径"
    )
    parser.add_argument(
        "--threads", type=int, required=False, default=16,
        help="plink2 线程数（默认：16）"
    )

    # —— 可视化分组参数 ——
    parser.add_argument(
        "--case_prefix", type=str, required=False, default="PHOM",
        help="病例样本 IID 前缀（用于识别病例）"
    )
    parser.add_argument(
        "--case_name", type=str, required=False, default="CTEPH",
        help="病例在图例中的名称"
    )
    parser.add_argument(
        "--control_name", type=str, required=False, default="AGP3K",
        help="对照在图例中的名称"
    )

    return parser.parse_args()


def main() -> None:
    """主流程：投影 → 可视化；异常时给出明确提示并退出非零。"""
    args = parse_args()

    # 1) 投影：生成 *.sscore
    try:
        sscore_path = run_bbj_projection(
            my_bed_prefix=args.my_bed_prefix,
            bbj_bed_prefix=args.bbj_bed_prefix,
            bbj_prune_in=args.bbj_prune_in,
            my_prefix_out=args.my_prefix_out,
            bbj_prefix_out=args.bbj_prefix_out,
            bbj_pca_acount=args.bbj_pca_acount,
            bbj_pca_eigenvec_allele=args.bbj_pca_eigenvec_allele,
            threads=args.threads,
        )
    except Exception as e:
        print(f"[错误] 投影阶段失败：{e}", file=sys.stderr)
        sys.exit(1)

    print(f"[信息] 投影结果：{sscore_path}")

    # 2) 可视化：PC 成对散点图（基于 sscore 的 PC*_AVG）
    try:
        pdf_path = plot_projection_pairwise_pdf(
            sscore_path=sscore_path,
            case_prefix=args.case_prefix,
            case_name=args.case_name,
            control_name=args.control_name,
            # 需要自定义输出文件名可在此传入 output_pdf=...
        )
    except Exception as e:
        print(f"[错误] 绘图阶段失败：{e}", file=sys.stderr)
        sys.exit(2)

    print(f"[信息] 图像 PDF：{pdf_path}")


if __name__ == "__main__":
    main()
