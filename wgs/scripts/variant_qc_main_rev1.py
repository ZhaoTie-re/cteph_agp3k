#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
变体质控主流程（CLI 版本）
================================
将原先的 Notebook/脚本版本改为命令行工具，支持参数传入并打印详细中文日志。

功能阶段：
1) 运行 plink2 变体级 QC，生成汇总表（variant_qc_summary）
2) 基于 VMISS：
   - mode='dp' 或 'case_ctrl'：绘制散点图（30X_VMISS vs 15X_VMISS / CTRL_VMISS vs CASE_VMISS），导出通过变体
   - mode='mix'：从 vmiss.json 读取全局 VMISS 阈值，绘制分布图/过滤，导出通过变体
3) 基于 HWE：从 hwe.json 读取阈值，绘制散点图并导出通过变体
4) 交集：输出 VMISS 与 HWE 双通过的变体，并导出对应 VCF/PLINK 子集（由 variant_qc_flags 实现）

用法示例（使用默认值）：
    python variant_qc_main_rev1.py \\
        --threads 16 \\
        --info_path "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/cteph_agp3k_jhrpv4.rev1.xlsx" \\
        --sample_col "ID" \\
        --target_dp_col "Target DP (JHRPv4)" \\
        --bed_prefix "cteph_agp3k.rand200x5000" \\
        --output_prefix "cteph_agp3k.sqc.vqc" \\
        --vmiss_json_path "vmiss.json" \\
        --vmiss_mode "dp" \\
        --hwe_json_path "hwe.json" \\
        --script_path "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/scripts"

注：
- 本脚本会在运行时将 --script_path 动态加入 sys.path，以加载本地模块。
- 日志输出为中文，便于在 HPC 日志里快速定位问题。

作者: ZHAO TIE
"""

import os
import sys
import json
import time
import logging
import argparse
import importlib
from datetime import datetime


def _setup_logger():
    """配置日志（中文），含时间、级别、信息。"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger("variant_qc_main")


def parse_args():
    """解析命令行参数（全部提供默认值，用户可覆盖）。"""
    parser = argparse.ArgumentParser(
        description="CTEPH 项目：变体质控主流程（VMISS + HWE + 交集）",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--threads", type=int, default=16, help="线程数")
    parser.add_argument(
        "--info_path",
        type=str,
        default="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/cteph_agp3k_jhrpv4.rev1.xlsx",
        help="样本 meta 信息 Excel 路径",
    )
    parser.add_argument("--sample_col", type=str, default="ID", help="样本 ID 列名（info_path 中）")
    parser.add_argument(
        "--target_dp_col",
        type=str,
        default="Target DP (JHRPv4)",
        help="目标测序深度列名（15x / 30x）",
    )
    parser.add_argument("--bed_prefix", type=str, default="cteph_agp3k.rand200x5000", help="PLINK 基因型前缀")
    parser.add_argument("--output_prefix", type=str, default="cteph_agp3k.sqc.vqc", help="输出文件前缀")
    parser.add_argument("--vmiss_json_path", type=str, default="vmiss.json", help="VMISS 阈值 JSON")
    parser.add_argument(
        "--vmiss_mode",
        type=str,
        default="dp",
        choices=["dp", "case_ctrl", "mix"],
        help="VMISS 评估模式：dp / case_ctrl / mix",
    )
    parser.add_argument("--hwe_json_path", type=str, default="hwe.json", help="HWE 阈值 JSON")
    parser.add_argument(
        "--script_path",
        type=str,
        default="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/scripts",
        help="本项目 Python 模块路径（会动态加入 sys.path）",
    )
    # 可按需扩展：--plink2_path, --tmpdir 等
    return parser.parse_args()


def main():
    t0 = time.time()
    logger = _setup_logger()
    args = parse_args()

    # ---------- 环境准备 ----------
    logger.info("初始化环境与参数 ...")
    logger.info("参数汇总：%s", vars(args))
    if args.script_path not in sys.path:
        sys.path.insert(0, args.script_path)
        logger.info("已将 script_path 加入 sys.path：%s", args.script_path)

    # 动态加载本地模块
    try:
        import variant_qc_calculator_robust
        import variant_qc_flags
        # 每次运行都强制重新加载模块（开发阶段）
        importlib.reload(variant_qc_calculator_robust)
        importlib.reload(variant_qc_flags)
        logger.info("开发模式：已强制重新加载模块 variant_qc_calculator_robust 与 variant_qc_flags")
    except Exception as e:
        logger.error("加载或重新加载模块失败，请检查 --script_path 是否正确：%s", e)
        sys.exit(1)

    # ---------- 阶段 1：运行 plink2 变体 QC ----------
    logger.info("阶段1：运行 plink2 变体 QC ...")
    try:
        from variant_qc_calculator_robust import run_plink2_variant_qc
        variant_qc_summary = run_plink2_variant_qc(
            bed_prefix=args.bed_prefix,
            tmpdir="/tmp/variant_qc",  # 如需参数化，可自行添加 CLI 选项
            plink2_path="/home/b/b37974/plink2",  # 如需参数化，可自行添加 CLI 选项
            threads=args.threads,
            output_prefix=args.output_prefix,
            verbose=True,
            info_path=args.info_path,
            sample_col=args.sample_col,
            target_dp_col=args.target_dp_col,
        )
        logger.info("QC 汇总表生成完毕：%s", variant_qc_summary)
    except Exception as e:
        logger.exception("运行变体 QC 失败：%s", e)
        sys.exit(1)

    # ---------- 阶段 2：VMISS 评估 ----------
    logger.info("阶段2：基于 VMISS 的绘图与过滤，模式：%s ...", args.vmiss_mode)
    pass_vmiss_path = None
    try:
        if args.vmiss_mode in ("dp", "case_ctrl"):
            from variant_qc_calculator_robust import plot_vmiss_scatter_by_maf_category
            pass_vmiss_path = plot_vmiss_scatter_by_maf_category(
                variant_qc_summary=variant_qc_summary,
                vmiss_json_path=args.vmiss_json_path,
                mode=args.vmiss_mode,
                plot_style="hex",  # 可选：'hex', 'hist2d', 'kde2d'
                density_norm="log",  # 可选：'log', 'linear'
                output_prefix=args.output_prefix,
            )
            logger.info("VMISS 通过变体（%s 模式）导出路径：%s", args.vmiss_mode, pass_vmiss_path)

        elif args.vmiss_mode == "mix":
            with open(args.vmiss_json_path, "r") as f:
                cfg = json.load(f)
            try:
                vmiss_threshold = float(cfg["mix"]["VMISS"])
            except Exception as e:
                raise ValueError(f"无法从 {args.vmiss_json_path} 读取 mix 的 VMISS 阈值：{e}")

            from variant_qc_flags import plot_vmiss_distribution_by_maf_category
            pass_vmiss_path = plot_vmiss_distribution_by_maf_category(
                variant_qc_summary=variant_qc_summary,
                vmiss_threshold=vmiss_threshold,
                output_prefix=args.output_prefix,
            )
            logger.info("VMISS 通过变体（mix 模式）导出路径：%s", pass_vmiss_path)
        else:
            raise ValueError(f"无效的 vmiss_mode：{args.vmiss_mode}")
    except Exception as e:
        logger.exception("VMISS 阶段失败：%s", e)
        sys.exit(1)

    # ---------- 阶段 3：HWE 评估 ----------
    logger.info("阶段3：基于 HWE 的绘图与过滤 ...")
    try:
        from variant_qc_flags import plot_hwe_scatter_by_maf_category
        with open(args.hwe_json_path, "r") as f:
            hwe_thresholds = json.load(f)
        pass_hwe_path = plot_hwe_scatter_by_maf_category(
            variant_qc_summary=variant_qc_summary,
            hwe_thresholds=hwe_thresholds,
            output_prefix=args.output_prefix,
        )
        logger.info("HWE 通过变体导出路径：%s", pass_hwe_path)
    except Exception as e:
        logger.exception("HWE 阶段失败：%s", e)
        sys.exit(1)

    # ---------- 阶段 4：交集 ----------
    logger.info("阶段4：导出 VMISS ∩ HWE 的交集变体，并生成子集数据 ...")
    try:
        from variant_qc_flags import extract_pass_variants_by_intersection
        pass_variants_path = extract_pass_variants_by_intersection(
            pass_vmiss_path=pass_vmiss_path,
            pass_hwe_path=pass_hwe_path,
            bed_prefix=args.bed_prefix,
            output_prefix=args.output_prefix,
            threads=args.threads,
        )
        logger.info("交集变体导出路径：%s", pass_variants_path)
    except Exception as e:
        logger.exception("交集阶段失败：%s", e)
        sys.exit(1)

    # ---------- 总结 ----------
    dt = time.time() - t0
    logger.info("✅ 全流程完成（用时 %.1f 秒）。", dt)
    logger.info("输出摘要：")
    logger.info("- QC 汇总表：%s", variant_qc_summary)
    logger.info("- VMISS 通过：%s", pass_vmiss_path)
    logger.info("- HWE 通过：%s", pass_hwe_path)
    logger.info("- 交集结果：%s", pass_variants_path)

    # 为了兼容外部调用，打印最终路径（最后一行便于 shell 获取）
    print(pass_variants_path)


if __name__ == "__main__":
    main()
