# -*- coding: utf-8 -*-
"""
脚本名称：panel_thr_main.py
=================================

【概述】
本脚本串联运行三个步骤，用于基于 ToMMo 对齐结果进行变体分组、
阈值扫描与拐点分析，并输出质量—数量权衡图及 knee 条件下的散点页：
  1) build_grouped_variant_tables：
     - 读取 `variant_qc_with_tommo.tsv`，按 CTRL_MAF → 主分组（rare/lowfreq/common），
       再按 IN_TOMMO/TOMMO_FILTER → 亚分组（a/b/c），
       对 c_in_pass 计算 ROBUST_Z（SciPy MAD）。
  2) summarize_c_in_pass_thresholds：
     - 仅针对三大主分组的 c_in_pass：阈值 |ROBUST_Z|=1..max_int（步长 0.2）
       计算每个阈值下：计数（count_variants）与 MSE(CTRL_AAF vs TOMMO_AAF)。
  3) plot_c_in_pass_threshold_tradeoff：
     - 在一页上绘制三联图（X=mse_ctrl_vs_tommo，Y=count_variants），
       并使用 KneeLocator 寻找拐点、添加星标与注释；
       额外增加一页：在拐点阈值下（|ROBUST_Z| < knee_thr）落入计数的变体，
       绘制 TOMMO_AAF (x) vs CTRL_AAF (y) 的散点图（PASS; 按 SNP/InDel 上色）。

【输入】
- variant_qc_with_tommo：run_plink2_variant_qc_with_tommo 的输出 TSV 路径。

【输出】
- manifest.json（由上述函数自动生成/更新，含统计与文件路径）
- tradeoff PDF（含 knee 标注页与 knee 条件散点页）

【使用示例】
python panel_thr_main.py \
  --variant_qc_with_tommo \
    "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/scripts/cteph_agp3k.variant_qc_summary.variant_qc_with_tommo.tsv" \
  --chunk_size 10000 \
  --knee_weight_y_map '{"rare": 1.0, "lowfreq": 2.3, "common": 2.3}'

注：knee_weight_y_map 为 JSON 字符串，键为主分组（rare/lowfreq/common）。

作者: ZHAO TIE
"""

from __future__ import annotations
import argparse
import json
import os
import sys
from typing import Dict, Optional

# 第三方库
import matplotlib.pyplot as plt  # noqa: F401  # 某些环境中后续函数内部使用

# 将工具模块加入路径
SCRIPT_DIR = os.path.abspath(
    "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/scripts"
)
if SCRIPT_DIR not in sys.path:
    sys.path.append(SCRIPT_DIR)

# 项目内工具函数
import importlib  # noqa: E402
import panel_compare_tools  # noqa: E402

# 开发调试：确保每次运行都加载到最新修改
importlib.reload(panel_compare_tools)
from panel_compare_tools import (  # noqa: E402
    build_grouped_variant_tables,
    summarize_c_in_pass_thresholds,
    plot_c_in_pass_threshold_tradeoff,
)


def parse_knee_weight_y_map(s: Optional[str]) -> Optional[Dict[str, float]]:
    """解析 knee_weight_y_map 的 JSON 字符串为字典。
    允许为 None；若提供则需能被 json.loads 正确解析。
    """
    if s is None:
        return None
    try:
        obj = json.loads(s)
        if not isinstance(obj, dict):
            raise ValueError("knee_weight_y_map 必须是字典 JSON（例如 '{\"rare\":2.0,...}' )")
        # 将 value 转为 float
        return {str(k): float(v) for k, v in obj.items()}
    except Exception as e:
        raise argparse.ArgumentTypeError(f"无法解析 knee_weight_y_map：{e}")


def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="panel_thr_main",
        description=(
            "基于 ToMMo 结果的 c_in_pass 变体阈值扫描与拐点分析，"
            "生成权衡曲线 PDF 以及 knee 条件下的 AAF 散点页。"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--variant_qc_with_tommo",
        type=str,
        required=True,
        help="variant_qc_with_tommo.tsv 的路径（run_plink2_variant_qc_with_tommo 的输出）",
    )
    parser.add_argument(
        "--chunk_size",
        type=int,
        default=10_000,
        help="分块读取大小（行数），用于大表内存友好处理",
    )
    parser.add_argument(
        "--knee_weight_y_map",
        type=parse_knee_weight_y_map,
        default={"rare": 1.0, "lowfreq": 2.3, "common": 2.3},
        help=(
            "KneeLocator 的按组权重（weight_y），JSON 字符串。"
            "键为主分组：rare/lowfreq/common；例如 '{\"rare\":1.0,\"lowfreq\":2.3,\"common\":2.3}'"
        ),
    )

    # 可选：图像渲染相关（保留专业化扩展接口）
    parser.add_argument(
        "--png_dpi",
        type=int,
        default=600,
        help="新增散点页（knee 条件）的 PNG 渲染分辨率 (DPI)",
    )
    parser.add_argument(
        "--no_keep_tmp",
        action="store_true",
        help="构建分组阶段完成后删除临时目录（默认保留）",
    )

    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_argparser().parse_args(argv)

    variant_qc_with_tommo: str = args.variant_qc_with_tommo
    chunk_size: int = args.chunk_size
    knee_weight_y_map: Optional[Dict[str, float]] = args.knee_weight_y_map

    # Step 1. 构建分组表与 manifest
    print("[Step1] 构建分组表并计算 c_in_pass 的 ROBUST_Z …", flush=True)
    manifest_path = build_grouped_variant_tables(
        variant_qc_with_tommo=variant_qc_with_tommo,
        output_dir=None,            # 最终输出放当前工作目录
        chunk_size=chunk_size,
        keep_tmp=not args.no_keep_tmp,
    )
    print(f"[Step1] 完成 → manifest: {manifest_path}")

    # Step 2. 概要阈值扫描（|ROBUST_Z|）
    print("[Step2] 汇总 c_in_pass：阈值 1..max_int（步长 0.2）下的 count 与 MSE …", flush=True)
    manifest_path = summarize_c_in_pass_thresholds(manifest_path)
    print(f"[Step2] 完成 → manifest 更新: {manifest_path}")

    # Step 3. 绘图（含 Knee 标注页 + Knee 条件散点页）
    print("[Step3] 绘制权衡曲线 PDF，并写回 knee 注释信息到 manifest.json …", flush=True)
    out_pdf = plot_c_in_pass_threshold_tradeoff(
        manifest_path=manifest_path,
        figsize=(14, 5),
        marker_size=18,
        line_width=1.5,
        y_min_zero=False,
        use_log_y=False,
        png_dpi=args.png_dpi,
        # KneeLocator 参数（如需统一设定权重，可在工具函数内提供全局 weight_y；
        # 这里按组指定 weight_y 覆盖）
        knee_weight_y_map=knee_weight_y_map,
    )
    print(f"[Step3] 完成 → PDF: {out_pdf}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
