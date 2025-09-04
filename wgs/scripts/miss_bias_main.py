#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
专业 CLI 脚本：miss_bias_main.py
目的：以命令行方式串联运行
  1) run_test_missing
  2) plot_raincloud_from_manifest
  3) remove_variants_by_fdr_from_manifest

特性：
  - 始终以“开发模式”导入 miss_bias_tools（每一步调用前都 reload）
  - 通过 argparse 传参（使用下划线参数名）
  - 对关键步骤打印简短的中文进度与产出路径，便于快速 debug

用法示例：
  python miss_bias_main.py \
    --bed_prefix "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/19.tommo_panel_filter/cteph_agp3k.lowfreq_common" \
    --out_prefix_run "cteph_agp3k.missing_bias" \
    --plink_path "/home/b/b37974/plink" \
    --use_midp \
    --color_by_variant \
    --fdr_threshold 0.05 \
    --threads 16 \
    --out_prefix_remove "cteph_agp3k.lowfreq_common.rm_q_lt_0.05"
"""

import os
import sys
import argparse
import importlib
from typing import Any, Dict

# 将脚本目录加入 sys.path，便于相对导入
SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.append(SCRIPT_DIR)

import miss_bias_tools  # 初次导入

def _reload_tools():
    """开发模式：每次调用前均 reload 工具模块，确保改动即时生效。"""
    importlib.reload(miss_bias_tools)
    return miss_bias_tools


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="基于 PLINK --test-missing 的缺失率偏倚分析一体化 CLI（开发模式 reload）。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # 1) run_test_missing 所需
    p.add_argument("--bed_prefix", type=str, required=True, help="PLINK/PLINK2 的 bfile 前缀（不含扩展名）")
    p.add_argument("--out_prefix_run", type=str, default="cteph_agp3k.missing_bias", help="run_test_missing 输出前缀")
    p.add_argument("--plink_path", type=str, default="/home/b/b37974/plink", help="plink 可执行文件路径")
    p.add_argument("--use_midp", action="store_true", help="是否使用 mid-p（Fisher mid-p）")
    p.add_argument("--no_use_midp", dest="use_midp", action="store_false", help="关闭 mid-p")
    p.set_defaults(use_midp=True)  # 按你的示例默认 True

    # 2) plot_raincloud_from_manifest 所需
    p.add_argument("--color_by_variant", action="store_true", help="雨滴散点是否按 SNP/InDel 上色")
    p.add_argument("--no_color_by_variant", dest="color_by_variant", action="store_false", help="关闭按变体类型上色")
    p.set_defaults(color_by_variant=True)  # 按你的示例默认 True

    # 3) remove_variants_by_fdr_from_manifest 所需
    p.add_argument("--fdr_threshold", type=float, default=0.05, help="BH FDR 阈值（提取 q < 阈值 的变体）")
    p.add_argument("--threads", type=int, default=16, help="plink2 线程数")
    p.add_argument("--plink2_path", type=str, default="/home/b/b37974/plink2", help="plink2 可执行文件路径")
    p.add_argument("--out_prefix_remove", type=str, default="cteph_agp3k.lowfreq_common.rm_q_lt_0.05",
                   help="去除变体后的 bfile 输出前缀")

    # 可选：只运行到某一步
    p.add_argument("--stop_after_run", action="store_true", help="仅运行 run_test_missing 后退出")
    p.add_argument("--stop_after_plot", action="store_true", help="运行 run + plot 后退出")

    return p.parse_args()


def main() -> int:
    args = parse_args()

    # ---------------- 1) run_test_missing ----------------
    print("==> [1/3] 运行 run_test_missing ...")
    tools = _reload_tools()
    try:
        man: Dict[str, Any] = tools.run_test_missing(
            bed_prefix=args.bed_prefix,
            out_prefix=args.out_prefix_run,
            plink_path=args.plink_path,
            use_midp=args.use_midp
        )
    except Exception as e:
        print(f"[ERROR] run_test_missing 失败：{e}", file=sys.stderr)
        return 1

    print(f"    - manifest 写出：{man.get('manifest_path')}")
    print(f"    - missing_path ：{man.get('missing_path')}")

    if args.stop_after_run:
        print("已按 --stop_after_run 要求在第 1 步后退出。")
        return 0

    # ---------------- 2) plot_raincloud_from_manifest ----------------
    print("==> [2/3] 绘制云雨图 plot_raincloud_from_manifest ...")
    tools = _reload_tools()
    try:
        plot_info: Dict[str, Any] = tools.plot_raincloud_from_manifest(
            manifest_path=man,
            color_by_variant=args.color_by_variant
        )
    except Exception as e:
        print(f"[ERROR] plot_raincloud_from_manifest 失败：{e}", file=sys.stderr)
        return 2

    print(f"    - 图像输出：{plot_info.get('output_path')}")
    print(f"    - 计数摘要：q<0.05={plot_info.get('q_le_0.05')}, q<0.01={plot_info.get('q_le_0.01')}, q<0.001={plot_info.get('q_le_0.001')}")

    if args.stop_after_plot:
        print("已按 --stop_after_plot 要求在第 2 步后退出。")
        return 0

    # ---------------- 3) remove_variants_by_fdr_from_manifest ----------------
    print("==> [3/3] 依据 FDR 阈值移除变体 remove_variants_by_fdr_from_manifest ...")
    tools = _reload_tools()
    try:
        rm_info: Dict[str, Any] = tools.remove_variants_by_fdr_from_manifest(
            manifest_path=man,
            fdr_threshold=args.fdr_threshold,
            plink2_path=args.plink2_path,
            threads=args.threads,
            out_prefix=args.out_prefix_remove
        )
    except Exception as e:
        print(f"[ERROR] remove_variants_by_fdr_from_manifest 失败：{e}", file=sys.stderr)
        return 3

    print(f"    - 变体列表：{rm_info.get('snplist_path')}（n={rm_info.get('exclude_n')}）")
    print(f"    - 新 bfile ：{rm_info.get('out_bed')}")
    print(f"    - 运行日志：{rm_info.get('log_path')}")

    print("==> 全流程完成。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
