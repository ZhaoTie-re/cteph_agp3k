# -*- coding: utf-8 -*-
"""
模块名称：bbj_sample_keep_main
=================================

【概述】
本脚本用于在 BBJ 投影（projection）坐标空间中选择目标样本，并基于所选样本从 PLINK 二进制基因型文件中导出子集。
流程包括两步：
1) 调用 `plot_projection_pc1_pc2_and_select` 绘制 PC1 vs PC2 散点图，并在给定矩形范围内自动选择样本；
2) 调用 `export_subset_bed` 使用 `--keep` 导出对应个体的 PLINK 子集（.bed/.bim/.fam）。

【适用场景】
- 使用 BBJ/参考面板的投影坐标进行样本聚类与筛选；
- 需要批量、可重复地生成“选择样本名单（keep.txt）”与对应的子集基因型数据。

【输入参数（CLI）】
- `--sscore_path` (str, 必填)：BBJ 投影后的 .sscore 文件路径。
- `--case_prefix` (str, 默认 "PHOM")：病例样本 IID 前缀，用于在图上分组着色。
- `--case_name` (str, 默认 "CTEPH")：病例组显示名称。
- `--control_name` (str, 默认 "AGP3K")：对照组显示名称。
- `--prefix_out` (str, 默认 "cteph_agp3k")：输出前缀（用于图片、keep 名单、以及导出的子集前缀）。
- `--rect_xlim` (str, 默认 "-0.028 0.016")：矩形选择的 X 轴范围，格式为 "xmin xmax"。
- `--rect_ylim` (str, 默认 "-0.024 0.033")：矩形选择的 Y 轴范围，格式为 "ymin ymax"。
- `--bed_prefix` (str, 必填)：输入的 PLINK 二进制前缀（不含扩展名）。
- `--threads` (int, 默认 16)：用于 PLINK 的线程数。
- `--reload_tools` (flag，可选)：开发调试用，强制重新加载 `bbj_projection_tools` 模块。

【输出】
- PNG 散点图：`{prefix_out}.bbj_projection.pc1_pc2.png`
- 保留样本名单：`{prefix_out}.bbj_projection.keep.txt`
- 导出的 PLINK 子集：`{prefix_out}.bbj_projected.subset.(bed|bim|fam)`

【使用示例】
python bbj_sample_keep_main.py \
  --sscore_path /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/bbj_projection/03.bbj_projection/cteph_agp3k.bbj.projection.sscore \
  --case_prefix PHOM \
  --case_name CTEPH \
  --control_name AGP3K \
  --prefix_out cteph_agp3k \
  --rect_xlim -0.028 0.016 \
  --rect_ylim -0.024 0.033 \
  --bed_prefix /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/13.run_variant_qc/cteph_agp3k.sqc.vqc \
  --threads 16

【注意事项】
- `--rect_xlim/--rect_ylim` 为空时将不使用矩形筛选（传入空字符串或 "none"/"None" 均视为 None）。
- 本脚本不会修改 `bbj_projection_tools` 内部逻辑，仅作为封装入口。
- 建议在 Conda/venv 环境下运行，确保依赖库（pandas、matplotlib、plink2 可执行程序等）可用。

作者: ZHAO TIE

"""

from __future__ import annotations
import argparse
import logging
import os
import sys
from typing import Optional, Tuple

# 确保可以导入同目录下的工具模块
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.append(SCRIPT_DIR)

import importlib
import bbj_projection_tools  # noqa: E402

# 需要的函数从工具模块引入
from bbj_projection_tools import (  # noqa: E402
    plot_projection_pc1_pc2_and_select,
    export_subset_bed,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bbj_sample_keep_main",
        description=(
            "在 BBJ 投影坐标上选样并导出 PLINK 子集的 CLI 工具。"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # 基础输入
    parser.add_argument("--sscore_path", required=True, help="BBJ 投影 .sscore 路径")
    parser.add_argument("--bed_prefix", required=True, help="输入 PLINK 二进制前缀")

    # 分组与命名
    parser.add_argument("--case_prefix", default="PHOM", help="病例样本 IID 前缀")
    parser.add_argument("--case_name", default="CTEPH", help="病例组显示名称")
    parser.add_argument("--control_name", default="AGP3K", help="对照组显示名称")
    parser.add_argument("--prefix_out", default="cteph_agp3k", help="输出前缀")

    # 选择区域
    parser.add_argument(
        "--rect_xlim",
        type=float,
        nargs=2,
        default=(-0.025, 0.016),
        help="X limits for the rectangle (default: (-0.025, 0.016))",
    )
    parser.add_argument(
        "--rect_ylim",
        type=float,
        nargs=2,
        default=(-0.025, 0.025),
        help="Y limits for the rectangle (default: (-0.025, 0.025))",
    )

    # 运行控制
    parser.add_argument("--threads", type=int, default=16, help="PLINK 线程数")
    parser.add_argument("--reload_tools", action="store_true", help="开发调试用：强制 reload 模块")

    return parser


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    # 日志设置
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
    )
    logging.info("参数解析完成：%s", vars(args))

    # 直接赋值矩形范围
    rect_xlim: Optional[Tuple[float, float]] = tuple(args.rect_xlim)
    rect_ylim: Optional[Tuple[float, float]] = tuple(args.rect_ylim)

    # 开发调试：可选强制 reload
    if args.reload_tools:
        logging.info("重新加载模块 bbj_projection_tools …")
        importlib.reload(bbj_projection_tools)  # noqa: F401
        # 重新导入函数以确保指向最新实现
        from bbj_projection_tools import (  # type: ignore
            plot_projection_pc1_pc2_and_select as _plot,
            export_subset_bed as _export,
        )
    else:
        _plot = plot_projection_pc1_pc2_and_select
        _export = export_subset_bed

    # Step 1: 投影散点图 + 选择
    logging.info("[Step 1] 绘制投影散点并选择矩形范围内样本 …")
    png_path, keep_path = _plot(
        sscore_path=args.sscore_path,
        case_prefix=args.case_prefix,
        case_name=args.case_name,
        control_name=args.control_name,
        prefix_out=args.prefix_out,
        rect_xlim=rect_xlim,
        rect_ylim=rect_ylim,
    )
    logging.info("生成图片: %s", png_path)
    logging.info("保留名单: %s", keep_path)

    # Step 2: 导出子集基因型
    logging.info("[Step 2] 使用 --keep 导出 PLINK 子集 …")
    subset_prefix = _export(
        bed_prefix=args.bed_prefix,
        keep_txt_path=keep_path,
        out_prefix=f"{args.prefix_out}.bbj_projected.subset",
        threads=args.threads,
    )

    # 控制台总结输出
    print("PNG:", png_path)
    print("KEEP:", keep_path)
    print("SUBSET_PREFIX:", subset_prefix)
    logging.info("处理完成 ✅")


if __name__ == "__main__":
    main()