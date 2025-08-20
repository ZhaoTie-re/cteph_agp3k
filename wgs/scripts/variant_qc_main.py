"""
变体质控主脚本（Variant QC Main）
=================================
用途：
  1) 调用 PLINK2 计算变体级 QC 指标（缺失率、HWE 等）。
  2) 按 MAF 分层绘制 VMISS 分布并基于阈值筛选。
  3) 绘制分层 HWE 散点并按阈值筛选。
  4) 取交集导出最终通过的变体列表，并可选将其从 bed 中导出。

命令行参数（可选项均有默认值）：
  --script_path        : 本工具集脚本所在目录，用于动态导入本地模块。
  --bed_prefix         : 输入 PLINK 二进制文件前缀（.bed/.bim/.fam）。
  --output_prefix      : 输出前缀（不带扩展名）。
  --threads            : 运行 PLINK2 的线程数。
  --vmiss_threshold    : 变体缺失率阈值（例如 0.05 表示 5%）。
  --hwe_json           : HWE 阈值字典的 JSON 文件路径（可选）。若未提供，使用内置默认值。

HWE 阈值 JSON 示例（传给 --hwe_json）：
{
  "Rare Variant (<0.01)": {"CTRL_HWE": null, "CASE_HWE": null},
  "Low Frequency Variant (0.01~0.05)": {"CTRL_HWE": 1e-6, "CASE_HWE": null},
  "Common Variant (>0.05)": {"CTRL_HWE": 1e-6, "CASE_HWE": 1e-10}
}

示例：
  python variant_qc_main.py \
    --script_path /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/scripts \
    --bed_prefix /LARGE0/.../cteph_agp3k.rm_maf0_vmiss1 \
    --output_prefix cteph_agp3k.sqc \
    --threads 16 \
    --vmiss_threshold 0.05

作者: ZHAO TIE

"""

from __future__ import annotations

import argparse
import json
import os
import sys

from importlib import import_module, reload

# 这些导入预留给可能的扩展（当前脚本内不直接使用）
import pandas as pd  # noqa: F401
import matplotlib  # noqa: F401
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: F401

# ----------------------------
# HWE 阈值配置校验函数
# ----------------------------
def validate_hwe_thresholds(data: dict) -> dict:
    """Validate and normalize HWE threshold mapping.
    Expected format:
    {
        "<Category>": {"CTRL_HWE": <float|None>, "CASE_HWE": <float|None>},
        ...
    }
    Returns the input dict if valid; exits with code 2 on error.
    """
    if not isinstance(data, dict):
        print("[错误] --hwe_json 解析后不是字典对象。")
        sys.exit(2)
    required_keys = {"CTRL_HWE", "CASE_HWE"}
    for cat, thr in data.items():
        if not isinstance(thr, dict):
            print(f"[错误] HWE 阈值类别 '{cat}' 的值应为字典。")
            sys.exit(2)
        if set(thr.keys()) != required_keys:
            print(f"[错误] HWE 阈值类别 '{cat}' 的键应为 {required_keys}，实际为 {set(thr.keys())}。")
            sys.exit(2)
        for k in required_keys:
            v = thr[k]
            if v is not None and not isinstance(v, (int, float)):
                print(f"[错误] HWE 阈值 '{cat}.{k}' 必须为数值或 null/None，实际为 {type(v)}。")
                sys.exit(2)
    return data


# ----------------------------
# 主程序
# ----------------------------

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
    "--script_path",
    default="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/scripts",
    help="本地工具脚本目录（用于导入 variant_qc_calculator / variant_qc_flags）",
)
parser.add_argument(
    "--bed_prefix",
    required=True,
    metavar="/path/to/prefix",
    help="输入 PLINK 二进制文件前缀",
)
parser.add_argument(
    "--output_prefix",
    default="cteph_agp3k.sqc",
    help="输出前缀",
)
parser.add_argument(
    "--threads",
    type=int,
    default=16,
    help="PLINK2 线程数",
)
parser.add_argument(
    "--vmiss_threshold",
    type=float,
    default=0.05,
    help="变体缺失率筛选阈值",
)
parser.add_argument(
    "--hwe_json",
    default=None,
    help="HWE 阈值配置（JSON 文件路径）；若未提供，使用内置默认",
)
args = parser.parse_args()

def main() -> None:
    # 基本参数健壮性检查
    if args.threads is None or args.threads < 1:
        print("[错误] --threads 必须为 >= 1 的整数。")
        sys.exit(2)

    # 确保输出目录存在（如果 output_prefix 包含目录）
    out_dir = os.path.dirname(os.path.abspath(args.output_prefix))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # 校验 PLINK 前缀文件是否存在
    required_exts = ["bed", "bim", "fam"]
    missing = [ext for ext in required_exts if not os.path.exists(f"{args.bed_prefix}.{ext}")]
    if missing:
        missing_paths = ", ".join(f"{args.bed_prefix}." + ext for ext in missing)
        print(f"[错误] 找不到以下 PLINK 文件：{missing_paths}")
        sys.exit(2)

    # 1) 动态导入并统一强制 reload
    script_path = os.path.abspath(args.script_path)
    if script_path not in sys.path:
        sys.path.append(script_path)

    variant_qc_calculator = import_module("variant_qc_calculator")  # type: ignore
    variant_qc_flags = import_module("variant_qc_flags")            # type: ignore

    # 强制重新加载（开发/调试阶段保持最新）
    variant_qc_calculator = reload(variant_qc_calculator)
    variant_qc_flags = reload(variant_qc_flags)

    from variant_qc_calculator import (
        init_globals_for_chunk,  # noqa: F401  # 保留兼容
        process_chunk,           # noqa: F401  # 保留兼容
        run_plink2_variant_qc,
    )
    from variant_qc_flags import (
        plot_vmiss_distribution_by_maf_category,
        plot_hwe_scatter_by_maf_category,
        extract_pass_variants_by_intersection,
    )

    # 2) 载入 HWE 阈值
    if args.hwe_json:
        try:
            with open(args.hwe_json, "r", encoding="utf-8") as f:
                hwe_thresholds = json.load(f)
        except json.JSONDecodeError as e:
            print(f"[错误] 无法解析 --hwe_json 文件：{e}")
            sys.exit(2)
        hwe_thresholds = validate_hwe_thresholds(hwe_thresholds)
    else:
        hwe_thresholds = {
            "Rare Variant (<0.01)": {"CTRL_HWE": None, "CASE_HWE": None},
            "Low Frequency Variant (0.01~0.05)": {"CTRL_HWE": 1e-6, "CASE_HWE": None},
            "Common Variant (>0.05)": {"CTRL_HWE": 1e-6, "CASE_HWE": 1e-10},
        }

    # 3) 运行变体级 QC 汇总
    variant_qc_summary = run_plink2_variant_qc(
        bed_prefix=args.bed_prefix,
        output_prefix=args.output_prefix,
        threads=args.threads,
    )

    # 4) VMISS 分布与筛选（使用参数中的阈值）
    pass_vmiss = plot_vmiss_distribution_by_maf_category(
        variant_qc_summary,
        vmiss_threshold=args.vmiss_threshold,
        output_prefix=args.output_prefix,
    )

    # 5) HWE 散点与筛选
    pass_hwe = plot_hwe_scatter_by_maf_category(
        variant_qc_summary,
        hwe_thresholds,
        output_prefix=args.output_prefix,
    )

    # 6) 交集导出最终通过列表（并可选筛出 bed 子集）
    pass_variants_path = extract_pass_variants_by_intersection(
        pass_vmiss_path=pass_vmiss,
        pass_hwe_path=pass_hwe,
        bed_prefix=args.bed_prefix,
        output_prefix=args.output_prefix,
        threads=args.threads,
    )

    print("\n[QC 完成]")
    print(f"- VMISS 阈值: {args.vmiss_threshold}")
    print("- HWE 阈值配置:")
    for category, thresholds in hwe_thresholds.items():
        ctrl = thresholds.get('CTRL_HWE')
        case = thresholds.get('CASE_HWE')
        ctrl_s = ctrl if ctrl is not None else 'skip'
        case_s = case if case is not None else 'skip'
        print(f"  {category}: CTRL_HWE={ctrl_s}, CASE_HWE={case_s}")
    print(f"- 输出前缀: {args.output_prefix}")
    print(f"- 通过变体列表: {pass_variants_path}")


if __name__ == "__main__":
    main()
