#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本名称：panel_filter_main.py
=================================
作者：ZHAO TIE
版本：v1.0

【用途概述】
本脚本将“面板对比与选择”的三个阶段串联为**一键式命令行流程**：
  1) 从 `manifest.json` 读取主表与阈值信息，生成**五列表** `out_summary`（含 `FILTER_STAT`）。
  2) 根据 `panel_select_config.json` 中的规则（逐组 rare/lowfreq/common，按四类 not_in_tommo / in_tommo_nonpass / in_tommo_pass_fail / in_tommo_pass_pass）
     对 `out_summary` 进行**灵活筛选**并输出每组的 `VARIANT_ID` 列表。
  3) 根据筛选出的变体列表，调用 **plink2** 对给定 `bed_prefix` 的基因型进行**子集提取**。可选择将
     lowfreq 与 common 合并提取，或分别提取三组。

【输入文件】
- `manifest.json`：包含字段
    - `input`: 大表路径（含 VARIANT_ID, CTRL_MAF, IN_TOMMO, TOMMO_FILTER 等）
    - `c_in_pass_knee`: 各组（rare/lowfreq/common）的 selected_variants_tsv 等信息
- `panel_select_config.json`：选择规则（每组四类布尔开关，详见示例）

【输出文件】
- `<input>.summary_filter.tsv`：五列结果表
  （`VARIANT_ID, GROUP, IN_TOMMO, PASS_TOMMO, PASS_GROUP_ROBUST_Z_FILTER, FILTER_STAT`）
- `<input>.summary_filter.counts.tsv`：统计（`GROUP × FILTER_STAT`）
- `*.selected_variants.tsv`：各组（或合并）用于 plink2 `--extract` 的变体列表
- `<out_prefix>.(rare|lowfreq|common|lowfreq_common).bed/bim/fam`：plink2 子集化结果

【FILTER_STAT 含义】
- `Stat_1`: IN_TOMMO==True 且 PASS_TOMMO==True 且 PASS_GROUP_ROBUST_Z_FILTER==True
- `Stat_2`: IN_TOMMO==True 且 PASS_TOMMO==True 且 PASS_GROUP_ROBUST_Z_FILTER==False
- `Stat_3`: IN_TOMMO==True 且 PASS_TOMMO==False
- `Stat_4`: IN_TOMMO==False

【用法示例】
```bash
python panel_filter_main.py \
  --manifest_path /LARGE0/.../manifest.json \
  --config_json   /LARGE0/.../panel_select_config.json \
  --bed_prefix    /LARGE0/.../cteph_agp3k.rand300x100000 \
  --out_prefix    cteph_agp3k \
  --threads 8 --chunk_size 1000 --max_workers 8 \
  --keep_tmp \
  --merge_low_common \
  --plink2_path /home/b/b37974/plink2 \
  --save_out_map --out_map_path /LARGE0/.../out_map.json
```

【专业提示】
- 确保 `.bim` 第二列的变体 ID 与各 `*.selected_variants.tsv` 中 `VARIANT_ID` 格式**完全一致**（通常推荐 `CHROM:POS:REF:ALT`）。
- 如出现 “0 variants remaining” 或者 plink2 报错，先用 `awk '{print $2}' *.bim | head` 检查 ID 格式；必要时在上游统一。
- 大文件下可增大 `--chunk_size` 并适度提高 `--max_workers`（I/O 成为瓶颈时不宜过高）。
"""

from __future__ import annotations
import argparse
import json
import os
import sys
import time

# 将当前脚本目录加入 sys.path，确保本地模块可导入
HERE = os.path.abspath(os.path.dirname(__file__))
if HERE not in sys.path:
    sys.path.append(HERE)

# 导入工具模块
import importlib
import variant_qc_calculator
import panel_compare_tools

# 可在开发阶段强制 reload（正式环境可关闭）
importlib.reload(variant_qc_calculator)
importlib.reload(panel_compare_tools)


def _bool_flag(parser: argparse.ArgumentParser, name: str, default: bool, help_true: str, help_false: str) -> None:
    """在 argparse 中添加互斥的布尔开关：--name / --name_off。
    仅使用下划线，不使用连字符，避免与变量命名风格不一致。
    """
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument(f"--{name}", dest=name, action="store_true", help=help_true)
    group.add_argument(f"--{name}_off", dest=name, action="store_false", help=help_false)
    parser.set_defaults(**{name: default})


def build_parser() -> argparse.ArgumentParser:
    """构建命令行参数解析器（全部中文帮助）。"""
    p = argparse.ArgumentParser(
        prog="panel_filter_main",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "基于 manifest → summary → 选择 → plink2 子集化 的一键流程。\n"
            "作者: ZHAO TIE"
        ),
    )
    # —— 阶段 1：summary ——
    p.add_argument("--manifest_path", required=True, help="manifest.json 的路径")
    p.add_argument("--chunk_size", type=int, default=1000, help="分块大小（pandas read_csv chunksize）")
    p.add_argument("--max_workers", type=int, default=8, help="并行 worker 数量（优先考虑 I/O 负载")
    _bool_flag(p, name="keep_tmp", default=True, help_true="保留可能的临时文件/调试中间件", help_false="不保留临时文件")

    # —— 阶段 2：选择 ——
    p.add_argument("--config_json", required=True, help="变体选择配置 JSON（per-group per-category 开关）")

    # —— 阶段 3：plink2 子集 ——
    p.add_argument("--bed_prefix", required=True, help="输入 plink 二进制前缀（.bed/.bim/.fam）")
    p.add_argument("--out_prefix", default="cteph_agp3k", help="plink2 输出前缀（将自动追加组名后缀）")
    p.add_argument("--threads", type=int, default=8, help="plink2 --threads，并用于部分内部并行")
    _bool_flag(p, name="merge_low_common", default=True,
               help_true="合并 lowfreq+common 为一个集合进行子集化",
               help_false="不合并；分别输出 lowfreq 与 common 子集")
    p.add_argument("--plink2_path", default="/home/b/b37974/plink2", help="plink2 可执行文件路径")

    # —— 额外输出控制 ——
    _bool_flag(p, name="save_out_map", default=False,
               help_true="保存阶段 2 选择结果 out_map 为 JSON 文件",
               help_false="不保存 out_map JSON")
    p.add_argument("--out_map_path", default=None, help="保存 out_map 的 JSON 路径（未提供则根据 out_prefix 推断）")

    return p


def _print_kv(title: str, d: dict) -> None:
    """小工具：以整齐的中文键值格式打印参数或路径字典。"""
    print(title)
    for k, v in d.items():
        print(f"  - {k}: {v}")


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    # —— 参数回显（中文）——
    print("========== 参数确认 ==========")
    _print_kv("[路径参数]", {
        "manifest_path": args.manifest_path,
        "config_json": args.config_json,
        "bed_prefix": args.bed_prefix,
        "out_prefix": args.out_prefix,
        "plink2_path": args.plink2_path,
    })
    _print_kv("[运行参数]", {
        "chunk_size": args.chunk_size,
        "max_workers": args.max_workers,
        "threads": args.threads,
        "keep_tmp": args.keep_tmp,
        "merge_low_common": args.merge_low_common,
        "save_out_map": args.save_out_map,
        "out_map_path": args.out_map_path or "(未指定，默认用 <out_prefix>.selected_map.json)",
    })

    # —— 计时器 ——
    t_all = time.time()

    # —— 阶段 1：summary ——
    print("\n========== 阶段 1/3：生成 summary（含 FILTER_STAT） ==========")
    try:
        out_summary = panel_compare_tools.summarize_variants_filter_from_manifest(
            manifest_path=args.manifest_path,
            chunk_size=args.chunk_size,
            max_workers=args.max_workers,
            keep_tmp=args.keep_tmp,
        )
        print(f"[完成] out_summary 生成 → {out_summary}")
    except Exception as e:
        print(f"[错误] 生成 summary 失败：{e}")
        return 2

    # —— 阶段 2：按 JSON 规则筛选 ——
    print("\n========== 阶段 2/3：按 JSON 规则筛选变体 ==========")
    try:
        out_map = panel_compare_tools.filter_variants_by_group_and_stat(
            out_summary=out_summary,
            config_json=args.config_json,
            chunk_size=args.chunk_size,
            max_workers=args.max_workers,
        )
        print("[完成] 选择结果（各组路径）：")
        print(json.dumps(out_map, indent=2, ensure_ascii=False))

        if args.save_out_map:
            # 推断保存路径：优先使用 --out_map_path；否则使用 <out_prefix>.selected_map.json 放在当前工作目录
            save_path = args.out_map_path
            if not save_path:
                base = os.path.basename(args.out_prefix) if args.out_prefix else "out"
                save_path = os.path.abspath(f"{base}.selected_map.json")
            try:
                with open(save_path, 'w', encoding='utf-8') as f:
                    json.dump(out_map, f, indent=2, ensure_ascii=False)
                print(f"[完成] 已保存 out_map JSON → {save_path}")
            except Exception as e:
                print(f"[警告] 无法保存 out_map JSON（已忽略）：{e}")
    except Exception as e:
        print(f"[错误] 变体筛选失败：{e}")
        return 3

    # —— 阶段 3：plink2 子集化 ——
    print("\n========== 阶段 3/3：plink2 子集化 ==========")
    try:
        outs = panel_compare_tools.subset_plink_by_selected_variants(
            bed_prefix=args.bed_prefix,
            out_prefix=args.out_prefix,
            selected_paths=out_map,
            threads=args.threads,
            merge_low_common=args.merge_low_common,
            plink2_path=args.plink2_path,
        )
        print("[完成] 子集化输出（前缀）如下：")
        print(json.dumps(outs, indent=2, ensure_ascii=False))
    except Exception as e:
        print(f"[错误] plink2 子集化失败：{e}")
        return 4

    # —— 总结 ——
    dt_all = time.time() - t_all
    print("\n========== 全流程完成 ==========")
    print(f"总耗时：{dt_all/60:.2f} 分钟")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
