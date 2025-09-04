#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
======================================================================
模块：miss_bias_tools.py
目的：围绕 PLINK 的 --test-missing 结果（.missing）进行BH校正、可视化、
      以及基于 FDR 阈值从 bfile 中移除可疑变体的实用工具集。

主要功能：
  1) run_test_missing：包装调用 PLINK/PLINK2 的 --test-missing，并将输出组织成 manifest；
     同时可选 mid-p / 置换检验（mperm）。结束后对 .missing 的 P 列进行 BH 矫正（SciPy）。
  2) plot_raincloud_from_manifest：基于 manifest 指向的 .missing + BH 列，绘制 -log10(q) 云雨图；
     仅按 BH(q) 进行计数；支持按 SNP/InDel 上色；阈值线默认绘制在 0.05/0.01/0.001 的位置。
  3) remove_variants_by_fdr_from_manifest：基于 BH(q) < 阈值 的变体列表，通过 plink2 从 bfile 中排除。

调用环境要求：
  - Python ≥ 3.8，SciPy 版本需支持 scipy.stats.false_discovery_control（若不支持将跳过 BH）。
  - 已安装 plink/plink2，对应可执行路径可在函数参数中指定。
  - Matplotlib 用于绘图；脚本在无显示（HPC）环境中也可运行（建议使用非交互式后端）。

使用示例（概要）：
  >>> mani = run_test_missing(bed_prefix="cohort", out_dir=".", use_midp=True)
  >>> summary = plot_raincloud_from_manifest(mani, color_by_variant=True)
  >>> res = remove_variants_by_fdr_from_manifest(mani, fdr_threshold=0.05, threads=16)

调试提示：
  - 本模块会把外部命令的完整 CMD、stdout/stderr 全量写入日志（*.run.log / *.remove_q.log）
    并把 BH 校正信息写入 manifest JSON，便于复现与排错。

作者：ZHAO TIE
======================================================================
"""

import json
import os
import shlex
import subprocess
from datetime import datetime
from typing import Dict, Any, Optional

import pandas as pd
import numpy as np
# SciPy BH (Benjamini–Hochberg) correction: controlled import to avoid hard crash on old SciPy
try:
    from scipy.stats import false_discovery_control as _scipy_fdc  # requires modern SciPy
except Exception as e:
    _scipy_fdc = None
    _SCIPY_IMPORT_ERROR = e
else:
    _SCIPY_IMPORT_ERROR = None
import scipy
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


def _check_bfile_exists(bed_prefix: str) -> None:
    """检查 {bed_prefix}.bed/.bim/.fam 是否存在。"""
    required = [f"{bed_prefix}.bed", f"{bed_prefix}.bim", f"{bed_prefix}.fam"]
    missing = [p for p in required if not os.path.isfile(p)]
    if missing:
        raise FileNotFoundError(f"缺少 bfile 组件: {missing}")


def _ensure_out_dir(out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)



def _open_log(path: str):
    return open(path, "w", encoding="utf-8")


# 日志时间戳辅助函数
def _log(log_fp, msg: str) -> None:
    """
    在给定日志文件句柄中写入一行带时间戳的调试信息。
    用途：在关键步骤增加可读的时间线，便于追踪问题。
    """
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    try:
        log_fp.write(f"[{ts}] {msg}\n")
    except Exception:
        # 兜底：即便日志写入失败也不影响主流程
        pass


def run_test_missing(
    bed_prefix: str,
    out_dir: str = ".",
    out_prefix: Optional[str] = None,
    plink_path: str = "/home/b/b37974/plink",
    use_permutation: bool = False,
    mperm: int = 100,
    use_midp: bool = True,
    perm_count: bool = False,
    extra_args: Optional[str] = None,
) -> Dict[str, Any]:
    """
    调用 PLINK 执行 --test-missing，返回输出文件路径清单。

    参数:
        bed_prefix     : bfile 前缀（必须存在 .bed/.bim/.fam）
        out_dir        : 输出目录（默认 "."）
        out_prefix     : 输出前缀（不含路径）；默认使用 basename(bed_prefix)+".missing_bias"
        plink_path     : plink 可执行文件路径
        use_permutation: 是否使用置换检验（默认 False）
        mperm          : 置换次数（仅当 use_permutation=True 时生效；典型值 100 或更高）
        use_midp       : 是否使用 Lancaster mid-p 修正
        perm_count     : 是否输出每次置换计数（仅当 use_permutation=True 且 PLINK 支持时生效）
        extra_args     : 透传给 plink 的其他参数（例如 "--threads 8 --allow-no-sex"）

    返回:
        dict 包含 out_prefix、生成文件路径等。
    """
    _check_bfile_exists(bed_prefix)
    _ensure_out_dir(out_dir)

    base_name = os.path.basename(bed_prefix)
    if out_prefix is None:
        out_prefix = f"{base_name}.missing_bias"

    out_prefix_path = os.path.join(out_dir, out_prefix)
    log_path = f"{out_prefix_path}.run.log"

    test_missing_tokens = ["--test-missing"]
    if use_midp:
        test_missing_tokens.append("midp")
    if use_permutation and mperm and mperm > 0:
        test_missing_tokens.append(f"mperm={mperm}")
        if perm_count:
            test_missing_tokens.append("perm-count")

    cmd = [
        plink_path,
        "--bfile", bed_prefix,
        *test_missing_tokens,
        "--out", out_prefix_path,
    ]

    if extra_args:
        # 将字符串拆分为 token，避免空格问题
        cmd.extend(shlex.split(extra_args))

    manifest: Dict[str, Any] = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "bed_prefix": bed_prefix,
        "plink_path": plink_path,
        "out_dir": out_dir,
        "out_prefix": out_prefix,
        "use_permutation": use_permutation,
        "mperm": mperm,
        "use_midp": use_midp,
        "perm_count": perm_count,
        "extra_args": extra_args,
        "cmd": cmd,
    }

    t0 = datetime.now()
    with _open_log(log_path) as log_fp:
        _log(log_fp, f"工作目录（cwd）: {os.getcwd()}")
        _log(log_fp, f"PLINK 路径: {plink_path}")
        _log(log_fp, f"bfile 前缀: {bed_prefix}")
        _log(log_fp, f"输出前缀路径: {out_prefix_path}")
        _log(log_fp, f"use_midp={use_midp}, use_permutation={use_permutation}, mperm={mperm}, perm_count={perm_count}")
        log_fp.write("[CMD]\n" + " ".join(shlex.quote(x) for x in cmd) + "\n\n")
        log_fp.flush()
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True)
        try:
            for line in proc.stdout:
                log_fp.write(line)
            ret = proc.wait()
            t1 = datetime.now()
            _log(log_fp, f"PLINK 运行结束；耗时(s): {(t1 - t0).total_seconds():.2f}")
            if ret != 0:
                log_fp.write(f"\n[ERROR] PLINK 运行失败，返回码: {ret}\n")
                raise RuntimeError(f"PLINK 执行失败，详见日志: {log_path}")
        finally:
            if proc.stdout:
                proc.stdout.close()

    missing_path = f"{out_prefix_path}.missing"
    if not os.path.isfile(missing_path):
        raise FileNotFoundError(f"未找到 PLINK 输出: {missing_path}")

    manifest.update({
        "missing_path": missing_path,
        "log_path": log_path,
        "manifest_path": f"{out_prefix_path}.manifest.json",
    })

    # 读取 .missing，进行 BH 矫正并回写新列（仅使用 SciPy）
    bh_info = {"bh_corrected": False, "bh_column": "P_BH", "bh_nonnull_n": 0, "bh_backend": None, "bh_version": getattr(scipy, "__version__", None)}
    try:
        if _scipy_fdc is None:
            raise RuntimeError(
                "SciPy 不支持 false_discovery_control（请升级 SciPy 到支持该函数的版本）。"
                f" 原始导入错误: {_SCIPY_IMPORT_ERROR}"
            )
        df = pd.read_csv(missing_path, delim_whitespace=True)
        if "P" in df.columns:
            p = pd.to_numeric(df["P"], errors="coerce")
            mask = p.notna() & np.isfinite(p)
            m = int(mask.sum())
            if m > 0:
                pv = p[mask].to_numpy()
                # 兼容不同 SciPy 版本的返回签名：
                # - 可能返回单个数组（调整后的 p 值）
                # - 也可能返回二元/三元组（通常第二个为调整后的 p 值）
                res_fdc = _scipy_fdc(pv, method="bh")
                bh_info["bh_return_type"] = str(type(res_fdc))
                try:
                    from collections.abc import Sequence as _Seq
                except Exception:
                    _Seq = (list, tuple)
                if isinstance(res_fdc, _Seq):
                    # 优先取第二个元素（常见为调整后的 p 值），否则取第一个
                    if len(res_fdc) >= 2:
                        adj_p = np.asarray(res_fdc[1], dtype=float)
                    else:
                        adj_p = np.asarray(res_fdc[0], dtype=float)
                    bh_info["bh_return_len"] = len(res_fdc)
                else:
                    adj_p = np.asarray(res_fdc, dtype=float)
                bh_info["bh_backend"] = "scipy.stats.false_discovery_control(bh)"
                adj_series = pd.Series(np.nan, index=p.index, dtype=float)
                adj_series.loc[mask] = adj_p
                df["P_BH"] = adj_series
                # 原子写入：先写到 .tmp 再替换
                tmp_path = missing_path + ".tmp"
                df.to_csv(tmp_path, sep="\t", index=False)
                os.replace(tmp_path, missing_path)
                bh_info.update({"bh_corrected": True, "bh_nonnull_n": m})
            else:
                bh_info["bh_reason"] = "无非缺失 P 值"
        else:
            bh_info["bh_reason"] = "未找到列 'P'"
    except Exception as e:
        bh_info["bh_error"] = str(e)
    manifest.update(bh_info)

    # 运行时长、工作目录等补充信息
    manifest["cwd"] = os.getcwd()
    try:
        manifest["duration_sec"] = (datetime.now() - t0).total_seconds()
    except Exception:
        pass

    # 将 BH 校正摘要追加写入运行日志，便于快速确认
    try:
        with open(log_path, "a", encoding="utf-8") as _ap:
            _log(_ap, f"BH 校正状态: {bh_info.get('bh_corrected')}, 列: {bh_info.get('bh_column')}, 非缺失数: {bh_info.get('bh_nonnull_n')}")
            if "bh_error" in bh_info:
                _log(_ap, f"BH 错误: {bh_info['bh_error']}")
            if "bh_reason" in bh_info:
                _log(_ap, f"BH 跳过原因: {bh_info['bh_reason']}")
            _log(_ap, f"manifest 已写出: {os.path.abspath(manifest['manifest_path'])}")
    except Exception:
        pass

    # 写出最终 manifest（包含 BH 信息）
    with open(manifest["manifest_path"], "w", encoding="utf-8") as fp:
        json.dump(manifest, fp, ensure_ascii=False, indent=2)

    return manifest


def _safe_read_missing_table(path: str) -> pd.DataFrame:
    """
    读取 .missing 表格（尽量鲁棒）：
    读取策略：
      1) 优先按制表符（\t）读取（因为我们写回时使用了 \t）
      2) 若只有 1 列，说明不是制表符分隔，回退为空白分隔（delim_whitespace=True）
      3) 如果仍失败，最终兜底为空白分隔读取
    返回：
      pandas.DataFrame
    调试：
      - 若解析异常，最终会回退到空白分隔读取；可在上层函数捕获并打印 path。
    """
    try:
        df = pd.read_csv(path, sep="\t")
        # 如果只有一列，可能是空白分隔，重读
        if df.shape[1] == 1:
            df = pd.read_csv(path, delim_whitespace=True)
        return df
    except Exception:
        # 最后兜底
        return pd.read_csv(path, delim_whitespace=True)

def plot_raincloud_from_manifest(
    manifest_path,
    title: Optional[str] = None,
    base_color: str = "#66c2a5",
    bw: float = 0.25,
    q_low: float = 0.00,
    q_high: float = 0.9999,
    max_points: int = 5000,
    output_path: Optional[str] = None,
    ax=None,
    color_by_variant: bool = False
) -> Dict[str, Any]:
    """
    根据 manifest 绘制云雨图（仅基于 P_BH/BH q 值）：
      - 读取 manifest 中的 missing_path 与 bh_column（默认 "P_BH"）
      - x 轴为 -log10(q)，仅对 q（BH）进行计数展示：<0.05、<0.01、<0.001
      - 云（半小提琴）右侧只画到“雨滴最大 -log10(q)”的位置
      - 若 color_by_variant=True，按 CHROM:POS:REF:ALT 判断 SNP（黑）/InDel（洋红）为雨滴着色
      - 阈值竖线仅绘制，不添加文字，具体数值显示在右侧信息框
    
    参数：
      manifest_path : dict 或 str
          run_test_missing 返回的 manifest 或其 JSON 路径
      title : str
          图标题；默认包含 -log10(bh_column) 与源文件名
      base_color : str
          基础配色（十六进制）
      bw : float
          KDE 带宽（传入 gaussian_kde 的 bw_method）
      q_low, q_high : float
          截取分位数范围以避免极端值影响 KDE（绘图稳健性）
      max_points : int
          雨滴最大采样点数（避免绘制极多点造成缓慢）
      output_path : str
          图保存路径；默认与 missing_path 同目录，命名为 "<basename>.raincloud.<bh_column>.png"
      ax : matplotlib.axes.Axes
          传入外部子图对象；若为 None 则内部创建
      color_by_variant : bool
          是否根据变体类型上色（SNP 黑、InDel 洋红）
    
    返回：
      dict，包含输出路径、样本量 N、阈值计数/比例等信息，便于在调用端记录与二次分析。
    """
    # 读取 manifest：既可接受 dict，也可接受 JSON 路径
    if isinstance(manifest_path, dict):
        mani = manifest_path
    else:
        with open(manifest_path, "r", encoding="utf-8") as fp:
            mani = json.load(fp)
    missing_path = mani.get("missing_path")
    if not missing_path or not os.path.isfile(missing_path):
        raise FileNotFoundError(f"manifest 中的 missing_path 不存在: {missing_path}")
    bh_col = mani.get("bh_column", "P_BH")

    # 强制使用默认样式
    plt.style.use('default')

    # 读取 missing 表
    miss_df = _safe_read_missing_table(missing_path)
    if bh_col not in miss_df.columns:
        raise KeyError(f"在 {missing_path} 中未找到列: {bh_col}")

    # 若需要按变体类型着色，定位变体列（优先 'SNP'，否则 'ID'）
    var_col = "SNP" if "SNP" in miss_df.columns else ("ID" if "ID" in miss_df.columns else None)
    if color_by_variant and var_col is None:
        raise KeyError("启用 color_by_variant 需要表中存在 'SNP' 或 'ID' 列以识别变体类型。")

    # --- 数据准备 ---
    col_p = bh_col
    p_all = pd.to_numeric(miss_df[col_p], errors="coerce").values
    q_all = p_all.copy()  # 这里即为 BH 调整后的 q 值（列名为 P_BH）

    # 作图数据：去掉 NaN
    mask = np.isfinite(p_all)
    p = p_all[mask]
    q = q_all[mask]

    eps = 1e-300
    x = -np.log10(np.clip(p, eps, None))   # -log10(q)

    # 若需要着色：构造与 x 对齐的 is_snp 数组（True=SNP, False=InDel）
    if color_by_variant:
        var_all = miss_df[var_col].astype(str).values
        var_masked = var_all[mask]

        def _is_snp_token(tok: str) -> bool:
            parts = tok.split(":")
            if len(parts) < 4:
                return False
            ref = parts[-2]
            alt = parts[-1]
            return (len(ref) == 1) and (len(alt) == 1) and (ref in "ACGT") and (alt in "ACGT")

        is_snp = np.array([_is_snp_token(t) for t in var_masked], dtype=bool)
    else:
        is_snp = None

    # robust 区间用于画半小提琴
    if x.size == 0:
        raise ValueError(f"{missing_path} 的列 {bh_col} 为空或无有效值，无法绘图。")
    lo, hi = np.quantile(x, [q_low, q_high])
    x_clip = x[(x >= lo) & (x <= hi)]
    if x_clip.size < 10:
        x_clip = x  # 太少则退化为全量

    # --- 画布 ---
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 4))
        created_fig = True

    # 默认标题
    if title is None:
        base_name = os.path.basename(missing_path)
        title = f"Raincloud of -log10({bh_col})\n{base_name}"

    # --- 先确定雨滴最大 x（用于限制云的右侧） ---
    n = x.size
    k = min(n, max_points)
    if k > 0:
        idx = np.random.choice(n, k, replace=False)
        xmax_drop = float(np.max(x[idx]))
    else:
        idx = np.array([], dtype=int)
        xmax_drop = float(np.max(x)) if n else 0.0

    # --- 半小提琴（云） --- 右侧仅绘制到雨滴最大 -log10(q)
    kde = gaussian_kde(x_clip, bw_method=bw)
    left  = x_clip.min()
    right = xmax_drop
    if right <= left:
        x_grid = np.linspace(left, left, 2)
        dens = np.zeros_like(x_grid)
    else:
        x_grid = np.linspace(left, right, 500)
        dens   = kde(x_grid)
        if np.max(dens) > 0:
            dens = dens / dens.max() * 0.8
        else:
            dens = np.zeros_like(dens)
    y0     = 0.0
    y_top  = y0 + 0.35
    ax.fill_between(x_grid, y_top, y_top + dens,
                    facecolor=base_color, edgecolor=base_color, alpha=0.35, linewidth=1.2)

    # --- 盒须（中位数/四分位/须） ---
    q1, q2, q3 = np.percentile(x, [25, 50, 75])
    iqr = q3 - q1
    lw = np.min(x[x >= q1 - 1.5 * iqr]) if np.any(x >= q1 - 1.5 * iqr) else x.min()
    uw = np.max(x[x <= q3 + 1.5 * iqr]) if np.any(x <= q3 + 1.5 * iqr) else x.max()
    box_h = 0.28
    ax.add_patch(plt.Rectangle((q1, y0 - box_h/2), q3 - q1, box_h,
                               facecolor=base_color, edgecolor='none', alpha=0.55))
    ax.plot([q2, q2], [y0 - box_h/2, y0 + box_h/2], color='black', lw=2)
    ax.plot([lw, q1], [y0, y0], color='black', lw=1.2)
    ax.plot([q3, uw], [y0, y0], color='black', lw=1.2)

    # --- 抖动散点（雨滴） ---
    if k > 0:
        jitter_y = (y0 - 0.55) + np.random.uniform(-0.12, 0.12, size=k)
        if color_by_variant and is_snp is not None and is_snp.size == x.size:
            idx_snp   = idx[is_snp[idx]]
            idx_indel = idx[~is_snp[idx]]
            if idx_snp.size > 0:
                ax.scatter(x[idx_snp], jitter_y[:idx_snp.size], s=6, alpha=0.35, color='#000000', linewidth=0)
            if idx_indel.size > 0:
                ax.scatter(x[idx_indel], jitter_y[idx_snp.size:idx_snp.size+idx_indel.size], s=6, alpha=0.35, color='#CC79A7', linewidth=0)
        else:
            ax.scatter(x[idx], jitter_y, s=6, alpha=0.35, color=base_color, linewidth=0)
        # 在右上角添加 SNP / InDel 图例（仅 color_by_variant=True 时）
        if color_by_variant and is_snp is not None and is_snp.size == x.size:
            from matplotlib.lines import Line2D
            handles = [
                Line2D([0], [0], marker='o', linestyle='None', markersize=6,
                       markerfacecolor='#000000', markeredgewidth=0, label='SNP'),
                Line2D([0], [0], marker='o', linestyle='None', markersize=6,
                       markerfacecolor='#CC79A7', markeredgewidth=0, label='InDel'),
            ]
            ax.legend(handles=handles, loc='upper right', frameon=True, framealpha=0.9, borderpad=0.6)

    # --- 仅基于 P_BH(q) 的阈值计数与可视化（<0.05, <0.01, <0.001） ---
    v05  = -np.log10(0.05)   # ≈ 1.30103
    v01  = -np.log10(0.01)   # = 2
    v001 = -np.log10(0.001)  # = 3

    cnt05  = int(np.sum(q < 0.05))
    cnt01  = int(np.sum(q < 0.01))
    cnt001 = int(np.sum(q < 0.001))
    pct05  = cnt05  / n if n else 0.0
    pct01  = cnt01  / n if n else 0.0
    pct001 = cnt001 / n if n else 0.0

    for v, _label, c, _a in [
        (v05,  None, 0.7, None),
        (v01,  None, 0.6, None),
        (v001, None, 0.5, None),
    ]:
        ax.axvline(v, ls="--", lw=1.0, color="gray", alpha=c)

    # --- 右侧统计框（仅给出基于 P_BH 的信息） ---
    q_sig_005  = cnt05
    q_sig_001  = cnt01
    q_sig_0001 = cnt001
    q_pct_005  = pct05
    q_pct_001  = pct01
    q_pct_0001 = pct001

    text = [
        f"Summary of -log10({bh_col})",
        f"N = {n:,}",
        f"mean   = {x.mean():.3f}",
        f"median = {np.median(x):.3f}",
        f"IQR    = [{q1:.3f}, {q3:.3f}] (w {iqr:.3f})",
        "",
        f"{bh_col} thresholds (q):",
        f"  {bh_col} < 0.05  : {cnt05:,} ({pct05:.3%})",
        f"  {bh_col} < 0.01  : {cnt01:,} ({pct01:.3%})",
        f"  {bh_col} < 0.001 : {cnt001:,} ({pct001:.3%})",
    ]
    inset = ax.inset_axes([1.02, 0.05, 0.48, 0.9], transform=ax.transAxes)
    inset.axis("off")
    inset.text(0, 1, "\n".join(text), ha="left", va="top",
               fontsize=9, family="monospace",
               bbox=dict(facecolor="white", edgecolor="gray", alpha=0.85))

    # --- 样式美化 ---
    ax.set_xlabel(f"-log10({bh_col})")
    ax.set_yticks([])
    ax.set_ylabel("Density")
    ax.set_title(title, fontsize=13, fontweight="bold", pad=10)
    ax.grid(axis="x", linestyle=":", alpha=0.9)
    for spine in ["right", "top", "left"]:
        ax.spines[spine].set_visible(False)

    # 保存图片
    if output_path is None:
        base = os.path.splitext(os.path.basename(missing_path))[0]
        output_path = os.path.join(os.path.dirname(missing_path), f"{base}.raincloud.{bh_col}.png")
    if created_fig:
        plt.tight_layout(rect=(0, 0, 0.82, 1))
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(ax.figure)

    return {
        "output_path": output_path,
        "N": n,
        "q_le_0.05": q_sig_005, "q_le_0.01": q_sig_001, "q_le_0.001": q_sig_0001,
        "q_pct_0.05": q_pct_005, "q_pct_0.01": q_pct_001, "q_pct_0.001": q_pct_0001,
        "bh_col": bh_col,
        "missing_path": missing_path,
        "manifest_path": mani.get("manifest_path", manifest_path),
    }



def remove_variants_by_fdr_from_manifest(
    manifest_path,
    fdr_threshold: float = 0.05,
    plink2_path: str = "/home/b/b37974/plink2",
    threads: int = 16,
    out_prefix: Optional[str] = None,
) -> Dict[str, Any]:
    """
    基于 manifest 的 .missing 表和 BH 校正列，提取 q < fdr_threshold 的变体列表，并用 plink2 从 bfile 中去除这些变体。

    参数
    ----
    manifest_path : str 或 dict
        run_test_missing 返回的 manifest 字典或其 JSON 路径。
    fdr_threshold : float, default=0.05
        FDR (BH) 阈值；提取 {bh_column} < fdr_threshold 的 SNP 列表。
    plink2_path : str, default="/home/b/b37974/plink2"
        plink2 可执行文件路径。
    threads : int, default=16
        plink2 使用的线程数。
    out_prefix : Optional[str]
        plink2 输出前缀（不含扩展名）。默认使用 "{bed_prefix}.rm_q_lt_{fdr_threshold:g}"。

    返回
    ----
    dict:
        {
          "snplist_path": <写出的变体列表文件路径>,
          "exclude_n": <被排除的变体个数>,
          "plink2_cmd": <运行的命令列表>,
          "out_prefix": <plink2 输出前缀>,
          "out_bed": <out_prefix + ".bed">,
          "out_bim": <out_prefix + ".bim">,
          "out_fam": <out_prefix + ".fam">,
          "log_path": <运行日志路径>
        }

    示例
    ----
    >>> res = remove_variants_by_fdr_from_manifest(
    ...     manifest_path="cohort.missing_bias.manifest.json",
    ...     fdr_threshold=0.05,
    ...     plink2_path="/home/b/b37974/plink2",
    ...     threads=16
    ... )

    调试提示
    --------
    - 该函数会将 CMD、标准输出与错误输出全量写入 `<out_prefix>.remove_q.log`
    - 若 `snps` 为空，仍会生成新的 bfile（不带 --exclude），日志中会提示 exclude_n=0
    - 若 plink2_path 不存在，仍尝试从 PATH 搜索 plink2，建议提前确认路径或在日志中查看错误输出
    """
    # 1) 读取 manifest
    if isinstance(manifest_path, dict):
        mani = manifest_path
    else:
        with open(manifest_path, "r", encoding="utf-8") as fp:
            mani = json.load(fp)

    missing_path = mani.get("missing_path")
    bh_col = mani.get("bh_column", "P_BH")
    bed_prefix = mani.get("bed_prefix")
    if not missing_path or not os.path.isfile(missing_path):
        raise FileNotFoundError(f"manifest['missing_path'] 不存在: {missing_path}")
    if not bed_prefix:
        raise KeyError("manifest 中缺少 'bed_prefix' 字段，无法定位 bfile。")

    # 2) 读取 .missing 表，基于 {bh_col} < fdr_threshold 提取 SNP/ID 列
    df = _safe_read_missing_table(missing_path)
    if bh_col not in df.columns:
        raise KeyError(f"在 {missing_path} 中未找到列: {bh_col}")
    var_col = "SNP" if "SNP" in df.columns else ("ID" if "ID" in df.columns else None)
    if var_col is None:
        raise KeyError(f"在 {missing_path} 中未找到用于变体 ID 的列: 'SNP' 或 'ID'")

    q = pd.to_numeric(df[bh_col], errors="coerce")
    mask = q.notna() & np.isfinite(q) & (q < fdr_threshold)
    snps = df.loc[mask, var_col].astype(str)

    # 3) 写出一列的变体列表（无表头）
    base = os.path.basename(bed_prefix)
    snplist_path = os.path.join(
        os.getcwd(),
        f"{base}.q_lt_{fdr_threshold:g}.snplist"
    )
    snps.to_csv(snplist_path, index=False, header=False)

    # 4) 调用 plink2 进行排除
    if out_prefix is None:
        out_prefix = os.path.join(os.getcwd(), f"{base}.rm_q_lt_{fdr_threshold:g}")
    out_dir = os.path.dirname(out_prefix) or "."
    os.makedirs(out_dir, exist_ok=True)

    cmd = [
        plink2_path,
        "--bfile", bed_prefix,
        "--make-bed",
        "--threads", str(threads),
        "--out", out_prefix,
    ]
    if len(snps) > 0:
        cmd[4:4] = ["--exclude", snplist_path]  # 在 --make-bed 之前插入 --exclude

    log_path = f"{out_prefix}.remove_q.log"
    with open(log_path, "w", encoding="utf-8") as log_fp:
        log_fp.write("[CMD]\n" + " ".join(shlex.quote(x) for x in cmd) + "\n\n")
        log_fp.flush()
        _log(log_fp, f"工作目录（cwd）: {os.getcwd()}")
        _log(log_fp, f"plink2 路径: {plink2_path}（存在: {os.path.isfile(plink2_path)}）")
        _log(log_fp, f"bed_prefix: {bed_prefix}")
        _log(log_fp, f"FDR 阈值: {fdr_threshold:g}；待排除变体数: {int(len(snps))}")
        if len(snps) > 0:
            # 打印前几条变体，便于快速核对
            preview = snps.head(min(5, len(snps))).tolist()
            _log(log_fp, f"变体列示例(最多5条): {preview}")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True)
        try:
            for line in proc.stdout:
                log_fp.write(line)
            ret = proc.wait()
            if ret != 0:
                log_fp.write(f"\n[ERROR] plink2 运行失败，返回码: {ret}\n")
                raise RuntimeError(f"plink2 执行失败，详见日志: {log_path}")
        finally:
            if proc.stdout:
                proc.stdout.close()

    return {
        "snplist_path": snplist_path,
        "exclude_n": int(len(snps)),
        "plink2_cmd": cmd,
        "out_prefix": out_prefix,
        "out_bed": f"{out_prefix}.bed",
        "out_bim": f"{out_prefix}.bim",
        "out_fam": f"{out_prefix}.fam",
        "log_path": log_path,
    }
