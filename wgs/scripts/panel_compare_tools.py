"""
模块名称：panel_compare_tools
=================================

【概述】
本模块提供了一系列针对 ToMMo 参考数据库进行变异位点交叉比对与统计可视化的工具函数，适用于全基因组测序（WGS）或芯片数据的质控与后续关联分析准备。

【功能组成】
1. run_plink2_variant_qc_with_tommo
   - 将 `variant_qc_summary` 与 ToMMo VCF 进行交叉比对，生成附带 ToMMo 信息的扩展表。
   - 输出包含 IN_TOMMO、TOMMO_AAF、TOMMO_FILTER 三列。

2. plot_tommo_panel_compare_pdf
   - 基于 `variant_qc_with_tommo.tsv` 绘制对比 PDF（共 5 页）。
   - 包括统计表、散点图（PASS/Non-PASS 与 SNP/InDel）、直方图等。

3. build_grouped_variant_tables
   - 按 CTRL_MAF 将变异划分为 Rare / LowFreq / Common 三类，并在各类内细分 a/b/c 三个亚组。
   - 针对 c_in_pass 组计算 DIFF 与 ROBUST_Z，输出分组化的 TSV.GZ 文件及 manifest.json。

4. summarize_c_in_pass_thresholds
   - 基于 manifest.json，扫描 c_in_pass 组，计算不同 |ROBUST_Z| 阈值下的 count 与 MSE。
   - 输出各分组的 summary TSV，并将结果更新写回 manifest.json。

5. plot_c_in_pass_threshold_tradeoff
   - 读取 summary TSV，在 PDF 中绘制 trade-off 曲线，并通过 KneeLocator 自动识别拐点。
   - 同时输出拐点阈值下的散点分布图和对应变体列表 TSV，更新至 manifest.json。

6. summarize_variants_filter_from_manifest
   - 结合 manifest.json 和输入表，生成带有 FILTER_STAT 标签的 summary 表。
   - 支持最终过滤状态的快速判定（Stat_1 ~ Stat_4），并输出分组统计结果。

【适用场景】
- 大规模变异数据的 ToMMo 参考面板比对与质量控制；
- 根据群体频率差异与稳健统计指标（ROBUST_Z）进行多阶段过滤；
- GWAS、pQTL 等关联研究前的数据预处理。

【实现特点】
- 全面采用 **分块读取** 与 **并行处理**，支持千万至上亿行规模的数据。
- 临时文件与中间结果均通过日志记录，便于排错与调试。
- 输出文件结构清晰：表格（TSV/TSV.GZ）、可视化（PDF）、摘要统计（JSON/TSV）。

作者: ZHAO TIE
"""
import subprocess
import os
import tempfile
import csv
import math
import uuid
from datetime import datetime
from typing import Optional
from typing import Dict
import concurrent.futures
import threading
import pandas as pd
import numpy as np
from itertools import islice
from collections import defaultdict
from typing import Dict, Tuple, Iterable
import shutil
import time

# ==== Unified MAF grouping (single source of truth) ====
# Definition (do not change without auditing all downstream consumers):
#   rare    : MAF < 0.01
#   lowfreq : 0.01 <= MAF <= 0.05
#   common  : MAF > 0.05
# NOTE: We intentionally do NOT clamp values to [0,1] here, to avoid silently
# changing legacy behavior; we only coerce to numeric with NaN on failures.

def assign_maf_group(maf_series):
    """
    Assign MAF groups with consistent boundaries, returning a pandas Series
    with values in {"rare", "lowfreq", "common"} or <NA> when not classifiable.
    Boundaries:
      - 'rare'    : value < 0.01
      - 'lowfreq' : 0.01 <= value <= 0.05
      - 'common'  : value > 0.05
    Only converts to numeric (errors='coerce'); does not clip negative/>1 values.
    """
    import pandas as _pd
    s = _pd.to_numeric(maf_series, errors='coerce')
    out = _pd.Series(_pd.NA, index=s.index, dtype='object')
    out = out.mask(s < 0.01, 'rare')
    out = out.mask((s >= 0.01) & (s <= 0.05), 'lowfreq')
    out = out.mask(s > 0.05, 'common')
    return out


def validate_grouping_counts(maf_series):
    """Lightweight diagnostic: return counts per group for quick sanity checks."""
    import pandas as _pd
    s = _pd.to_numeric(maf_series, errors='coerce')
    g = assign_maf_group(s)
    total = int(s.notna().sum())
    rare = int((g == 'rare').sum())
    lowf = int((g == 'lowfreq').sum())
    comm = int((g == 'common').sum())
    return {
        'total_non_nan': total,
        'rare': rare,
        'lowfreq': lowf,
        'common': comm,
        'ungrouped_non_nan': total - (rare + lowf + comm),
    }



def run_plink2_variant_qc_with_tommo(
    variant_qc_summary: str,
    tommo_vcf_path: str,
    output_path: Optional[str] = None,
    bcftools_path: str = "bcftools",
    threads: int = 8,
    chunk_size: int = 500_000,
    max_workers: Optional[int] = None,
    regions_chunk_lines: int = 50_000,
    keep_tmp: bool = False,
) -> str:
    """
    模块函数：run_plink2_variant_qc_with_tommo
    ========================================
    【功能】
    - 针对 *非常长* 的 `variant_qc_summary`（列含 VARIANT_ID=CHROM:POS:REF:ALT），
      先根据 CHROM:POS 生成 bcftools 可读的 regions 文件（按染色体拆分并去重）；
    - 并行调用 `bcftools query -R` 从 ToMMo VCF 中抽取位点信息，
      使用 per-allele 展开格式确保按 REF/ALT 精确匹配；
    - 以**流式分块**方式读取原表并合并 3 列：
        * IN_TOMMO: bool（是否存在完全匹配的 CHROM:POS:REF:ALT）
        * TOMMO_AAF: float（ToMMo 的 INFO/AF，对应 ALT 等位）
        * TOMMO_FILTER: str（该记录的 FILTER）
    - 最终写出 TSV（默认后缀 `.variant_qc_with_tommo.tsv`）。

    【重要实现要点】
    - regions 文件采用两列 1-based 的 `CHROM\tPOS`（**不要**混用 BED 坐标）。
    - 使用 `bcftools query` 的 per-allele 展开：
        格式串：`%CHROM\t%POS[\t%REF\t%ALT\t%FILTER\t%INFO/AF]\n`
      方括号 `[]` 会对 ALT 逐等位展开，保证 REF/ALT 一一对应。
    - 合并阶段不会把整张 ToMMo 或整张 summary 全部载入内存：
        * 第一步仅生成每条染色体的去重位置列表（磁盘中转 + `sort -u` 去重）。
        * 第二步对每条染色体独立 `bcftools query` 并写出中间映射表。
        * 第三步**分块**读取 summary，分组到染色体后仅按需要的键子集
          从映射表中"按需加载"对应的少量行，映射完成即丢弃。
    - 染色体名需与 VCF 保持完全一致（例如 `chr20` ≠ `20`）。

    参数
    ----
    variant_qc_summary : str
        `run_plink2_variant_qc` 产出的 `*.variant_qc_summary.tsv` 路径。
    tommo_vcf_path : str
        ToMMo 的 bgzip 压缩并建立 `.tbi` 索引的 VCF 路径。
    output_path : Optional[str]
        输出路径；默认与输入同名，后缀改为 `.variant_qc_with_tommo.tsv`。
    bcftools_path : str
        `bcftools` 可执行程序路径（默认走环境中的 `bcftools`）。
    threads : int
        传给 bcftools 的线程数（读写解压线程）。
    chunk_size : int
        读取 `variant_qc_summary` 的 pandas 分块大小。
    max_workers : Optional[int]
        并行运行 `bcftools query` 的最大并发数；默认等于 `min(4, 可用CPU)`。
    regions_chunk_lines : int
        将每条染色体的 region 列表再按行数切分为小块，逐块执行 bcftools（默认 1,000,000 行/块），便于打印进度与控制单次调用规模。
    keep_tmp : bool
        是否保留临时目录以便排错。

    返回
    ----
    str
        生成的 `variant_qc_summary_with_tommo` 的文件路径。
    """
    import sys
    import shlex
    import tempfile
    import pandas as pd
    import numpy as np
    import concurrent.futures
    from collections import OrderedDict

    def _progress(msg: str):
        print(f"[run_plink2_variant_qc_with_tommo] {msg}", file=sys.stderr, flush=True)

    def _parse_vid(vid: str) -> Tuple[str, int, str, str]:
        # 期望形如 CHROM:POS:REF:ALT
        parts = vid.split(":", 3)
        if len(parts) != 4:
            raise ValueError(f"VARIANT_ID 不是 CHROM:POS:REF:ALT 格式: {vid}")
        chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
        return chrom, int(pos), ref, alt

    def _run_cmd(cmd_list: Iterable[str], stdout_path: str = None) -> int: # type: ignore
        if stdout_path is None:
            proc = subprocess.run(cmd_list, check=False)
        else:
            with open(stdout_path, "wb") as fo:
                proc = subprocess.run(cmd_list, check=False, stdout=fo) # type: ignore
        return proc.returncode

    # 1) 临时工作目录 & 输出路径
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    workdir = os.path.abspath(os.path.join(os.getcwd(), f"tommo_merge_{ts}_{uuid.uuid4().hex[:8]}"))
    os.makedirs(workdir, exist_ok=True)
    if output_path is None:
        base, _ = os.path.splitext(variant_qc_summary)
        output_path = base + ".variant_qc_with_tommo.tsv"
    _progress(f"临时目录: {workdir}")

    # 进度日志文件（便于在 Nextflow /.command.log 之外进行追踪）
    progress_log_path = os.path.join(workdir, "progress.log")
    def _log(msg: str):
        _progress(msg)
        try:
            with open(progress_log_path, 'a') as pf:
                pf.write(msg + "\n")
        except Exception:
            pass

    # 2) 第一步：扫描 summary，按染色体写出 CHROM\tPOS 的 region 原始文件
    #    为避免高内存，占位文件按染色体拆分，并最终使用 `sort -u` 去重。
    region_tmp_files: Dict[str, str] = {}
    # 以分块方式读取，只需要 VARIANT_ID 一列
    reader = pd.read_csv(
        variant_qc_summary, sep='\t', usecols=["VARIANT_ID"], dtype={"VARIANT_ID": "string"},
        chunksize=chunk_size, engine='c'
    )
    total_rows = 0
    for chunk in reader:
        total_rows += len(chunk)
        # 向量化解析 VARIANT_ID -> CHROM, POS
        sp = chunk["VARIANT_ID"].str.split(":", n=3, expand=True)
        sp.columns = ["CHROM", "POS", "REF", "ALT"]
        sp = sp[["CHROM", "POS"]]
        # 逐染色体写入（允许重复，后续 sort -u 去重）
        for chrom, sub in sp.groupby("CHROM"):
            path = region_tmp_files.get(chrom)
            if path is None:
                path = os.path.join(workdir, f"regions.{chrom}.tsv")
                region_tmp_files[chrom] = path
            # 只写两列，POS 按原始字符串即可（1-based）
            sub[["CHROM", "POS"]].to_csv(
                path, sep='\t', header=False, index=False, mode='a'
            )
    _progress(f"已扫描 {total_rows:,} 行，生成 {len(region_tmp_files)} 个染色体 region 文件（未去重）")

    # 3) 对每个染色体的 region 文件做 sort -u 去重，得到 .uniq 文件
    uniq_region_files: Dict[str, str] = {}
    for chrom, raw_path in region_tmp_files.items():
        uniq_path = os.path.join(workdir, f"regions.{chrom}.uniq.tsv")
        # 使用系统 sort -u，按 (CHROM, POS) 去重并保证数值排序
        # 注：LC_ALL=C 可显著加速排序
        cmd = [
            "bash", "-lc",
            f"LC_ALL=C sort -u -t$'\t' -k1,1 -k2,2n {shlex.quote(raw_path)} > {shlex.quote(uniq_path)}"
        ]
        ret = _run_cmd(cmd)
        if ret != 0:
            raise RuntimeError(f"sort -u 去重失败: {raw_path}")
        uniq_region_files[chrom] = uniq_path
    _progress("已完成每条染色体的 region 去重")

    # 4) 并行运行 bcftools query 生成每条染色体的等位基因级映射表
    #    输出格式：CHROM POS REF ALT FILTER AF（每行一条 ALT 等位）
    mapping_files: Dict[str, str] = {}
    # 注意：bcftools `[]` 迭代的是 FORMAT/样本字段，不是 ALT/INFO 数组；
    # 因此这里打印 ALT 与 INFO/AF 的逗号分隔列表，后续在 Python 侧按等位一一展开。
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AF\n"

    def _run_bcftools_for_chrom(chrom: str) -> Tuple[str, str]:
        region_file = uniq_region_files[chrom]
        out_map = os.path.join(workdir, f"tommo.map.{chrom}.tsv")

        # 将 region 文件按行数切分为小块，便于打印进度并缩短单次 bcftools 调用
        parts_dir = os.path.join(workdir, f"regions.{chrom}.parts")
        os.makedirs(parts_dir, exist_ok=True)
        part_paths = []
        with open(region_file, 'r') as fin:
            part_idx = 0
            buf = []
            for ln, line in enumerate(fin, start=1):
                buf.append(line)
                if ln % regions_chunk_lines == 0:
                    part_idx += 1
                    p = os.path.join(parts_dir, f"part_{part_idx:05d}.tsv")
                    with open(p, 'w') as fout:
                        fout.writelines(buf)
                    part_paths.append(p)
                    buf.clear()
            if buf:
                part_idx += 1
                p = os.path.join(parts_dir, f"part_{part_idx:05d}.tsv")
                with open(p, 'w') as fout:
                    fout.writelines(buf)
                part_paths.append(p)
        _log(f"[{chrom}] region 切分为 {len(part_paths)} 块（每块 ≤ {regions_chunk_lines:,} 行）")

        # 清空/创建输出文件
        open(out_map, 'wb').close()

        # 注意：部分环境中的 `bcftools query` 不支持 --threads。
        # 线程>1：`view --threads` 管道 + query；否则：直接 query。
        t0 = time.time()
        for i, p in enumerate(part_paths, start=1):
            _log(f"[{chrom}] 处理 chunk {i}/{len(part_paths)}: {os.path.basename(p)}")
            if isinstance(threads, int) and threads > 1:
                fmt_q = shlex.quote(fmt)
                cmdline = (
                    f"{shlex.quote(bcftools_path)} view --threads {threads} -R {shlex.quote(p)} "
                    f"-Ou {shlex.quote(tommo_vcf_path)} | "
                    f"{shlex.quote(bcftools_path)} query -f {fmt_q} >> {shlex.quote(out_map)}"
                )
                cmd = ["bash", "-lc", cmdline]
                ret = _run_cmd(cmd)
            else:
                cmd = [
                    bcftools_path, "query",
                    "-R", p,
                    "-f", fmt,
                    tommo_vcf_path,
                ]
                # 以追加方式写入
                with open(out_map, 'ab') as fout:
                    proc = subprocess.run(cmd, check=False, stdout=fout)
                    ret = proc.returncode
            if ret != 0:
                raise RuntimeError(f"bcftools query 失败: 染色体 {chrom} 块 {i}/{len(part_paths)}")
            # 进度条打印（文本版，适配日志文件）：
            done = i
            total = len(part_paths)
            pct = done / total if total else 1.0
            bar_len = 28
            filled = int(bar_len * pct)
            bar = '█' * filled + '-' * (bar_len - filled)
            elapsed = time.time() - t0
            rate = done / max(elapsed, 1e-6)
            eta = int((total - done) / max(rate, 1e-9))
            eta_m, eta_s = divmod(eta, 60)
            _log(f"[{chrom}] [{bar}] {done}/{total} ({pct*100:.1f}%) 速率 {rate:.2f} chunk/s  预计剩余 {eta_m:02d}:{eta_s:02d}")

        return chrom, out_map

    if max_workers is None:
        try:
            import multiprocessing as _mp
            max_workers = max(1, min(4, _mp.cpu_count()))
        except Exception:
            max_workers = 2

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = [ex.submit(_run_bcftools_for_chrom, chrom) for chrom in uniq_region_files.keys()]
        for fut in concurrent.futures.as_completed(futs):
            chrom, out_map = fut.result()
            mapping_files[chrom] = out_map
            _progress(f"bcftools 完成: {chrom}")
    _progress(f"已生成 {len(mapping_files)} 条染色体映射表")

    # 诊断：若所有映射表均为空，提示可能的染色体命名不一致问题
    empty_maps = 0
    for _c, _path in mapping_files.items():
        try:
            _size = os.path.getsize(_path)
        except OSError:
            _size = 0
        if _size == 0:
            empty_maps += 1
    if empty_maps == len(mapping_files):
        _progress("[warn] 所有染色体映射表均为空。请检查：1) VARIANT_ID 中的染色体前缀是否与 ToMMo VCF 一致（例如 chr1 vs 1）；2) -R 区域文件是否为 1-based 两列格式；3) VCF 是否有索引 .tbi 且路径正确。")

    # 5) 合并阶段：按需加载映射子集，分块写出结果
    #    - 输出列为原表所有列 + [IN_TOMMO, TOMMO_AAF, TOMMO_FILTER]
    #    - TOMMO_AAF 用 float，不能解析时为 NaN；IN_TOMMO 为 True/False。

    # 获取原表列名（防止顺序变化）
    with open(variant_qc_summary, 'r') as fi:
        header_line = fi.readline().rstrip('\n')
    base_cols = header_line.split('\t')
    out_cols = base_cols + ["IN_TOMMO", "TOMMO_AAF", "TOMMO_FILTER"]

    # 输出文件：先写表头
    with open(output_path, 'w') as fo:
        fo.write('\t'.join(out_cols) + '\n')

    def _load_mapping_subset_for_chrom(chrom: str, needed_keys: set) -> Dict[str, Tuple[str, str]]:
        """仅加载该染色体映射表中 *需要* 的键，返回 {VID: (FILTER, AF)}。
        支持 ALT/AF 逗号分隔的多等位展开。
        """
        out: Dict[str, Tuple[str, str]] = {}
        map_path = mapping_files.get(chrom)
        if (map_path is None) or (not os.path.exists(map_path)):
            return out
        with open(map_path, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                if not line:
                    continue
                # 兼容某些 shell 传递下，fmt 未被转义为真实制表符导致输出含字面量 "\\t" 的情况
                if "\\t" in line and "\t" not in line:
                    line = line.replace("\\t", "\t")
                cols = line.split('\t')
                # 期望：固定6列：CHROM POS REF ALT(s) FILTER AF(s)
                if len(cols) < 6:
                    continue
                c, p, r, alts_s, flt, afs_s = cols[0], cols[1], cols[2], cols[3], cols[4], cols[5]
                # ALT 与 AF 都可能是逗号分隔的多等位数组；需要一一配对
                alts = alts_s.split(",") if alts_s != "." else []
                afs = afs_s.split(",") if afs_s not in (".", "") else []
                # 对齐长度：若 AF 缺失或长度与 ALT 不同，仅在可配对位置输出
                n = min(len(alts), len(afs)) if afs else len(alts)
                for j in range(n):
                    a = alts[j]
                    af = afs[j] if j < len(afs) else "."
                    vid = f"{c}:{p}:{r}:{a}"
                    if (not needed_keys) or (vid in needed_keys):
                        out[vid] = (flt, af)
        # 轻量诊断：如需要可打印命中数（仅当存在需要键时）
        if needed_keys:
            hit_cnt = sum(1 for k in needed_keys if k in out)
            _log(f"[diag] {chrom}: 映射子集载入 {len(out):,} 条，命中 {hit_cnt:,} / {len(needed_keys):,}")
        return out

    # 分块读取原表，合并并写出
    reader2 = pd.read_csv(
        variant_qc_summary, sep='\t', dtype="string", chunksize=chunk_size, engine='c'
    )

    processed = 0
    for chunk in reader2:
        processed += len(chunk)
        # 默认值
        chunk["IN_TOMMO"] = False
        chunk["TOMMO_AAF"] = pd.Series([pd.NA] * len(chunk), dtype="string")
        chunk["TOMMO_FILTER"] = pd.Series([pd.NA] * len(chunk), dtype="string")

        # 解析 VARIANT_ID -> CHROM
        sp = chunk["VARIANT_ID"].str.split(":", n=3, expand=True)
        sp.columns = ["CHROM", "POS", "REF", "ALT"]
        chunk["__CHROM__"] = sp["CHROM"].astype("string")

        # 按染色体处理，减少一次读取的映射量
        for chrom, idx in chunk.groupby("__CHROM__").groups.items():
            sub = chunk.loc[idx]
            need_keys = set(sub["VARIANT_ID"].tolist())
            mapping = _load_mapping_subset_for_chrom(chrom, need_keys)
            if not mapping:
                continue
            # 命中掩码
            hit_mask = sub["VARIANT_ID"].isin(mapping.keys())
            if not hit_mask.any():
                continue
            vids_hit = sub.loc[hit_mask, "VARIANT_ID"]
            # 单独构造映射字典以便矢量化 map
            to_filter = {k: v[0] for k, v in mapping.items()}
            to_af = {k: v[1] for k, v in mapping.items()}

            chunk.loc[idx[hit_mask], "IN_TOMMO"] = True
            # TOMMO_FILTER 直接映射到字符串
            chunk.loc[idx[hit_mask], "TOMMO_FILTER"] = vids_hit.map(to_filter).astype("string").values
            # TOMMO_AAF 先作为字符串接收，再安全转为 float（'.' -> NaN）
            af_str = vids_hit.map(to_af).astype("string")
            # 将 '.' 或无法解析的转为 NaN
            af_num = pd.to_numeric(af_str, errors='coerce')
            chunk.loc[idx[hit_mask], "TOMMO_AAF"] = af_num.astype("Float32").astype("string")

        # 移除临时列
        chunk = chunk.drop(columns=["__CHROM__"])

        # 写出（使用字符串 dtype，缺失值显示为 'nan' 以与既有代码风格一致）
        chunk.to_csv(
            output_path, sep='\t', header=False, index=False, mode='a', na_rep='nan'
        )
        _progress(f"已处理 {processed:,} 行")

    # 6) 清理或保留临时目录
    if keep_tmp:
        _log(f"保留临时目录: {workdir}")
    else:
        try:
            shutil.rmtree(workdir)
            _log(f"已清理临时目录: {workdir}")
        except Exception as e:
            _log(f"[warn] 清理临时目录失败: {workdir} -> {e}")

    _log(f"完成。输出: {output_path}")
    return output_path




import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator, FuncFormatter

plt.rcParams['pdf.compression'] = 0

# Publication-oriented defaults
plt.rcParams.update({
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans'],
    'font.size': 10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.titlesize': 12,
    'axes.linewidth': 0.8,
    'grid.linewidth': 0.4,
    'grid.color': '#CCCCCC',
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})

plt.style.use('default')  # 使用默认样式

def plot_tommo_panel_compare_pdf(
    variant_qc_with_tommo: str,
    output_pdf: Optional[str] = None,
    max_points: Optional[int] = None,
    png_dpi: int = 600,
    page23_figsize: tuple = (14, 5),
    base_fontsize: int = 10,
    theme: str = 'okabe_ito',
    page4_figsize: tuple = (14, 5),
    page4_wspace: float = 0.30,
    snp_color: Optional[str] = None,
    indel_color: Optional[str] = None,
):
    """
    函数名称：plot_tommo_panel_compare_pdf
    =====================================
    【功能】
    - 输入 `run_plink2_variant_qc_with_tommo` 生成的结果表（*.variant_qc_with_tommo.tsv），
      输出一个包含 **5 页** 的 PDF 文件，每页含 3 个分图（Rare/LowFreq/Common 三组）。

    【页面内容】
    1. **第 1 页**：三组统计表，列出
        - Total Variants：该分组的变体总数
        - In ToMMo：在 ToMMo 中存在的变体数
        - Pass ToMMo：ToMMo FILTER=PASS 的变体数
        - Percent：相对于组内总数的百分比（1 位小数），Count 使用千分位逗号格式。
    2. **第 2 页**：三组散点图，X=TOMMO_AAF，Y=CTRL_AAF。颜色区分 ToMMo FILTER 是否 PASS。
    3. **第 3 页**：三组散点图，X=TOMMO_AAF，Y=CTRL_AAF。颜色区分变异类型（SNP vs InDel），不筛选 PASS。
    4. **第 4 页**：三组散点图，X=TOMMO_AAF，Y=CTRL_AAF。仅保留 ToMMo FILTER=PASS 的点，颜色区分 SNP vs InDel。
    5. **第 5 页**：三组直方图，X=CTRL_AAF-TOMMO_AAF 的分布，仅 PASS 变体；标题中显示 $\\mu \\pm \\sigma$。

    【参数说明】
    - variant_qc_with_tommo : str
        输入表路径。
    - output_pdf : Optional[str]
        输出 PDF 路径（默认与输入同名，后缀改为 `.panel_compare.pdf`）。
    - max_points : Optional[int]
        每组最大采样点数（None 表示不抽样）。
    - png_dpi : int
        第 2~4 页散点图先渲染为 PNG 再插入 PDF，此参数为分辨率。
    - page23_figsize : tuple
        第 2~3 页整体画布尺寸（英寸）。
    - page4_figsize : tuple
        第 4 页整体画布尺寸。
    - page4_wspace : float
        第 4 页直方图之间的水平间距。
    - snp_color / indel_color : Optional[str]
        第 3~4 页 SNP 与 InDel 的点颜色；默认使用 Okabe–Ito 或 Matplotlib 配色。

    【实现要点】
    - 自动将必要列转换为数值类型，丢弃 NaN。
    - 基于 CTRL_MAF 将变体划分为三组：Rare (<0.01)、LowFreq (0.01~0.05)、Common (>0.05)。
    - 为提升性能，第 2~4 页的散点图均先渲染为 PNG 再插入 PDF，避免 PDF 过大或渲染缓慢。
    - 直方图保持矢量绘制，便于放大。

    返回
    ----
    str
        输出 PDF 文件路径。
    """
    import os
    import numpy as np
    import pandas as pd
    import tempfile

    # --- Styling knobs (publication-ready) ---
    plt.rcParams['font.size'] = base_fontsize
    if theme == 'okabe_ito':
        # Okabe–Ito colorblind-safe
        pass_color = '#0072B2'     # blue
        nonpass_color = '#ff6f01'  # vermillion
        refline_color = '#CC0000'
        # SNP/InDel 颜色（与 PASS/Non-PASS 区分开）
        if snp_color is None:
            snp_color = '#000000'   # black (max contrast on white)
        if indel_color is None:
            indel_color = '#CC79A7' # magenta (Okabe–Ito)
    else:
        pass_color = '#1f77b4'
        nonpass_color = '#ff7f0e'
        refline_color = 'red'
        if snp_color is None:
            snp_color = '#000000'   # black (max contrast)
        if indel_color is None:
            indel_color = '#9467bd'

    # marker and line sizes
    ms_pass = 16
    ms_nonpass = 18
    alpha_pass = 0.7
    alpha_nonpass = 0.8
    refline_ls = (0, (4, 2))  # dashed pattern
    refline_lw = 1.0

    # 逐块累积，避免一次性占用大量内存
    needed_cols = [
        'VARIANT_ID', 'CTRL_MAF', 'CTRL_AAF', 'TOMMO_AAF', 'IN_TOMMO', 'TOMMO_FILTER'
    ]
    chunk_size = 500_000  # 可根据内存情况调整
    dfs = []
    for chunk in pd.read_csv(variant_qc_with_tommo, sep='\t', usecols=needed_cols,
                             dtype={'VARIANT_ID': 'string',
                                    'CTRL_MAF': 'float64',
                                    'CTRL_AAF': 'float64',
                                    'TOMMO_AAF': 'float64',
                                    'IN_TOMMO': 'object',
                                    'TOMMO_FILTER': 'string'},
                             chunksize=chunk_size, engine='c'):
        if 'IN_TOMMO' in chunk.columns:
            chunk['IN_TOMMO'] = chunk['IN_TOMMO'].map(
                lambda x: True if x in (True, 1, '1', 'True', 'TRUE') else (False if x in (False, 0, '0', 'False', 'FALSE') else pd.NA)
            )
        dfs.append(chunk)
    df = pd.concat(dfs, ignore_index=True, copy=False)


    # 统一列名期望，并进行类型转换
    needed_cols = [
        'VARIANT_ID', 'CTRL_MAF', 'CTRL_AAF', 'TOMMO_AAF', 'IN_TOMMO', 'TOMMO_FILTER'
    ]
    for col in needed_cols:
        if col not in df.columns:
            raise ValueError(f"输入缺少必要列: {col}")

    # 类型转换
    for col in ['CTRL_MAF', 'CTRL_AAF', 'TOMMO_AAF']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df['IN_TOMMO'] = df['IN_TOMMO'].astype('boolean')
    # 规范 TOMMO_FILTER，大写并将缺失视为非 PASS
    df['TOMMO_FILTER'] = df['TOMMO_FILTER'].astype(str).str.upper()

    # 基于 CTRL_MAF 的三组划分（统一助手，NaN 不计入）
    df['__GROUP__'] = assign_maf_group(df['CTRL_MAF'])
    rare = df[df['__GROUP__'] == 'rare'].copy()
    lowf = df[df['__GROUP__'] == 'lowfreq'].copy()
    comm = df[df['__GROUP__'] == 'common'].copy()

    groups = [
        ("Rare Variant (<0.01)", rare),
        ("Low Frequency Variant (0.01~0.05)", lowf),
        ("Common Variant (>0.05)", comm),
    ]

    def _subset_for_scatter(g):
        g2 = g.copy()
        g2 = g2[['CTRL_AAF', 'TOMMO_AAF', 'TOMMO_FILTER']].dropna()
        if max_points is not None and len(g2) > max_points:
            g2 = g2.sample(n=max_points, random_state=42)
        return g2

    def _subset_for_hist(g):
        g2 = g.copy()
        g2 = g2[['CTRL_AAF', 'TOMMO_AAF']].dropna()
        if max_points is not None and len(g2) > max_points:
            g2 = g2.sample(n=max_points, random_state=42)
        g2['DIFF'] = g2['CTRL_AAF'] - g2['TOMMO_AAF']
        return g2

    if output_pdf is None:
        base = os.path.splitext(variant_qc_with_tommo)[0]
        output_pdf = base + ".panel_compare.pdf"

    with PdfPages(output_pdf) as pdf, tempfile.TemporaryDirectory(prefix="panel_png_") as pngdir:
        # 第 1 页：三组统计表（含百分比列 + Count 使用三位逗号分隔）
        fig, axes = plt.subplots(1, 3, figsize=(14, 5))
        for ax, (title, g) in zip(axes, groups):
            g_stats = g.copy()
            # 只对定义组内的有效行计数（CTRL_MAF 非空）
            g_stats = g_stats[~g_stats['CTRL_MAF'].isna()]
            total = int(len(g_stats))
            in_tommo = int((g_stats['IN_TOMMO'] == True).sum())
            pass_tommo = int((g_stats['TOMMO_FILTER'] == 'PASS').sum())

            def _pct(n):
                return (n / total * 100.0) if total > 0 else 0.0

            table_df = pd.DataFrame({
                'Metric': ['Total Variants', 'In ToMMo', 'Pass ToMMo'],
                'Count': [f"{total:,}", f"{in_tommo:,}", f"{pass_tommo:,}"],
                'Percent': [f"{100.0:.1f}%", f"{_pct(in_tommo):.1f}%", f"{_pct(pass_tommo):.1f}%"],
            })
            ax.axis('off')
            ax.set_title(title, fontsize=12)
            tbl = ax.table(cellText=table_df.values,
                           colLabels=table_df.columns,
                           cellLoc='center', loc='center')
            tbl.scale(1, 1.3)
            for (row, col), cell in tbl.get_celld().items():
                if row == 0:  # header row
                    cell.set_text_props(weight='bold')
            # 轻量脚注（只在第一个子图放一次）
            if ax is axes[0]:
                ax.text(0.0, -0.15,
                        'Counts use comma separators; Percent is relative to each group\'s Total.',
                        transform=ax.transAxes, fontsize=8, ha='left', va='top', color='#555555')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # ===== 帮助函数：将三张PNG以网格方式插入到一个PDF页面 =====
        def _compose_three_pngs_to_pdf_page(png_paths, titles, page_title_suffix=None):
            # 更紧凑的布局，三图之间几乎无间距
            fig, axes = plt.subplots(1, 3, figsize=page23_figsize, gridspec_kw={'wspace': 0.0})

            # 先绘图，不在 axes 上放标题，避免标题与图像错位
            full_titles = []
            for ax, path, t in zip(axes, png_paths, titles):
                full_t = t + ('' if page_title_suffix is None else page_title_suffix)
                full_titles.append(full_t)
                if path is None:
                    ax.text(0.5, 0.5, 'No Data', ha='center', va='center')
                    ax.axis('off')
                else:
                    img = plt.imread(path)
                    # 使用 aspect='auto' 使位图填满轴域，避免留空导致“标题看起来偏移”
                    ax.imshow(img, interpolation='none', aspect='auto')
                    ax.axis('off')
                # 去掉任何默认边距
                ax.margins(0)

            # 极小页边距；标题统一用 fig.text 精准居中到各轴上方
            fig.subplots_adjust(left=0.01, right=0.99, bottom=0.06, top=0.94)
            for ax, full_t in zip(axes, full_titles):
                bbox = ax.get_position()
                cx = (bbox.x0 + bbox.x1) / 2.0 + 0.025  # 中心位置 + 少许偏移
                ty = bbox.y1 + 0.001
                fig.text(cx, ty, full_t, ha='center', va='bottom', fontsize=12)

            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # ===== 页面 2：三组散点（颜色=是否 PASS），每个分图先渲染为PNG =====
        png_paths_page2 = []
        titles_page2 = []
        for (title, g) in groups:
            g2 = _subset_for_scatter(g)
            if g2.empty:
                png_paths_page2.append(None)
                titles_page2.append(title)
                continue
            # 单独渲染一个小图为PNG
            f_sc, ax_sc = plt.subplots(figsize=(6, 6))  # 方形画布，减少后续缩放插值
            is_pass = (g2['TOMMO_FILTER'] == 'PASS')
            # 先画 PASS，再画 Non-PASS，让 Non-PASS 叠在上层
            ax_sc.scatter(
                g2.loc[is_pass, 'TOMMO_AAF'], g2.loc[is_pass, 'CTRL_AAF'],
                s=ms_pass, alpha=alpha_pass, linewidths=0, label='PASS in ToMMo', color=pass_color, zorder=9
            )
            ax_sc.scatter(
                g2.loc[~is_pass, 'TOMMO_AAF'], g2.loc[~is_pass, 'CTRL_AAF'],
                s=ms_nonpass, alpha=alpha_nonpass, linewidths=0, label='Non-PASS in ToMMo', color=nonpass_color, zorder=11
            )
            # y=x 参考线
            ax_sc.plot([0, 1], [0, 1], linestyle=refline_ls, linewidth=refline_lw, color=refline_color, zorder=12)
            ax_sc.minorticks_on()
            ax_sc.legend(frameon=False, fontsize=9)
            ax_sc.set_xlabel('TOMMO_AAF')
            ax_sc.set_ylabel('CTRL_AAF')
            ax_sc.set_xlim(0, 1)
            ax_sc.set_ylim(0, 1)
            ax_sc.set_box_aspect(1)
            ax_sc.xaxis.set_major_locator(MaxNLocator(nbins=5))
            ax_sc.yaxis.set_major_locator(MaxNLocator(nbins=5))
            ax_sc.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
            f_sc.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.94)
            out_png = os.path.join(pngdir, f"page2_{title.replace(' ', '_').replace('>', 'gt').replace('<', 'lt')}.png")
            f_sc.savefig(out_png, dpi=png_dpi)  # 避免 tight 导致的再次缩放
            plt.close(f_sc)
            png_paths_page2.append(out_png)
            titles_page2.append(title)
        _compose_three_pngs_to_pdf_page(png_paths_page2, titles_page2, page_title_suffix='')

        # ===== 页面 3：三组散点（颜色=变异类型：SNP vs InDel；不筛 PASS），每个分图先渲染为PNG =====
        png_paths_page3_all = []
        titles_page3_all = []

        # 复用已有的 SNP/InDel 判别函数
        def _classify_variant_type_from_vid(series_vid): # type: ignore
            """根据 VARIANT_ID( CHROM:POS:REF:ALT ) 判断变异类型。
            - 若 REF 和 ALT 均为单碱基且属于 {A,T,C,G} → 'SNP'
            - 否则 → 'InDel'
            返回同长度的 Series，值为 'SNP' 或 'InDel'。
            """
            atcg = {"A", "T", "C", "G"}
            def _one(vid):
                if not isinstance(vid, str):
                    return 'InDel'
                parts = vid.split(':', 3)
                if len(parts) != 4:
                    return 'InDel'
                ref, alt = parts[2], parts[3]
                if ref in atcg and alt in atcg and len(ref) == 1 and len(alt) == 1:
                    return 'SNP'
                return 'InDel'
            return series_vid.apply(_one)

        for (title, g) in groups:
            cols_needed = ['VARIANT_ID', 'CTRL_AAF', 'TOMMO_AAF']
            if not all(c in g.columns for c in cols_needed):
                png_paths_page3_all.append(None)
                titles_page3_all.append(title + ' (All variants)')
                continue
            g2 = g[cols_needed].dropna()
            if max_points is not None and len(g2) > max_points:
                g2 = g2.sample(n=max_points, random_state=42)
            if g2.empty:
                png_paths_page3_all.append(None)
                titles_page3_all.append(title + ' (All variants)')
                continue

            # 变异类型标注
            g2 = g2.assign(TYPE=_classify_variant_type_from_vid(g2['VARIANT_ID']))
            is_snp = (g2['TYPE'] == 'SNP')
            is_indel = (g2['TYPE'] == 'InDel')

            # 绘制（SNP 用 snp_color，InDel 用 indel_color）
            f_sc, ax_sc = plt.subplots(figsize=(6, 6))  # 方形画布
            if is_snp.any():
                ax_sc.scatter(
                    g2.loc[is_snp, 'TOMMO_AAF'], g2.loc[is_snp, 'CTRL_AAF'],
                    s=ms_pass, alpha=alpha_pass, linewidths=0,
                    color=snp_color, label='SNP', zorder=3
                )
            if is_indel.any():
                ax_sc.scatter(
                    g2.loc[is_indel, 'TOMMO_AAF'], g2.loc[is_indel, 'CTRL_AAF'],
                    s=ms_nonpass, alpha=alpha_nonpass, linewidths=0,
                    color=indel_color, label='InDel', zorder=4
                )

            # 参考线 y=x（置于最上）
            ax_sc.plot([0, 1], [0, 1], linestyle=refline_ls, linewidth=refline_lw, color=refline_color, zorder=10)
            ax_sc.minorticks_on()
            ax_sc.legend(frameon=False, fontsize=9)
            ax_sc.set_xlabel('TOMMO_AAF')
            ax_sc.set_ylabel('CTRL_AAF')
            ax_sc.set_xlim(0, 1)
            ax_sc.set_ylim(0, 1)
            ax_sc.set_box_aspect(1)
            ax_sc.xaxis.set_major_locator(MaxNLocator(nbins=5))
            ax_sc.yaxis.set_major_locator(MaxNLocator(nbins=5))
            ax_sc.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
            f_sc.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.94)

            out_png = os.path.join(pngdir, f"page3_all_{title.replace(' ', '_').replace('>', 'gt').replace('<', 'lt')}_SNP_InDel.png")
            f_sc.savefig(out_png, dpi=png_dpi)
            plt.close(f_sc)
            png_paths_page3_all.append(out_png)
            titles_page3_all.append(title + ' (All variants)')

        _compose_three_pngs_to_pdf_page(png_paths_page3_all, titles_page3_all, page_title_suffix='')

        # ===== 页面 4：三组散点（仅 PASS；颜色=变异类型：SNP vs InDel），每个分图先渲染为PNG =====
        png_paths_page3 = []
        titles_page3 = []

        def _classify_variant_type_from_vid(series_vid):
            """根据 VARIANT_ID( CHROM:POS:REF:ALT ) 判断变异类型。
            - 若 REF 和 ALT 均为单碱基且属于 {A,T,C,G} → 'SNP'
            - 否则 → 'InDel'
            返回同长度的 Series，值为 'SNP' 或 'InDel'。
            """
            atcg = {"A", "T", "C", "G"}
            def _one(vid):
                if not isinstance(vid, str):
                    return 'InDel'
                parts = vid.split(':', 3)
                if len(parts) != 4:
                    return 'InDel'
                ref, alt = parts[2], parts[3]
                # 严格按单碱基界定 SNP
                if ref in atcg and alt in atcg and len(ref) == 1 and len(alt) == 1:
                    return 'SNP'
                return 'InDel'
            return series_vid.apply(_one)

        for (title, g) in groups:
            # 仅 PASS
            g_pass = g[g['TOMMO_FILTER'] == 'PASS'].copy()
            # 仅保留绘图所需列，并移除缺失
            cols_needed = ['VARIANT_ID', 'CTRL_AAF', 'TOMMO_AAF']
            if not all(c in g_pass.columns for c in cols_needed):
                png_paths_page3.append(None)
                titles_page3.append(title + ' (PASS ToMMo)')
                continue
            g2 = g_pass[cols_needed].dropna()
            if max_points is not None and len(g2) > max_points:
                g2 = g2.sample(n=max_points, random_state=42)
            if g2.empty:
                png_paths_page3.append(None)
                titles_page3.append(title + ' (PASS ToMMo)')
                continue

            # 变异类型标注
            g2 = g2.assign(TYPE=_classify_variant_type_from_vid(g2['VARIANT_ID']))
            is_snp = (g2['TYPE'] == 'SNP')
            is_indel = (g2['TYPE'] == 'InDel')

            # 绘制（SNP 用 snp_color，InDel 用 indel_color）
            f_sc, ax_sc = plt.subplots(figsize=(6, 6))  # 方形画布
            if is_snp.any():
                ax_sc.scatter(
                    g2.loc[is_snp, 'TOMMO_AAF'], g2.loc[is_snp, 'CTRL_AAF'],
                    s=ms_pass, alpha=alpha_pass, linewidths=0,
                    color=snp_color, label='SNP', zorder=3
                )
            if is_indel.any():
                ax_sc.scatter(
                    g2.loc[is_indel, 'TOMMO_AAF'], g2.loc[is_indel, 'CTRL_AAF'],
                    s=ms_nonpass, alpha=alpha_nonpass, linewidths=0,
                    color=indel_color, label='InDel', zorder=4
                )

            # 参考线 y=x（置于最上）
            ax_sc.plot([0, 1], [0, 1], linestyle=refline_ls, linewidth=refline_lw, color=refline_color, zorder=10)
            ax_sc.minorticks_on()
            ax_sc.legend(frameon=False, fontsize=9)
            ax_sc.set_xlabel('TOMMO_AAF')
            ax_sc.set_ylabel('CTRL_AAF')
            ax_sc.set_xlim(0, 1)
            ax_sc.set_ylim(0, 1)
            ax_sc.set_box_aspect(1)
            ax_sc.xaxis.set_major_locator(MaxNLocator(nbins=5))
            ax_sc.yaxis.set_major_locator(MaxNLocator(nbins=5))
            ax_sc.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
            f_sc.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.94)

            out_png = os.path.join(pngdir, f"page3_{title.replace(' ', '_').replace('>', 'gt').replace('<', 'lt')}_SNP_InDel.png")
            f_sc.savefig(out_png, dpi=png_dpi)
            plt.close(f_sc)
            png_paths_page3.append(out_png)
            titles_page3.append(title + ' (PASS ToMMo)')

        _compose_three_pngs_to_pdf_page(png_paths_page3, titles_page3, page_title_suffix='')

        # ===== 页面 5：三组直方图（仅 PASS，CTRL_AAF - TOMMO_AAF），保留为矢量 =====
        fig, axes = plt.subplots(1, 3, figsize=page4_figsize, gridspec_kw={'wspace': page4_wspace})
        for ax, (title, g) in zip(axes, groups):
            g_pass = g[g['TOMMO_FILTER'] == 'PASS']
            g3 = _subset_for_hist(g_pass)
            if g3.empty:
                ax.text(0.5, 0.5, 'No Data (PASS ToMMo)', ha='center', va='center')
            else:
                ax.hist(g3['DIFF'], bins=50, alpha=0.8)
                # format y-axis with thousands separators
                ax.yaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{int(x):,}'))
                # mean line
                mean_diff = float(np.nanmean(g3['DIFF'])) if len(g3) else float('nan')
                std_diff = float(np.nanstd(g3['DIFF'])) if len(g3) else float('nan')
                ax.axvline(mean_diff, linestyle='-', linewidth=1.2, color='#555555', alpha=0.9, zorder=11)
                # emphasize zero line
                ax.axvline(0, linestyle=refline_ls, linewidth=refline_lw, color=refline_color, zorder=12)
                ax.set_title(
                    f"{title} (PASS ToMMo)\n$\\mu\\,\\pm\\,\\sigma = {mean_diff:.4f}\\,\\pm\\,{std_diff:.4f}$",
                    fontsize=12
                )
            ax.set_xlabel('CTRL_AAF - TOMMO_AAF')
            ax.set_ylabel('Count')
            ax.set_box_aspect(1)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
            ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
        fig.subplots_adjust(left=0.01, right=0.99, bottom=0.06, top=0.94)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    return output_pdf


def build_grouped_variant_tables(
    variant_qc_with_tommo: str,
    output_dir: Optional[str] = None,
    chunk_size: int = 500_000,
    robust_target: str = 'CTRL_AAF_minus_TOMMO_AAF',
    keep_tmp: bool = False,
) -> str:
    """
    函数名称：build_grouped_variant_tables
    ====================================
    【功能】
    以 `run_plink2_variant_qc_with_tommo` 的输出（*.variant_qc_with_tommo.tsv）为输入，
    按 CTRL_MAF 使用统一 assign_maf_group() 将变体划分为三个**主分组**（rare<0.01；lowfreq 0.01–0.05（含）；common>0.05），并在每个主分组内
    进一步细分 **3 个亚分组**：
      a) TOMMO 中无记录（IN_TOMMO==False）
      b) TOMMO 中有记录，但 TOMMO_FILTER != 'PASS'
      c) TOMMO 中有记录，且 TOMMO_FILTER == 'PASS'

    对每个主分组的 c 亚分组，计算每个变体的 **ROBUST_Z**，作为新列添加：
      - 默认以差值 DIFF = CTRL_AAF - TOMMO_AAF 为目标变量；
      - ROBUST_Z = (DIFF - median(DIFF)) / (1.4826 * MAD(DIFF))，当 MAD==0 时置为 NaN。

    【实现要点】
    - 采用**两遍扫描**与**分块**读取，避免一次性加载大表：
      第一遍仅收集各主分组(c)的 DIFF 值，分别写入临时文件；随后计算每组的 median/MAD；
      第二遍再分块读取并将行路由到 3×3 输出文件，同时为 (c) 组计算并写入 ROBUST_Z。
    - 输出采用按分组拆分的 TSV.GZ 文件，并附带一个 JSON manifest，便于后续按条件过滤。

    参数
    ----
    variant_qc_with_tommo : str
        输入表路径（列至少包含 VARIANT_ID, CTRL_MAF, IN_TOMMO, TOMMO_FILTER, CTRL_AAF, TOMMO_AAF）。
    output_dir : Optional[str]
        输出目录（默认与输入同目录）。
    chunk_size : int
        Pandas 分块大小（默认 500k 行）。
    robust_target : str
        ROBUST_Z 目标度量，当前仅支持 'CTRL_AAF_minus_TOMMO_AAF'。
    keep_tmp : bool
        是否保留临时中间文件（调试用）。

    返回
    ----
    str
        生成的 manifest JSON 路径；其中包含每个输出文件的位置与行数统计、
        以及各主分组(c)用于 ROBUST_Z 的 median/MAD。
    """
    import json
    import gzip
    import shutil
    import numpy as np
    import pandas as pd
    from datetime import datetime

    # ---- 并行参数：合理设置线程数（不改动函数入参与返回） ----
    try:
        import multiprocessing as _mp
        _MAX_WORKERS = max(1, min(8, _mp.cpu_count()))
    except Exception:
        _MAX_WORKERS = 4

    # ---- 路径与输出命名 ----
    in_path = os.path.abspath(variant_qc_with_tommo)
    # 将所有最终输出放在“当下运行文件夹”（当前工作目录）下；
    # 仅临时文件放在带时间标签的子目录中。
    if output_dir is None:
        output_dir = os.getcwd()
    os.makedirs(output_dir, exist_ok=True)

    ts = datetime.now().strftime('%Y%m%d_%H%M%S')
    workdir = os.path.join(os.getcwd(), f"tmp_group_{ts}_{uuid.uuid4().hex[:8]}")
    os.makedirs(workdir, exist_ok=True)

    # ---- 工具函数：解析/归一化列 ----
    def _normalize_bool_series(s):
        return s.map(lambda x: True if x in (True, 1, '1', 'True', 'TRUE')
                     else (False if x in (False, 0, '0', 'False', 'FALSE') else pd.NA))

    def _assign_subgroup(in_tommo: pd.Series, tommo_filter: pd.Series) -> pd.Series:
        tf = tommo_filter.astype(str).str.upper()
        it = _normalize_bool_series(in_tommo)
        sub = pd.Series(pd.NA, index=tf.index, dtype='object')
        sub = sub.mask(it == False, 'a_not_in_tommo')
        sub = sub.mask((it == True) & (tf != 'PASS'), 'b_in_nonpass')
        sub = sub.mask((it == True) & (tf == 'PASS'), 'c_in_pass')
        return sub

    # ---- 需要的列 ----
    usecols = [
        'VARIANT_ID', 'CTRL_MAF', 'IN_TOMMO', 'TOMMO_FILTER', 'CTRL_AAF', 'TOMMO_AAF'
    ]

    # ---- 第一遍（并行计算、串行落盘）：收集各主分组(c)的 DIFF ----
    diff_tmp_paths = {
        'rare':   os.path.join(workdir, 'diff_rare.txt'),
        'lowfreq':os.path.join(workdir, 'diff_lowfreq.txt'),
        'common': os.path.join(workdir, 'diff_common.txt'),
    }
    # 清空文件
    for p in diff_tmp_paths.values():
        open(p, 'w').close()

    def _first_pass_chunk(chunk: pd.DataFrame) -> Dict[str, np.ndarray]:
        """工作线程：从一个分块中提取 (c_in_pass) 的 DIFF，按主分组返回数组。磁盘写入在主线程完成。"""
        mg = assign_maf_group(chunk['CTRL_MAF'])
        sg = _assign_subgroup(chunk['IN_TOMMO'], chunk['TOMMO_FILTER'])
        # 仅 (c) 组且两列可数值化
        mask_c = (sg == 'c_in_pass') & chunk['CTRL_AAF'].notna() & chunk['TOMMO_AAF'].notna()
        if not mask_c.any():
            return {'rare': np.array([], dtype='float64'),
                    'lowfreq': np.array([], dtype='float64'),
                    'common': np.array([], dtype='float64')}
        csub = chunk.loc[mask_c, ['CTRL_AAF', 'TOMMO_AAF']].apply(pd.to_numeric, errors='coerce')
        csub['DIFF'] = csub['CTRL_AAF'] - csub['TOMMO_AAF']
        out = {}
        for gname in ('rare', 'lowfreq', 'common'):
            idx = (mg[mask_c] == gname)
            if idx.any():
                vals = csub.loc[idx.values, 'DIFF'].dropna().to_numpy(dtype='float64') # type: ignore
            else:
                vals = np.array([], dtype='float64')
            out[gname] = vals
        return out

    total_rows = 0
    total_chunks = 0
    print(f"[第一遍] 开始扫描 {in_path}，并行处理分块以收集 DIFF ...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=_MAX_WORKERS) as ex:
        futures = []
        for chunk in pd.read_csv(in_path, sep='\t', usecols=usecols, dtype='string', chunksize=chunk_size, engine='c'):
            total_rows += len(chunk)
            total_chunks += 1
            futures.append(ex.submit(_first_pass_chunk, chunk))
            # 控制队列长度，避免内存占用过高
            if len(futures) >= _MAX_WORKERS * 4:
                done, futures = futures[:], []
                for fut in concurrent.futures.as_completed(done):
                    res = fut.result()
                    for gname, arr in res.items():
                        if arr.size:
                            with open(diff_tmp_paths[gname], 'a') as f:
                                f.write('\n'.join(f"{v:.10g}" for v in arr) + '\n')
                print(f"[第一遍] 并行已完成 {total_chunks} 个分块，共 {total_rows:,} 行 ...")
        # flush remaining
        for fut in concurrent.futures.as_completed(futures):
            res = fut.result()
            for gname, arr in res.items():
                if arr.size:
                    with open(diff_tmp_paths[gname], 'a') as f:
                        f.write('\n'.join(f"{v:.10g}" for v in arr) + '\n')
        print(f"[第一遍] 并行已完成 {total_chunks} 个分块，共 {total_rows:,} 行 ...")
    print(f"[第一遍] 完成。总分块数: {total_chunks}, 总行数: {total_rows:,}.")

    # 计算各主分组的 median 与 MAD
    from scipy.stats import median_abs_deviation
    robust_stats = {}
    for gname, path in diff_tmp_paths.items():
        if os.path.getsize(path) == 0:
            robust_stats[gname] = {'median': None, 'mad': None}
            continue
        arr = np.loadtxt(path, dtype=float, ndmin=1)
        if arr.size == 0:
            robust_stats[gname] = {'median': None, 'mad': None}
            continue
        med = float(np.median(arr))
        mad = float(median_abs_deviation(arr, scale=1, nan_policy="omit"))
        robust_stats[gname] = {'median': med, 'mad': mad}

    # ---- 准备输出文件句柄（按 3×3 拆分） ----
    def _outfile(main_g: str, sub_g: str) -> str:
        return os.path.join(output_dir, f"{os.path.basename(in_path)}.{main_g}.{sub_g}.tsv.gz")

    out_paths = {g: {h: _outfile(g, h) for h in ('a_not_in_tommo','b_in_nonpass','c_in_pass')} for g in ('rare','lowfreq','common')}
    out_files = {g: {h: gzip.open(out_paths[g][h], 'wt') for h in out_paths[g]} for g in out_paths}

    # 写表头（追加两列：GROUP_MAIN, GROUP_SUB；并在 (c) 组包含 DIFF 与 ROBUST_Z）
    header_cols = ['VARIANT_ID','MAF','VMISS','CASE_AAF','CTRL_AAF','CTRL_MAF','CASE_HWE','CTRL_HWE','IN_TOMMO','TOMMO_AAF','TOMMO_FILTER','GROUP_MAIN','GROUP_SUB','DIFF','ROBUST_Z']
    header_line = '\t'.join(header_cols) + '\n'
    for g in out_files:
        for h in out_files[g]:
            out_files[g][h].write(header_line)

    # ---- 第二遍（并行计算、串行落盘）：路由写出，并计算 ROBUST_Z（仅 c 组） ----
    counts = {g: {h: 0 for h in ('a_not_in_tommo','b_in_nonpass','c_in_pass')} for g in ('rare','lowfreq','common')}

    # 为每个输出文件准备一个锁，保证写入的原子性（顺序无关但避免交叉）
    _locks = {g: {h: threading.Lock() for h in out_files[g]} for g in out_files}

    def _second_pass_chunk(chunk: pd.DataFrame) -> Dict[Tuple[str, str], pd.DataFrame]:
        # 统一列
        for col in ['CTRL_AAF','TOMMO_AAF','CTRL_MAF']:
            if col in chunk.columns:
                chunk[col] = pd.to_numeric(chunk[col], errors='coerce')
        mg = assign_maf_group(chunk['CTRL_MAF'])
        sg = _assign_subgroup(chunk['IN_TOMMO'], chunk['TOMMO_FILTER'])
        diff = (chunk['CTRL_AAF'] - chunk['TOMMO_AAF']).astype('float64')

        outdfs: Dict[Tuple[str, str], pd.DataFrame] = {}
        for main_g in ('rare','lowfreq','common'):
            mask_main = (mg == main_g)
            if not mask_main.any():
                continue
            for sub_g in ('a_not_in_tommo','b_in_nonpass','c_in_pass'):
                mask = mask_main & (sg == sub_g)
                if not mask.any():
                    continue
                subdf = chunk.loc[mask, :].copy()
                subdf['GROUP_MAIN'] = main_g
                subdf['GROUP_SUB'] = sub_g
                subdf['DIFF'] = pd.NA
                subdf['ROBUST_Z'] = pd.NA
                if sub_g == 'c_in_pass':
                    med = robust_stats.get(main_g, {}).get('median', None)
                    mad = robust_stats.get(main_g, {}).get('mad', None)
                    if med is not None and mad is not None and mad != 0:
                        subdf['DIFF'] = diff.loc[mask].values
                        subdf['ROBUST_Z'] = ((subdf['DIFF'] - med) / (1.4826 * mad)).astype('float64')
                    elif med is not None and (mad == 0):
                        subdf['DIFF'] = diff.loc[mask].values
                        subdf['ROBUST_Z'] = 0.0
                    else:
                        subdf['DIFF'] = diff.loc[mask].values
                subdf = subdf[['VARIANT_ID','MAF','VMISS','CASE_AAF','CTRL_AAF','CTRL_MAF','CASE_HWE','CTRL_HWE','IN_TOMMO','TOMMO_AAF','TOMMO_FILTER','GROUP_MAIN','GROUP_SUB','DIFF','ROBUST_Z']]
                outdfs[(main_g, sub_g)] = subdf
        return outdfs

    # 流式读取 → 并行转换 → 主线程串行写出
    total_rows2 = 0
    total_chunks2 = 0
    print(f"[第二遍] 开始扫描 {in_path}，并行处理分块并路由写出 ...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=_MAX_WORKERS) as ex:
        futures = []
        for chunk in pd.read_csv(in_path, sep='\t', dtype='string', chunksize=chunk_size, engine='c'):
            futures.append(ex.submit(_second_pass_chunk, chunk))
            total_rows2 += len(chunk)
            total_chunks2 += 1
            if len(futures) >= _MAX_WORKERS * 2:
                done, futures = futures[:], []
                for fut in concurrent.futures.as_completed(done):
                    res = fut.result()
                    for (main_g, sub_g), subdf in res.items():
                        if len(subdf) == 0:
                            continue
                        out = out_files[main_g][sub_g]
                        with _locks[main_g][sub_g]:
                            subdf.to_csv(out, sep='\t', header=False, index=False, na_rep='nan')
                        counts[main_g][sub_g] += len(subdf)
                print(f"[第二遍] 并行已完成 {total_chunks2} 个分块，共 {total_rows2:,} 行 ...")
        for fut in concurrent.futures.as_completed(futures):
            res = fut.result()
            for (main_g, sub_g), subdf in res.items():
                if len(subdf) == 0:
                    continue
                out = out_files[main_g][sub_g]
                with _locks[main_g][sub_g]:
                    subdf.to_csv(out, sep='\t', header=False, index=False, na_rep='nan')
                counts[main_g][sub_g] += len(subdf)
        print(f"[第二遍] 并行已完成 {total_chunks2} 个分块，共 {total_rows2:,} 行 ...")
    print(f"[第二遍] 完成。总分块数: {total_chunks2}, 总行数: {total_rows2:,}.")

    # 关闭文件
    for g in out_files:
        for h in out_files[g]:
            out_files[g][h].close()

    # 清理临时目录
    if keep_tmp:
        tmp_note = os.path.join(workdir, 'KEEP_TMP.txt')
        with open(tmp_note, 'w') as f:
            f.write('临时文件保留以便排错。')
    else:
        try:
            shutil.rmtree(workdir)
        except Exception:
            pass

    # 写 manifest.json，便于后续引用
    manifest = {
        'input': in_path,
        'output_dir': output_dir,
        'robust_target': robust_target,
        'robust_stats': robust_stats,
        'counts': counts,
        'files': out_paths,
    }
    manifest_path = os.path.join(output_dir, 'manifest.json')
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    return manifest_path


def summarize_c_in_pass_thresholds(
    manifest_path: str,
    output_dir: Optional[str] = None,
    chunk_size: int = 1_000_000,
) -> str:
    """
    函数名称：summarize_c_in_pass_thresholds
    ======================================
    【功能】
    读取 `build_grouped_variant_tables` 生成的 manifest.json，仅针对三个主分组
    ('rare', 'lowfreq', 'common') 的 **c_in_pass** 亚分组，基于该亚分组中的 `ROBUST_Z`
    计算“|ROBUST_Z| 阈值”从 1 到 `max`，**步长为 0.2**（例如 t=1..3 → 阈值为 1.0, 1.2, 1.4, ..., 2.8, 3.0；
    其中 `max=ceil(max(|ROBUST_Z|))`，若 <1 则取 1）的统计摘要：
      - `count`：满足 |ROBUST_Z| < 阈值 的 `VARIANT_ID` 数量（同时要求 CTRL_AAF/TOMMO_AAF 可用）；
      - `mse`：在该条件下的 MSE(CTRL_AAF vs TOMMO_AAF)。

    【实现】
    - 对每个主分组的 c_in_pass 文件，先进行**第一遍**分块扫描，获得全组的 `max_abs_robust_z`；
    - 计算 `max_int = max(1, ceil(max_abs_robust_z))`；
    - 再进行**第二遍**分块扫描：对 t=1..max_int，步长为 0.2，累计计数与平方误差和（SSE），最终得到 MSE=SSE/count。
    - 将汇总结果写出为每组一个 TSV 文件，并将路径与关键统计更新写回 manifest.json。

    参数
    ----
    manifest_path : str
        `build_grouped_variant_tables` 输出的 manifest.json 路径。
    output_dir : Optional[str]
        摘要 TSV 输出目录；默认与 manifest.json 同目录。
    chunk_size : int
        分块大小（默认 1,000,000 行）。

    返回
    ----
    str
        更新后的 manifest.json 路径（与输入相同）。
    """
    import json
    import os
    import numpy as np
    import pandas as pd

    # 读取 manifest
    manifest_path = os.path.abspath(manifest_path)
    with open(manifest_path, 'r') as f:
        manifest = json.load(f)

    if output_dir is None:
        output_dir = os.path.dirname(manifest_path) or os.getcwd()
    os.makedirs(output_dir, exist_ok=True)

    # —— 并行设置（不改变输入输出）：合理控制并发线程数 ——
    try:
        import multiprocessing as _mp
        _MAX_WORKERS = max(1, min(8, _mp.cpu_count()))
    except Exception:
        _MAX_WORKERS = 4

    files_map = manifest.get('files', {})
    if not files_map:
        raise ValueError("manifest.json 缺少 'files' 字段，或为空。")

    # 仅处理三个主分组
    main_groups = ['rare', 'lowfreq', 'common']
    subgroup = 'c_in_pass'

    # 结果写回区域
    out_key = 'c_in_pass_thresholds'
    manifest[out_key] = manifest.get(out_key, {})

    def _first_pass_max_abs(path: str) -> float:
        """第一遍：并行扫描 ROBUST_Z 的绝对值最大值。忽略 NaN。
        采用线程池对每个 chunk 的局部最大值并行计算，再在主线程归约为全局最大值。
        同时打印中文进度：分块数与并行处理提示。
        """
        import concurrent.futures
        max_abs = 0.0
        usecols = ['ROBUST_Z']
        total_chunks = 0
        total_rows = 0
        print(f"[第一遍] 开始扫描 {path}，并行处理分块以计算 max(|ROBUST_Z|) ...")

        def _local_max(chunk) -> float:
            if 'ROBUST_Z' not in chunk:
                return 0.0
            z = chunk['ROBUST_Z'].to_numpy(dtype='float64', copy=False)
            if z.size == 0:
                return 0.0
            z = z[~np.isnan(z)]  # 去 NaN
            if z.size == 0:
                return 0.0
            return float(np.max(np.abs(z)))

        futures = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=_MAX_WORKERS) as ex:
            for chunk in pd.read_csv(path, sep='\t', usecols=usecols, dtype='float64',
                                     chunksize=chunk_size, engine='c', na_values=['nan', 'NaN', 'NA', '.']):
                total_chunks += 1
                total_rows += len(chunk)
                futures.append(ex.submit(_local_max, chunk))
                # 背压：控制在飞任务，防止内存上涨
                if len(futures) >= _MAX_WORKERS * 4:
                    done, futures = futures[:], []
                    for fut in concurrent.futures.as_completed(done):
                        loc = fut.result()
                        if loc > max_abs:
                            max_abs = loc
                    print(f"[第一遍] 并行已完成 {total_chunks} 个分块，共 {total_rows:,} 行 ...")
            # flush remaining
            for fut in concurrent.futures.as_completed(futures):
                loc = fut.result()
                if loc > max_abs:
                    max_abs = loc
        print(f"[第一遍] 完成。总分块数: {total_chunks}, 总行数: {total_rows:,}, 最大绝对值: {max_abs}")
        return max_abs

    def _second_pass_summary(path: str, max_int: int) -> pd.DataFrame:
        """第二遍：并行累计 t=1..max_int（步长0.2）下的 count 与 SSE，最终计算 MSE。
        返回包含列：threshold_abs_robust_z, count_variants, mse_ctrl_vs_tommo。
        同时打印中文进度：分块数与并行处理提示。
        """
        import concurrent.futures
        thresholds = np.arange(1.0, max_int + 0.2, 0.2, dtype='float64')  # 1.0, 1.2, ..., max_int
        k = thresholds.size
        sse_total = np.zeros(k, dtype='float64')
        cnt_total = np.zeros(k, dtype='int64')
        usecols = ['VARIANT_ID', 'CTRL_AAF', 'TOMMO_AAF', 'ROBUST_Z']
        total_chunks = 0
        total_rows = 0
        print(f"[第二遍] 开始统计 {path}，并行处理阈值 1..{max_int} (步长=0.2) ...")

        def _acc_chunk(chunk) -> tuple:
            # 转换类型并构建掩码
            z = pd.to_numeric(chunk['ROBUST_Z'], errors='coerce').to_numpy(dtype='float64')
            ca = pd.to_numeric(chunk['CTRL_AAF'], errors='coerce').to_numpy(dtype='float64')
            ta = pd.to_numeric(chunk['TOMMO_AAF'], errors='coerce').to_numpy(dtype='float64')
            valid = (~np.isnan(z)) & (~np.isnan(ca)) & (~np.isnan(ta))
            if not np.any(valid):
                return (np.zeros(k, dtype='float64'), np.zeros(k, dtype='int64'), True)
            z_abs = np.abs(z[valid])
            err2 = (ca[valid] - ta[valid]) ** 2
            sse = np.zeros(k, dtype='float64')
            cnt = np.zeros(k, dtype='int64')
            for i, t in enumerate(thresholds):
                mask = z_abs < t
                if np.any(mask):
                    cnt[i] = int(mask.sum())
                    sse[i] = float(err2[mask].sum())
            return (sse, cnt, False)

        futures = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=_MAX_WORKERS) as ex:
            for chunk in pd.read_csv(path, sep='\t', usecols=usecols, dtype='string',
                                     chunksize=chunk_size, engine='c', na_values=['nan', 'NaN', 'NA', '.']):
                total_chunks += 1
                total_rows += len(chunk)
                futures.append(ex.submit(_acc_chunk, chunk))
                if len(futures) >= _MAX_WORKERS * 3:
                    done, futures = futures[:], []
                    for fut in concurrent.futures.as_completed(done):
                        sse, cnt, empty = fut.result()
                        sse_total += sse
                        cnt_total += cnt
                    print(f"[第二遍] 并行已完成 {total_chunks} 个分块，共 {total_rows:,} 行 ...")
            for fut in concurrent.futures.as_completed(futures):
                sse, cnt, empty = fut.result()
                sse_total += sse
                cnt_total += cnt
        print(f"[第二遍] 完成。总分块数: {total_chunks}, 总行数: {total_rows:,}.")

        mse = np.full_like(sse_total, fill_value=np.nan, dtype='float64')
        nz = cnt_total > 0
        mse[nz] = sse_total[nz] / cnt_total[nz]
        out = pd.DataFrame({
            'threshold_abs_robust_z': thresholds,
            'count_variants': cnt_total,
            'mse_ctrl_vs_tommo': mse,
        })
        return out

    for main in main_groups:
        path = files_map.get(main, {}).get(subgroup)
        if not path:
            continue
        if not os.path.exists(path):
            # 若路径为相对路径，尝试相对 manifest 的目录
            cand = os.path.join(os.path.dirname(manifest_path), os.path.basename(path))
            if os.path.exists(cand):
                path = cand
            else:
                # 记录缺失
                manifest[out_key][main] = {
                    'summary_tsv': None,
                    'max_abs_robust_z': None,
                    'max_int': None,
                    'note': f"missing file: {path}"
                }
                continue

        # 第一遍：找最大绝对值
        max_abs = _first_pass_max_abs(path)
        max_int = int(np.ceil(max_abs))
        if max_int < 1:
            max_int = 1

        # 第二遍：统计 1..max_int 的 count 与 MSE
        df_sum = _second_pass_summary(path, max_int)

        # 写出该主分组的摘要 TSV
        base = os.path.basename(path)
        out_tsv = os.path.join(
            output_dir,
            f"{base}.c_in_pass.robustz_threshold_summary.tsv"
        )
        df_sum.to_csv(out_tsv, sep='\t', index=False)

        # 回写 manifest
        manifest[out_key][main] = {
            'summary_tsv': out_tsv,
            'max_abs_robust_z': float(max_abs) if np.isfinite(max_abs) else None,
            'max_int': int(max_int),
            'note': 'thresholds are 1..max_int with 0.2 step; count excludes rows without CTRL_AAF/TOMMO_AAF/ROBUST_Z.'
        }

    # 保存更新后的 manifest
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    return manifest_path



# new function based the updated manifest_path structure for plotting
def plot_c_in_pass_threshold_tradeoff(
    manifest_path: str,
    output_pdf: Optional[str] = None,
    figsize: tuple = (14, 5),
    marker_size: int = 18,
    line_width: float = 1.5,
    y_min_zero: bool = True,
    use_log_y: bool = False,
    png_dpi: int = 600,
    # KneeLocator 可调参数（默认按你的要求）
    knee_curve: str = 'concave',
    knee_direction: str = 'increasing',
    knee_S: float = 1.0,
    knee_weight_x: float = 1.0,
    knee_weight_y: float = 1.0,
    knee_weight_y_map: Optional[Dict[str, float]] = None,
) -> str:
    """
    函数名称：plot_c_in_pass_threshold_tradeoff
    =========================================
    【功能】
    读取 `summarize_c_in_pass_thresholds` 更新后的 manifest.json，仅针对 `c_in_pass_thresholds`
    的 3 个主分组（rare/lowfreq/common），分别加载其 `summary_tsv`，
    在一页 PDF 中绘制 3 个子图：
      - X 轴：`mse_ctrl_vs_tommo`
      - Y 轴：`count_variants`
      - 点按 `threshold_abs_robust_z` 升序连接（scatter+line）
    输出 PDF 文件路径。
    并将每个子图标注（阈值/计数与比例/MSE）写回 manifest.json 的 `c_in_pass_knee` 字段。

    现在支持为不同主分组（rare, lowfreq, common）分别指定 knee_weight_y，覆盖全局默认值。

    另新增一页：针对每个主分组，在 `|ROBUST_Z| < knee_threshold` 条件下，
    绘制这些**被计入 count 的变体**的散点图（X=TOMMO_AAF, Y=CTRL_AAF），
    采用 PNG (dpi=600) 先渲染后插入 PDF，提高性能。

    参数
    ----
    manifest_path : str
        `build_grouped_variant_tables` → `summarize_c_in_pass_thresholds` 后的 manifest.json 路径。
    output_pdf : Optional[str]
        输出 PDF 路径；默认与 manifest 同目录，文件名为 `c_in_pass_threshold_tradeoff.pdf`。
    figsize : tuple
        单页画布尺寸（英寸）。
    marker_size : int
        散点大小。
    line_width : float
        折线宽度。
    y_min_zero : bool
        若为 True，则 y 轴下界以 0 起（当数据允许）。
    use_log_y : bool
        若为 True，则使用对数 y 轴（适合跨度很大时）。
    png_dpi : int
        新增散点页的位图渲染分辨率（DPI），默认 600。
    knee_curve, knee_direction, knee_S :
        传给 KneeLocator 的 `curve`, `direction`, `S`，默认分别为 `'concave'`, `'increasing'`, `1.0`。
    knee_weight_x, knee_weight_y :
        KneeLocator 的权重（若安装的 `kneed` 版本不支持，将自动回退为不带权重调用）。默认 `1.0` 与 `2.3`。
    knee_weight_y_map : dict, optional
        针对不同主分组设置的 weight_y，例如 {'rare':2.0,'lowfreq':2.3,'common':1.5}。
        若未提供，则使用全局 knee_weight_y。

    返回
    ----
    str
        输出 PDF 文件路径。
    """
    import os
    import json
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.ticker import MaxNLocator, FuncFormatter
    from kneed import KneeLocator

    manifest_path = os.path.abspath(manifest_path)
    with open(manifest_path, 'r') as f:
        manifest = json.load(f)

    ckey = 'c_in_pass_thresholds'
    if ckey not in manifest or not manifest[ckey]:
        raise ValueError("manifest.json 中未找到 'c_in_pass_thresholds' 字段或其为空。请先运行 summarize_c_in_pass_thresholds().")

    # 3 个主分组
    order = [('Rare Variant', 'rare'), ('Low Frequency Variant', 'lowfreq'), ('Common Variant', 'common')]
    key_map = { 'Rare Variant': 'rare', 'Low Frequency Variant': 'lowfreq', 'Common Variant': 'common' }

    # 保存拐点注释，稍后写回 manifest
    knee_notes = {}

    # 读取各组 TSV
    plots = []
    for title, key in order:
        info = manifest[ckey].get(key, {})
        tsv = info.get('summary_tsv')
        if not tsv or not os.path.exists(tsv):
            plots.append((title, None))
            continue
        df = pd.read_csv(tsv, sep='\t')
        # 期望列
        need = ['threshold_abs_robust_z', 'count_variants', 'mse_ctrl_vs_tommo']
        if not all(c in df.columns for c in need):
            plots.append((title, None))
            continue
        # 清理 & 排序
        df = df[need].dropna()
        df = df.sort_values('threshold_abs_robust_z', kind='mergesort').reset_index(drop=True)
        # 确保类型
        df['count_variants'] = pd.to_numeric(df['count_variants'], errors='coerce')
        df['mse_ctrl_vs_tommo'] = pd.to_numeric(df['mse_ctrl_vs_tommo'], errors='coerce')
        df = df.dropna()
        plots.append((title, df))

    # 确定输出路径
    if output_pdf is None:
        out_dir = os.path.dirname(manifest_path) or os.getcwd()
        output_pdf = os.path.join(out_dir, 'c_in_pass_threshold_tradeoff.pdf')

    import tempfile
    with PdfPages(output_pdf) as pdf:
        fig, axes = plt.subplots(1, 3, figsize=figsize, gridspec_kw={'wspace': 0.20})
        # 统一样式
        plt.rcParams.update({
            'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans'],
            'font.size': 10,
            'axes.titlesize': 11,
            'axes.labelsize': 10,
            'xtick.labelsize': 9,
            'ytick.labelsize': 9,
            'axes.linewidth': 0.8,
        })

        title_map = {
            'Rare Variant': 'Rare Variant (<0.01) (PASS ToMMo)',
            'Low Frequency Variant': 'Low Frequency Variant (0.01~0.05) (PASS ToMMo)',
            'Common Variant': 'Common Variant (>0.05) (PASS ToMMo)'
        }

        for ax, (title, df) in zip(axes, plots):
            ax.set_title(title_map.get(title, title))
            if df is None or df.empty:
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center')
                ax.set_xlabel('mse_ctrl_vs_tommo')
                ax.set_ylabel('count_variants')
                ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
                continue

            # 从 manifest 里取该主分组 c_in_pass 的总数，用于显示比例
            mg_key = key_map.get(title, None)
            total_c_in_pass = None
            if mg_key is not None:
                total_c_in_pass = manifest.get('counts', {}).get(mg_key, {}).get('c_in_pass', None)

            x = df['mse_ctrl_vs_tommo'].to_numpy(dtype='float64')
            y = df['count_variants'].to_numpy(dtype='float64')

            # 按阈值顺序连接（df 已按 threshold_abs_robust_z 升序）
            ax.plot(x, y, linewidth=line_width, alpha=0.9, zorder=2)
            ax.scatter(x, y, s=marker_size, alpha=0.8, zorder=3)

            # === 使用 KneeLocator 寻找拐点 ===
            # 允许为不同主分组指定 weight_y
            wy = knee_weight_y_map.get(mg_key, knee_weight_y) if knee_weight_y_map else knee_weight_y # type: ignore
            knee_idx = None
            knee_x = knee_y = None
            knee_thr = None
            try:
                if len(x) >= 3:
                    order_x = np.argsort(x)
                    x_sorted = x[order_x]
                    y_sorted = y[order_x]
                    try:
                        kl = KneeLocator(
                            x_sorted, y_sorted,
                            curve=knee_curve, direction=knee_direction, S=knee_S,
                            weight_x=knee_weight_x, weight_y=wy,
                        )
                    except TypeError:
                        kl = KneeLocator(
                            x_sorted, y_sorted,
                            curve=knee_curve, direction=knee_direction, S=knee_S,
                        )
                    knee_x = kl.knee
                    if knee_x is None:
                        knee_x = kl.elbow
                    if knee_x is not None:
                        orig_idx = int(np.nanargmin(np.abs(x - knee_x)))
                        knee_idx = orig_idx
                        knee_y = float(y[knee_idx])
                        knee_thr = float(df.loc[df.index[knee_idx], 'threshold_abs_robust_z'])
                        # 标记拐点：红色星星
                        ax.scatter([knee_x], [knee_y], marker='*', s=220, color='red', zorder=10)
                        # 注释文本：threshold / count / mse
                        if isinstance(total_c_in_pass, (int, float)) and total_c_in_pass and total_c_in_pass > 0:
                            pct = 100.0 * (knee_y / float(total_c_in_pass))
                            count_str = f"{int(knee_y):,}/{int(total_c_in_pass):,} ({pct:.1f}%)"
                        else:
                            count_str = f"{int(knee_y):,}"
                        ann = (
                            f"thr={knee_thr:.2f}\n"
                            f"count={count_str}\n"
                            f"mse={knee_x:.6g}"
                        )
                        # 自适应注释偏移，尽量避开边缘与曲线
                        xlim = ax.get_xlim(); ylim = ax.get_ylim()
                        xn = 0.0 if xlim[1] == xlim[0] else (knee_x - xlim[0]) / (xlim[1] - xlim[0])
                        yn = 0.0 if ylim[1] == ylim[0] else (knee_y - ylim[0]) / (ylim[1] - ylim[0])
                        if xn >= 0.5 and yn >= 0.5:
                            offset = (-12, -12); ha, va = 'right', 'top'
                        elif xn < 0.5 and yn >= 0.5:
                            offset = (12, -12); ha, va = 'left', 'top'
                        elif xn >= 0.5 and yn < 0.5:
                            offset = (-12, 12); ha, va = 'right', 'bottom'
                        else:
                            offset = (12, 12); ha, va = 'left', 'bottom'
                        ax.annotate(
                            ann,
                            xy=(knee_x, knee_y),
                            xytext=offset, textcoords='offset points',
                            fontsize=9, color='red', ha=ha, va=va,
                            bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='red', lw=0.8, alpha=0.95),
                            arrowprops=dict(arrowstyle='->', lw=0.8, color='red')
                        )
                        # 记录注释内容以写回 manifest
                        knee_notes[mg_key] = {
                            'threshold_abs_robust_z': float(knee_thr),
                            'count_variants_at_knee': int(knee_y),
                            'total_c_in_pass': int(total_c_in_pass) if isinstance(total_c_in_pass, (int, float)) and total_c_in_pass is not None else None,
                            'percent_of_c_in_pass': (float(knee_y) / float(total_c_in_pass) * 100.0) if isinstance(total_c_in_pass, (int, float)) and total_c_in_pass and total_c_in_pass > 0 else None,
                            'mse_ctrl_vs_tommo_at_knee': float(knee_x),
                            'annotation_text': ann,
                            'knee_params': {
                                'curve': knee_curve,
                                'direction': knee_direction,
                                'S': float(knee_S),
                                'weight_x': float(knee_weight_x),
                                'weight_y': float(wy),
                            }
                        }
            except Exception:
                # 避免因个别数据异常导致绘图失败；静默跳过拐点标注
                pass

            ax.set_xlabel('mse_ctrl_vs_tommo')
            ax.set_ylabel('count_variants')
            if y_min_zero:
                try:
                    xmin = min(0.0, float(np.nanmin(x)))
                except ValueError:
                    xmin = 0.0
                ax.set_xlim(left=xmin)
            if use_log_y:
                ax.set_xscale('log')

            # 让坐标轴整数刻度更友好
            ax.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=False))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=6, integer=False))
            ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)

        fig.subplots_adjust(left=0.05, right=0.98, bottom=0.10, top=0.90)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # ===== 新增页面：按 knee 阈值选中的变体散点（X=TOMMO_AAF, Y=CTRL_AAF），每个分组一张PNG再合成一页 =====
        # 需要：各分组的 c_in_pass 路径 + knee 阈值
        files_map = manifest.get('files', {})
        title_map2 = {
            'rare': 'Rare Variant (<0.01) (PASS ToMMo)',
            'lowfreq': 'Low Frequency Variant (0.01~0.05) (PASS ToMMo)',
            'common': 'Common Variant (>0.05) (PASS ToMMo)'
        }

        # 小工具：组合三张 PNG 到单页
        def _compose_three_pngs_to_pdf_page(png_paths, titles, subtitle_suffixes=None):
            fig2, axes2 = plt.subplots(1, 3, figsize=figsize, gridspec_kw={'wspace': 0.0})
            for ax2, path2, t2, sub2 in zip(axes2, png_paths, titles, subtitle_suffixes or ['']*3):
                full_t = t2 + ('' if not sub2 else f"\n{sub2}")
                if path2 is None:
                    ax2.text(0.5, 0.5, 'No Data', ha='center', va='center')
                else:
                    img2 = plt.imread(path2)
                    ax2.imshow(img2, interpolation='none', aspect='auto')
                # 使用轴自身标题，避免 fig.text 布局偏移导致错位
                ax2.set_title(full_t, fontsize=12, pad=6)
                ax2.title.set_x(0.57)
                ax2.set_xticks([])
                ax2.set_yticks([])
                for spine in ax2.spines.values():
                    spine.set_visible(False)
            fig2.subplots_adjust(left=0.01, right=0.99, bottom=0.06, top=0.94)
            pdf.savefig(fig2, bbox_inches='tight')
            plt.close(fig2)

        # ===== 新增页面并行渲染：按 knee 阈值选中的变体散点（X=TOMMO_AAF, Y=CTRL_AAF）=====
        # 使用线程池并行绘制每个分组的 PNG，全部完成后再合成一页 PDF
        import concurrent.futures

        png_paths = []
        titles2 = []
        subtitles = []

        # 需要：各分组的 c_in_pass 路径 + knee 阈值
        files_map = manifest.get('files', {})
        title_map2 = {
            'rare': 'Rare Variant (<0.01) (PASS ToMMo)',
            'lowfreq': 'Low Frequency Variant (0.01~0.05) (PASS ToMMo)',
            'common': 'Common Variant (>0.05) (PASS ToMMo)'
        }

        # 并行渲染的工作函数（每个线程独立创建 Figure，线程安全）
        def _render_one_group_png(mg_key: str, disp_title: str):
            kn = knee_notes.get(mg_key)
            c_path = files_map.get(mg_key, {}).get('c_in_pass')
            disp = title_map2.get(mg_key, disp_title)
            if not kn or not c_path or not os.path.exists(c_path):
                return (mg_key, None, disp, '', None, 0)

            thr = float(kn.get('threshold_abs_robust_z')) if kn.get('threshold_abs_robust_z') is not None else None
            if thr is None or not np.isfinite(thr):
                return (mg_key, None, disp, '', None, 0)

            cols_need = ['VARIANT_ID','CTRL_AAF','TOMMO_AAF','ROBUST_Z']
            pts = []  # (TOMMO_AAF, CTRL_AAF, TYPE)
            atcg = {"A","T","C","G"}
            # Prepare output TSV for selected VARIANT_IDs
            out_tsv_path = os.path.join(os.path.dirname(output_pdf), f"knee_variants.{mg_key}.tsv")
            n_written = 0
            wrote_header = False

            for chunk in pd.read_csv(c_path, sep='\t', usecols=cols_need, dtype='string',
                                     chunksize=500_000, engine='c', na_values=['nan','NaN','NA','.']):
                z = pd.to_numeric(chunk['ROBUST_Z'], errors='coerce')
                ca = pd.to_numeric(chunk['CTRL_AAF'], errors='coerce')
                ta = pd.to_numeric(chunk['TOMMO_AAF'], errors='coerce')
                mask = z.notna() & ca.notna() & ta.notna() & (z.abs() < thr)
                # Stream-write VARIANT_IDs to per-group TSV
                if mask.any():
                    vids = chunk.loc[mask, 'VARIANT_ID'].astype('string')
                    # Write header lazily
                    if not wrote_header:
                        with open(out_tsv_path, 'w') as fo:
                            fo.write('VARIANT_ID\n')
                        wrote_header = True
                    with open(out_tsv_path, 'a') as fo:
                        fo.write('\n'.join(vids.tolist()) + '\n')
                    n_written += int(mask.sum())
                if not mask.any():
                    continue
                sub = chunk.loc[mask, ['VARIANT_ID']].copy()
                sub['CTRL_AAF'] = ca[mask].to_numpy()
                sub['TOMMO_AAF'] = ta[mask].to_numpy()

                def _typ(vid):
                    if not isinstance(vid, str):
                        return 'InDel'
                    parts = vid.split(':', 3)
                    if len(parts) != 4:
                        return 'InDel'
                    ref, alt = parts[2], parts[3]
                    if ref in atcg and alt in atcg and len(ref)==1 and len(alt)==1:
                        return 'SNP'
                    return 'InDel'
                sub['TYPE'] = sub['VARIANT_ID'].apply(_typ)
                pts.append(sub[['TOMMO_AAF','CTRL_AAF','TYPE']])

            if pts:
                df_pts = pd.concat(pts, ignore_index=True)
            else:
                df_pts = pd.DataFrame(columns=['TOMMO_AAF','CTRL_AAF','TYPE'])

            # 绘制各自 PNG（单独 Figure，避免线程共享状态）
            fpng, axpng = plt.subplots(figsize=(6,6))
            if not df_pts.empty:
                is_snp = df_pts['TYPE'].eq('SNP')
                is_indel = df_pts['TYPE'].eq('InDel')
                if is_snp.any():
                    axpng.scatter(df_pts.loc[is_snp,'TOMMO_AAF'], df_pts.loc[is_snp,'CTRL_AAF'],
                                  s=16, alpha=0.7, linewidths=0, label='SNP', color='#000000', zorder=3)
                if is_indel.any():
                    axpng.scatter(df_pts.loc[is_indel,'TOMMO_AAF'], df_pts.loc[is_indel,'CTRL_AAF'],
                                  s=18, alpha=0.8, linewidths=0, label='InDel', color='#CC79A7', zorder=4)
            else:
                axpng.text(0.5, 0.5, 'No Data', ha='center', va='center')
            axpng.plot([0,1],[0,1], linestyle=(0,(4,2)), linewidth=1.0, color='#CC0000', zorder=10)
            axpng.minorticks_on()
            axpng.legend(frameon=False, fontsize=9)
            axpng.set_xlabel('TOMMO_AAF')
            axpng.set_ylabel('CTRL_AAF')
            axpng.set_xlim(0,1)
            axpng.set_ylim(0,1)
            axpng.set_box_aspect(1)
            axpng.xaxis.set_major_locator(MaxNLocator(nbins=5))
            axpng.yaxis.set_major_locator(MaxNLocator(nbins=5))
            axpng.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
            fpng.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.94)
            out_png2 = os.path.join(os.path.dirname(output_pdf), f"knee_scatter_{mg_key}.png")
            fpng.savefig(out_png2, dpi=png_dpi)
            plt.close(fpng)
            return (mg_key, out_png2, disp, f"|ROBUST_Z| < {thr:.2f}", out_tsv_path, n_written)

        # 提交三个任务并行执行
        groups_for_png = [('Rare Variant','rare'), ('Low Frequency Variant','lowfreq'), ('Common Variant','common')]
        results = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
            fut_map = {ex.submit(_render_one_group_png, mg_key, disp_title): mg_key
                       for (disp_title, mg_key) in groups_for_png}
            for fut in concurrent.futures.as_completed(fut_map):
                mg_key = fut_map[fut]
                try:
                    key, path, disp, sub, tsv_path, nsel = fut.result()
                except Exception:
                    key, path, disp, sub, tsv_path, nsel = mg_key, None, title_map2.get(mg_key, mg_key), '', None, 0
                results[key] = (path, disp, sub, tsv_path, nsel)

        # 按固定顺序收集
        ordered_pngs = []
        ordered_titles = []
        ordered_subs = []
        for _, mg_key in groups_for_png:
            path, disp, sub, tsv_path, nsel = results.get(mg_key, (None, title_map2.get(mg_key, mg_key), '', None, 0))
            ordered_pngs.append(path)
            ordered_titles.append(disp)
            ordered_subs.append(sub)
            # 把 TSV 信息并入 knee_notes，以便稍后写回 manifest
            if mg_key in knee_notes:
                knee_notes[mg_key]['selected_variants_tsv'] = tsv_path
                knee_notes[mg_key]['selected_variants_count'] = int(nsel) if nsel is not None else 0

        _compose_three_pngs_to_pdf_page(ordered_pngs, ordered_titles, ordered_subs)

    # 将拐点注释写回 manifest.json
    try:
        if knee_notes:
            manifest.setdefault('c_in_pass_knee', {})
            # 更新而不是覆盖其它组的数据
            for k, v in knee_notes.items():
                manifest['c_in_pass_knee'][k] = v
            with open(manifest_path, 'w') as f:
                json.dump(manifest, f, indent=2)
    except Exception:
        # 写回失败不影响绘图输出
        pass

    return output_pdf


def summarize_variants_filter_from_manifest(
    manifest_path: str,
    output_dir: Optional[str] = None,
    chunk_size: int = 500_000,
    max_workers: Optional[int] = None,
    keep_tmp: bool = False,
) -> dict:
    """
    函数名称：summarize_variants_filter_from_manifest
    ==============================================
    【功能】
    读取 `manifest.json`（由 build_grouped_variant_tables/summarize_c_in_pass_thresholds/plot_c_in_pass_threshold_tradeoff 产生）
    中的 `input` 表路径，按用户指定逻辑生成一个 **summary 表**：
      - 列：VARIANT_ID, GROUP, IN_TOMMO, PASS_TOMMO, PASS_GROUP_ROBUST_Z_FILTER, FILTER_STAT
      - GROUP：基于 CTRL_MAF，使用 assign_maf_group()：rare(<0.01)、lowfreq(0.01–0.05，含 0.05)、common(>0.05)
      - IN_TOMMO：来自 input 的 IN_TOMMO
      - PASS_TOMMO：来自 input 的 TOMMO_FILTER（NaN 保持 NaN；'PASS'→True；其它→False）
      - PASS_GROUP_ROBUST_Z_FILTER：仅当 IN_TOMMO 和 PASS_TOMMO 同时为 True 时才评估；
        根据该行 GROUP 选择 manifest['c_in_pass_knee'][GROUP]['selected_variants_tsv']，若 VARIANT_ID 出现在该 TSV 中 → True；
        若不在 → False；否则（前置条件不满足）填 NaN。
      - FILTER_STAT：基于三步判定生成的最终状态标签（便于后续一键过滤）：
        * Stat_1：IN_TOMMO==True 且 PASS_TOMMO==True 且 PASS_GROUP_ROBUST_Z_FILTER==True
        * Stat_2：IN_TOMMO==True 且 PASS_TOMMO==True 且 PASS_GROUP_ROBUST_Z_FILTER==False
        * Stat_3：IN_TOMMO==True 且 PASS_TOMMO==False（此时 PASS_GROUP_ROBUST_Z_FILTER 为 NaN）
        * Stat_4：IN_TOMMO==False（此时 PASS_TOMMO 与 PASS_GROUP_ROBUST_Z_FILTER 分别为 NaN）

    【实现】
    - 使用 pandas 分块读取大表（chunk_size 可配），并**并行**处理每个分块（ThreadPoolExecutor）。
    - 结果以流式追加写出，避免占用大量内存。
    - 最终额外输出一个 2 维统计：按 (GROUP, FILTER_STAT) 计数，并在统计表中新增两列：
      `GROUP_TOTAL`（该 GROUP 的总数）与 `PERCENT`（COUNT/该 GROUP 总数，百分比，保留 4 位小数）。

    Stat_1~Stat_4 说明：
      - FILTER_STAT：基于三步判定生成的最终状态标签（便于后续一键过滤）：
        * Stat_1：IN_TOMMO==True 且 PASS_TOMMO==True 且 PASS_GROUP_ROBUST_Z_FILTER==True
        * Stat_2：IN_TOMMO==True 且 PASS_TOMMO==True 且 PASS_GROUP_ROBUST_Z_FILTER==False
        * Stat_3：IN_TOMMO==True 且 PASS_TOMMO==False（此时 PASS_GROUP_ROBUST_Z_FILTER 为 NaN）
        * Stat_4：IN_TOMMO==False（此时 PASS_TOMMO 与 PASS_GROUP_ROBUST_Z_FILTER 分别为 NaN）

    参数
    ----
    manifest_path : str
        manifest.json 路径。
    output_dir : Optional[str]
        输出目录；默认与 manifest 同目录。
    chunk_size : int
        pandas 分块大小（默认 1,000,000 行）。
    max_workers : Optional[int]
        并行工作线程数；默认 min(8, CPU)。
    keep_tmp : bool
        是否保留可能的临时文件（当前函数仅打印日志，不创建大临时文件）。

    返回
    ----
    dict
        包含输出路径与计数字典的简要信息。
    """
    import os
    import json
    import math
    import time
    import threading
    import concurrent.futures
    from collections import defaultdict

    import pandas as pd
    import numpy as np

    t0 = time.time()
    manifest_path = os.path.abspath(manifest_path)
    with open(manifest_path, 'r') as f:
        manifest = json.load(f)

    # —— 路径解析 ——
    in_path = manifest.get('input')
    if not in_path:
        raise ValueError("manifest.json 缺少 'input' 字段")
    if output_dir is None:
        output_dir = os.getcwd()
    os.makedirs(output_dir, exist_ok=True)

    # —— 读取每组 selected_variants_tsv → set ——
    knee = manifest.get('c_in_pass_knee', {}) or {}
    sel_sets = {}
    for gkey in ('rare', 'lowfreq', 'common'):
        path = (knee.get(gkey) or {}).get('selected_variants_tsv')
        s = set()
        if path and os.path.exists(path):
            try:
                for chunk in pd.read_csv(path, sep='\t', usecols=['VARIANT_ID'], dtype='string', chunksize=500_000):
                    s.update(v for v in chunk['VARIANT_ID'].dropna().tolist())
            except Exception:
                # 尝试无表头读取
                for chunk in pd.read_csv(path, sep='\t', header=0, names=['VARIANT_ID'], dtype='string', chunksize=500_000):
                    s.update(v for v in chunk['VARIANT_ID'].dropna().tolist())
        sel_sets[gkey] = s

    # —— 输出文件 ——
    base_name = os.path.basename(in_path)
    out_summary = os.path.join(output_dir, f"{base_name}.summary_filter.tsv")
    out_stats = os.path.join(output_dir, f"{base_name}.summary_filter.counts.tsv")

    # 写表头
    with open(out_summary, 'w') as fo:
        fo.write('\t'.join(['VARIANT_ID', 'GROUP', 'IN_TOMMO', 'PASS_TOMMO', 'PASS_GROUP_ROBUST_Z_FILTER', 'FILTER_STAT']) + '\n')

    # 统计计数（分块累积）
    counts = defaultdict(int)  # key: (GROUP, FILTER_STAT)
    counts_lock = threading.Lock()

    def _pass_tommo_from_filter(s: pd.Series) -> pd.Series:
        # NaN 保持 NaN；PASS→True；其它→False
        tf = s.astype('string')
        is_na = tf.isna() | tf.str.lower().isin(['nan', 'na', '.'])
        out = pd.Series(pd.NA, index=tf.index, dtype='object')
        out.loc[~is_na & (tf.str.upper() == 'PASS')] = True
        out.loc[~is_na & (tf.str.upper() != 'PASS')] = False
        return out

    # 单块处理函数（返回处理后的 DataFrame 以及局部计数）
    def _process_chunk(chunk: pd.DataFrame):
        sub = pd.DataFrame(index=chunk.index)
        sub['VARIANT_ID'] = chunk['VARIANT_ID'].astype('string')
        sub['GROUP'] = assign_maf_group(chunk['CTRL_MAF'])
        # IN_TOMMO → 归一化为布尔/NA
        s_in = chunk['IN_TOMMO']
        sub['IN_TOMMO'] = s_in.map(lambda x: True if x in (True, 1, '1', 'True', 'TRUE') else (False if x in (False, 0, '0', 'False', 'FALSE') else pd.NA))
        # PASS_TOMMO
        sub['PASS_TOMMO'] = _pass_tommo_from_filter(chunk['TOMMO_FILTER'])
        # 先置 NaN
        sub['PASS_GROUP_ROBUST_Z_FILTER'] = pd.Series([pd.NA]*len(sub), index=sub.index, dtype='object')

        # 仅对 (IN_TOMMO==True & PASS_TOMMO==True) 的行进行 membership 检查
        mask_eval = (sub['IN_TOMMO'] == True) & (sub['PASS_TOMMO'] == True) & sub['GROUP'].notna()
        if mask_eval.any():
            for gkey in ('rare', 'lowfreq', 'common'):
                idx = mask_eval & (sub['GROUP'] == gkey)
                if idx.any():
                    vids = sub.loc[idx, 'VARIANT_ID']
                    hit = vids.isin(sel_sets.get(gkey, set()))
                    sub.loc[idx, 'PASS_GROUP_ROBUST_Z_FILTER'] = hit.astype('bool').astype('object')

        # 生成 FILTER_STAT（Stat_1~4）
        stat = pd.Series(pd.NA, index=sub.index, dtype='string')
        it = sub['IN_TOMMO']
        pt = sub['PASS_TOMMO']
        pr = sub['PASS_GROUP_ROBUST_Z_FILTER']
        # Stat_1: in_tommo & pass_tommo & pass_robust
        m1 = (it == True) & (pt == True) & (pr == True)
        # Stat_2: in_tommo & pass_tommo & fail_robust
        m2 = (it == True) & (pt == True) & (pr == False)
        # Stat_3: in_tommo & nonpass_tommo
        m3 = (it == True) & (pt == False)
        # Stat_4: not in tommo
        m4 = (it == False)
        stat.loc[m1] = 'Stat_1'
        stat.loc[m2] = 'Stat_2'
        stat.loc[m3] = 'Stat_3'
        stat.loc[m4] = 'Stat_4'
        sub['FILTER_STAT'] = stat

        # 统计：仅按 (GROUP, FILTER_STAT) 计数
        local_counts = defaultdict(int)
        g_vals = sub['GROUP'].astype('string').where(sub['GROUP'].notna(), other='NA')
        fs_vals = sub['FILTER_STAT'].astype('string').where(sub['FILTER_STAT'].notna(), other='NA')
        for g, fs in zip(g_vals.tolist(), fs_vals.tolist()):
            local_counts[(g, fs)] += 1

        return sub[['VARIANT_ID','GROUP','IN_TOMMO','PASS_TOMMO','PASS_GROUP_ROBUST_Z_FILTER','FILTER_STAT']], local_counts

    # 并行设置
    if max_workers is None:
        try:
            import multiprocessing as _mp
            max_workers = max(1, min(8, _mp.cpu_count()))
        except Exception:
            max_workers = 4

    total_rows = 0
    total_chunks = 0

    # 写入锁，保证多线程安全地追加到同一文件
    write_lock = threading.Lock()

    print(f"[summary] 读取: {in_path}")
    print(f"[summary] 输出: {out_summary}")
    print(f"[summary] 统计输出: {out_stats}")

    usecols = ['VARIANT_ID','CTRL_MAF','IN_TOMMO','TOMMO_FILTER']
    print(f"[summary] 并行 worker 数: {max_workers}")

    # --- 新并行流水线：producer (reader) -> N workers (process & write part files) -> merge ---
    import queue
    import uuid
    try:
        import shutil
    except ImportError:
        shutil = None

    q: "queue.Queue[pd.DataFrame]" = queue.Queue(maxsize=max(2, max_workers * 3))
    worker_counts: dict = {}
    part_paths: dict = {}

    def _worker_loop(wid: int):
        part_path = os.path.join(output_dir, f"{base_name}.summary_filter.part{wid:02d}.tsv")
        part_paths[wid] = part_path
        # 确保空文件存在（无表头；主文件已写表头）
        open(part_path, 'wb').close()
        local = defaultdict(int)
        n_batches = 0
        while True:
            chunk = q.get()
            if chunk is None:  # 哨兵
                q.task_done()
                break
            df_part, lc = _process_chunk(chunk)
            # 直接写入该 worker 的 part 文件，避免全局写锁争用
            df_part.to_csv(part_path, sep='\t', header=False, index=False, mode='a', na_rep='nan')
            # 合并计数
            for k, v in lc.items():
                local[k] += v
            n_batches += 1
            q.task_done()
        worker_counts[wid] = local
        _log(f"[worker-{wid}] 退出；处理分块 {n_batches} 个 → {os.path.basename(part_path)}")

    # 启动 worker 线程
    threads = []
    for wid in range(max_workers):
        t = threading.Thread(target=_worker_loop, args=(wid,), daemon=True)
        t.start()
        threads.append(t)
    def _log(msg): print(msg)
    _log(f"[summary] 已启动 {len(threads)} 个 worker 线程进行并行处理")

    # 生产者：读取分块并投递到队列
    reader = pd.read_csv(
        in_path, sep='\t', usecols=usecols, dtype='string',
        chunksize=chunk_size, engine='c', na_values=['nan','NaN','NA','.']
    )
    for chunk in reader:
        total_chunks += 1
        total_rows += len(chunk)
        q.put(chunk)  # 阻塞式，结合 queue 大小形成背压
        if total_chunks % 10 == 0:
            print(f"[summary] 已投递 {total_chunks} 个分块，共 {total_rows:,} 行 … 队列积压 {q.qsize()} / {q.maxsize}")

    # 投递哨兵，通知所有 worker 退出
    for _ in range(max_workers):
        q.put(None) # type: ignore
    q.join()  # 等待所有任务完成

    # 汇总计数
    for wid, lc in worker_counts.items():
        with counts_lock:
            for k, v in lc.items():
                counts[k] += v

    print(f"[summary] 完成。总分块数: {total_chunks}, 总行数: {total_rows:,}。")

    # 合并所有 part 文件到最终 out_summary（已写过表头，这里仅追加数据行）
    if shutil is None:
        import shutil
    with open(out_summary, 'a') as fout:
        for wid in sorted(part_paths.keys()):
            p = part_paths[wid]
            if not os.path.exists(p) or os.path.getsize(p) == 0:
                continue
            with open(p, 'r') as fin:
                shutil.copyfileobj(fin, fout)
    # 清理 part 文件
    for p in part_paths.values():
        try:
            os.remove(p)
        except Exception:
            pass

    # 计算每个 GROUP 的总数（GROUP_TOTAL），作为 COUNT 的分母
    import csv as _csv
    stat_notes = {
        'Stat_1': 'IN_TOMMO==True & PASS_TOMMO==True & PASS_GROUP_ROBUST_Z_FILTER==True',
        'Stat_2': 'IN_TOMMO==True & PASS_TOMMO==True & PASS_GROUP_ROBUST_Z_FILTER==False',
        'Stat_3': 'IN_TOMMO==True & PASS_TOMMO==False',
        'Stat_4': 'IN_TOMMO==False',
        'NA':     '状态不适用/未知（例如缺失值）'
    }
    # 计算每个 GROUP 的总数（GROUP_TOTAL），作为 COUNT 的分母
    group_totals = {}
    for (g, fs), c in counts.items():
        if g not in group_totals:
            group_totals[g] = 0
        group_totals[g] += c
    with open(out_stats, 'w', newline='') as fo:
        # 注释行（以 # 开头）
        fo.write('# FILTER_STAT 含义:\n')
        for k in ('Stat_1','Stat_2','Stat_3','Stat_4'):
            fo.write(f"#   {k}: {stat_notes[k]}\n")
        fo.write('#   NA: 未能归类的记录或缺失\n')
        # 表头（新增 GROUP_TOTAL 与 PERCENT 两列）
        w = _csv.writer(fo, delimiter='\t')
        w.writerow(['GROUP','FILTER_STAT','COUNT','GROUP_TOTAL','PERCENT'])
        # 数据行：附带各 GROUP 的总数与百分比（小数点后4位）
        for (g, fs), c in sorted(counts.items()):
            gt = group_totals.get(g, 0)
            if gt:
                pct_str = f"{(c / gt) * 100:.4f}%"
            else:
                pct_str = 'nan'
            w.writerow([g, fs, c, gt, pct_str])

    dt = time.time() - t0
    print(f"[summary] 耗时: {dt/60:.2f} min; 输出: {out_summary}; 统计: {out_stats}")

    return out_summary # type: ignore


def filter_variants_by_group_and_stat(
    out_summary: str,
    config_json: str,
    output_dir: Optional[str] = None,
    chunk_size: int = 1_000_000,
    max_workers: Optional[int] = None,
) -> Dict[str, str]:
    """
    函数名称：filter_variants_by_group_and_stat
    ========================================
    【功能】
    基于 summarize_variants_filter_from_manifest 产出的 `out_summary` 表，按照 **配置 JSON** 中为
    各 GROUP 指定的 FILTER_STAT 白名单，筛选出对应的 VARIANT_ID，并分别生成 3 个文件：
      - rare 组：<basename>.rare.selected_variants.tsv
      - lowfreq 组：<basename>.lowfreq.selected_variants.tsv
      - common 组：<basename>.common.selected_variants.tsv
    文件第一行写入表头 `VARIANT_ID`（此表头行不计入变体列表，满足“第一行这个 VARIANT_ID 不用纳入进来”的要求）。

    【输入】
    - out_summary : str
        `*.summary_filter.tsv`（至少包含列：VARIANT_ID, GROUP, FILTER_STAT）。
    - config_json : str
        JSON 配置文件路径，指定每个 GROUP 允许保留的 FILTER_STAT 集合。
        示例（保存为 config.json）：
        {
          "rare":    ["Stat_1", "Stat_2", "Stat_3", "Stat_4"],
          "lowfreq": ["Stat_1"],
          "common":  ["Stat_1"]
        }
    - output_dir : Optional[str]
        输出目录；默认与 out_summary 同目录。
    - chunk_size : int
        Pandas 分块大小（默认 1,000,000 行）。
    - max_workers : Optional[int]
        并行 worker 数（默认 min(8, CPU)）。

    【输出】
    - 返回 dict：{'rare': path_rare, 'lowfreq': path_low, 'common': path_common}
      每个文件均只包含一列 `VARIANT_ID`，已按 **CHROM: chr1..chr22 主序，POS 次序** 排好。

    【排序规则】
    - VARIANT_ID 形如 CHROM:POS:REF:ALT，按 CHROM 的自然顺序（chr1..chr22）分桶，并在每个桶内按 POS 升序。
    - 若 CHROM 不在 {chr1..chr22}，将被忽略（仅输出常染 1–22）。

    【实现要点（面向大文件内存安全）】
    1) 读入分块（chunksize），仅保留所需列，并**并行**处理过滤逻辑；
    2) 对每个 GROUP，采用“每条染色体一个临时文件”的外排序策略：临时写入 `POS\\tVARIANT_ID`；
    3) 末尾再对每个染色体临时文件进行排序并合并写回最终文件；
    4) 最终文件第一行写入 `VARIANT_ID` 表头。
    """
    import os, json, shutil, threading, tempfile, queue, concurrent.futures
    from typing import Dict, Tuple, List
    import pandas as pd

    # ---- 路径与输出 ----
    out_summary = os.path.abspath(out_summary)
    if output_dir is None:
        output_dir = os.path.dirname(out_summary) or os.getcwd()
    os.makedirs(output_dir, exist_ok=True)

    base = os.path.basename(out_summary)
    out_paths = {
        'rare':    os.path.join(output_dir, f"{base}.rare.selected_variants.tsv"),
        'lowfreq': os.path.join(output_dir, f"{base}.lowfreq.selected_variants.tsv"),
        'common':  os.path.join(output_dir, f"{base}.common.selected_variants.tsv"),
    }

    # ---- 读取 JSON 配置 ----
    with open(config_json, 'r') as f:
        cfg = json.load(f)
    allow_stats = {
        'rare':    set(map(str, cfg.get('rare',    []))),
        'lowfreq': set(map(str, cfg.get('lowfreq', []))),
        'common':  set(map(str, cfg.get('common',  []))),
    }

    # ---- 染色体桶（chr1..chr22）与临时文件结构 ----
    chrom_buckets = [f"chr{i}" for i in range(1, 23)]
    workdir = tempfile.mkdtemp(prefix="filter_gstat_")
    tmp_paths: Dict[str, Dict[str, str]] = {g: {} for g in out_paths}
    locks: Dict[str, Dict[str, threading.Lock]] = {g: {} for g in out_paths}
    for g in out_paths:
        gdir = os.path.join(workdir, g)
        os.makedirs(gdir, exist_ok=True)
        for ch in chrom_buckets:
            p = os.path.join(gdir, f"{ch}.pos_vid.tmp")
            open(p, 'wb').close()
            tmp_paths[g][ch] = p
            locks[g][ch] = threading.Lock()

    # ---- 并行设置 ----
    if max_workers is None:
        try:
            import multiprocessing as _mp
            max_workers = max(1, min(8, _mp.cpu_count()))
        except Exception:
            max_workers = 4

    def _parse_vid(vid: str) -> Tuple[str, int]:
        parts = vid.split(':', 3)
        if len(parts) != 4:
            raise ValueError("bad vid")
        chrom, pos = parts[0], int(parts[1])
        return chrom, pos

    def _process_chunk(chunk: pd.DataFrame) -> int:
        # 仅保留必要列
        c = chunk[['VARIANT_ID', 'GROUP', 'FILTER_STAT']].copy()
        c['GROUP'] = c['GROUP'].astype('string')
        c['FILTER_STAT'] = c['FILTER_STAT'].astype('string')
        c = c[c['GROUP'].isin(['rare','lowfreq','common'])]
        if c.empty:
            return 0

        out_frames = []
        for g in ('rare','lowfreq','common'):
            allow = allow_stats.get(g, set())
            if not allow:
                continue
            sub = c[(c['GROUP'] == g) & (c['FILTER_STAT'].isin(allow))][['VARIANT_ID']]
            if not sub.empty:
                sub['GROUP'] = g
                out_frames.append(sub)
        if not out_frames:
            return 0

        cc = pd.concat(out_frames, ignore_index=True)
        vids = cc['VARIANT_ID'].astype('string').tolist()
        groups = cc['GROUP'].tolist()

        kept = 0
        for vid, g in zip(vids, groups):
            try:
                chrom, pos = _parse_vid(vid)
            except Exception:
                continue
            if chrom not in chrom_buckets:
                continue
            with locks[g][chrom]:
                with open(tmp_paths[g][chrom], 'a') as fo:
                    fo.write(f"{pos}\t{vid}\n")
            kept += 1
        return kept

    # ---- 生产者-消费者并发框架（队列） ----
    q: "queue.Queue[pd.DataFrame]" = queue.Queue(maxsize=max(2, max_workers*3))
    stop = object()
    stats = {'chunks': 0, 'rows': 0, 'kept': 0}

    def _worker():
        while True:
            obj = q.get()
            if obj is stop:
                q.task_done()
                break
            kept_local = _process_chunk(obj)
            stats['kept'] += kept_local
            q.task_done()

    workers = []
    for _ in range(max_workers):
        t = threading.Thread(target=_worker, daemon=True)
        t.start()
        workers.append(t)

    usecols = ['VARIANT_ID','GROUP','FILTER_STAT']
    for chunk in pd.read_csv(out_summary, sep='\t', usecols=usecols, dtype='string',
                             chunksize=chunk_size, engine='c',
                             na_values=['nan','NaN','NA','.']):
        stats['chunks'] += 1
        stats['rows'] += len(chunk)
        q.put(chunk)

    for _ in workers:
        q.put(stop) # type: ignore
    q.join()
    for t in workers:
        t.join()

    # ---- 合并：按 chr1..chr22 + POS 排序并写最终文件（无表头，仅 VARIANT_ID，每行一个） ----
    for g, outp in out_paths.items():
        # 初始化空文件，不写入表头
        open(outp, 'w').close()
        total_written = 0
        for ch in chrom_buckets:
            tmp = tmp_paths[g][ch]
            if not os.path.exists(tmp) or os.path.getsize(tmp) == 0:
                continue
            pos_vids: List[Tuple[int, str]] = []
            with open(tmp, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        pos_s, vid = line.split('\t', 1)
                        pos = int(pos_s)
                    except Exception:
                        continue
                    pos_vids.append((pos, vid))
            if pos_vids:
                pos_vids.sort(key=lambda x: x[0])
                with open(outp, 'a') as fo:
                    fo.write('\n'.join(v for _, v in pos_vids) + '\n')
                total_written += len(pos_vids)
        print(f"[filter] {g}: wrote {total_written:,} variants -> {outp}")

    # ---- 清理临时目录 ----
    try:
        shutil.rmtree(workdir)
    except Exception:
        pass

    return out_paths


from typing import Dict, Optional, Tuple, List
import subprocess

def subset_plink_by_selected_variants(
    bed_prefix: str,
    out_prefix: str,
    selected_paths: Dict[str, str],
    threads: int = 8,
    merge_low_common: bool = True,
    plink2_path: str = "/home/b/b37974/plink2",
) -> Dict[str, str]:
    """
    函数名称：subset_plink_by_selected_variants
    ========================================
    【功能】
    基于 `filter_variants_by_group_and_stat()` 的输出（各 GROUP 的 `VARIANT_ID` 列表），
    使用 plink2 对给定 bed_prefix 的基因型数据进行子集提取（subset）。

    【输入】
    - bed_prefix : str
        plink 二进制基因型前缀（.bed/.bim/.fam）。
    - out_prefix : str
        plink2 输出前缀（函数将基于此追加 `.rare`、`.lowfreq_common` 或 `.lowfreq`/`.common`）。
    - selected_paths : Dict[str, str]
        `filter_variants_by_group_and_stat()` 的返回字典，形如：
        {
          'rare': '/path/to/...rare.selected_variants.tsv',
          'lowfreq': '/path/to/...lowfreq.selected_variants.tsv',
          'common': '/path/to/...common.selected_variants.tsv'
        }
        每个文件仅一列 `VARIANT_ID`，且可能非常大。
    - threads : int = 8
        plink2 的 `--threads` 数量，同时也用于本函数的并发处理数量。
    - merge_low_common : bool = True
        若 True：合并 lowfreq 与 common 的变体列表，进行一次 subset（输出后缀 `.lowfreq_common`）。
        若 False：分别对 rare/lowfreq/common 三组各做一次 subset。
    - plink2_path : str = '/home/b/b37974/plink2'
        plink2 可执行文件路径。

    【实现要点】
    1) 读取 lowfreq 与 common 的 `VARIANT_ID`（当 merge_low_common=True），进行**外排序**：
       - 兼容 'chr1' 与 '1' 的染色体标记；仅保留 1..22 常染色体；
       - 采用“每条染色体一个临时文件”的分桶策略，并行处理/写入；
       - 最终在主线程按 POS 升序合并各桶，写成单列 `VARIANT_ID` 文件；
    2) 调用 plink2：`plink2 --bfile <bed_prefix> --extract <list.tsv> --make-bed --out <out_prefix.suffix> --threads <threads>`；
    3) 打印清晰日志（输入/输出、每步耗时、行数、可能的空列表提示）。

    【返回】
    - 返回字典，键为输出数据集的逻辑名（'rare'、'lowfreq_common' 或 'lowfreq'、'common'），
      值为对应的 plink2 输出前缀路径（即 `.bed/.bim/.fam` 的公共前缀）。
    """
    import os
    import time
    import tempfile
    import math
    import threading
    import queue

    t0 = time.time()
    bed_prefix = os.path.abspath(bed_prefix)
    out_prefix = os.path.abspath(out_prefix)
    plink2_path = os.path.abspath(plink2_path)

    if not os.path.exists(bed_prefix + '.bed'):
        raise FileNotFoundError(f"找不到输入 .bed: {bed_prefix}.bed")
    if not os.path.exists(bed_prefix + '.bim'):
        raise FileNotFoundError(f"找不到输入 .bim: {bed_prefix}.bim")
    if not os.path.exists(bed_prefix + '.fam'):
        raise FileNotFoundError(f"找不到输入 .fam: {bed_prefix}.fam")

    # --- 工具函数：解析 VARIANT_ID 为 (chrom, pos) 并做标准化 ---
    def _parse_vid(vid: str) -> Optional[Tuple[str, int]]:
        if not isinstance(vid, str) or not vid:
            return None
        p = vid.split(':', 3)
        if len(p) != 4:
            return None
        chrom = p[0]
        if chrom.lower().startswith('chr'):
            chrom = chrom[3:]
        try:
            pos = int(p[1])
        except Exception:
            return None
        # 仅保留 1..22
        if chrom.isdigit():
            cnum = int(chrom)
            if 1 <= cnum <= 22:
                return (f"chr{cnum}", pos)
        return None

    # --- 内部：将一个变体列表文件分桶到临时文件（每条染色体一个 tmp） ---
    def _bucketize_variants(list_path: str, tmp_dir: str, label: str) -> Dict[str, str]:
        chroms = [f"chr{i}" for i in range(1, 23)]
        paths = {ch: os.path.join(tmp_dir, f"{label}.{ch}.pos_vid.tmp") for ch in chroms}
        for p in paths.values():
            open(p, 'wb').close()

        # 生产者-消费者：分块读取 & 并发写桶
        q: "queue.Queue[List[str]]" = queue.Queue(maxsize=max(2, threads*3))
        stop = object()
        locks = {ch: threading.Lock() for ch in chroms}

        def _worker():
            while True:
                obj = q.get()
                if obj is stop:
                    q.task_done(); break
                for line in obj:
                    vid = line.strip()
                    if not vid or vid == 'VARIANT_ID':
                        continue
                    parsed = _parse_vid(vid)
                    if parsed is None:
                        continue
                    ch, pos = parsed
                    with locks[ch]:
                        with open(paths[ch], 'a') as fo:
                            fo.write(f"{pos}\t{vid}\n")
                q.task_done()

        workers = []
        for _ in range(max(1, threads)):
            t = threading.Thread(target=_worker, daemon=True)
            t.start(); workers.append(t)

        total = 0
        with open(list_path, 'r') as f:
            buf: List[str] = []
            for line in f:
                buf.append(line)
                if len(buf) >= 200000:  # ~200k 行一批，避免内存过大
                    q.put(buf); buf = []
            if buf:
                q.put(buf)
        for _ in workers:
            q.put(stop) # type: ignore
        q.join()
        for t in workers:
            t.join()
        return paths

    # --- 合并 lowfreq + common 列表（如需） ---
    outputs: Dict[str, str] = {}
    tmp_root = tempfile.mkdtemp(prefix="plink_subset_")
    try:
        if merge_low_common:
            low_path = selected_paths.get('lowfreq')
            com_path = selected_paths.get('common')
            if not low_path and not com_path:
                print("[subset] ⚠️ 未提供 lowfreq/common 变体列表，跳过合并集。")
            else:
                print("[subset] 合并 lowfreq + common 变体列表，并进行排序（chr1..22，POS 升序）...")
                # 先各自分桶
                merge_dir = os.path.join(tmp_root, 'merge')
                os.makedirs(merge_dir, exist_ok=True)
                paths_low = _bucketize_variants(low_path, merge_dir, 'lowfreq') if low_path else {}
                paths_com = _bucketize_variants(com_path, merge_dir, 'common') if com_path else {}

                # 合并每条染色体桶 & 排序 & 写入最终合并文件
                merged_list = os.path.join(tmp_root, 'lowfreq_common.merged.sorted.tsv')
                open(merged_list, 'w').close()
                written = 0
                for i in range(1, 23):
                    ch = f"chr{i}"
                    candidates = [p for p in [paths_low.get(ch), paths_com.get(ch)] if p and os.path.exists(p)]
                    pos_vids: List[Tuple[int, str]] = []
                    for p in candidates:
                        with open(p, 'r') as f:
                            for line in f:
                                line = line.strip()
                                if not line:
                                    continue
                                pos_s, vid = line.split('\t', 1)
                                try:
                                    pos = int(pos_s)
                                except Exception:
                                    continue
                                pos_vids.append((pos, vid))
                    if pos_vids:
                        pos_vids.sort(key=lambda x: x[0])
                        with open(merged_list, 'a') as fo:
                            fo.write('\n'.join(v for _, v in pos_vids) + '\n')
                        written += len(pos_vids)
                print(f"[subset] 合并后总变体数: {written:,}")

                # 调用 plink2 进行合并集的 subset
                out_lc = f"{out_prefix}.lowfreq_common"
                cmd = [
                    plink2_path,
                    '--bfile', bed_prefix,
                    '--extract', merged_list,
                    '--make-bed',
                    '--threads', str(int(max(1, threads))),
                    '--out', out_lc,
                ]
                print("[subset] 运行:", ' '.join(cmd))
                try:
                    subprocess.run(cmd, check=True)
                except subprocess.CalledProcessError as e:
                    raise RuntimeError(f"plink2 子集（lowfreq_common）失败: {e}")
                outputs['lowfreq_common'] = out_lc

        else:
            # 分别对子三组进行 subset
            for g in ('lowfreq','common'):
                path = selected_paths.get(g)
                if not path:
                    print(f"[subset] ⚠️ 未提供 {g} 变体列表，跳过该组。")
                    continue
                outg = f"{out_prefix}.{g}"
                cmd = [
                    plink2_path,
                    '--bfile', bed_prefix,
                    '--extract', os.path.abspath(path),
                    '--make-bed',
                    '--threads', str(int(max(1, threads))),
                    '--out', outg,
                ]
                print("[subset] 运行:", ' '.join(cmd))
                try:
                    subprocess.run(cmd, check=True)
                except subprocess.CalledProcessError as e:
                    raise RuntimeError(f"plink2 子集（{g}）失败: {e}")
                outputs[g] = outg

        # 无论是否合并，都要处理 rare
        rare_path = selected_paths.get('rare')
        if not rare_path:
            print("[subset] ⚠️ 未提供 rare 变体列表，跳过 rare 子集。")
        else:
            out_r = f"{out_prefix}.rare"
            cmd = [
                plink2_path,
                '--bfile', bed_prefix,
                '--extract', os.path.abspath(rare_path),
                '--make-bed',
                '--threads', str(int(max(1, threads))),
                '--out', out_r,
            ]
            print("[subset] 运行:", ' '.join(cmd))
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"plink2 子集（rare）失败: {e}")
            outputs['rare'] = out_r

    finally:
        # 清理临时目录
        try:
            import shutil
            shutil.rmtree(tmp_root)
        except Exception:
            pass

    dt = time.time() - t0
    print(f"[subset] 完成。耗时 {dt/60:.2f} 分钟；输出前缀：{outputs}")
    return outputs

