import subprocess
import os
import tempfile
import csv
import math
import uuid
from typing import Optional
import concurrent.futures
import pandas as pd
import numpy as np
from itertools import islice
from collections import defaultdict
from typing import Dict, Tuple, Iterable


def run_plink2_variant_qc_with_tommo(
    variant_qc_summary: str,
    tommo_vcf_path: str,
    output_path: Optional[str] = None,
    bcftools_path: str = "bcftools",
    threads: int = 8,
    chunk_size: int = 500_000,
    max_workers: Optional[int] = None,
    keep_temp: bool = False,
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
    keep_temp : bool
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

    def _run_cmd(cmd_list: Iterable[str], stdout_path: str = None) -> int:
        if stdout_path is None:
            proc = subprocess.run(cmd_list, check=False)
        else:
            with open(stdout_path, "wb") as fo:
                proc = subprocess.run(cmd_list, check=False, stdout=fo)
        return proc.returncode

    # 1) 临时工作目录 & 输出路径
    workdir_obj = tempfile.TemporaryDirectory(prefix="tommo_merge_")
    workdir = workdir_obj.name
    if output_path is None:
        base, _ = os.path.splitext(variant_qc_summary)
        output_path = base + ".variant_qc_with_tommo.tsv"
    _progress(f"临时目录: {workdir}")

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
        # 注意：部分环境中的 `bcftools query` 不支持 --threads。
        # 为了同时获得多线程解压能力与兼容性，这里：
        #   threads>1 时：使用管道 `bcftools view --threads N -R ... -Ou VCF | bcftools query -f FMT`
        #   否则：直接 `bcftools query -R ... -f FMT VCF`。
        if isinstance(threads, int) and threads > 1:
            fmt_q = shlex.quote(fmt)
            cmdline = (
                f"{shlex.quote(bcftools_path)} view --threads {threads} -R {shlex.quote(region_file)} "
                f"-Ou {shlex.quote(tommo_vcf_path)} | "
                f"{shlex.quote(bcftools_path)} query -f {fmt_q} > {shlex.quote(out_map)}"
            )
            cmd = ["bash", "-lc", cmdline]
            ret = _run_cmd(cmd)
        else:
            cmd = [
                bcftools_path, "query",
                "-R", region_file,
                "-f", fmt,
                tommo_vcf_path,
            ]
            ret = _run_cmd(cmd, stdout_path=out_map)
        if ret != 0:
            raise RuntimeError(f"bcftools query 失败: 染色体 {chrom}")
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
            _progress(f"[diag] {chrom}: 映射子集载入 {len(out):,} 条，命中 {hit_cnt:,} / {len(needed_keys):,}")
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
    if keep_temp:
        _progress(f"保留临时目录: {workdir}")
    else:
        try:
            workdir_obj.cleanup()
        except Exception:
            pass

    _progress(f"完成。输出: {output_path}")
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
):
    """
    函数名称：plot_tommo_panel_compare_pdf
    =====================================
    【功能】
    读取 `run_plink2_variant_qc_with_tommo` 生成的结果（*.variant_qc_with_tommo.tsv），
    并输出一个包含 4 页内容的 PDF。每一页包含 3 个分图（Rare/LowFreq/Common 三组）。

    【输入】
    - variant_qc_with_tommo: str
        由 run_plink2_variant_qc_with_tommo 产出的结果表路径。
    - output_pdf: Optional[str]
        输出 PDF 文件路径；默认与输入同名但后缀为 `.panel_compare.pdf`。
    - max_points: Optional[int]
        为了绘图速度，可对每个分组随机抽样此数量的点（None 表示不抽样）。
    - png_dpi: int = 600
        第 2～3 页的每个分图先以高清 PNG（位图）保存，再插入到 PDF 中；此参数用于调节 PNG 的分辨率（DPI）。
    - page23_figsize: tuple = (14, 5)
        第 2～3 页组合页面的整体尺寸（英寸），可调大以放大每个分图在 PDF 页内的可视大小。
    - page4_figsize: tuple = (14, 5)
        第 4 页（直方图页）的整体画布尺寸（英寸）。建议与第 2～3 页一致以保持版式统一。
    - page4_wspace: float = 0.10
        第 4 页三幅直方图之间的水平间距（0~1，越大间距越大）。

    【分组定义（基于 CTRL_MAF）】
    - Rare Variant (<0.01)
    - Low Frequency Variant (0.01 ~ 0.05)
    - Common Variant (>0.05)

    【页面说明】
    1. 第 1 页：3 个小表（每个分组一个表），列：
        - Total Variants（该分组内总计：按定义可用的 CTRL_MAF）
        - In ToMMo（IN_TOMMO 为 True 的个数）
        - Pass ToMMo（TOMMO_FILTER == 'PASS' 的个数）
        - Percent（相对于该组 Total 的百分比，保留 1 位小数）
       注：Count 列使用 3 位逗号分隔格式。
    2. 第 2 页：3 个散点图（每个分组一个图）。X=TOMMO_AAF，Y=CTRL_AAF。
        颜色按 TOMMO_FILTER 是否 PASS（PASS / 非 PASS）。
        **实现细节**：为了避免载入过多 artist 导致 PDF 体积大与渲染缓慢，
        每个分图先渲染为单独的高清 PNG（由 `png_dpi` 控制），再以位图插入 PDF。
    3. 第 3 页：与第 2 页相同，但**仅保留 TOMMO_FILTER=='PASS'** 的点（同样以 PNG 先渲染后插入）。
    4. 第 4 页：3 个直方图（仅 PASS），绘制 CTRL_AAF - TOMMO_AAF 的分布（仍为矢量）。

    【注意】
    - 自动将相关列转为数值类型并丢弃无法转换的记录（NaN）。
    - 频率轴范围设为 [0,1]（直方图除外）。
    - 支持可选抽样以提升绘图速度。
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
    else:
        pass_color = '#1f77b4'
        nonpass_color = '#ff7f0e'
        refline_color = 'red'

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
                                    'CTRL_MAF': 'float32',
                                    'CTRL_AAF': 'float32',
                                    'TOMMO_AAF': 'float32',
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

    # 基于 CTRL_MAF 的三组划分（NaN 将被排除出三组统计）
    rare = df[df['CTRL_MAF'] < 0.01].copy()
    lowf = df[(df['CTRL_MAF'] >= 0.01) & (df['CTRL_MAF'] < 0.05)].copy()
    comm = df[df['CTRL_MAF'] >= 0.05].copy()

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

        # ===== 页面 3：三组散点（仅 PASS），每个分图先渲染为PNG =====
        png_paths_page3 = []
        titles_page3 = []
        for (title, g) in groups:
            g_pass = g[g['TOMMO_FILTER'] == 'PASS']
            g2 = _subset_for_scatter(g_pass)
            if g2.empty:
                png_paths_page3.append(None)
                titles_page3.append(title + ' (PASS ToMMo)')
                continue
            f_sc, ax_sc = plt.subplots(figsize=(6, 6))  # 方形画布
            ax_sc.scatter(
                g2['TOMMO_AAF'], g2['CTRL_AAF'], s=ms_pass, alpha=alpha_pass, linewidths=0,
                color=pass_color, label='PASS in ToMMo', zorder=2
            )
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
            out_png = os.path.join(pngdir, f"page3_{title.replace(' ', '_').replace('>', 'gt').replace('<', 'lt')}.png")
            f_sc.savefig(out_png, dpi=png_dpi)
            plt.close(f_sc)
            png_paths_page3.append(out_png)
            titles_page3.append(title + ' (PASS ToMMo)')
        _compose_three_pngs_to_pdf_page(png_paths_page3, titles_page3, page_title_suffix='')

        # ===== 页面 4：三组直方图（仅 PASS，CTRL_AAF - TOMMO_AAF），保留为矢量 =====
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
