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
import pysam

_global_lookup = {}

def _parse_variant_id(vid: str):
    """支持以下格式：
    1) chr1:968384:G:A
    2) 1:968384:G:A
    3) chr1:968384_A_T
    4) 1:968384_A_T
    返回 (chrom, pos(int), ref, alt) 或 None
    """
    try:
        if '_' in vid and ':' in vid:
            # e.g. chr1:968384_A_T
            chrom_pos, ref, alt = vid.split('_', 2)
            chrom, pos_str = chrom_pos.split(':', 1)
            return chrom, int(pos_str), ref, alt
        # colon-delimited form
        parts = vid.split(':')
        if len(parts) >= 4:
            chrom, pos_str, ref, alt = parts[0], parts[1], parts[2], parts[3]
            return chrom, int(pos_str), ref, alt
    except Exception:
        return None
    return None

def _query_tommo(tabix: pysam.TabixFile, chrom: str, pos: int, ref: str, alt: str):
    """查询ToMMo，返回 (IN_TOMMO(bool), TOMMO_AAF(float|nan), TOMMO_FILTER(str|nan))。
    处理chr/非chr，以及多等位/多AF的对应关系。
    """
    if tabix is None:
        return False, float("nan"), float("nan")
    chrom_candidates = [chrom]
    if chrom.startswith('chr'):
        chrom_candidates.append(chrom[3:])
    else:
        chrom_candidates.append('chr' + chrom)
    for c in chrom_candidates:
        try:
            for record in tabix.fetch(c, max(0, pos - 1), pos):
                fields = record.split('\t')
                vcf_pos = int(fields[1])
                if vcf_pos != pos:
                    continue
                vcf_ref = fields[3]
                if vcf_ref != ref:
                    continue
                vcf_alts = fields[4].split(',')
                if alt not in vcf_alts:
                    continue
                alt_idx = vcf_alts.index(alt)
                filter_field = fields[6]
                info_field = fields[7]
                af_val = float('nan')
                for part in info_field.split(';'):
                    if part.startswith('AF='):
                        af_str = part.split('=', 1)[1]
                        # AF 可能是逗号分隔，对应各ALT
                        try:
                            af_items = [float(x) if x not in ('', '.') else float('nan') for x in af_str.split(',')]
                            if alt_idx < len(af_items):
                                af_val = af_items[alt_idx]
                        except Exception:
                            af_val = float('nan')
                        break
                return True, af_val, filter_field
        except Exception:
            continue
    return False, float('nan'), float('nan')

def init_globals_for_chunk(vmiss_dict, case_aaf_dict, ctrl_aaf_dict, case_hwe_dict, ctrl_hwe_dict, tommo_vcf_path=None):
    global _global_lookup
    _global_lookup['vmiss_dict'] = vmiss_dict
    _global_lookup['case_aaf_dict'] = case_aaf_dict
    _global_lookup['ctrl_aaf_dict'] = ctrl_aaf_dict
    _global_lookup['case_hwe_dict'] = case_hwe_dict
    _global_lookup['ctrl_hwe_dict'] = ctrl_hwe_dict
    if tommo_vcf_path is not None:
        _global_lookup['tommo_tabix'] = pysam.TabixFile(tommo_vcf_path)
    else:
        _global_lookup['tommo_tabix'] = None

def process_chunk(chunk, idx, tmpdir):
    """
    对变异数据块进行处理，将每个变异的QC指标提取并写入临时文件。

    参数:
        chunk (DataFrame): 当前数据块，包含多个变异的统计信息
        idx (int): 数据块编号，用于命名输出文件
        tmpdir (str): 临时目录路径，用于存放中间结果

    返回:
        str: 写入结果的临时文件路径
    """
    vmiss_dict = _global_lookup.get('vmiss_dict', {})
    case_aaf_dict = _global_lookup.get('case_aaf_dict', {})
    ctrl_aaf_dict = _global_lookup.get('ctrl_aaf_dict', {})
    case_hwe_dict = _global_lookup.get('case_hwe_dict', {})
    ctrl_hwe_dict = _global_lookup.get('ctrl_hwe_dict', {})
    tommo_tabix = _global_lookup.get('tommo_tabix', None)

    tmp_output = os.path.join(tmpdir, f"chunk_{idx}_{uuid.uuid4().hex}.tsv")
    with open(tmp_output, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        for _, row in chunk.iterrows():
            vid = row["ID"]  # 变异ID
            aaf = row["ALT_FREQS"]  # 等位基因频率
            maf = min(aaf, 1 - aaf) if pd.notnull(aaf) else float("nan")  # 计算MAF
            vmiss = vmiss_dict.get(vid, float("nan"))  # 缺失率VMISS
            case_aaf = case_aaf_dict.get(vid, float("nan"))  # 病例组等位基因频率
            ctrl_aaf = ctrl_aaf_dict.get(vid, float("nan"))  # 对照组等位基因频率
            ctrl_maf = min(ctrl_aaf, 1 - ctrl_aaf) if pd.notnull(ctrl_aaf) else float("nan")  # 对照组MAF
            case_hwe = case_hwe_dict.get(vid, float("nan"))  # 病例组HWE p值
            ctrl_hwe = ctrl_hwe_dict.get(vid, float("nan"))  # 对照组HWE p值

            IN_TOMMO = False
            TOMMO_AAF = float("nan")
            TOMMO_FILTER = float("nan")

            if tommo_tabix is not None:
                parsed = _parse_variant_id(vid)
                if parsed is not None:
                    chrom, pos, ref, alt = parsed
                    IN_TOMMO, TOMMO_AAF, TOMMO_FILTER = _query_tommo(tommo_tabix, chrom, pos, ref, alt)

            # 写入列顺序: 变异ID, MAF, VMISS, 病例组AAF, 对照组AAF, 对照组MAF, 病例组HWE p值, 对照组HWE p值, IN_TOMMO, TOMMO_AAF, TOMMO_FILTER
            writer.writerow([vid, maf, vmiss, case_aaf, ctrl_aaf, ctrl_maf, case_hwe, ctrl_hwe, IN_TOMMO, TOMMO_AAF, TOMMO_FILTER])
    return tmp_output

def run_plink2_variant_qc_with_tommo(
    bed_prefix: str,
    tmpdir: str = "/tmp/variant_qc",
    plink2_path: str = "/home/b/b37974/plink2",
    threads: int = 8,
    output_prefix: str = "cteph_agp3k",
    verbose: bool = True,
    tommo_vcf_path: str = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz"
) -> str:
    """
    使用 plink2 对 Plink 格式基因型文件进行变异层面的QC计算。

    输出字段包括：
        - MAF（全部样本）
        - VMISS（全部样本）
        - CASE/CONTROL AAF
        - CONTROL MAF
        - CASE/CONTROL HWE P值
        - ToMMo VCF注释: 是否存在, 等位基因频率, FILTER状态

    参数:
        bed_prefix (str): 输入文件的 Plink 数据前缀（.bed/.bim/.fam）
        tmpdir (str): 临时目录路径，用于存放中间结果
        plink2_path (str): plink2 执行路径
        threads (int): 并行使用的线程数
        output_prefix (str): 输出 QC 结果文件的前缀
        verbose (bool): 是否打印进度信息
        tommo_vcf_path (str): ToMMo VCF文件路径，用于注释变异频率和过滤状态

    返回:
        str: 输出 QC 汇总结果文件的路径（.variant_qc_with_tommo.tsv）
    """

    # 确保 tmpdir 目录存在
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp()
    else:
        os.makedirs(tmpdir, exist_ok=True)
    if verbose:
        print(f"[INFO] 使用临时目录: {tmpdir}")

    def run_cmd(cmd, desc: Optional[str] = None):
        if verbose and desc:
            print(f"[INFO] 运行: {desc}")
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] 命令运行失败: {' '.join(cmd)}")
            raise e

    # Step 1: freq + vmiss
    out_all = os.path.join(tmpdir, "all_samples")
    # 计算全体样本的等位基因频率和缺失率
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--threads", str(threads),
        "--freq", "--missing", "--out", out_all
    ], "全体样本 AAF + VMISS")

    # Step 2: 分组 IID
    fam_df = pd.read_csv(f"{bed_prefix}.fam", sep=r"\s+", header=None)
    fam_df.columns = ["FID", "IID", "PID", "MID", "SEX", "PHENO"]
    # 根据PHENO字段区分病例和对照样本ID
    case_iids = fam_df[fam_df["PHENO"] == 2][["FID", "IID"]]
    ctrl_iids = fam_df[fam_df["PHENO"] == 1][["FID", "IID"]]

    case_iid_path = os.path.join(tmpdir, "case_iids.txt")
    ctrl_iid_path = os.path.join(tmpdir, "ctrl_iids.txt")
    # 保存病例和对照样本ID，用于plink2的--keep参数
    case_iids.to_csv(case_iid_path, sep="\t", index=False, header=False)
    ctrl_iids.to_csv(ctrl_iid_path, sep="\t", index=False, header=False)

    # Step 3: 分组 freq + HWE
    out_case = os.path.join(tmpdir, "case")
    out_ctrl = os.path.join(tmpdir, "ctrl")
    # 计算病例组的等位基因频率和HWE检验
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--keep", case_iid_path,
        "--threads", str(threads), "--freq", "--hardy", "--out", out_case
    ], "Case AAF + HWE")

    # 计算对照组的等位基因频率和HWE检验
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--keep", ctrl_iid_path,
        "--threads", str(threads), "--freq", "--hardy", "--out", out_ctrl
    ], "Control AAF + HWE")

    # Step 4: load lookup tables
    if verbose:
        print("[INFO] 读取辅助统计表...")

    # 读取全体样本缺失率字典
    vmiss_dict = dict(pd.read_csv(out_all + ".vmiss", sep=r"\s+")[["ID", "F_MISS"]].values)
    # 读取病例组等位基因频率字典
    case_aaf_dict = dict(pd.read_csv(out_case + ".afreq", sep=r"\s+")[["ID", "ALT_FREQS"]].values)
    # 读取对照组等位基因频率字典
    ctrl_aaf_dict = dict(pd.read_csv(out_ctrl + ".afreq", sep=r"\s+")[["ID", "ALT_FREQS"]].values)

    # 读取病例组和对照组HWE检验p值字典
    hwe_case_df = pd.read_csv(out_case + ".hardy", sep=r"\s+")
    hwe_ctrl_df = pd.read_csv(out_ctrl + ".hardy", sep=r"\s+")
    case_hwe_dict = dict(hwe_case_df[["ID", "P"]].values)
    ctrl_hwe_dict = dict(hwe_ctrl_df[["ID", "P"]].values)

    # Step 5: parallel chunk processing
    output_file = output_prefix + ".variant_qc_with_tommo.tsv"
    if verbose:
        test_chunk = pd.read_csv(out_all + ".afreq", sep=r"\s+", nrows=5)
        print("[DEBUG] .afreq 字段名:", list(test_chunk.columns))


    chunk_files = []
    reader = pd.read_csv(out_all + ".afreq", sep=r"\s+", chunksize=10000)
    with concurrent.futures.ProcessPoolExecutor(max_workers=10, initializer=init_globals_for_chunk,
                                                initargs=(vmiss_dict, case_aaf_dict, ctrl_aaf_dict, case_hwe_dict, ctrl_hwe_dict, tommo_vcf_path)) as executor:
        futures = []
        for i, chunk in enumerate(reader):
            if verbose:
                print(f"[INFO] 正在提交 chunk {i + 1} 任务...")
            futures.append(
                executor.submit(
                    process_chunk,
                    chunk, i, tmpdir
                )
            )
        for i, future in enumerate(futures):
            chunk_file = future.result()
            if verbose:
                print(f"[INFO] chunk {i + 1} 处理完成，结果文件: {chunk_file}")
            chunk_files.append(chunk_file)

    # merge chunk files
    with open(output_file, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        # 写入表头
        writer.writerow([
            "VARIANT_ID", "MAF", "VMISS", "CASE_AAF",
            "CTRL_AAF", "CTRL_MAF", "CASE_HWE", "CTRL_HWE",
            "IN_TOMMO", "TOMMO_AAF", "TOMMO_FILTER"
        ])
        # 合并所有chunk的结果文件
        for chunk_file in chunk_files:
            with open(chunk_file, "r") as fin:
                for line in fin:
                    fout.write(line)

    if verbose:
        print(f"[INFO] 输出完成: {output_file}")
        print(f"[INFO] 结果文件路径（可用于后续加载）: {output_file}")
    return output_file

# new funtion to plot panel comapre results
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

    # 以分块方式读取，降低内存占用（仅加载必要列）
    needed_cols = [
        'VARIANT_ID', 'CTRL_MAF', 'CTRL_AAF', 'TOMMO_AAF', 'IN_TOMMO', 'TOMMO_FILTER'
    ]
    chunk_size = 500_000  # 可根据内存情况调整
    chunks = []
    for chunk in pd.read_csv(variant_qc_with_tommo, sep='\t', usecols=needed_cols,
                             dtype={'VARIANT_ID': 'string',
                                    'CTRL_MAF': 'float32',
                                    'CTRL_AAF': 'float32',
                                    'TOMMO_AAF': 'float32',
                                    'IN_TOMMO': 'object',  # 先以 object 读取，后续统一转为 boolean
                                    'TOMMO_FILTER': 'string'},
                             chunksize=chunk_size):
        # 逐块轻量清洗，减少最终拼接的开销
        # 将 IN_TOMMO 规范为 pandas NA/True/False 的字符串或布尔值
        if 'IN_TOMMO' in chunk.columns:
            # 兼容 True/False/1/0/"True"/"False"
            chunk['IN_TOMMO'] = chunk['IN_TOMMO'].map(
                lambda x: True if x in (True, 1, '1', 'True', 'TRUE') else (False if x in (False, 0, '0', 'False', 'FALSE') else pd.NA)
            )
        chunks.append(chunk)
    df = pd.concat(chunks, ignore_index=True)


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
