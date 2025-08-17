"""
variant_qc_flags.py
===================

脚本用途概述
-----------
本脚本用于大规模基因型变体质控（QC）过滤、统计与可视化，适用于WGS/WES等大样本量的群体遗传学分析场景。
主要面向GWAS、罕见变异关联分析等场景下的变体筛选与QC流程。脚本支持对变体缺失率（VMISS）、Hardy-Weinberg平衡（HWE）等指标进行分布绘图、分类型统计，并输出通过QC的变体ID列表，支持与PLINK二进制文件联动提取最终变体集合。

输入文件要求
------------
主输入为一个变体质控汇总文件（variant_qc_summary），通常为TSV格式，须包含如下字段：
  - VARIANT_ID（或指定的变体ID列名）：如chr:pos:ref:alt格式
  - CTRL_MAF：对照组等位基因频率（float, 0~1）
  - VMISS：变体缺失率（float, 0~1）
  - CTRL_HWE：对照组HWE检验p值（float, 0~1）
  - CASE_HWE：病例组HWE检验p值（float, 0~1）
文件需包含上述字段中的至少部分，具体依赖于所调用的函数。

主要函数功能简介
----------------
1. plot_vmiss_distribution_by_maf_category
   - 功能：按CTRL_MAF分为稀有、低频、常见三类，绘制VMISS分布直方图，统计每类通过VMISS阈值的变体数，并输出通过变体ID列表。
   - 主要参数：
     - variant_qc_summary (str): 变体QC汇总文件路径
     - vmiss_threshold (float): VMISS阈值，低于此值视为通过
     - variant_id_col (str): 变体ID列名
     - output_tsv (str): 输出通过VMISS变体ID的TSV文件路径
     - output_prefix (str): 输出图表文件名前缀

2. plot_hwe_scatter_by_maf_category
   - 功能：按CTRL_MAF分为三类，绘制CTRL_HWE vs CASE_HWE散点图（含边缘KDE分布），支持为不同类别设定HWE阈值线，仅用于可视化。输出通过HWE阈值的变体ID列表。
   - 主要参数：
     - variant_qc_summary (str): 变体QC汇总文件路径
     - hwe_thresholds (dict): 每个类别的HWE阈值字典（如{"Rare Variant (<0.01)": {"CTRL_HWE":1e-6,"CASE_HWE":1e-6}, ...}）
     - variant_id_col (str): 变体ID列名
     - output_tsv (str): 输出通过HWE变体ID的TSV文件路径
     - output_prefix (str): 输出图表文件名前缀

3. extract_pass_variants_by_intersection
   - 功能：对两个通过QC的变体ID列表（如VMISS和HWE），取交集并按染色体/位置排序输出。可选调用PLINK2从给定bed_prefix的基因型数据中提取这些变体。
   - 主要参数：
     - pass_vmiss_path (str): VMISS通过变体ID文件路径
     - pass_hwe_path (str): HWE通过变体ID文件路径
     - output_path (str): 输出交集变体ID文件路径
     - bed_prefix (str): 输入PLINK二进制文件前缀（可选）
     - output_prefix (str): 输出PLINK提取结果前缀
     - threads (int): PLINK2线程数
     - plink2 (str): plink2可执行文件路径

输出内容说明
-----------
1. PDF图表：每个主函数均会生成对应的PDF格式QC分布图，文件名由output_prefix指定。
2. TSV文件：每个QC步骤输出通过筛选的变体ID列表（无表头，每行一个ID）。
3. PLINK提取结果（可选）：若extract_pass_variants_by_intersection指定bed_prefix，则自动调用plink2生成新的.bed/.bim/.fam文件，仅包含最终通过变体。

使用建议
--------
- 本脚本适合大规模变体数据（百万~千万SNP/INDEL），采用pandas分块读取，节省内存。
- 推荐在内存≥16GB的服务器上运行，尤其是全基因组大队列数据。
- 若输入文件极大，建议先用bgzip/tabix等工具按需预筛选或拆分。
- 调用PLINK2提取时，请确保plink2在环境PATH或指定其全路径，且磁盘空间充足。

作者与修改记录
--------------
- 作者：ZHAO TIE
- 修改记录：2025-08-15
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

plt.style.use('default')

def plot_vmiss_distribution_by_maf_category(
    variant_qc_summary: str,
    vmiss_threshold: float = 0.05,
    variant_id_col: str = "VARIANT_ID",
    output_tsv: str = "vmiss_pass_variants.tsv",
    output_prefix: str = "cteph_agp3k"
) -> str:
    """
    从大型variant_qc_summary中读取并绘制VMISS分布图，分类展示，统计信息纳入legend框，保存通过阈值的变体ID。

    参数:
    variant_qc_summary (str): 输入的变体质控汇总文件路径，TSV格式。
    vmiss_threshold (float): VMISS过滤阈值，低于该值的变体被视为通过。
    variant_id_col (str): 变体ID所在列名，默认"VARIANT_ID"。
    output_tsv (str): 输出通过VMISS阈值的变体ID文件路径。
    output_prefix (str): 输出图表文件名前缀。

    返回:
    str: 输出的通过VMISS阈值变体ID文件路径。
    """

    # 初始化列表用于存储分块读取的数据，统计NA行数
    data_chunks = []
    total_na_rows = 0

    # 指定需要读取的列，分块读取大文件以节省内存
    cols_needed = ["CTRL_MAF", "VMISS", variant_id_col]
    reader = pd.read_csv(variant_qc_summary, sep="\t", usecols=cols_needed, chunksize=500000)

    for i, chunk in enumerate(reader):
        # 过滤掉CTRL_MAF或VMISS为NA的行，统计丢弃数量
        before = len(chunk)
        chunk = chunk[chunk["CTRL_MAF"].notna() & chunk["VMISS"].notna()]
        na_dropped = before - len(chunk)
        total_na_rows += na_dropped
        if na_dropped > 0:
            print(f"[Warning] Chunk {i}: Dropped {na_dropped} rows with NA in CTRL_MAF or VMISS")

        # 基于CTRL_MAF划分类别
        chunk["Category"] = chunk["CTRL_MAF"].apply(lambda maf: (
            "Rare Variant (<0.01)" if 0 <= maf < 0.01 else
            "Low Frequency Variant (0.01~0.05)" if 0.01 <= maf <= 0.05 else
            "Common Variant (>0.05)" if maf > 0.05 else None
        ))
        chunk = chunk[chunk["Category"].notna()]
        data_chunks.append(chunk)

    # 合并所有分块数据
    df = pd.concat(data_chunks, ignore_index=True)

    # 构造VMISS直方图的bin边界
    x_min, x_max = 0.0, 1.0
    bin_edges = np.arange(x_min, x_max + 0.02, 0.02)

    # 设置绘图样式和子图布局
    sns.set(style="whitegrid")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=False)

    categories = [
        "Rare Variant (<0.01)",
        "Low Frequency Variant (0.01~0.05)",
        "Common Variant (>0.05)"
    ]
    colors = ["#1f77b4", "#2ca02c", "#ff7f0e"]

    # 计算全局统计信息
    total_variant_count = len(df)
    total_below_threshold = (df["VMISS"] < vmiss_threshold).sum()

    # 保存通过VMISS阈值的变体ID到TSV文件
    passed_variants = df[df["VMISS"] < vmiss_threshold][variant_id_col].dropna().unique()
    pd.Series(passed_variants).to_csv(output_tsv, sep="\t", index=False, header=False)

    # 逐类别绘制VMISS分布直方图
    for ax, category, color in zip(axes, categories, colors):
        sub_df = df[df["Category"] == category]
        count_total = len(sub_df)
        count_pass = (sub_df["VMISS"] < vmiss_threshold).sum()
        percent_pass = (100 * count_pass / count_total) if count_total > 0 else 0

        # 绘制带核密度估计的直方图
        sns.histplot(sub_df["VMISS"], bins=bin_edges, kde=True, color=color,
                     edgecolor="black", ax=ax)
        # 添加阈值线
        ax.axvline(vmiss_threshold, linestyle="--", color="red", label=f"Threshold = {vmiss_threshold}")
        ax.set_xlim(x_min, x_max)
        ax.set_title(category)
        ax.set_xlabel("VMISS")

        # 添加统计信息到legend
        summary_label = (
            f"Total: {count_total:,}\n"
            f"Pass: {count_pass:,} ({percent_pass:.1f}%)"
        )
        dummy_patch = Patch(color='none', label=summary_label)

        # 合并图例元素，显示阈值线和统计信息
        handles, labels = ax.get_legend_handles_labels()
        handles.append(dummy_patch)
        ax.legend(handles=handles, loc='upper right', frameon=True, fontsize=9)

    # 设置主标题，包含全局统计信息
    global_percent = 100 * total_below_threshold / total_variant_count
    fig.suptitle(
        f"VMISS Distribution by CTRL_MAF Category\n"
        f"Total Variants: {total_variant_count:,} Pass Threshold (< {vmiss_threshold}): "
        f"{total_below_threshold:,} ({global_percent:.1f}%)",
        fontsize=16
    )

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plt.savefig(f"{output_prefix}.vmiss.pdf")
    plt.show()

    # 打印NA值丢弃总计
    if total_na_rows > 0:
        print(f"\n[Summary] Total rows dropped due to NA: {total_na_rows:,}")

    return output_tsv
    

def plot_hwe_scatter_by_maf_category(
    variant_qc_summary: str,
    hwe_thresholds: dict = None,
    variant_id_col: str = "VARIANT_ID",
    output_tsv: str = "hwe_pass_variants.tsv",
    output_prefix: str = "cteph_agp3k"
) -> str:
    """
    根据CTRL_MAF分类，并绘制CTRL_HWE vs CASE_HWE散点图（3个分图，2x2布局，右下角空白）。
    hwe_thresholds是一个可选字典，指定每个类别的HWE阈值线（不进行过滤，仅用于可视化）。

    参数:
    variant_qc_summary (str): 输入的变体质控汇总文件路径，TSV格式。
    hwe_thresholds (dict): 每个类别对应的HWE阈值字典，格式示例：
        {
            "Rare Variant (<0.01)": {"CTRL_HWE": 1e-6, "CASE_HWE": 1e-6},
            "Low Frequency Variant (0.01~0.05)": {"CTRL_HWE": 1e-5, "CASE_HWE": 1e-5},
            "Common Variant (>0.05)": {"CTRL_HWE": 1e-4, "CASE_HWE": 1e-4}
        }
    variant_id_col (str): 变体ID所在列名，默认"VARIANT_ID"。
    output_tsv (str): 输出通过HWE阈值的变体ID文件路径。
    output_prefix (str): 输出图表文件名前缀。

    返回:
    str: 输出的通过HWE阈值变体ID文件路径。
    """
    import matplotlib.gridspec as gridspec
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    import gc

    # 只读取必要列，分块读取以节省内存
    usecols = [variant_id_col, "CTRL_MAF", "CTRL_HWE", "CASE_HWE"]
    def classify_maf(maf):
        if 0 <= maf < 0.01:
            return "Rare Variant (<0.01)"
        elif 0.01 <= maf <= 0.05:
            return "Low Frequency Variant (0.01~0.05)"
        elif maf > 0.05:
            return "Common Variant (>0.05)"
        else:
            return None

    data_chunks = []
    reader = pd.read_csv(variant_qc_summary, sep="\t", usecols=usecols, chunksize=500000)

    for i, chunk in enumerate(reader):
        # 不丢弃NA，保留所有行，分类MAF
        chunk["Category"] = chunk["CTRL_MAF"].apply(classify_maf)
        chunk = chunk[chunk["Category"].notna()]
        data_chunks.append(chunk)
        del chunk
        gc.collect()

    df = pd.concat(data_chunks, ignore_index=True)

    pass_variant_ids = []

    categories = [
        "Rare Variant (<0.01)",
        "Low Frequency Variant (0.01~0.05)",
        "Common Variant (>0.05)"
    ]
    colors = ["#1f77b4", "#2ca02c", "#ff7f0e"]

    def plot_one_category(ax_main, ax_top, ax_right, sub_df, label, color, threshold_control, threshold_case):
        # 使用原始数据进行阈值判断；为绘图准备清洗后的数据（log 轴下非正数不可见）
        x = sub_df["CTRL_HWE"]
        y = sub_df["CASE_HWE"]

        # 计算通过阈值的变体数量及四象限计数（基于原始值，不做替换）
        if threshold_control is not None and threshold_case is not None:
            q1 = ((x >= threshold_control) & (y >= threshold_case)).sum()
            pass_ids = sub_df[(x >= threshold_control) & (y >= threshold_case)].index
            pass_variant_ids.extend(pass_ids)
            q2 = ((x < threshold_control) & (y >= threshold_case)).sum()
            q3 = ((x < threshold_control) & (y < threshold_case)).sum()
            q4 = ((x >= threshold_control) & (y < threshold_case)).sum()
        elif threshold_control is not None and threshold_case is None:
            q1 = (x >= threshold_control).sum()
            pass_ids = sub_df[(x >= threshold_control)].index
            pass_variant_ids.extend(pass_ids)
            q2 = (x < threshold_control).sum()
            q3 = 0
            q4 = 0
        else:
            q1 = len(sub_df)
            pass_ids = sub_df.index
            pass_variant_ids.extend(pass_ids)
            q2 = q3 = q4 = 0

        # ---- 散点图（仅用于显示，过滤掉非正数以适配对数坐标） ----
        x_scatter = pd.to_numeric(x, errors='coerce')
        y_scatter = pd.to_numeric(y, errors='coerce')
        valid_scatter = np.isfinite(x_scatter) & np.isfinite(y_scatter) & (x_scatter > 0) & (y_scatter > 0)
        ax_main.scatter(x_scatter[valid_scatter], y_scatter[valid_scatter], alpha=0.3, c=color, s=20)
        if threshold_control is not None:
            ax_main.axvline(x=threshold_control, color='red', linestyle='--')
        if threshold_case is not None:
            ax_main.axhline(y=threshold_case, color='blue', linestyle='--')
        ax_main.set_xscale('log')
        ax_main.set_yscale('log')
        ax_main.set_xlabel('CTRL_HWE')
        ax_main.set_ylabel('CASE_HWE')
        ax_main.grid(True)

        # ---- 顶部 KDE（CTRL_HWE） ----
        x_plot = pd.to_numeric(x, errors='coerce')
        x_plot = x_plot[np.isfinite(x_plot) & (x_plot > 0)]
        if x_plot.nunique() >= 2:
            sns.kdeplot(x_plot, ax=ax_top, fill=True, color='gray', linewidth=1.5, cut=0)
            ax_top.set_xscale('log')
            if threshold_control is not None:
                ax_top.axvline(x=threshold_control, color='red', linestyle='--')
        else:
            ax_top.text(0.5, 0.5, 'KDE skipped\n(non-positive or low variance)',
                        ha='center', va='center', transform=ax_top.transAxes, fontsize=8)
        ax_top.set_xlabel('')
        ax_top.set_ylabel('')
        ax_top.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                           labelbottom=False, labelleft=False)
        ax_top.grid(False)

        # ---- 右侧 KDE（CASE_HWE） ----
        y_plot = pd.to_numeric(y, errors='coerce')
        y_plot = y_plot[np.isfinite(y_plot) & (y_plot > 0)]
        if y_plot.nunique() >= 2:
            sns.kdeplot(y_plot, ax=ax_right, fill=True, color='gray', linewidth=1.5, vertical=True, cut=0)
            ax_right.set_yscale('log')
            if threshold_case is not None:
                ax_right.axhline(y=threshold_case, color='blue', linestyle='--')
        else:
            ax_right.text(0.5, 0.5, 'KDE skipped\n(non-positive or low variance)',
                          ha='center', va='center', transform=ax_right.transAxes, fontsize=8)
        ax_right.set_xlabel('')
        ax_right.set_ylabel('')
        ax_right.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                             labelbottom=False, labelleft=False)
        ax_right.grid(False)
        # 插图显示四象限计数或总数
        top_ylim = ax_top.get_ylim()
        right_xlim = ax_right.get_xlim()
        inset_ax = inset_axes(ax_main, width="25%", height="25%", loc='upper right',
                              bbox_to_anchor=(0.28, 0.28, 1, 1), bbox_transform=ax_main.transAxes, borderpad=0)
        inset_ax.set_xlim(right_xlim)
        inset_ax.set_ylim(top_ylim)
        for spine in inset_ax.spines.values():
            spine.set_visible(True)
        inset_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        inset_ax.grid(False)
        xmid = (right_xlim[0] + right_xlim[1]) / 2
        ymid = (top_ylim[0] + top_ylim[1]) / 2

        if threshold_control is not None and threshold_case is not None:
            # 绘制四象限分割线及计数
            inset_ax.axvline(x=xmid, linestyle='--', color='red', linewidth=1.5)
            inset_ax.axhline(y=ymid, linestyle='--', color='blue', linewidth=1.5)
            inset_ax.text((xmid + right_xlim[1]) / 2, (ymid + top_ylim[1]) / 2, f'{q1:,}', ha='center', va='center', fontsize=7, fontweight='bold')
            inset_ax.text((right_xlim[0] + xmid) / 2, (ymid + top_ylim[1]) / 2, f'{q2:,}', ha='center', va='center', fontsize=7, fontweight='bold')
            inset_ax.text((right_xlim[0] + xmid) / 2, (top_ylim[0] + ymid) / 2, f'{q3:,}', ha='center', va='center', fontsize=7, fontweight='bold')
            inset_ax.text((xmid + right_xlim[1]) / 2, (top_ylim[0] + ymid) / 2, f'{q4:,}', ha='center', va='center', fontsize=7, fontweight='bold')
        elif threshold_control is not None and threshold_case is None:
            # 仅绘制CTRL_HWE阈值线及左右计数
            inset_ax.axvline(x=xmid, linestyle='--', color='red', linewidth=1.5)
            left = q2
            right = q1
            inset_ax.text((right_xlim[0] + xmid) / 2, ymid, f'{left:,}', ha='center', va='center', fontsize=7, fontweight='bold')
            inset_ax.text((xmid + right_xlim[1]) / 2, ymid, f'{right:,}', ha='center', va='center', fontsize=7, fontweight='bold')
        else:
            # 无阈值，展示总数
            inset_ax.text(xmid, ymid, f'{q1:,}', ha='center', va='center', fontsize=7, fontweight='bold')

    # 创建2x2网格布局，右下角空白
    fig = plt.figure(figsize=(14, 14))
    outer_gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
    subplots = [outer_gs[0, 0], outer_gs[0, 1], outer_gs[1, 0]]

    summary_rows = []

    # 遍历三个类别，绘制对应子图
    for subplot_spec, category, color in zip(subplots, categories, colors):
        sub_df = df[df["Category"] == category]
        threshold = hwe_thresholds.get(category, {}) if hwe_thresholds else {}
        threshold_control = threshold.get("CTRL_HWE")
        threshold_case = threshold.get("CASE_HWE")

        # 创建内嵌子网格，主图和两个KDE图
        inner_gs = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=subplot_spec,
                                                    width_ratios=[6, 1.5], height_ratios=[1.5, 6],
                                                    wspace=0.05, hspace=0.05)
        ax_main = plt.Subplot(fig, inner_gs[1, 0])
        ax_top = plt.Subplot(fig, inner_gs[0, 0], sharex=ax_main)
        ax_right = plt.Subplot(fig, inner_gs[1, 1], sharey=ax_main)
        fig.add_subplot(ax_main)
        fig.add_subplot(ax_top)
        fig.add_subplot(ax_right)

        # 汇总信息，用于右下角表格
        summary_rows.append({
            "": category.replace(" (", "\n("),
            "raw_category": category,
            "CTRL_HWE Threshold": f"{threshold_control:.1e}" if threshold_control is not None else "None",
            "CASE_HWE Threshold": f"{threshold_case:.1e}" if threshold_case is not None else "None",
        })

        plot_one_category(ax_main, ax_top, ax_right, sub_df, category, color, threshold_control, threshold_case)

    # 右下角添加汇总表格
    ax_legend = fig.add_subplot(outer_gs[1, 1])
    ax_legend.axis('off')
    summary_df = pd.DataFrame(summary_rows)
    summary_df_display = summary_df.drop(columns=["raw_category"])
    table = ax_legend.table(cellText=summary_df_display.values,
                            colLabels=summary_df_display.columns,
                            cellLoc='center',
                            colWidths=[0.35, 0.325, 0.325],
                            loc='center')
    table.scale(1.2, 1.6)
    table.auto_set_font_size(False)
    table.set_fontsize(9)

    # 设置表头样式，灰底黑字加粗
    for col_idx in range(len(summary_df_display.columns)):
        cell = table[(0, col_idx)]
        cell.set_facecolor('#f0f0f0')
        cell.set_text_props(color='black', weight='bold')

    # 设置行名背景色（带透明度），黑色加粗字体
    for row_idx, raw_category in enumerate(summary_df["raw_category"]):
        color = colors[categories.index(raw_category)]
        facecolor_rgba = plt.matplotlib.colors.to_rgba(color, alpha=0.35)
        cell = table[(row_idx + 1, 0)]
        cell.set_facecolor(facecolor_rgba)
        cell.set_text_props(color='black', weight='bold')

    fig.suptitle("CTRL_HWE vs CASE_HWE by CTRL_MAF Category", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(f"{output_prefix}.hwe.png", dpi=300) # 保存为PNG格式（PDF太大）
    plt.show()

    # 导出通过阈值的变体ID，若无该列则警告
    if variant_id_col in df.columns:
        pass_variants = df.loc[pass_variant_ids, variant_id_col].dropna().unique()
        pd.Series(pass_variants).to_csv(output_tsv, sep="\t", index=False, header=False)
    else:
        print(f"[Warning] {variant_id_col} column not found; skipping export of passed variants.")
    
    return output_tsv

def plot_hwe_scatter_by_maf_category(
    variant_qc_summary: str,
    hwe_thresholds: dict = None,
    variant_id_col: str = "VARIANT_ID",
    output_tsv: str = "hwe_pass_variants.tsv",
    output_prefix: str = "cteph_agp3k"
) -> str:
    """
    根据CTRL_MAF分类，并绘制CTRL_HWE vs CASE_HWE散点图（3个分图，2x2布局，右下角空白）。
    hwe_thresholds是一个可选字典，指定每个类别的HWE阈值线（不进行过滤，仅用于可视化）。

    参数:
    variant_qc_summary (str): 输入的变体质控汇总文件路径，TSV格式。
    hwe_thresholds (dict): 每个类别对应的HWE阈值字典，格式示例：
        {
            "Rare Variant (<0.01)": {"CTRL_HWE": 1e-6, "CASE_HWE": 1e-6},
            "Low Frequency Variant (0.01~0.05)": {"CTRL_HWE": 1e-5, "CASE_HWE": 1e-5},
            "Common Variant (>0.05)": {"CTRL_HWE": 1e-4, "CASE_HWE": 1e-4}
        }
    variant_id_col (str): 变体ID所在列名，默认"VARIANT_ID"。
    output_tsv (str): 输出通过HWE阈值的变体ID文件路径。
    output_prefix (str): 输出图表文件名前缀。

    返回:
    str: 输出的通过HWE阈值变体ID文件路径。
    """
    import matplotlib.gridspec as gridspec
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    import gc

    # 只读取必要列，分块读取以节省内存
    usecols = [variant_id_col, "CTRL_MAF", "CTRL_HWE", "CASE_HWE"]
    def classify_maf(maf):
        if 0 <= maf < 0.01:
            return "Rare Variant (<0.01)"
        elif 0.01 <= maf <= 0.05:
            return "Low Frequency Variant (0.01~0.05)"
        elif maf > 0.05:
            return "Common Variant (>0.05)"
        else:
            return None

    data_chunks = []
    reader = pd.read_csv(variant_qc_summary, sep="\t", usecols=usecols, chunksize=500000)

    for i, chunk in enumerate(reader):
        # 不丢弃NA，保留所有行，分类MAF
        chunk["Category"] = chunk["CTRL_MAF"].apply(classify_maf)
        chunk = chunk[chunk["Category"].notna()]
        data_chunks.append(chunk)
        del chunk
        gc.collect()

    df = pd.concat(data_chunks, ignore_index=True)

    pass_variant_ids = []

    categories = [
        "Rare Variant (<0.01)",
        "Low Frequency Variant (0.01~0.05)",
        "Common Variant (>0.05)"
    ]
    colors = ["#1f77b4", "#2ca02c", "#ff7f0e"]

    def plot_one_category(ax_main, ax_top, ax_right, sub_df, label, color, threshold_control, threshold_case):
        # 使用原始数据进行阈值判断；为绘图准备清洗后的数据（log 轴下非正数不可见）
        x = sub_df["CTRL_HWE"]
        y = sub_df["CASE_HWE"]

        # 计算通过阈值的变体数量及四象限计数（基于原始值，不做替换）
        if threshold_control is not None and threshold_case is not None:
            q1 = ((x >= threshold_control) & (y >= threshold_case)).sum()
            pass_ids = sub_df[(x >= threshold_control) & (y >= threshold_case)].index
            pass_variant_ids.extend(pass_ids)
            q2 = ((x < threshold_control) & (y >= threshold_case)).sum()
            q3 = ((x < threshold_control) & (y < threshold_case)).sum()
            q4 = ((x >= threshold_control) & (y < threshold_case)).sum()
        elif threshold_control is not None and threshold_case is None:
            q1 = (x >= threshold_control).sum()
            pass_ids = sub_df[(x >= threshold_control)].index
            pass_variant_ids.extend(pass_ids)
            q2 = (x < threshold_control).sum()
            q3 = 0
            q4 = 0
        else:
            q1 = len(sub_df)
            pass_ids = sub_df.index
            pass_variant_ids.extend(pass_ids)
            q2 = q3 = q4 = 0

        # ---- 散点图（仅用于显示，过滤掉非正数以适配对数坐标） ----
        x_scatter = pd.to_numeric(x, errors='coerce')
        y_scatter = pd.to_numeric(y, errors='coerce')
        valid_scatter = np.isfinite(x_scatter) & np.isfinite(y_scatter) & (x_scatter > 0) & (y_scatter > 0)
        ax_main.scatter(x_scatter[valid_scatter], y_scatter[valid_scatter], alpha=0.3, c=color, s=20)
        if threshold_control is not None:
            ax_main.axvline(x=threshold_control, color='red', linestyle='--')
        if threshold_case is not None:
            ax_main.axhline(y=threshold_case, color='blue', linestyle='--')
        ax_main.set_xscale('log')
        ax_main.set_yscale('log')
        ax_main.set_xlabel('CTRL_HWE')
        ax_main.set_ylabel('CASE_HWE')
        ax_main.grid(True)

        # ---- 顶部 KDE（CTRL_HWE） ----
        x_plot = pd.to_numeric(x, errors='coerce')
        x_plot = x_plot[np.isfinite(x_plot) & (x_plot > 0)]
        if x_plot.nunique() >= 2:
            sns.kdeplot(x_plot, ax=ax_top, fill=True, color='gray', linewidth=1.5, cut=0)
            ax_top.set_xscale('log')
            if threshold_control is not None:
                ax_top.axvline(x=threshold_control, color='red', linestyle='--')
        else:
            ax_top.text(0.5, 0.5, 'KDE skipped\n(non-positive or low variance)',
                        ha='center', va='center', transform=ax_top.transAxes, fontsize=8)
        ax_top.set_xlabel('')
        ax_top.set_ylabel('')
        ax_top.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                           labelbottom=False, labelleft=False)
        ax_top.grid(False)

        # ---- 右侧 KDE（CASE_HWE） ----
        y_plot = pd.to_numeric(y, errors='coerce')
        y_plot = y_plot[np.isfinite(y_plot) & (y_plot > 0)]
        if y_plot.nunique() >= 2:
            sns.kdeplot(y_plot, ax=ax_right, fill=True, color='gray', linewidth=1.5, vertical=True, cut=0)
            ax_right.set_yscale('log')
            if threshold_case is not None:
                ax_right.axhline(y=threshold_case, color='blue', linestyle='--')
        else:
            ax_right.text(0.5, 0.5, 'KDE skipped\n(non-positive or low variance)',
                          ha='center', va='center', transform=ax_right.transAxes, fontsize=8)
        ax_right.set_xlabel('')
        ax_right.set_ylabel('')
        ax_right.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                             labelbottom=False, labelleft=False)
        ax_right.grid(False)
        # 插图显示四象限计数或总数
        top_ylim = ax_top.get_ylim()
        right_xlim = ax_right.get_xlim()
        inset_ax = inset_axes(ax_main, width="25%", height="25%", loc='upper right',
                              bbox_to_anchor=(0.28, 0.28, 1, 1), bbox_transform=ax_main.transAxes, borderpad=0)
        inset_ax.set_xlim(right_xlim)
        inset_ax.set_ylim(top_ylim)
        for spine in inset_ax.spines.values():
            spine.set_visible(True)
        inset_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        inset_ax.grid(False)
        xmid = (right_xlim[0] + right_xlim[1]) / 2
        ymid = (top_ylim[0] + top_ylim[1]) / 2

        if threshold_control is not None and threshold_case is not None:
            # 绘制四象限分割线及计数
            inset_ax.axvline(x=xmid, linestyle='--', color='red', linewidth=1.5)
            inset_ax.axhline(y=ymid, linestyle='--', color='blue', linewidth=1.5)
            inset_ax.text((xmid + right_xlim[1]) / 2, (ymid + top_ylim[1]) / 2, f'{q1:,}', ha='center', va='center', fontsize=7, fontweight='bold')
            inset_ax.text((right_xlim[0] + xmid) / 2, (ymid + top_ylim[1]) / 2, f'{q2:,}', ha='center', va='center', fontsize=7, fontweight='bold')
            inset_ax.text((right_xlim[0] + xmid) / 2, (top_ylim[0] + ymid) / 2, f'{q3:,}', ha='center', va='center', fontsize=7, fontweight='bold')
            inset_ax.text((xmid + right_xlim[1]) / 2, (top_ylim[0] + ymid) / 2, f'{q4:,}', ha='center', va='center', fontsize=7, fontweight='bold')
        elif threshold_control is not None and threshold_case is None:
            # 仅绘制CTRL_HWE阈值线及左右计数
            inset_ax.axvline(x=xmid, linestyle='--', color='red', linewidth=1.5)
            left = q2
            right = q1
            inset_ax.text((right_xlim[0] + xmid) / 2, ymid, f'{left:,}', ha='center', va='center', fontsize=7, fontweight='bold')
            inset_ax.text((xmid + right_xlim[1]) / 2, ymid, f'{right:,}', ha='center', va='center', fontsize=7, fontweight='bold')
        else:
            # 无阈值，展示总数
            inset_ax.text(xmid, ymid, f'{q1:,}', ha='center', va='center', fontsize=7, fontweight='bold')

    # 创建2x2网格布局，右下角空白
    fig = plt.figure(figsize=(14, 14))
    outer_gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
    subplots = [outer_gs[0, 0], outer_gs[0, 1], outer_gs[1, 0]]

    summary_rows = []

    # 遍历三个类别，绘制对应子图
    for subplot_spec, category, color in zip(subplots, categories, colors):
        sub_df = df[df["Category"] == category]
        threshold = hwe_thresholds.get(category, {}) if hwe_thresholds else {}
        threshold_control = threshold.get("CTRL_HWE")
        threshold_case = threshold.get("CASE_HWE")

        # 创建内嵌子网格，主图和两个KDE图
        inner_gs = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=subplot_spec,
                                                    width_ratios=[6, 1.5], height_ratios=[1.5, 6],
                                                    wspace=0.05, hspace=0.05)
        ax_main = plt.Subplot(fig, inner_gs[1, 0])
        ax_top = plt.Subplot(fig, inner_gs[0, 0], sharex=ax_main)
        ax_right = plt.Subplot(fig, inner_gs[1, 1], sharey=ax_main)
        fig.add_subplot(ax_main)
        fig.add_subplot(ax_top)
        fig.add_subplot(ax_right)

        # 汇总信息，用于右下角表格
        summary_rows.append({
            "": category.replace(" (", "\n("),
            "raw_category": category,
            "CTRL_HWE Threshold": f"{threshold_control:.1e}" if threshold_control is not None else "None",
            "CASE_HWE Threshold": f"{threshold_case:.1e}" if threshold_case is not None else "None",
        })

        plot_one_category(ax_main, ax_top, ax_right, sub_df, category, color, threshold_control, threshold_case)

    # 右下角添加汇总表格
    ax_legend = fig.add_subplot(outer_gs[1, 1])
    ax_legend.axis('off')
    summary_df = pd.DataFrame(summary_rows)
    summary_df_display = summary_df.drop(columns=["raw_category"])
    table = ax_legend.table(cellText=summary_df_display.values,
                            colLabels=summary_df_display.columns,
                            cellLoc='center',
                            colWidths=[0.35, 0.325, 0.325],
                            loc='center')
    table.scale(1.2, 1.6)
    table.auto_set_font_size(False)
    table.set_fontsize(9)

    # 设置表头样式，灰底黑字加粗
    for col_idx in range(len(summary_df_display.columns)):
        cell = table[(0, col_idx)]
        cell.set_facecolor('#f0f0f0')
        cell.set_text_props(color='black', weight='bold')

    # 设置行名背景色（带透明度），黑色加粗字体
    for row_idx, raw_category in enumerate(summary_df["raw_category"]):
        color = colors[categories.index(raw_category)]
        facecolor_rgba = plt.matplotlib.colors.to_rgba(color, alpha=0.35)
        cell = table[(row_idx + 1, 0)]
        cell.set_facecolor(facecolor_rgba)
        cell.set_text_props(color='black', weight='bold')

    fig.suptitle("CTRL_HWE vs CASE_HWE by CTRL_MAF Category", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(f"{output_prefix}.hwe.png", dpi=300) # 保存为PNG格式（PDF太大）
    plt.show()

    # 导出通过阈值的变体ID，若无该列则警告
    if variant_id_col in df.columns:
        pass_variants = df.loc[pass_variant_ids, variant_id_col].dropna().unique()
        pd.Series(pass_variants).to_csv(output_tsv, sep="\t", index=False, header=False)
    else:
        print(f"[Warning] {variant_id_col} column not found; skipping export of passed variants.")
    
    return output_tsv


def extract_pass_variants_by_intersection(
    pass_vmiss_path,
    pass_hwe_path,
    output_path="pass_variants.tsv",
    bed_prefix=None,
    output_prefix="filtered",
    threads=4,
    plink2="/home/b/b37974/plink2"
) -> str:
    """
    提取同时通过VMISS和HWE筛选的变体ID，并根据其染色体与位置排序后输出；
    可选地调用plink2从指定的bed_prefix基因型数据中提取这些变体对应的数据。

    参数:
    - pass_vmiss_path (str): VMISS通过的变体ID文件路径，每行一个变体ID
    - pass_hwe_path (str): HWE通过的变体ID文件路径，每行一个变体ID
    - output_path (str): 输出交集并排序后的变体ID文件路径（默认: pass_variants.tsv）
    - bed_prefix (str): 输入plink二进制文件前缀（.bed/.bim/.fam），用于提取变体
    - output_prefix (str): 输出plink提取结果的前缀名
    - threads (int): plink2运行时使用的线程数
    - plink2 (str): plink2可执行文件路径

    返回:
    - str: 排序后的交集变体ID文件路径
    """

    import tempfile
    import os
    import subprocess

    def parse_variant_id(vid) -> tuple[str, int, str]:
        fields = str(vid).split(":")
        if len(fields) < 2:
            return ("", -1, vid)
        chr_str = fields[0].replace("chr", "")
        try:
            pos = int(fields[1])
        except Exception:
            pos = -1
        return (chr_str, pos, vid)

    chr_order = [str(i) for i in range(1, 23)]

    # 第一步：将一个文件（较小者）载入内存建索引（假设HWE更小）
    hwe_ids = set()
    with open(pass_hwe_path) as f:
        for line in f:
            line = line.strip()
            if line:
                hwe_ids.add(line)

    # 第二步：逐行读取另一个文件（VMISS），筛选交集写入临时文件
    temp_fd, temp_path = tempfile.mkstemp()
    with os.fdopen(temp_fd, "w") as temp_out:
        with open(pass_vmiss_path) as f:
            for line in f:
                vid = line.strip()
                if vid in hwe_ids:
                    temp_out.write(f"{vid}\n")

    # 第三步：读取交集临时文件，排序并写入最终输出
    parsed = []
    with open(temp_path) as f:
        for line in f:
            vid = line.strip()
            chr_str, pos, _ = parse_variant_id(vid)
            if chr_str in chr_order:
                parsed.append((chr_order.index(chr_str), pos, vid))

    parsed.sort()
    with open(output_path, "w") as f_out:
        for _, _, vid in parsed:
            f_out.write(f"{vid}\n")

    # 删除中间临时文件
    os.remove(temp_path)

    # 可选：调用plink2提取基因型数据
    if bed_prefix is not None:
        plink_cmd = [
            plink2,
            "--bfile", bed_prefix,
            "--extract", output_path,
            "--make-bed",
            "--out", output_prefix,
            "--threads", str(threads)
        ]
        subprocess.run(plink_cmd, check=True)

    return output_path