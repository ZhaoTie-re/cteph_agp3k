"""
sample_qc_flags.py

样本质量控制（Sample QC）工具集

本脚本包含用于遗传学数据分析中样本级别质量控制（QC）的多个实用函数，主要功能包括：
  
1. 绘制 Heterozygosity F 与 Sample Missing Rate 的分布图，辅助识别异常样本；
2. 绘制 Mean DP 的 robust Z-score 与 Heterozygosity F 的分布图；
3. 绘制 PI_HAT 分布图并进行亲缘网络可视化，识别冗余样本；
4. 根据设定阈值生成每个样本是否通过各项 QC 检查的布尔标记表。

每个函数都提供独立的输入参数与 PDF 输出结果，适用于 QC 自动化流程。

作者: ZHAO TIE
"""
# ===== 标准库导入 =====
import re
import logging
from pathlib import Path

# ===== 第三方库导入 =====
import numpy as np
import pandas as pd
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch, Circle, Rectangle
from networkx.drawing.nx_agraph import graphviz_layout
from networkx.algorithms.approximation import min_weighted_vertex_cover

# ===== 日志配置 =====
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')


def plot_pi_hat_distribution(pihat_df, pi_threshold, case_prefix, out_prefix):
    """
    绘制 PI_HAT 值的分布图，区分样本对的类型（case-case, control-control, case-control）。
    
    参数：
    pihat_df (pd.DataFrame): 包含样本对及其 PI_HAT 值的数据框，必须包含列 'IID1', 'IID2', 'PI_HAT'。
                             该数据框应已过滤，确保所有 PI_HAT 值均大于 pi_threshold，无需在函数内再次筛选。
    pi_threshold (float): PI_HAT 的阈值，绘图时只考虑大于该阈值的样本对。
    case_prefix (str): 用于判定样本是否为 case 的前缀字符串。
    out_prefix (str): 输出文件名前缀。
    
    功能：
    根据样本ID判断样本身份，划分样本对类型，绘制不同类型样本对的 PI_HAT 分布直方图，并保存成 PDF 文件。
    """
    logging.info("开始绘制 PI_HAT 分布图...")

    # ==== 判断每个样本是 case 还是 control ====
    def get_case_status(iid):
        return "case" if iid.startswith(case_prefix) else "control"

    pihat_df = pihat_df.copy()
    pihat_df["IID1_status"] = pihat_df["IID1"].apply(get_case_status)
    pihat_df["IID2_status"] = pihat_df["IID2"].apply(get_case_status)

    # ==== 判断每对样本属于哪种组合类型 ====
    def kinship_type(row):
        if row["IID1_status"] == "case" and row["IID2_status"] == "case":
            return "case-case"
        elif row["IID1_status"] == "control" and row["IID2_status"] == "control":
            return "control-control"
        else:
            return "case-control"

    pihat_df["pair_type"] = pihat_df.apply(kinship_type, axis=1)

    # ==== 设置学术配色与统一 bin 边界 ====
    colors = {
        "case-case": "#D55E00",        # 橙红
        "control-control": "#0072B2",  # 蓝色
        "case-control": "#009E73"      # 绿色
    }
    bins = np.linspace(pi_threshold, 1.0, num=51)

    # ==== 绘图 ====
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharex=True)

    for i, pair_type in enumerate(["case-case", "control-control", "case-control"]):
        ax = axes[i]
        subset = pihat_df[
            (pihat_df["pair_type"] == pair_type)
        ]
        iids = pd.concat([
            subset[["IID1", "IID1_status"]].rename(columns={"IID1": "IID", "IID1_status": "status"}),
            subset[["IID2", "IID2_status"]].rename(columns={"IID2": "IID", "IID2_status": "status"})
        ])
        unique_iids = iids.drop_duplicates("IID")
        counts = unique_iids["status"].value_counts()
        n_case = counts.get("case", 0)
        n_ctrl = counts.get("control", 0)
        ax.hist(subset["PI_HAT"], bins=bins, color=colors[pair_type], edgecolor='black')
        ax.set_title(
            f"{pair_type} pairs\nn_case = {n_case}, n_ctrl = {n_ctrl}",
            fontsize=12,
        )
        ax.set_xlabel("PI_HAT", fontsize=11)
        ax.set_ylabel("Number of Pairs" if i == 0 else "", fontsize=11)
        ax.set_xlim(pi_threshold, 1)
        ax.tick_params(axis="x", labelrotation=45)

    # ==== 添加主标题 ====
    fig.suptitle(
        f"Distribution of PI_HAT (Filtered: PI_HAT > {pi_threshold:.2f})",
        fontsize=15,
        y=1.08,
        weight="bold"
    )
    plt.tight_layout()
    plt.style.use('default')

    # ==== 保存图像 ====
    out_path = Path(f"{out_prefix}.pi_hat{pi_threshold:.2f}.hist.pdf")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_kinship_network_and_prune(pihat_df, sample_qc_summary, pi_threshold, case_prefix, out_prefix):
    """
    构建并可视化基于 PI_HAT 阈值的亲缘关系网络，利用最小加权顶点覆盖算法进行节点剪枝。
    
    参数：
    pihat_df (pd.DataFrame): 包含样本对及其 PI_HAT 值的数据框，必须包含列 'IID1', 'IID2', 'PI_HAT'。
                             该数据框应已过滤，确保所有 PI_HAT 值均大于 pi_threshold，无需在函数内再次筛选。
    sample_qc_summary (pd.DataFrame): 样本质量控制汇总表，包含至少 'IID' 和 'SMISS' 列。
    pi_threshold (float): PI_HAT 的阈值，网络构建时仅考虑大于该阈值的样本对。
    case_prefix (str): 用于判定样本是否为 case 的前缀字符串。
    out_prefix (str): 输出文件名前缀。
    
    返回：
    list: 被标记为删除的节点列表（IID），按字母序排序。
    
    功能：
    1. 根据样本身份构建带权无向图，权重结合样本缺失率与身份（case/control）。
    2. 计算最小加权顶点覆盖，确定需删除的节点以消除高亲缘关系冲突。
    3. 绘制剪枝前后网络图，节点颜色区分身份，删除节点用红圈标记。
    4. 保存网络图至 PDF 文件。
    """
    logging.info("开始构建并绘制亲缘关系网络...")

    # ==== 设置绘图风格 ====
    plt.style.use("default")

    # ---- 构建图并添加节点，节点标记身份（case/control） ----
    G_all = nx.Graph()
    for iid in sorted(set(pihat_df["IID1"]).union(set(pihat_df["IID2"]))):
        status = "case" if iid.startswith(case_prefix) else "control"
        G_all.add_node(iid, status=status)

    # ---- 添加边，确保边按字母序排序以保证一致性 ----
    edges_sorted = sorted((min(a, b), max(a, b)) for a, b in zip(pihat_df["IID1"], pihat_df["IID2"]))
    G_all.add_edges_from(edges_sorted)

    # ==== 构建样本缺失率字典，用于节点权重计算 ====
    smiss_dict = dict(zip(sample_qc_summary['IID'], sample_qc_summary['SMISS']))

    # ==== 赋予节点权重，case 权重显著高于 control，且与缺失率成反比 ====
    for n in G_all.nodes:
        fmiss = smiss_dict.get(n, 0.001)  # 避免除以零
        if G_all.nodes[n]["status"] == "case":
            G_all.nodes[n]["weight"] = 1e6 / fmiss
        else:
            G_all.nodes[n]["weight"] = 1.0 / fmiss

    # ==== 计算最小加权顶点覆盖，确定需删除的节点集合 ====
    vc = min_weighted_vertex_cover(G_all, weight="weight")
    to_remove = set(vc)
    remaining_nodes = set(G_all.nodes) - to_remove

    # ==== 计算节点布局及颜色映射 ====
    pos = graphviz_layout(G_all, prog="neato")
    color_map = {"case": "#E64B35", "control": "#4DBBD5"}
    node_colors_before = [color_map[G_all.nodes[n]["status"]] for n in G_all.nodes]
    node_colors_after = [color_map[G_all.nodes[n]["status"]] for n in remaining_nodes]

    # ==== 构建图例元素 ====
    legend_elements = [
        Patch(facecolor='#E64B35', label='Case'),
        Patch(facecolor='#4DBBD5', label='Control'),
        Patch(edgecolor='red', facecolor='none', label='Deleted', linewidth=1.5)
    ]

    # ==== 计算绘图边界，添加适当边距 ====
    x_vals = [p[0] for p in pos.values()]
    y_vals = [p[1] for p in pos.values()]
    pad = 20
    x_min, x_max = min(x_vals) - pad, max(x_vals) + pad
    y_min, y_max = min(y_vals) - pad, max(y_vals) + pad

    # ==== 绘制剪枝前后网络图 ====
    fig, axes = plt.subplots(1, 2, figsize=(18, 10))

    # 统计节点数量信息
    total_nodes = len(G_all.nodes)
    deleted_nodes = len(to_remove)
    remaining_nodes_count = len(remaining_nodes)

    # --- 左图：剪枝前网络 ---
    ax = axes[0]
    ax.set_title("Before Pruning", fontsize=16, fontweight='bold')
    nx.draw_networkx_nodes(G_all, pos, node_color=node_colors_before, node_size=60, alpha=0.95, ax=ax)
    nx.draw_networkx_edges(G_all, pos, edge_color="#444444", width=1.2, alpha=0.5, ax=ax)
    for node in to_remove:
        if node in pos:
            x, y = pos[node]
            ax.add_patch(Circle((x, y), radius=12, fill=False, edgecolor='red', linewidth=1.5))
    ax.add_patch(Rectangle((x_min, y_min), x_max - x_min, y_max - y_min,
                           fill=False, edgecolor='gray', linewidth=1.5))
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.axis('off')
    ax.text(
        0.01, 0.01,
        f"Total Nodes: {total_nodes}\nNodes Marked for Deletion: {deleted_nodes}",
        transform=ax.transAxes,
        ha='left', va='bottom',
        fontsize=12,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=1)
    )

    # --- 右图：剪枝后网络 ---
    ax = axes[1]
    ax.set_title("After Pruning", fontsize=16, fontweight='bold')
    G_final = nx.Graph()
    for n in remaining_nodes:
        G_final.add_node(n, status=G_all.nodes[n]["status"])
    nx.draw_networkx_nodes(G_final, pos, node_color=node_colors_after, node_size=60, alpha=0.95, ax=ax)
    ax.add_patch(Rectangle((x_min, y_min), x_max - x_min, y_max - y_min,
                           fill=False, edgecolor='gray', linewidth=1.5))
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.axis('off')
    ax.text(
        0.99, 0.01,
        f"Remaining Nodes: {remaining_nodes_count}",
        transform=ax.transAxes,
        ha='right', va='bottom',
        fontsize=12,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=1)
    )

    # --- 总标题与图例 ---
    fig.suptitle(
        "Kinship Network Pruning via Minimum Weighted Vertex Cover\n"
        f"Filtered: PI_HAT > {pi_threshold:.2f}",
        fontsize=20, fontweight='bold', y=1.0
    )
    fig.legend(
        handles=legend_elements,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.93),
        ncol=3,
        fontsize=16,
        frameon=False
    )
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    out_path = Path(f"{out_prefix}.pi_hat{pi_threshold:.2f}.net.pdf")
    fig.savefig(out_path, bbox_inches='tight', dpi=300)
    plt.close()
    return sorted(to_remove)


def plot_hetf_vs_meandp_table(sample_qc_summary, dp_robust_z_threshold, het_threshold, case_prefix, out_prefix):
    """
    可视化 Heterozygosity F 与 Mean DP 的关系，并按平台和离群类型进行分类展示，附带样本数量统计表格。

    参数：
    sample_qc_summary (pd.DataFrame): 样本质量控制汇总数据，必须包含列 ['IID', 'HET_F', 'ROBUST_Z_DP', 'TARGET_DP']。
    dp_robust_z_threshold (float): 判断 DP 异常的 robust Z-score 阈值（例如 -3.0）。
    het_threshold (str): HET_F 的判定标准，格式为 "{n}sd"，如 "5sd" 表示均值±5倍标准差。
    case_prefix (str): 用于识别 Case 样本的 ID 前缀。
    out_prefix (str): 输出 PDF 文件的路径前缀。

    输出：
    生成 PDF 文件：{out_prefix}.hetF_vs_DP.dp_z{阈值}_f{het_threshold}.pdf

    功能：
    1. 按平台区分绘制 Mean DP 的 KDE 分布；
    2. 按 HET_F 和 DP 异常标记样本分类，并绘制散点图；
    3. 构建各分类与平台组合的 Case / Ctrl / Total 样本数表格；
    4. 可视化所有信息并导出 PDF 图像。
    """
    logging.info("开始绘制 Heterozygosity F 与 Mean DP 关系图及表格...")

    # === 色系映射（QC_Class + 平台）===
    custom_palette = {
        'Pass (15x)': '#5B9BD5', 'Pass (30x)': '#A9CDEB',
        'Only F Outlier (15x)': '#228B22', 'Only F Outlier (30x)': '#7FC97F',
        'Only DP Outlier (15x)': '#FF8C00', 'Only DP Outlier (30x)': '#FFD180',
        'Both Outlier (15x)': '#B22222', 'Both Outlier (30x)': '#E57373',
    }

    # === 点形状映射（按平台）===
    marker_dict = {'15x': 'o', '30x': 's'}

    # === HET_F 阈值计算 ===
    # 解析 het_threshold 格式（如 "5sd"），提取数值
    match = re.match(r"(\d+(?:\.\d+)?)sd", het_threshold)
    if match:
        het_threshold_value = int(match.group(1))
    else:
        raise ValueError(f"Invalid het_threshold format: {het_threshold}")

    mean_f = sample_qc_summary['HET_F'].mean()
    std_f = sample_qc_summary['HET_F'].std()
    f_upper = mean_f + het_threshold_value * std_f
    f_lower = mean_f - het_threshold_value * std_f

    # === 标记 Case / Control ===
    sample_qc_summary['GROUP'] = sample_qc_summary['IID'].apply(
        lambda x: 'Case' if str(x).startswith(case_prefix) else 'Control'
    )

    # === 分类 QC Class ===
    # 判断 HET_F 是否为离群值
    cond_f_outlier = (sample_qc_summary['HET_F'] > f_upper) | (sample_qc_summary['HET_F'] < f_lower)
    # 判断 DP 是否为离群值
    cond_dp_outlier = sample_qc_summary['ROBUST_Z_DP'] < dp_robust_z_threshold
    only_f_outlier = cond_f_outlier & (~cond_dp_outlier)
    only_dp_outlier = (~cond_f_outlier) & cond_dp_outlier
    both_outlier = cond_f_outlier & cond_dp_outlier

    def classify(row):
        # 按优先级分类
        if both_outlier.loc[row.name]:
            return "Both Outlier"
        elif only_f_outlier.loc[row.name]:
            return "Only F Outlier"
        elif only_dp_outlier.loc[row.name]:
            return "Only DP Outlier"
        else:
            return "Pass"

    sample_qc_summary['QC_Class'] = sample_qc_summary.apply(classify, axis=1)
    sample_qc_summary['TARGET_DP'] = sample_qc_summary['TARGET_DP'].astype(str)
    sample_qc_summary['QC_Class_Plat'] = sample_qc_summary['QC_Class'] + " (" + sample_qc_summary['TARGET_DP'] + ")"

    # === 构建 Case / Control 数表格 ===
    qc_classes = ['Pass', 'Only F Outlier', 'Only DP Outlier', 'Both Outlier']
    target_dps = sorted(sample_qc_summary['TARGET_DP'].dropna().unique())
    all_combinations = [f"{qc} ({dp})" for qc in qc_classes for dp in target_dps]

    table_data = []
    row_labels = []

    for label in all_combinations:
        qc_class, dp = re.match(r"(.+?) \((.+)\)", label).groups()
        subset = sample_qc_summary[
            (sample_qc_summary['QC_Class'] == qc_class) &
            (sample_qc_summary['TARGET_DP'] == dp)
        ]
        case_count = (subset['GROUP'] == 'Case').sum()
        ctrl_count = (subset['GROUP'] == 'Control').sum()
        total = len(subset)
        row_labels.append(label)
        table_data.append([case_count, ctrl_count, total])

    table_df = pd.DataFrame(table_data, columns=['Case', 'Ctrl', 'Total'], index=row_labels)

    # === 主图绘制 ===
    plt.style.use('default')
    fig = plt.figure(figsize=(14, 6))
    gs = gridspec.GridSpec(2, 2, width_ratios=[4.5, 2.0], height_ratios=[1, 4], wspace=0.6, hspace=0.1)

    # KDE 图：按平台展示 MEAN_DP 分布
    ax0 = fig.add_subplot(gs[0, 0])
    for target_dp, group_df in sample_qc_summary.groupby("TARGET_DP"):
        valid_dp = group_df["MEAN_DP"].dropna()
        if len(valid_dp) > 0:
            linestyle = '--' if target_dp == '15x' else '-' if target_dp == '30x' else '-.'
            sns.kdeplot(
                valid_dp,
                ax=ax0,
                label=str(target_dp),
                color='black',
                linestyle=linestyle,
                linewidth=2.0
            )
    ax0.set_ylabel('')
    ax0.set_xlim(left=0)
    ax0.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax0.set_xlabel('')
    ax0.set_yticks([])
    ax0.legend(title='')

    # 散点图：HET_F vs MEAN_DP，按分类和平台着色
    ax1 = fig.add_subplot(gs[1, 0])
    for label in all_combinations:
        qc_class, dp = re.match(r"(.+?) \((.+)\)", label).groups()
        subset = sample_qc_summary[sample_qc_summary['QC_Class_Plat'] == label]
        marker = marker_dict.get(dp, 'o')
        ax1.scatter(
            subset['MEAN_DP'],
            subset['HET_F'],
            color=custom_palette.get(label, 'gray'),
            edgecolor='black',
            marker=marker,
            s=20,
            linewidth=0.5,
            alpha=0.8
        )
    ax1.axhline(f_upper, color='darkred', linestyle='--', lw=1)
    ax1.axhline(f_lower, color='darkred', linestyle='--', lw=1)
    ax1.set_xlabel('Mean DP', fontsize=11)
    ax1.set_ylabel('Heterozygosity F', fontsize=11)
    ax1.grid(True, linestyle='--', linewidth=0.5)

    # 表格：Case / Control 数量统计
    ax2 = fig.add_subplot(gs[1, 1])
    ax2.axis('off')
    table = ax2.table(
        cellText=table_df.values,
        rowLabels=table_df.index,
        colLabels=table_df.columns,
        cellLoc='center',
        loc='center'
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.3, 1.3)

    # 设置统一列宽
    max_text_length = max(
        max(table_df.astype(str).applymap(len).max()),
        max([len(str(col)) for col in table_df.columns])
    )
    def calc_width(char_len, scale=0.07):
        return char_len * scale
    uniform_width = calc_width(max_text_length)
    for (row, col), cell in table.get_celld().items():
        if col < len(table_df.columns):
            cell.set_width(uniform_width)

    # 设置表头和行标签样式
    for col in range(len(table_df.columns)):
        header_cell = table[0, col]
        header_cell.set_facecolor('#D3D3D3')
        header_cell.set_text_props(weight='bold', color='black')
    for i, label in enumerate(table_df.index):
        row_label_cell = table[i + 1, -1]
        row_label_cell.set_facecolor(to_rgba(custom_palette.get(label, 'gray'), alpha=0.5))
        row_label_cell.set_text_props(weight='bold', color='black')
        table[i + 1, 0].set_x(0.45)
    for row in range(1, len(table_df) + 1):
        for col in range(len(table_df.columns)):
            cell = table[row, col]
            cell.set_facecolor('#F5F5F5')
            cell.set_text_props(color='black')

    # 图标题与导出
    fig.suptitle(
        f'Heterozygosity F vs. Mean DP\n(DP Outlier: Robust Z < {dp_robust_z_threshold}, '
        f'F threshold: ±{het_threshold_value}SD)',
        fontsize=13, fontweight='bold'
    )
    suffix = f".hetF_vs_DP.dp_z{'neg' if dp_robust_z_threshold < 0 else ''}{abs(dp_robust_z_threshold):.1f}_f{het_threshold}.pdf"
    out_path = Path(f"{out_prefix}{suffix}")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    


def plot_hetf_vs_smiss_table(sample_qc_summary, smiss_threshold, het_threshold, case_prefix, out_prefix):
    """
    可视化 Heterozygosity F 与 Sample Missing Rate (SMISS) 的关系，并按平台和离群类型进行分类展示，附带样本数量统计表格。

    参数：
    sample_qc_summary (pd.DataFrame): 样本质量控制汇总数据，必须包含列 ['IID', 'HET_F', 'SMISS', 'TARGET_DP']。
    smiss_threshold (float): SMISS 判定为异常的阈值。
    het_threshold (str): HET_F 的判定标准，格式为 "{n}sd"，如 "5sd" 表示均值±5倍标准差。
    case_prefix (str): 用于识别 Case 样本的 ID 前缀。
    out_prefix (str): 输出 PDF 文件的路径前缀。

    功能：
    1. 按平台区分绘制 log10(SMISS) 的 KDE 分布；
    2. 按 HET_F 和 SMISS 离群标记样本分类，并绘制散点图；
    3. 构建各分类与平台组合的 Case / Ctrl / Total 样本数表格；
    4. 可视化所有信息并导出 PDF 图像。
    """
    logging.info("开始绘制 Heterozygosity F 与 SMISS 散点图及表格...")

    # === 色系映射（QC_Class + 平台）===
    custom_palette = {
        'Pass (15x)': '#5B9BD5', 'Pass (30x)': '#A9CDEB',
        'Only F Outlier (15x)': '#228B22', 'Only F Outlier (30x)': '#7FC97F',
        'Only SMISS Outlier (15x)': '#FF8C00', 'Only SMISS Outlier (30x)': '#FFD180',
        'Both Outlier (15x)': '#B22222', 'Both Outlier (30x)': '#E57373',
    }

    # === 点形状映射（按平台）===
    marker_dict = {
        '15x': 'o',
        '30x': 's',
    }

    # === 提取 het_threshold 数值 ===
    match = re.match(r"(\d+(?:\.\d+)?)sd", het_threshold)
    if match:
        het_threshold_value = int(match.group(1))
    else:
        raise ValueError(f"Invalid het_threshold format: {het_threshold}")

    # === F 值阈值计算 ===
    mean_f = sample_qc_summary['HET_F'].mean()
    std_f = sample_qc_summary['HET_F'].std()
    f_upper = mean_f + het_threshold_value * std_f
    f_lower = mean_f - het_threshold_value * std_f

    # === 添加组别标签（Case / Control）===
    sample_qc_summary['GROUP'] = sample_qc_summary['IID'].apply(
        lambda x: 'Case' if str(x).startswith(case_prefix) else 'Control'
    )

    sample_qc_summary['TARGET_DP'] = sample_qc_summary['TARGET_DP'].astype(str)

    # === 分类条件定义 ===
    cond_f_outlier = (sample_qc_summary['HET_F'] > f_upper) | (sample_qc_summary['HET_F'] < f_lower)
    cond_smiss_outlier = sample_qc_summary['SMISS'] > smiss_threshold
    only_f_outlier = cond_f_outlier & (~cond_smiss_outlier)
    only_smiss_outlier = (~cond_f_outlier) & cond_smiss_outlier
    both_outlier = cond_f_outlier & cond_smiss_outlier

    # === 分类函数 ===
    def classify(row):
        if both_outlier.loc[row.name]:
            return "Both Outlier"
        elif only_f_outlier.loc[row.name]:
            return "Only F Outlier"
        elif only_smiss_outlier.loc[row.name]:
            return "Only SMISS Outlier"
        else:
            return "Pass"

    sample_qc_summary['QC_Class'] = sample_qc_summary.apply(classify, axis=1)
    sample_qc_summary['QC_Class_Plat'] = sample_qc_summary['QC_Class'] + " (" + sample_qc_summary['TARGET_DP'] + ")"

    # === 构造表格数据（确保所有组合都包括）===
    qc_classes = ['Pass', 'Only F Outlier', 'Only SMISS Outlier', 'Both Outlier']
    platforms = sorted(sample_qc_summary['TARGET_DP'].unique())
    all_combinations = [f"{qc} ({dp})" for qc in qc_classes for dp in platforms]

    table_data = []
    row_labels = []

    for label in all_combinations:
        qc_class, dp = re.match(r"(.+?) \((.+)\)", label).groups()
        subset = sample_qc_summary[
            (sample_qc_summary['QC_Class'] == qc_class) &
            (sample_qc_summary['TARGET_DP'] == dp)
        ]
        case_count = (subset['GROUP'] == 'Case').sum()
        ctrl_count = (subset['GROUP'] == 'Control').sum()
        total = len(subset)
        row_labels.append(label)
        table_data.append([case_count, ctrl_count, total])

    table_df = pd.DataFrame(table_data, columns=['Case', 'Ctrl', 'Total'], index=row_labels)

    # === 图像绘制 ===
    plt.style.use('default')
    fig = plt.figure(figsize=(14, 6))
    gs = gridspec.GridSpec(2, 2, width_ratios=[4.5, 2.0], height_ratios=[1, 4], wspace=0.6, hspace=0.1)

    # === KDE subplot (按平台分组绘制 log10 SMISS) ===
    ax0 = fig.add_subplot(gs[0, 0])
    for dp, group_df in sample_qc_summary.groupby("TARGET_DP"):
        smiss_vals = group_df['SMISS']
        smiss_vals_log = np.log10(smiss_vals[smiss_vals > 0])
        if len(smiss_vals_log) > 0:
            # 设置线型
            linestyle = '--' if dp == '15x' else '-' if dp == '30x' else '-.'
            sns.kdeplot(
                smiss_vals_log,
                ax=ax0,
                label=dp,
                color='black',
                linestyle=linestyle,
                linewidth=2.0
            )

    ax0.set_ylabel('')
    # 若所有平台都没 SMISS>0，则 smiss_vals_log 为空，避免报错
    if not sample_qc_summary['SMISS'][sample_qc_summary['SMISS'] > 0].empty:
        ax0.set_xlim(left=np.log10(sample_qc_summary['SMISS'][sample_qc_summary['SMISS'] > 0].min()))
    ax0.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax0.set_xlabel('')
    ax0.set_yticks([])
    ax0.legend(title='')

    # === 主图散点图（颜色 + 形状）===
    ax1 = fig.add_subplot(gs[1, 0])
    for label in all_combinations:
        qc_class, dp = re.match(r"(.+?) \((.+)\)", label).groups()
        subset = sample_qc_summary[sample_qc_summary['QC_Class_Plat'] == label]
        marker = marker_dict.get(dp, 'o')
        ax1.scatter(
            subset['SMISS'],
            subset['HET_F'],
            color=custom_palette.get(label, 'gray'),
            edgecolor='black',
            marker=marker,
            s=20,
            linewidth=0.5,
            alpha=0.8
        )

    ax1.axvline(smiss_threshold, color='darkred', linestyle='--', lw=1)
    ax1.axhline(f_upper, color='darkred', linestyle='--', lw=1)
    ax1.axhline(f_lower, color='darkred', linestyle='--', lw=1)
    ax1.set_xscale('log')
    # ax1.set_xticks([0.01, 0.1])
    ax1.set_xlabel('Sample Missing Rate (log scale)', fontsize=11)
    ax1.set_ylabel('Heterozygosity F', fontsize=11)
    ax1.grid(True, linestyle='--', linewidth=0.5)

    # === 表格 subplot ===
    ax2 = fig.add_subplot(gs[1, 1])
    ax2.axis('off')

    table = ax2.table(
        cellText=table_df.values,
        rowLabels=table_df.index,
        colLabels=table_df.columns,
        cellLoc='center',
        loc='center'
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.3, 1.3)

    # 列宽调整
    max_text_length = max(
        max(table_df.astype(str).applymap(len).max()),
        max([len(str(col)) for col in table_df.columns])
    )
    def calc_width(char_len, scale=0.07):
        return char_len * scale
    uniform_width = calc_width(max_text_length)
    for (row, col), cell in table.get_celld().items():
        if col < len(table_df.columns):
            cell.set_width(uniform_width)

    # 表头样式
    for col in range(len(table_df.columns)):
        header_cell = table[0, col]
        header_cell.set_facecolor('#D3D3D3')
        header_cell.set_text_props(weight='bold', color='black')

    # 行标签染色 + 对齐
    for i, label in enumerate(table_df.index):
        row_label_cell = table[i + 1, -1]
        row_label_cell.set_facecolor(to_rgba(custom_palette.get(label, 'gray'), alpha=0.5))
        row_label_cell.set_text_props(weight='bold', color='black')
        table[i + 1, 0].set_x(0.45)

    # 内容单元格底纹
    for row in range(1, len(table_df) + 1):
        for col in range(len(table_df.columns)):
            cell = table[row, col]
            cell.set_facecolor('#F5F5F5')
            cell.set_text_props(color='black')

    # === 图标题 ===
    fig.suptitle(
        f'Heterozygosity F vs. Sample Missing Rate\n'
        f'(SMISS threshold: {smiss_threshold}, F threshold: ±{het_threshold_value}SD)',
        fontsize=13, fontweight='bold'
    )
    # === 保存 PDF 文件 ===
    out_path = Path(f"{out_prefix}.hetF_vs_SMISS.smiss{smiss_threshold}_f{het_threshold}.pdf")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    

def generate_sample_qc_flags(sample_qc_summary, smiss_threshold, dp_robust_z_threshold, het_threshold, pi_outlier, out_prefix):
    """
    生成样本级别的 QC 标志（True/False），并保存为 CSV 文件。

    参数:
    - sample_qc_summary: pd.DataFrame, 包含样本信息，需包含列 ['#FID', 'IID', 'SMISS', 'ROBUST_Z_DP', 'HET_F']
    - smiss_threshold: float, SMISS 的阈值（如 0.1），大于该值的样本设为 False
    - dp_robust_z_threshold: float, ROBUST_Z_DP 的下限阈值（如 -3.0），低于该值设为 False
    - het_threshold: str, 形如 "5sd"，表示 HET_F 的标准差范围阈值
    - pi_outlier: set[str], PI_HAT 离群样本集合（已通过亲缘网络分析生成的样本 ID 集合）
    - out_prefix: str, 输出文件名前缀，将保存为 {out_prefix}.sample_qc_flags.csv

    返回:
    - pd.DataFrame: 包含 ['#FID', 'IID', 'PASS_SMISS', 'PASS_MEAN_DP', 'PASS_HET_F', 'PASS_PI_HAT']
    """
    logging.info("开始生成样本 QC 标志并保存为 CSV ...")

    # 提取 HET_F 阈值范围
    match = re.match(r"(\d+(?:\.\d+)?)sd", het_threshold)
    if not match:
        raise ValueError(f"无效的 het_threshold 格式: {het_threshold}，应形如 '5sd'")
    sd_value = float(match.group(1))

    het_f_mean = sample_qc_summary["HET_F"].mean()
    het_f_std = sample_qc_summary["HET_F"].std()
    het_f_upper = het_f_mean + sd_value * het_f_std
    het_f_lower = het_f_mean - sd_value * het_f_std

    # 创建输出 DataFrame
    flags_df = pd.DataFrame()
    flags_df["#FID"] = sample_qc_summary["#FID"]
    flags_df["IID"] = sample_qc_summary["IID"]
    flags_df["PASS_SMISS"] = sample_qc_summary["SMISS"] <= smiss_threshold
    flags_df["PASS_MEAN_DP"] = sample_qc_summary["ROBUST_Z_DP"] >= dp_robust_z_threshold
    flags_df["PASS_HET_F"] = sample_qc_summary["HET_F"].between(het_f_lower, het_f_upper)
    flags_df["PASS_PI_HAT"] = ~sample_qc_summary["IID"].isin(pi_outlier)

    # 保存文件
    out_path = Path(f"{out_prefix}.sample_qc_flags.csv")
    flags_df.to_csv(out_path, index=False)
    return flags_df

# ===== main 流程模板（可选入口）=====
if __name__ == "__main__":
    logging.info("此脚本为样本QC工具函数集，请在主脚本中导入使用。")