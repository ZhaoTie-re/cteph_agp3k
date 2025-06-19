# %%
"""
脚本名称：基因型一致性与缺失率评估脚本

功能简介：
本脚本用于从给定的原始评估字典（pkl文件）中提取变体层面的缺失率表（vmiss_df）与基因型混淆矩阵（confusion_df），
并对不同平台（ALL, 15X, 30X）分别计算一系列基因型一致性评估指标，最终保存为结构化的 .summary.pkl 文件。

输入：
- raw_pkl: 每个染色体、过滤条件下生成的 .pkl 文件，包含：
    {
      frozenset([...]): (vmiss_df, confusion_df)
    }

输出：
- .summary.pkl: 每个平台分别输出：
    {
      frozenset([...]): (vmiss_df, confusion_df, summary_df)
    }

每个 summary_df 包含的统计指标包括：
- TOTAL_GENOTYPE, GENOTYPE_CONCORDANCE
- GENOTYPE_MISS_RATE, FALSE_POSITIVE_RATE, FALSE_NEGATIVE_RATE
- 各类错误类型的 COUNT 与 RATE（如 HET>HOMVAR）
- VMISS_FREQ 的均值、标准差、中位数
"""
import pandas as pd
import numpy as np
import pickle
import argparse

# ---参数设置---
parser = argparse.ArgumentParser(description="Genotype Concordance and VMISS Summary Script") 
parser.add_argument('--chr', type=str, required=True, help='Chromosome name (e.g., chr22)')
parser.add_argument('--global_dp', type=int, default=20, help='Global DP threshold')
parser.add_argument('--global_gq', type=int, default=30, help='Global GQ threshold')
parser.add_argument('--global_laf', type=float, default=0.25, help='Global LAF threshold')
parser.add_argument('--global_haf', type=float, default=0.75, help='Global HAF threshold')
parser.add_argument('--raw_pkl', type=str, required=True, help='Path to the raw .pkl file containing vmiss and confusion data')
args = parser.parse_args()  

# --- 参数设置 ---
chr = args.chr
global_dp = args.global_dp
global_gq = args.global_gq
global_laf = args.global_laf
global_haf = args.global_haf

raw_pkl = args.raw_pkl
raw_dict = pd.read_pickle(raw_pkl)

# --- 解包函数 ---
def unpack_concordance_vmiss_dict(raw_dict):
    all_vmiss, all_confusion = [], []

    for key, (vmiss_df, confusion_df) in raw_dict.items():
        key_strs = list(key)
        dp = int(next(s[2:] for s in key_strs if s.startswith('DP')))
        gq = int(next(s[2:] for s in key_strs if s.startswith('GQ')))
        laf = float(next(s[3:] for s in key_strs if s.startswith('LAF')))
        haf = float(next(s[3:] for s in key_strs if s.startswith('HAF')))
        platform = next(s for s in key_strs if s in ['ALL', '15X', '30X'])

        for df in (vmiss_df, confusion_df):
            df["DP"], df["GQ"], df["LAF"], df["HAF"], df["PLATFORM"] = dp, gq, laf, haf, platform

        all_vmiss.append(vmiss_df)
        all_confusion.append(confusion_df)

    return pd.concat(all_vmiss, ignore_index=True), pd.concat(all_confusion, ignore_index=True)

vmiss_df, confusion_df = unpack_concordance_vmiss_dict(raw_dict)

# --- 矩阵转换函数 ---
def confusion_to_matrix(confusion_df):
    ordered = ["0", "1", "2", "NA"]
    return (
        confusion_df.pivot_table(index="TRUE_GENOTYPE", columns="CALL_GENOTYPE", values="COUNT", aggfunc="sum", fill_value=0)
        .reindex(index=ordered, columns=ordered, fill_value=0)
    )

# --- 汇总函数 ---
import pandas as pd

def summarize_genotype_concordance(matrix, vmiss_df, global_dp, global_gq, global_laf, global_haf, platform):
    """
    根据给定的混淆矩阵和缺失率表格，汇总基因型一致性统计信息。

    参数：
        matrix: pd.DataFrame, 行为 TRUE_GENOTYPE，列为 CALL_GENOTYPE，元素为 COUNT
        global_dp, global_gq, global_laf, global_haf: 全局阈值参数
        platform: str, 平台名，如 "15X", "30X", "ALL"
        vmiss_df: pd.DataFrame or None, 包含 "VMISS_FREQ" 列

    返回：
        summary_df: 含统计指标的 DataFrame（1 行）
    """
    # 防止字符串型 index/columns
    matrix = matrix.copy()
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)

    # --- 总体数值 ---
    true_values = ["0", "1", "2"]
    call_values = ["0", "1", "2"]

    total = matrix.loc[true_values, call_values].values.sum()
    total_with_empty = matrix.values.sum() - matrix.loc["NA"].sum() if "NA" in matrix.index else matrix.values.sum()

    # --- 一致性 ---
    correct = sum(matrix.at[gt, gt] for gt in ["0", "1", "2"] if gt in matrix.index and gt in matrix.columns)
    concordance = correct / total if total > 0 else None
    concordance_with_empty = correct / total_with_empty if total_with_empty > 0 else None

    # --- 缺失率 ---
    na_calls = matrix.loc[true_values, "NA"].sum() if "NA" in matrix.columns else 0
    miss_rate = na_calls / total_with_empty if total_with_empty > 0 else None

    # --- 错误分类 ---
    het_to_homvar = matrix.at["1", "2"] if "1" in matrix.index and "2" in matrix.columns else 0
    homref_to_homvar = matrix.at["0", "2"] if "0" in matrix.index and "2" in matrix.columns else 0
    homref_to_het = matrix.at["0", "1"] if "0" in matrix.index and "1" in matrix.columns else 0

    het_to_homvar_rate = het_to_homvar / total if total > 0 else None
    homref_to_homvar_rate = homref_to_homvar / total if total > 0 else None
    homref_to_het_rate = homref_to_het / total if total > 0 else None

    # --- 假阳性与假阴性 ---
    # 假阳性: CALL 为变异 (1 or 2)，TRUE 为非对应的低等位基因（0 或 1）
    false_positive = (
        matrix.at["0", "1"] if ("0" in matrix.index and "1" in matrix.columns) else 0
    ) + (
        matrix.at["0", "2"] if ("0" in matrix.index and "2" in matrix.columns) else 0
    ) + (
        matrix.at["1", "2"] if ("1" in matrix.index and "2" in matrix.columns) else 0
    )

    # 假阴性: CALL 为 REF 或较低等级的等位基因，TRUE 为变异 (1 or 2)
    false_negative = (
        matrix.at["1", "0"] if ("1" in matrix.index and "0" in matrix.columns) else 0
    ) + (
        matrix.at["2", "0"] if ("2" in matrix.index and "0" in matrix.columns) else 0
    ) + (
        matrix.at["2", "1"] if ("2" in matrix.index and "1" in matrix.columns) else 0
    )

    false_positive_rate = false_positive / total if total > 0 else None
    false_negative_rate = false_negative / total if total > 0 else None
    
    # --- VMISS 统计 ---
    vmiss_mean = vmiss_df["VMISS_FREQ"].mean() if vmiss_df is not None else None
    vmiss_sd = vmiss_df["VMISS_FREQ"].std() if vmiss_df is not None else None
    vmiss_median = vmiss_df["VMISS_FREQ"].median() if vmiss_df is not None else None

    # --- 汇总为 DataFrame ---
    summary_df = pd.DataFrame([{
        'DP': global_dp,
        'GQ': global_gq,
        'LAF': global_laf,
        'HAF': global_haf,
        'PLATFORM': platform,
        'TOTAL_GENOTYPE': int(total),
        'TOTAL_GENOTYPE(WITH_EMPTY)': int(total_with_empty),
        'TOTAL_GENOTYPE(WITH_ALL_NA)': int(matrix.values.sum()),
        'GENOTYPE_CONCORDANCE': concordance,
        'GENOTYPE_CONCORDANCE(WITH_EMPTY)': concordance_with_empty,
        'GENOTYPE_MISS_RATE': miss_rate,
        'FALSE_POSITIVE_RATE': false_positive_rate,
        'FALSE_NEGATIVE_RATE': false_negative_rate,
        'HET>HOMVAR_COUNT': int(het_to_homvar),
        'HOMREF>HOMVAR_COUNT': int(homref_to_homvar),
        'HOMREF>HET_COUNT': int(homref_to_het),
        'HET>HOMVAR_RATE': het_to_homvar_rate,
        'HOMREF>HOMVAR_RATE': homref_to_homvar_rate,
        'HOMREF>HET_RATE': homref_to_het_rate,
        'VMISS_FREQ_MEAN': vmiss_mean,
        'VMISS_FREQ_SD': vmiss_sd,
        'VMISS_FREQ_MEDIAN': vmiss_median,
    }])

    return summary_df

# --- 分平台处理 ---
platforms = ['ALL', '15X', '30X']
results_dict = {}

for platform in platforms:
    vmiss_sub = vmiss_df[vmiss_df["PLATFORM"] == platform]
    confusion_sub = confusion_df[confusion_df["PLATFORM"] == platform]
    matrix = confusion_to_matrix(confusion_sub)
    summary = summarize_genotype_concordance(matrix, vmiss_sub, global_dp, global_gq, global_laf, global_haf, platform)

    key = frozenset([f'DP{global_dp}', f'GQ{global_gq}', f'LAF{global_laf}', f'HAF{global_haf}', platform])
    results_dict[key] = (vmiss_sub, confusion_sub, summary)

# --- 写出结果 ---
summary_path = f'{chr}._DP{global_dp}_GQ{global_gq}_LAF{global_laf}_HAF{global_haf}_.summary.pkl'
with open(summary_path, 'wb') as f:
    pickle.dump(results_dict, f)


