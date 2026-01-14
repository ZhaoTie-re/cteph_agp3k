# %%
"""
脚本名称：基因型一致性与缺失率汇总脚本

功能简介：
本脚本用于从每个染色体在特定过滤条件下生成的评估结果字典（pkl文件）中，
解包并整理三类表格信息：
  - 变体缺失率表（vmiss_df）
  - 样本缺失率表（smiss_df）
  - 基因型混淆矩阵（confusion_df）

并分别对三个测序平台（ALL、15X、30X）计算统计指标，包括：
  - 基因型总数
  - 一致性、缺失率、假阳性率、假阴性率
  - 错配错误分类统计与比例
  - VMISS 与 SMISS 的频率均值、标准差、中位数

输入参数：
  --chr: 染色体名称（如 chr22）
  --global_dp: 全局DP阈值
  --global_gq: 全局GQ阈值
  --global_laf: 全局LAF阈值
  --global_haf: 全局HAF阈值
  --raw_pkl: 输入的原始pkl文件，内容结构为：
    {
      frozenset([...]): (vmiss_df, smiss_df, confusion_df)
    }

输出：
  - .summary.pkl 文件，结构为：
    {
      frozenset([...]): (vmiss_sub, smiss_sub, confusion_sub, summary_df)
    }
  每个平台分别汇总并保存。

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
# 解包评估结果字典，将平台信息、过滤参数添加为每个表的列，并拼接成总表。
def unpack_concordance_vmiss_dict(raw_dict):
    all_vmiss, all_smiss, all_confusion = [], [], []

    for key, (vmiss_df, smiss_df, confusion_df) in raw_dict.items():
        key_strs = list(key)
        dp = int(next(s[2:] for s in key_strs if s.startswith('DP')))
        gq = int(next(s[2:] for s in key_strs if s.startswith('GQ')))
        laf = float(next(s[3:] for s in key_strs if s.startswith('LAF')))
        haf = float(next(s[3:] for s in key_strs if s.startswith('HAF')))
        platform = next(s for s in key_strs if s in ['ALL', '15X', '30X'])

        for df in (vmiss_df, smiss_df, confusion_df):
            df["DP"], df["GQ"], df["LAF"], df["HAF"], df["PLATFORM"] = dp, gq, laf, haf, platform

        all_vmiss.append(vmiss_df)
        all_smiss.append(smiss_df)
        all_confusion.append(confusion_df)

    return pd.concat(all_vmiss, ignore_index=True), pd.concat(all_smiss, ignore_index=True), pd.concat(all_confusion, ignore_index=True)

vmiss_df, smiss_df, confusion_df = unpack_concordance_vmiss_dict(raw_dict)

# --- 矩阵转换函数 ---
# 将 confusion_df 转换为标准混淆矩阵（行：真实基因型，列：调用基因型），
# 并补齐 0/1/2/NA 四类组合。
def confusion_to_matrix(confusion_df):
    ordered = ["0", "1", "2", "NA"]
    return (
        confusion_df.pivot_table(index="TRUE_GENOTYPE", columns="CALL_GENOTYPE", values="COUNT", aggfunc="sum", fill_value=0)
        .reindex(index=ordered, columns=ordered, fill_value=0)
    )

# --- 汇总函数 ---
# 输入单个平台的混淆矩阵与缺失率表，计算一系列统计指标。
import pandas as pd

def summarize_genotype_concordance(matrix, vmiss_df, smiss_df, global_dp, global_gq, global_laf, global_haf, platform):
    """
    根据给定的混淆矩阵和缺失率表格，汇总基因型一致性统计信息。

    参数：
        matrix: pd.DataFrame, 行为 TRUE_GENOTYPE，列为 CALL_GENOTYPE，元素为 COUNT
        global_dp, global_gq, global_laf, global_haf: 全局阈值参数
        platform: str, 平台名，如 "15X", "30X", "ALL"
        vmiss_df: pd.DataFrame or None, 包含 "VMISS_FREQ" 列
        smiss_df: pd.DataFrame or None, 包含 "SMISS_FREQ" 列

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

    # --- 错误分类（假阳性） ---
    het_to_homvar = matrix.at["1", "2"] if "1" in matrix.index and "2" in matrix.columns else 0
    homref_to_homvar = matrix.at["0", "2"] if "0" in matrix.index and "2" in matrix.columns else 0
    homref_to_het = matrix.at["0", "1"] if "0" in matrix.index and "1" in matrix.columns else 0

    het_to_homvar_rate = het_to_homvar / total if total > 0 else None
    homref_to_homvar_rate = homref_to_homvar / total if total > 0 else None
    homref_to_het_rate = homref_to_het / total if total > 0 else None
    
    # ---错误分类（假阴性）---
    het_to_homref = matrix.at["1", "0"] if "1" in matrix.index and "0" in matrix.columns else 0
    homvar_to_homref = matrix.at["2", "0"] if "2" in matrix.index and "0" in matrix.columns else 0
    homvar_to_het = matrix.at["2", "1"] if "2" in matrix.index and "1" in matrix.columns else 0
    
    het_to_homref_rate = het_to_homref / total if total > 0 else None
    homvar_to_homref_rate = homvar_to_homref / total if total > 0 else None
    homvar_to_het_rate = homvar_to_het / total if total > 0 else None
    
    # ---正确分类---
    homref_to_homref = matrix.at["0", "0"] if "0" in matrix.index and "0" in matrix.columns else 0
    het_to_het = matrix.at["1", "1"] if "1" in matrix.index and "1" in matrix.columns else 0
    homvar_to_homvar = matrix.at["2", "2"] if "2" in matrix.index and "2" in matrix.columns else 0
    
    homref_to_homref_rate = homref_to_homref / total if total > 0 else None
    het_to_het_rate = het_to_het / total if total > 0 else None
    homvar_to_homvar_rate = homvar_to_homvar / total if total > 0 else None

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
    
    # ---SMISS 统计---
    smiss_mean = smiss_df["SMISS_FREQ"].mean() if smiss_df is not None else None
    smiss_sd = smiss_df["SMISS_FREQ"].std() if smiss_df is not None else None
    smiss_median = smiss_df["SMISS_FREQ"].median() if smiss_df is not None else None

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
        'HET>HOMREF_COUNT': int(het_to_homref),
        'HOMVAR>HOMREF_COUNT': int(homvar_to_homref),
        'HOMVAR>HET_COUNT': int(homvar_to_het),
        'HOMREF>HOMREF_COUNT': int(homref_to_homref),
        'HET>HET_COUNT': int(het_to_het),
        'HOMVAR>HOMVAR_COUNT': int(homvar_to_homvar),
        'HET>HOMVAR_RATE': het_to_homvar_rate,
        'HOMREF>HOMVAR_RATE': homref_to_homvar_rate,
        'HOMREF>HET_RATE': homref_to_het_rate,
        'HET>HOMREF_RATE': het_to_homref_rate,
        'HOMVAR>HOMREF_RATE': homvar_to_homref_rate,
        'HOMVAR>HET_RATE': homvar_to_het_rate,
        'HOMREF>HOMREF_RATE': homref_to_homref_rate,
        'HET>HET_RATE': het_to_het_rate,
        'HOMVAR>HOMVAR_RATE': homvar_to_homvar_rate,
        'VMISS_FREQ_MEAN': vmiss_mean,
        'VMISS_FREQ_SD': vmiss_sd,
        'VMISS_FREQ_MEDIAN': vmiss_median,
        'SMISS_FREQ_MEAN': smiss_mean,
        'SMISS_FREQ_SD': smiss_sd,
        'SMISS_FREQ_MEDIAN': smiss_median
    }])

    return summary_df

# --- 分平台处理 ---
# 按照平台（ALL、15X、30X）分别提取对应子集，计算 summary，并重新封装为输出字典。
platforms = ['ALL', '15X', '30X']
results_dict = {}

for platform in platforms:
    vmiss_sub = vmiss_df[vmiss_df["PLATFORM"] == platform]
    smiss_sub = smiss_df[smiss_df["PLATFORM"] == platform]
    confusion_sub = confusion_df[confusion_df["PLATFORM"] == platform]
    matrix = confusion_to_matrix(confusion_sub)
    summary = summarize_genotype_concordance(matrix, vmiss_sub, smiss_sub, global_dp, global_gq, global_laf, global_haf, platform)

    key = frozenset([f'DP{global_dp}', f'GQ{global_gq}', f'LAF{global_laf}', f'HAF{global_haf}', platform])
    results_dict[key] = (vmiss_sub, smiss_sub, confusion_sub, summary)

# --- 写出结果 ---
# 保存汇总结果到指定命名的 .summary.pkl 文件中。
summary_path = f'{chr}._DP{global_dp}_GQ{global_gq}_LAF{global_laf}_HAF{global_haf}_.summary.pkl'
with open(summary_path, 'wb') as f:
    pickle.dump(results_dict, f)
