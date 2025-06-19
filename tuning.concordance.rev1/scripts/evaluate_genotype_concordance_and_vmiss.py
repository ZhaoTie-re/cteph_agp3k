# %%
"""
脚本名称：evaluate_genotype_concordance_and_vmiss.py
功能概述：评估WGS与array基因型一致性，并计算变体层面的缺失率（VMISS）

描述：
    读取WGS和array数据的基因型、测序深度（DP）、基因型质量（GQ）、等位基因比例（Allele Fraction, AF）等指标，
    通过设定阈值过滤低质量基因型，计算基因型一致性矩阵，
    同时计算每个位点的VMISS（变体缺失率）。

输入文件：
    - {chr}.wgs.GT.tsv.gz: WGS基因型
    - {chr}.wgs.DP.tsv.gz: WGS测序深度
    - {chr}.wgs.GQ.tsv.gz: WGS基因型质量
    - {chr}.wgs.AF.tsv.gz: WGS等位基因比例（AF）
    - {chr}.array.GT.tsv.gz: array基因型真值
    - {chr}.shared_samples.txt: 共享样本列表
    - {chr}.shared_variant_ids.txt: 共享变体列表
    - wgs_array_dp.csv: 样本测序深度信息（区分15x和30x）

输出内容：
    - results_dict：字典，键为参数组合，值为(vmiss_df, confusion_df)
    - 保存文件名格式：
      DP{dp}_GQ{gq}_LAF{low_af}_HAF{high_af}.pkl

注意事项：
    - 代码结构支持多进程加速
    - 可扩展至更多质控指标（如PL, AD、GQmiss、sample-level missing等）
"""
# %%
import multiprocessing as mp
from functools import partial
import pandas as pd
import argparse

# === 参数设置 ===
parser = argparse.ArgumentParser(description="Evaluate genotype concordance and VMISS for WGS vs array data.")
parser.add_argument("--chr", type=str, required=True, help="Chromosome to process (e.g., chr22)")
parser.add_argument("--global_dp", type=int, default=8, help="Global DP threshold for filtering")
parser.add_argument("--global_gq", type=int, default=20, help="Global GQ threshold for filtering")
parser.add_argument("--global_low_af", type=float, default=0.25, help="Low allele fraction threshold for filtering")
parser.add_argument("--global_high_af", type=float, default=0.75, help="High allele fraction threshold for filtering")
parser.add_argument("--batch_size", type=int, default=1000, help="Batch size for processing variants")
parser.add_argument("--num_threads", type=int, default=4, help="Number of threads for multiprocessing")
parser.add_argument("--call_gt", type=str, required=True, help="Path to WGS genotype file (GT)")
parser.add_argument("--call_dp", type=str, required=True, help="Path to WGS depth file (DP)")
parser.add_argument("--call_gq", type=str, required=True, help="Path to WGS genotype quality file (GQ)")
parser.add_argument("--call_af", type=str, required=True, help="Path to WGS allele fraction file (AF)")
parser.add_argument("--true_gt", type=str, required=True, help="Path to true genotype file (GT) for array data")
parser.add_argument("--samples_txt", type=str, required=True, help="Path to shared samples text file")
parser.add_argument("--variants_txt", type=str, required=True, help="Path to shared variant IDs text file")
parser.add_argument("--sample_info_csv", type=str, required=True, help="Path to sample info CSV file with DP information")
args = parser.parse_args()

# === 解析参数 ===
chr = args.chr
global_dp = args.global_dp  # Global DP threshold
global_gq = args.global_gq  # Global GQ threshold
global_low_af = args.global_low_af  # Allele Fraction (low confidence)
global_high_af = args.global_high_af  # Allele Fraction (high confidence)
batch_size = args.batch_size  # Batch size for processing variants
num_threads = args.num_threads  # Number of threads for multiprocessing

call_gt = args.call_gt  # Path to WGS genotype file (GT)
call_dp = args.call_dp  # Path to WGS depth file (DP)
call_gq = args.call_gq  # Path to WGS genotype quality file (GQ)
call_af = args.call_af  # Path to WGS allele fraction file (AF)
true_gt = args.true_gt  # Path to true genotype file (GT) for array data

samples_txt = args.samples_txt  # Path to shared samples text file
variants_txt = args.variants_txt  # Path to shared variant IDs text file
sample_info_csv = args.sample_info_csv  # Path to sample info CSV file with DP information

# %%
# === 读取样本和变体 ID ===
with open(samples_txt) as f:
    sample_ids = [line.strip() for line in f if line.strip()]
with open(variants_txt) as f:
    variant_ids = [line.strip() for line in f if line.strip()]

# === 读取样本测序深度信息 ===
sample_info = pd.read_csv(sample_info_csv, usecols=["ID", "Target DP (JHRPv4)"])

# === 提取测序深度为 15x 和 30x 的样本 ID 并与 sample_ids 交集 ===
sample_ids_set = set(sample_ids)

sample_ids_15x = sample_info.loc[sample_info["Target DP (JHRPv4)"] == "15x", "ID"]
sample_ids_30x = sample_info.loc[sample_info["Target DP (JHRPv4)"] == "30x", "ID"]

sample_ids_15x = list(sample_ids_set & set(sample_ids_15x))
sample_ids_30x = list(sample_ids_set & set(sample_ids_30x))

# %%
# === 顶层函数：可被多进程调用 ===
def process_variant_block(cols, sample_ids_subset):
    try:
        usecols = ["sample_id"] + cols

        gt = pd.read_csv(call_gt, sep="\t", compression="gzip", usecols=usecols, index_col=0).loc[sample_ids_subset].astype("Int64")
        dp = pd.read_csv(call_dp, sep="\t", compression="gzip", usecols=usecols, index_col=0).loc[sample_ids_subset]
        gq = pd.read_csv(call_gq, sep="\t", compression="gzip", usecols=usecols, index_col=0).loc[sample_ids_subset]
        af = pd.read_csv(call_af, sep="\t", compression="gzip", usecols=usecols, index_col=0).loc[sample_ids_subset]
        tru = pd.read_csv(true_gt, sep="\t", compression="gzip", usecols=usecols, index_col=0).loc[sample_ids_subset].astype("Int64")

        assert list(gt.columns) == list(tru.columns), "Variant ID columns not aligned"
        assert list(gt.index) == list(tru.index), "Sample ID rows not aligned"

        qc_gt = gt.copy()
        mask = (dp < global_dp) | (gq < global_gq)
        mask |= ((gt == 1) & (af < global_low_af))
        mask |= ((gt == 1) & (af > global_high_af))
        qc_gt[mask] = pd.NA

        merged = pd.DataFrame({
            "CALL_GENOTYPE": qc_gt.values.ravel(),
            "TRUE_GENOTYPE": tru.values.ravel()
        })
        counts = merged.groupby(["CALL_GENOTYPE", "TRUE_GENOTYPE"], dropna=False).size().reset_index(name="COUNT")

        vmiss_count = qc_gt.isna().sum()
        total_sample = qc_gt.shape[0]
        vmiss_freq = (vmiss_count / total_sample).round(4)
        vmiss = pd.DataFrame({
            "VARIANT_ID": vmiss_count.index,
            "VMISS_COUNT": vmiss_count.values,
            "TOTAL_SAMPLE": total_sample,
            "VMISS_FREQ": vmiss_freq.values
        })

        return counts, vmiss

    except Exception as e:
        print(f"[ERROR] Failed on columns: {cols[:5]} → {e}")
        empty_conf = pd.DataFrame(columns=["CALL_GENOTYPE", "TRUE_GENOTYPE", "COUNT"])
        empty_vmiss = pd.DataFrame(columns=["VARIANT_ID", "VMISS_COUNT", "TOTAL_SAMPLE", "VMISS_FREQ"])
        return empty_conf, empty_vmiss


# === 执行函数 ===
results_dict = {}

def run_and_save(name, sample_ids):
    block_func = partial(process_variant_block, sample_ids_subset=sample_ids)
    with mp.Pool(processes=num_threads) as pool:
        results = pool.map(block_func, blocks)

    # === 合并结果 ===
    confusion_tables, vmiss_tables = zip(*results)
    confusion_df = pd.concat(confusion_tables, ignore_index=True)
    confusion_df = confusion_df.groupby(["CALL_GENOTYPE", "TRUE_GENOTYPE"], dropna=False).sum().reset_index()

    def format_value(x):
        return "NA" if pd.isna(x) else str(int(x))
    confusion_df["CALL_GENOTYPE"] = confusion_df["CALL_GENOTYPE"].map(format_value)
    confusion_df["TRUE_GENOTYPE"] = confusion_df["TRUE_GENOTYPE"].map(format_value)

    vmiss_df = pd.concat(vmiss_tables, ignore_index=True)

    # === 构建键 ===
    param_key = frozenset([
        f'DP{global_dp}',
        f'GQ{global_gq}',
        f'LAF{global_low_af}',
        f'HAF{global_high_af}',
        name.upper()  # ALL / 15X / 30X
    ])

    # === 存入全局字典 ===
    results_dict[param_key] = (vmiss_df, confusion_df)
    print(f"[✓] Results stored in results_dict under key → {param_key}")


# === 分批变体 ===
blocks = [variant_ids[i:i + batch_size] for i in range(0, len(variant_ids), batch_size)]

# === 分别运行 all, 15x, 30x ===
run_and_save("all", sample_ids)
run_and_save("15x", sample_ids_15x)
run_and_save("30x", sample_ids_30x)

# %%
import pickle

with open(f'{chr}._DP{global_dp}_GQ{global_gq}_LAF{global_low_af}_HAF{global_high_af}_.pkl', 'wb') as f:
    pickle.dump(results_dict, f)
