
#!/usr/bin/env python3

"""
脚本说明：
该脚本用于合并多个染色体级别的.pkl文件中记录的混淆矩阵（confusion）和缺失率表（vmiss）。
每个.pkl文件包含多个不同过滤参数组合（DP, GQ, LAF, HAF）下的结果。脚本会将其中
参数值与命令行参数匹配的结果合并，并按照预设标签（ALL, 15X, 30X）分别保存。

合并逻辑：
- vmiss：按照“VARIANT_ID”列进行升序排序，排序规则为先按染色体编号（chr1–chr22），再按位点位置。
- confusion：将相同的 (CALL_GENOTYPE, TRUE_GENOTYPE) 键的计数累加。

输入参数：
--dp            最小测序深度（DP）过滤阈值（int）
--gq            最小基因型质量（GQ）过滤阈值（int）
--laf           最小低等位基因频率（LAF）过滤阈值（float）
--haf           最大高等位基因频率（HAF）过滤阈值（float）
--input_chr_pkl 多个.pkl文件路径，每个文件包含一个染色体的结果

输出：
- 合并后的结果保存在一个.pkl文件中，文件名格式为：
  merged._DP{dp}_GQ{gq}_LAF{laf}_HAF{haf}_.pkl
- 其中为每个标签（ALL, 15X, 30X）保存一个合并后的 vmss 表和 confusion 表。

使用示例：
python3 merge_concordance_results.py --dp 10 --gq 20 --laf 0.1 --haf 0.9 --input_chr_pkl chr1.pkl chr2.pkl ...
"""

import argparse
import pickle
import pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

parser = argparse.ArgumentParser(description="Merge confusion and vmiss data from multiple .pkl files")
parser.add_argument('--dp', type=int, required=True)
parser.add_argument('--gq', type=int, required=True)
parser.add_argument('--laf', type=float, required=True)
parser.add_argument('--haf', type=float, required=True)
parser.add_argument('--input_chr_pkl', nargs='+', required=True, help='List of input .pkl files (one per chromosome)')
args = parser.parse_args()

def merge_tag(tag, args, all_vmiss, all_confusion):
    if all_vmiss[tag]:
        merged_vmiss_df = pd.concat(all_vmiss[tag], ignore_index=True)
        # 插入排序逻辑
        if not merged_vmiss_df.empty and "VARIANT_ID" in merged_vmiss_df.columns:
            merged_vmiss_df["__CHROM_INT"] = merged_vmiss_df["VARIANT_ID"].str.split(":").str[0].str.replace("chr", "").astype(int)
            merged_vmiss_df["__POS_INT"] = merged_vmiss_df["VARIANT_ID"].str.split(":").str[1].astype(int)
            merged_vmiss_df.sort_values(["__CHROM_INT", "__POS_INT"], inplace=True)
            merged_vmiss_df.drop(columns=["__CHROM_INT", "__POS_INT"], inplace=True)
    else:
        merged_vmiss_df = pd.DataFrame()

    confusion_df_combined = pd.DataFrame([
        {"CALL_GENOTYPE": k[0], "TRUE_GENOTYPE": k[1], "COUNT": v}
        for k, v in all_confusion[tag].items()
    ])

    output_key = frozenset([
        f'DP{args.dp}',
        f'GQ{args.gq}',
        f'LAF{args.laf}',
        f'HAF{args.haf}',
        tag
    ])
    return output_key, (merged_vmiss_df, confusion_df_combined)



tags = ["ALL", "15X", "30X"]
param_keys = {
    tag: frozenset([
        f'DP{args.dp}',
        f'GQ{args.gq}',
        f'LAF{args.laf}',
        f'HAF{args.haf}',
        tag
    ]) for tag in tags
}

# === 初始化容器 ===
all_vmiss = {tag: [] for tag in tags}
all_confusion = {tag: defaultdict(int) for tag in tags}

# === 遍历pkl文件 ===
for pkl_file in args.input_chr_pkl:
    with open(pkl_file, "rb") as f:
        data = pickle.load(f)

    for key, (vmiss_df, confusion_df) in data.items():
        for tag, param_key in param_keys.items():
            if key == param_key:
                if vmiss_df is not None:
                    all_vmiss[tag].append(vmiss_df)
                if confusion_df is not None:
                    for _, row in confusion_df.iterrows():
                        tup = (row["CALL_GENOTYPE"], row["TRUE_GENOTYPE"])
                        all_confusion[tag][tup] += row["COUNT"]

with ProcessPoolExecutor() as executor:
    futures = [executor.submit(merge_tag, tag, args, all_vmiss, all_confusion) for tag in tags]
    merged_results_dict = dict(f.result() for f in futures)

output_path = f"merged._DP{args.dp}_GQ{args.gq}_LAF{args.laf}_HAF{args.haf}_.pkl"
with open(output_path, "wb") as out_f:
    pickle.dump(merged_results_dict, out_f)
print(f"[✓] Merged results for ALL/15X/30X saved → {output_path}")