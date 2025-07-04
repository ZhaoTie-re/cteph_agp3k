# %%
#usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GWAS显著位点多模型整合与系统注释分析脚本

本脚本面向全基因组关联分析（GWAS）中显著位点（SNP）的多模型（加性、显性、隐性）结果整合与深度注释，旨在为遗传学研究提供高效、自动化的数据处理与注释流程。通过集成外部数据库（如ToMMo、JHRPv4）及主流生物信息学工具（plink2、bcftools、tabix），实现对GWAS显著变异的多维度特征提取与功能解读。

核心功能包括：
1. 自动读取并筛选加性、显性、隐性三种遗传模型下的显著SNP（支持自定义P值阈值），并对结果进行标准化处理。
2. 按位点整合多模型统计信息，实现跨模型的综合比较与展示。
3. 检查显著位点在JHRPv4 VCF数据库中的多等位点（multi-allelic site）状态，输出等位基因组合及相关统计。
4. 基于plink2，分别统计case/control组的基因型分布、等位基因频率（AAF）、Hardy-Weinberg平衡检验（HWE）等群体遗传学指标。
5. 调用bcftools/tabix接口，批量查询ToMMo数据库，获取变异的rsID、基因注释、过滤状态、等位基因频率等核心信息。
6. 自动提取并结构化整理ANN功能注释，便于快速评估变异的潜在生物学影响。
7. 所有注释与统计信息最终合并为结构化表格，支持下游生信分析与可视化。

适用场景：
- GWAS显著位点的多模型整合与系统注释
- 变异功能解读与群体遗传学特征分析
- 生物信息学下游分析的数据准备与结果汇总

依赖环境：
- Python3
- pandas, numpy
- plink2, bcftools, tabix
"""
# %%
import pandas as pd
import numpy as np
import subprocess
import os
import tempfile
import argparse

# ==== 解析命令行参数 ====
parser = argparse.ArgumentParser(description="Process & summary GWAS results.")
parser.add_argument("--bed_prefix", type=str, required=True, help="Prefix for the bed file.")
parser.add_argument("--tommo_vcf_file", type=str, required=True, help="Path to ToMMo VCF file.")
parser.add_argument("--jhrp4_vcf_path", type=str, required=True, help="Path to JHRPv4 VCF files (folder).")
parser.add_argument("--jhrp4_vcf_prefix", type=str, required=True, help="Prefix for JHRPv4 VCF files.")
parser.add_argument("--p_threshold", type=float, default=5e-8, help="P-value threshold for significant SNPs.")
parser.add_argument("--add_path", type=str, default=None, help="Path to additive model results.")
parser.add_argument("--dom_path", type=str, default=None, help="Path to dominant model results.")
parser.add_argument("--rec_path", type=str, default=None, help="Path to recessive model results.")
args = parser.parse_args()

# ==== Path and common variable definitions ====
bed_prefix = args.bed_prefix
tommo_vcf_file = args.tommo_vcf_file
jhrp4_vcf_path = args.jhrp4_vcf_path
jhrp4_vcf_prefix = args.jhrp4_vcf_prefix
p_threshold = args.p_threshold

add_path = args.add_path
dom_path = args.dom_path
rec_path = args.rec_path

# ==== Function modules ====

plink2_path = "/home/b/b37974/plink2" # 请根据实际情况修改plink2的路径

# --- 读取GWAS结果 ---
# 定义所需的列名
desired_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'TEST', 'OR', 'CI95', 'P']

def read_significant_snps(plink2_result_path, p_threshold=5e-8):
    """
    从plink2结果文件中读取P值小于阈值的记录，仅保留指定列，
    并将L95和U95合并为CI95字段。
    """
    with open(plink2_result_path, 'r') as f:
        header = f.readline().strip().split()
        idx_chrom = header.index('#CHROM')
        idx_pos   = header.index('POS')
        idx_id    = header.index('ID')
        idx_ref   = header.index('REF')
        idx_alt   = header.index('ALT')
        idx_test  = header.index('TEST')
        idx_or    = header.index('OR')
        idx_l95   = header.index('L95')
        idx_u95   = header.index('U95')
        idx_p     = header.index('P')
        significant_rows = []
        for line in f:
            fields = line.strip().split()
            try:
                p_value = float(fields[idx_p])
                if p_value < p_threshold:
                    row = [
                        fields[idx_chrom],
                        fields[idx_pos],
                        fields[idx_id],
                        fields[idx_ref],
                        fields[idx_alt],
                        fields[idx_test],
                        fields[idx_or],
                        (fields[idx_l95], fields[idx_u95]),
                        fields[idx_p]
                    ]
                    significant_rows.append(row)
            except ValueError:
                continue
    return pd.DataFrame(significant_rows, columns=desired_columns) if significant_rows else pd.DataFrame(columns=desired_columns)

# --- 合并GWAS模型 ---
def merge_gwas_models(df_add, df_dom, df_rec):
    """
    合并 ADD / DOM / REC 三种模型的结果，按位点合并，字段顺序为 ADD | DOM | REC。
    支持输入为空的DataFrame。
    """
    df_add = df_add.copy()
    df_dom = df_dom.copy()
    df_rec = df_rec.copy()
    if not df_add.empty:
        df_add['MODEL'] = 'ADD'
    if not df_dom.empty:
        df_dom['MODEL'] = 'DOM'
    if not df_rec.empty:
        df_rec['MODEL'] = 'REC'
    frames = [df for df in [df_add, df_dom, df_rec] if not df.empty]
    if not frames:
        return pd.DataFrame(columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'TEST', 'OR', 'CI95', 'P'])
    df_all = pd.concat(frames, ignore_index=True)
    key_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT']
    grouped = df_all.groupby(key_cols)
    merged_rows = []
    for key, group in grouped:
        model_to_row = {row['MODEL']: row for _, row in group.iterrows()}
        row_out = list(key)
        tests, ors, cis, ps = [], [], [], []
        for model in ['ADD', 'DOM', 'REC']:
            if model in model_to_row:
                row = model_to_row[model]
                tests.append(row['TEST'])
                ors.append(row['OR'])
                cis.append(row['CI95'])
                ps.append(row['P'])
        merged_rows.append(row_out + [tests, ors, cis, ps])
    merged_df = pd.DataFrame(
        merged_rows,
        columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'TEST', 'OR', 'CI95', 'P']
    )
    chrom_col = merged_df['#CHROM']
    all_chr_prefixed = chrom_col.str.startswith('chr').all()
    none_chr_prefixed = (~chrom_col.str.startswith('chr')).all()
    strip_chr = all_chr_prefixed or none_chr_prefixed
    def chrom_sort_key(chrom):
        if strip_chr:
            chrom = chrom.lower().replace('chr', '')
        try:
            return (int(chrom),)
        except ValueError:
            return (float('inf'), chrom)
    merged_df = merged_df.sort_values(
        by=['#CHROM', 'POS'],
        key=lambda col: col.map(chrom_sort_key) if col.name == '#CHROM' else pd.to_numeric(col)
    ).reset_index(drop=True)
    return merged_df

# --- 多等位点检测 ---
def check_multiallelic_sites(focus_loci, vcf_path, vcf_prefix):
    """
    检查给定关注位点（focus_loci）在指定VCF文件中是否为多等位点（multi-allelic site）。
    """
    multiallelic_tag = {}
    multiallelic_num = {}
    multiallelic_records_dict = {}
    for locus in focus_loci:
        chrom, pos, ref, alt = locus.split(":")
        vcf_file = f"{vcf_path}/{vcf_prefix}.{chrom}.vcf.gz"
        region = f"{chrom}:{pos}-{pos}"
        try:
            result = subprocess.run(
                ["tabix", vcf_file, region],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            lines = [line for line in result.stdout.strip().split("\n") if line and not line.startswith("#")]
            alts = set()
            for line in lines:
                fields = line.split("\t")
                if fields[0] == chrom and fields[1] == pos and fields[3] == ref:
                    for alt_candidate in fields[4].split(","):
                        alts.add(alt_candidate)
            multiallelic_tag[locus] = len(alts) > 1
            multiallelic_num[locus] = len(alts) + 1  # 包括参考等位基因
            multiallelic_records_dict[locus] = [f"{chrom}:{pos}:{ref}:{alt_candidate}" for alt_candidate in alts]
        except Exception:
            multiallelic_tag[locus] = None
            multiallelic_num[locus] = None
            multiallelic_records_dict[locus] = None
    multiallelic_df = pd.DataFrame({
        "ID": focus_loci,
        "Multi-Allelic TAG": [multiallelic_tag.get(locus) for locus in focus_loci],
        "Allelic NUM": [multiallelic_num.get(locus) for locus in focus_loci],
        "Allelic Records": [multiallelic_records_dict.get(locus) for locus in focus_loci]
    })
    return multiallelic_df

# --- Plink2子集运行及统计 ---
def run_plink_subset(plink2_path, bed_prefix, sample_df, out_suffix, plink_args, focus_loci=None):
    """
    通用的plink2子集运行函数，支持--keep, --extract, 以及自定义参数。
    """
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as keep_file:
        sample_df.to_csv(keep_file.name, sep="\t", index=False, header=False)
    extract_path = None
    if focus_loci:
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as extract_file:
            for snp in focus_loci:
                extract_file.write(f"{snp}\n")
            extract_path = extract_file.name
    with tempfile.NamedTemporaryFile(suffix=out_suffix, delete=False) as out_file:
        out_prefix = out_file.name.rsplit(".", 1)[0]
    cmd = [plink2_path, "--bfile", bed_prefix, "--keep", keep_file.name] + plink_args + ["--out", out_prefix]
    if extract_path:
        cmd += ["--extract", extract_path]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"PLINK2 error: {e.stderr if hasattr(e,'stderr') else str(e)}")
        raise
    result_path = f"{out_prefix}{out_suffix}"
    df = pd.read_csv(result_path, delim_whitespace=True)
    os.remove(keep_file.name)
    if extract_path:
        os.remove(extract_path)
    os.remove(result_path)
    return df

def count_genotypes_by_group(bed_prefix, focus_loci=None):
    """
    使用plink2统计case和control中感兴趣位点的基因型分布。
    """
    fam_path = f"{bed_prefix}.fam"
    fam = pd.read_csv(fam_path, sep="\s+", header=None)
    fam.columns = ['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO']
    case_iids = fam[fam['PHENO'] == 2][['FID', 'IID']]
    ctrl_iids = fam[fam['PHENO'] == 1][['FID', 'IID']]
    plink_args = ["--geno-counts"]
    df_case = run_plink_subset(plink2_path, bed_prefix, case_iids, ".gcount", plink_args, focus_loci)
    df_ctrl = run_plink_subset(plink2_path, bed_prefix, ctrl_iids, ".gcount", plink_args, focus_loci)
    for df in [df_case, df_ctrl]:
        df["TOTAL_CT"] = df[["HOM_REF_CT", "HET_REF_ALT_CTS", "TWO_ALT_GENO_CTS", "MISSING_CT"]].sum(axis=1)
        df["HOM_REF_FREQ"] = np.round(np.where(df["TOTAL_CT"] != 0, df["HOM_REF_CT"] / df["TOTAL_CT"], np.nan), 4)
        df["HET_REF_ALT_FREQ"] = np.round(np.where(df["TOTAL_CT"] != 0, df["HET_REF_ALT_CTS"] / df["TOTAL_CT"], np.nan), 4)
        df["TWO_ALT_GENO_FREQ"] = np.round(np.where(df["TOTAL_CT"] != 0, df["TWO_ALT_GENO_CTS"] / df["TOTAL_CT"], np.nan), 4)
        df["MISSING_FREQ"] = np.round(np.where(df["TOTAL_CT"] != 0, df["MISSING_CT"] / df["TOTAL_CT"], np.nan), 4)
    return df_case, df_ctrl

def compute_hwe_by_group(bed_prefix, focus_loci=None):
    """
    使用plink2计算case和control样本中感兴趣位点的Hardy-Weinberg平衡检验（HWE）。
    """
    fam_path = f"{bed_prefix}.fam"
    fam = pd.read_csv(fam_path, sep="\s+", header=None)
    fam.columns = ['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO']
    case_iids = fam[fam['PHENO'] == 2][['FID', 'IID']]
    ctrl_iids = fam[fam['PHENO'] == 1][['FID', 'IID']]
    plink_args = ["--hardy"]
    df_hwe_case = run_plink_subset(plink2_path, bed_prefix, case_iids, ".hardy", plink_args, focus_loci)
    df_hwe_ctrl = run_plink_subset(plink2_path, bed_prefix, ctrl_iids, ".hardy", plink_args, focus_loci)
    df_hwe_case = df_hwe_case[["ID", "P"]]
    df_hwe_ctrl = df_hwe_ctrl[["ID", "P"]]
    return df_hwe_case, df_hwe_ctrl

def compute_aaf_by_group(bed_prefix, focus_loci=None):
    """
    使用plink2分别计算case和control样本的ALT allele frequency（AAF）。
    """
    fam_path = f"{bed_prefix}.fam"
    fam = pd.read_csv(fam_path, sep="\s+", header=None)
    fam.columns = ['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO']
    case_iids = fam[fam['PHENO'] == 2][['FID', 'IID']]
    ctrl_iids = fam[fam['PHENO'] == 1][['FID', 'IID']]
    plink_args = ["--freq"]
    df_aaf_case = run_plink_subset(plink2_path, bed_prefix, case_iids, ".afreq", plink_args, focus_loci)
    df_aaf_ctrl = run_plink_subset(plink2_path, bed_prefix, ctrl_iids, ".afreq", plink_args, focus_loci)
    return df_aaf_case[["ID", "ALT_FREQS"]], df_aaf_ctrl[["ID", "ALT_FREQS"]]

# --- 变异注释 ---
def query_variant_bcftools(vcf_path, variant_key):
    """
    查询VCF中指定变异的信息。
    输入:
        variant_key: 字符串形式 "CHROM:POS:REF:ALT"
    返回:
        rsid: ID列
        gene_symbols: ANN注释中的gene symbol集合
        filter_status: FILTER列
        allele_freq: INFO字段中的AF值（float or list）
    """
    try:
        chrom, pos, ref, alt = variant_key.strip().split(":")
        region = f"{chrom}:{pos}-{pos}"
    except ValueError:
        raise ValueError("variant_key格式应为 CHROM:POS:REF:ALT")
    cmd = [
        "bcftools", "view", "-r", region, vcf_path
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[bcftools ERROR]\n{e.stderr}")
        return None, None, None, None
    rsid = None
    gene_symbols = set()
    filter_status = None
    allele_freq = None
    for line in result.stdout.splitlines():
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        vcf_ref = fields[3]
        vcf_alts = fields[4].split(",")
        info = fields[7]
        filt = fields[6]
        if vcf_ref == ref and alt in vcf_alts:
            rsid = fields[2]
            filter_status = filt
            af = None
            for entry in info.split(";"):
                if entry.startswith("AF="):
                    af = entry.split("=", 1)[1]
                    allele_freq = af.split(",") if "," in af else [af]
                    break
            for entry in info.split(";"):
                if entry.startswith("ANN="):
                    ann_data = entry[4:]
                    for ann in ann_data.split(","):
                        parts = ann.split("|")
                        if len(parts) > 3 and parts[3]:
                            gene_symbols.add(parts[3])
    return rsid, gene_symbols, filter_status, allele_freq

def batch_query_variants(vcf_path, variant_list):
    """
    批量查询VCF中多个变异的信息，返回DataFrame
    参数:
        vcf_path: VCF文件路径
        variant_list: List[str]，每项为"CHROM:POS:REF:ALT"
    返回:
        DataFrame，列为 ["ID", "rsID", "Gene", "ToMMo Filter", "ToMMo AAF"]
    """
    records = []
    for variant_key in variant_list:
        rsid, gene_symbols, filt, af = query_variant_bcftools(vcf_path, variant_key)
        af_str = ",".join(af) if af else None
        records.append({
            "ID": variant_key,
            "rsID": rsid,
            "Gene": gene_symbols if gene_symbols else pd.NA,
            "ToMMo Filter": filt,
            "ToMMo AAF": af_str
        })
    return pd.DataFrame(records)

import gzip
def extract_ann_records_from_variant(vcf_path, variant_str):
    """
    从 VCF 文件中提取一个变异字符串（如 chr17:43092418:T:C）的 ANN 注释记录。
    返回：List[Dict]，每个 Dict 包含一个 ANN 注释字段。
    """
    def vcf_uses_chr_prefix(path):
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                chrom = line.split("\t")[0]
                return chrom.startswith("chr")
        return False
    def parse_ann_to_dicts(info_field, alt_filter=None):
        ann_key = "ANN="
        ann_string = None
        for entry in info_field.split(";"):
            if entry.startswith(ann_key):
                ann_string = entry[len(ann_key):]
                break
        if not ann_string:
            return []
        ann_fields = [
            "Allele", "Consequence", "Impact", "Gene_symbol", "Gene_ID",
            "Feature_type", "Transcript_ID", "Biotype", "Rank",
            "HGVS_c", "HGVS_p", "cDNA_pos", "CDS_pos", "AA_pos",
            "Distance", "Errors"
        ]
        annotations = []
        for ann_entry in ann_string.split(","):
            values = ann_entry.split("|")
            ann_dict = {ann_fields[i]: values[i] if i < len(values) else "" for i in range(len(ann_fields))}
            if alt_filter and ann_dict["Allele"] != alt_filter:
                continue
            annotations.append(ann_dict)
        return annotations
    chrom, pos, ref, alt = variant_str.strip().split(":")
    pos = int(pos)
    use_chr = vcf_uses_chr_prefix(vcf_path)
    if chrom.startswith("chr") and not use_chr:
        chrom = chrom[3:]
    if not chrom.startswith("chr") and use_chr:
        chrom = "chr" + chrom
    region = f"{chrom}:{pos}-{pos}"
    cmd = ["bcftools", "view", "-r", region, vcf_path]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"bcftools error: {e.stderr}")
        return []
    for line in result.stdout.splitlines():
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        vcf_ref = fields[3]
        vcf_alts = fields[4].split(",")
        info = fields[7]
        if vcf_ref == ref and alt in vcf_alts:
            return parse_ann_to_dicts(info, alt_filter=alt)
    return []

# ==== Main execution ====

df_add = read_significant_snps(add_path, p_threshold=p_threshold)
df_dom = read_significant_snps(dom_path, p_threshold=p_threshold)
df_rec = read_significant_snps(rec_path, p_threshold=p_threshold)
merged_df = merge_gwas_models(df_add, df_dom, df_rec)

# 获取所有关注的位点
focus_loci = merged_df['ID'].tolist()

# 多等位点检测
multiallelic_df = check_multiallelic_sites(
    focus_loci=focus_loci,
    vcf_path=jhrp4_vcf_path,
    vcf_prefix=jhrp4_vcf_prefix
)

# 基因型分布
df_geno_case, df_geno_ctrl = count_genotypes_by_group(bed_prefix, focus_loci=focus_loci)

# HWE
df_hwe_case, df_hwe_ctrl = compute_hwe_by_group(bed_prefix, focus_loci=focus_loci)

# AAF
df_aaf_case, df_aaf_ctrl = compute_aaf_by_group(bed_prefix, focus_loci=focus_loci)

# ToMMo注释
df_tommo_info = batch_query_variants(tommo_vcf_file, focus_loci)

# ANN注释
records = []
for variant_str in focus_loci:
    ann_recs = extract_ann_records_from_variant(tommo_vcf_file, variant_str)
    if not ann_recs:
        summary = pd.NA
        ann_recs = pd.NA
    else:
        summary_lines = [
            f"{r['Gene_symbol']} ({r['Transcript_ID']}): {r['Consequence']} [{r['Biotype']}]"
            for r in ann_recs
        ]
        summary = "\n".join([f"[{i+1}] {s}" for i, s in enumerate(summary_lines)])
    records.append({
        "ID": variant_str,
        "ANN Summary": summary,
        "ANN Records": ann_recs
    })
df_ann = pd.DataFrame(records)

# ==== 列重命名函数 ====
def rename_columns(df, group, col_map):
    return df.rename(columns={k: f"{group} {v}" for k, v in col_map.items()})

# ==== 列名映射表 ====
col_maps = {
    "geno": {
        "HOM_REF_CT": "HomRef Count",
        "HET_REF_ALT_CTS": "Het Count",
        "TWO_ALT_GENO_CTS": "HomAlt Count",
        "MISSING_CT": "Miss Count",
        "TOTAL_CT": "Total Count",
        "HOM_REF_FREQ": "HomRef Freq",
        "HET_REF_ALT_FREQ": "Het Freq",
        "TWO_ALT_GENO_FREQ": "HomAlt Freq",
        "MISSING_FREQ": "Miss Freq"
    },
    "aaf": {
        "ALT_FREQS": "AAF"
    },
    "hwe": {
        "P": "HWE P-value"
    }
}

# ==== 合并所有DataFrame ====
merge_sources = [
    (multiallelic_df, "ID"),
    (rename_columns(df_geno_case, "Case", col_maps["geno"]), "ID"),
    (rename_columns(df_geno_ctrl, "Ctrl", col_maps["geno"]), "ID"),
    (rename_columns(df_aaf_case, "Case", col_maps["aaf"]), "ID"),
    (rename_columns(df_aaf_ctrl, "Ctrl", col_maps["aaf"]), "ID"),
    (rename_columns(df_hwe_case, "Case", col_maps["hwe"]), "ID"),
    (rename_columns(df_hwe_ctrl, "Ctrl", col_maps["hwe"]), "ID"),
    (df_tommo_info, "ID"),
    (df_ann, "ID"),
]
for df_source, key in merge_sources:
    merged_df = merged_df.merge(df_source, on=key, how="left")
merged_df["Case-Ctrl Total Miss Freq"] = np.round(
    (merged_df["Case Miss Count"] + merged_df["Ctrl Miss Count"]) / (merged_df["Case Total Count"] + merged_df["Ctrl Total Count"]), 4
)
# ==== 列顺序调整 ====
col_order = [
    '#CHROM', 'POS', 'ID', 'rsID', 'Gene', 'REF', 'ALT',
    'TEST', 'OR', 'CI95', 'P',
    'Multi-Allelic TAG', 'Allelic NUM', 'Allelic Records',
    'Case HomRef Count', 'Case Het Count', 'Case HomAlt Count', 'Case Miss Count', 'Case Total Count',
    'Ctrl HomRef Count', 'Ctrl Het Count', 'Ctrl HomAlt Count', 'Ctrl Miss Count', 'Ctrl Total Count',
    'Case HomRef Freq', 'Case Het Freq', 'Case HomAlt Freq', 'Case Miss Freq',
    'Ctrl HomRef Freq', 'Ctrl Het Freq', 'Ctrl HomAlt Freq', 'Ctrl Miss Freq',
    'Case-Ctrl Total Miss Freq',
    'Case AAF', 'Ctrl AAF', 'ToMMo AAF', 'ToMMo Filter',
    'Case HWE P-value', 'Ctrl HWE P-value',
    'ANN Summary', 'ANN Records'
]
merged_df = merged_df[col_order]

# ==== 输出结果 ====
merged_df.to_csv("gwas_summary.plink2.csv", index=False)
