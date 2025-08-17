"""
sample_qc_calculator.py

功能简介：
    本模块为人类全基因组关联分析（GWAS）中的样本质量控制（Sample QC）提供一系列标准化计算工具，
    包括性别/表型更新、杂合性（Het F）、样本缺失率（SMISS）、亲缘关系（PI_HAT）、Robust Z 分数等，
    特别适用于使用 PLINK 生成的中间文件（.bed/.bim/.fam）格式。

依赖环境：
    - pandas, numpy
    - scipy, statsmodels
    - PLINK1.9, PLINK2（命令行调用）

模块接口：
    - update_fam_and_pheno: 基于 info 表更新 sex 字段，生成 phenotype 文件，并创建新的 .bed 数据。
    - run_sample_qc_summary: 基于 BED 数据自动提取 QC 指标，生成 sample_qc_summary.csv 和 pi_hat.csv。
    - compute_robust_z_by_group: 在样本分组内对指定值列计算 robust Z 分数，输出显著性与 FDR。

作者: ZHAO TIE
"""
import pandas as pd
import numpy as np
import subprocess
import os
import tempfile
from scipy.stats import median_abs_deviation, norm
from statsmodels.stats.multitest import fdrcorrection
from typing import Literal

def update_fam_and_pheno(
    info_file: str,
    bed_prefix: str,
    case_prefix: str,
    out_prefix: str,
    plink_path: str = "/home/b/b37974/plink2"
) -> str:
    """
    使用样本信息表更新 FAM 文件中的性别字段，并生成用于 PLINK2 分析的表型文件。

    Parameters
    ----------
    info_file : str
        Excel 文件路径，需包含 'ID' 和 'Sex' 两列。
    bed_prefix : str
        原始 PLINK 文件前缀（不含扩展名 .bed/.bim/.fam）。
    case_prefix : str
        用于标记 case 样本的 ID 前缀，例如 "PHOM"。
    out_prefix : str
        输出文件的命名前缀，生成 {out_prefix}.bed/.bim/.fam。
    plink_path : str, default="/home/b/b37974/plink2"
        plink2 执行路径。

    Returns
    -------
    str
        新生成的 .bed 数据文件前缀路径。

    Raises
    ------
    RuntimeError
        若 Plink 命令执行失败。
    """
    
    # Step 1: 读取包含样本信息（ID 和 性别）的 Excel 文件
    info_df = pd.read_excel(info_file, header=0)

    # Step 2: 构建 sex 文件 DataFrame
    def _generate_sex_file(df: pd.DataFrame) -> pd.DataFrame:
        """构建用于 --update-sex 的 DataFrame。"""
        return pd.DataFrame({
            '#FID': df['ID'],
            'IID': df['ID'],
            'SEX': df['Sex']
        })

    # Step 3: 构建 phenotype 文件 DataFrame
    def _generate_pheno_file(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
        """
        如果样本 ID 以 prefix 开头，标记为 case (2)，否则标记为 control (1)。
        PLINK2 默认读取 PHENO1 列作为表型分类。
        """
        return pd.DataFrame({
            '#FID': df['ID'],
            'IID': df['ID'],
            'PHENO1': [2 if id.startswith(prefix) else 1 for id in df['ID']]
        })

    sex_df = _generate_sex_file(info_df)
    pheno_df = _generate_pheno_file(info_df, case_prefix)

    # Step 4: 保存中间文本文件
    sex_file = f"{out_prefix}.sex.txt"
    pheno_file = f"{out_prefix}.pheno.txt"
    sex_df.to_csv(sex_file, sep='\t', index=False)
    pheno_df.to_csv(pheno_file, sep='\t', index=False)

    # Step 5: 构造并运行 plink2 命令
    plink_command = [
        plink_path,
        "--bfile", bed_prefix,
        "--pheno", pheno_file,
        "--update-sex", sex_file,
        "--make-bed",
        "--out", out_prefix
    ]

    try:
        subprocess.run(plink_command, check=True)
        print(f"[成功] Plink2 执行成功，输出文件前缀为: {out_prefix}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"[错误] Plink2 执行失败：{e}")

    # Step 6: 返回最终生成文件的前缀路径
    return out_prefix

def run_sample_qc_summary(
    info_file: str, 
    bed_prefix: str, 
    high_ld_file: str,
    threads: int = 8
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    基于给定的 .fam 文件样本列表和 info 表格，使用 PLINK 工具生成样本质量控制指标。
    包括：平台信息清洗、样本缺失率（SMISS）、杂合性F统计（HET_F）、亲缘关系（PI_HAT）。

    Parameters
    ----------
    info_file : str
        Excel 信息表路径，包含样本相关信息。
    bed_prefix : str
        PLINK 文件前缀路径（不含扩展名）。
    high_ld_file : str
        高 LD 区域 range 文件路径。
    threads : int, optional
        运行线程数，默认 8。

    Returns
    -------
    tuple of pd.DataFrame
        sample_qc_summary : pd.DataFrame
            合并后的样本信息、缺失率、杂合统计。
        pihat_df : pd.DataFrame
            PI_HAT 亲缘关系数据表。

    Raises
    ------
    RuntimeError
        PLINK 命令执行失败时抛出。
    FileNotFoundError
        预期的 PLINK 输出文件未找到时抛出。
    """

    plink1_path = "/home/b/b37974/plink"
    plink2_path = "/home/b/b37974/plink2"

    def clean_info_file(info_file: str, bed_prefix: str) -> pd.DataFrame:
        """
        清理 info 表格，统一测序深度字段，并仅保留 FAM 文件中样本。

        Parameters
        ----------
        info_file : str
            Excel 信息表路径。
        bed_prefix : str
            PLINK 文件前缀。

        Returns
        -------
        pd.DataFrame
            清洗并格式化后的样本信息表。
        """
        df = pd.read_excel(info_file)

        def process_row(row):
            target_dp = row['Target DP (JHRPv4)']
            wgs = row['WGS']
            dp_val = row['DP (JHRPv4)']

            if target_dp == '15x | 30x':
                row['Target DP (JHRPv4)'] = '15x'
                row['WGS'] = wgs.split(' | ')[0]
                row['DP (JHRPv4)'] = dp_val.split(' | ')[0]
            elif target_dp == '15x | 15x':
                row['Target DP (JHRPv4)'] = '15x'
                row['WGS'] = 'HiSeqX 15x'
                row['DP (JHRPv4)'] = dp_val.split(' | ')[1]
            return row

        df = df.apply(process_row, axis=1)
        df['DP (JHRPv4)'] = df['DP (JHRPv4)'].astype(float)

        fam_file = f"{bed_prefix}.fam"
        fam_ids = pd.read_csv(fam_file, delim_whitespace=True, header=None, usecols=[1])[1].tolist()
        df = df[df['ID'].isin(fam_ids)].reset_index(drop=True)

        df = df[["ID", "WGS", "Target DP (JHRPv4)", "DP (JHRPv4)"]].rename(columns={
            "ID": "IID",
            "WGS": "PLATFORM",
            "Target DP (JHRPv4)": "TARGET_DP",
            "DP (JHRPv4)": "MEAN_DP"
        })
        df.insert(0, "#FID", df["IID"])
        return df

    def get_sample_missingness(plink_path: str, bed_prefix: str) -> pd.DataFrame:
        """
        计算每个样本的缺失率（Sample missingness），使用 plink2 的 --missing 选项。

        Parameters
        ----------
        plink_path : str
            plink2 执行路径。
        bed_prefix : str
            PLINK 文件前缀。

        Returns
        -------
        pd.DataFrame
            包含样本缺失率的 DataFrame，列包括 '#FID', 'IID', 'SMISS'。

        Raises
        ------
        RuntimeError
            plink2 执行失败。
        FileNotFoundError
            预期 .smiss 文件未找到。
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            out_prefix = os.path.join(tmpdir, "temp_output")
            cmd = [plink_path, "--bfile", bed_prefix, "--missing", "--maf", "0.05", "--out", out_prefix]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                raise RuntimeError(f"PLINK2 failed:\n{result.stderr}")
            smiss_file = f"{out_prefix}.smiss"
            if not os.path.exists(smiss_file):
                raise FileNotFoundError("Expected .smiss file not found.")
            smiss_df = pd.read_csv(smiss_file, delim_whitespace=True)
        return smiss_df[['#FID', 'IID', 'F_MISS']].rename(columns={"F_MISS": "SMISS"})

    def compute_het_and_pi_hat_clean(
        plink2_path: str,
        plink1_path: str,
        input_prefix: str,
        high_ld_file: str,
        threads: int = 8,
        mode: str = "qc"
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        计算每个样本的杂合性统计（Het F）和亲缘关系（PI_HAT）。
        - 先对变异进行 LD-prune，并过滤高LD区域
        - 用 plink2 计算 Het F 值，用 plink1.9 计算 PI_HAT

        Parameters
        ----------
        plink2_path : str
            plink2 执行路径。
        plink1_path : str
            plink1.9 执行路径。
        input_prefix : str
            输入 PLINK 文件前缀。
        high_ld_file : str
            高 LD 区域 range 文件路径。
        threads : int, optional
            线程数，默认 8。
        mode : str, optional
            运行模式，默认 "qc"。

        Returns
        -------
        tuple of pd.DataFrame
            het_df: 包含 '#FID', 'IID', 'HET_F' 列。
            pihat_df: PI_HAT 亲缘关系数据表。
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            temp_ld_prefix = os.path.join(tmpdir, "ld_pruned")
            temp_pruned_prefix = os.path.join(tmpdir, "pruned_data")
            het_out = os.path.join(tmpdir, "het")
            genome_out = os.path.join(tmpdir, "genome")

            cmd_prune = [
                plink2_path, "--bfile", input_prefix,
                "--snps-only", "just-acgt",
                "--exclude", "range", high_ld_file,
                "--indep-pairwise", "50", "5", "0.2",
                "--threads", str(threads),
                "--out", temp_ld_prefix
            ]
            if mode == "qc":
                cmd_prune.insert(cmd_prune.index("--exclude"), "--maf")
                cmd_prune.insert(cmd_prune.index("--maf") + 1, "0.05")

            subprocess.run(cmd_prune, check=True)

            subprocess.run([
                plink2_path, "--bfile", input_prefix,
                "--extract", temp_ld_prefix + ".prune.in",
                "--make-bed", "--out", temp_pruned_prefix,
                "--threads", str(threads)
            ], check=True)

            subprocess.run([
                plink2_path, "--bfile", temp_pruned_prefix,
                "--het", "--out", het_out, "--threads", str(threads)
            ], check=True)

            subprocess.run([
                plink1_path, "--bfile", temp_pruned_prefix,
                "--genome", "--out", genome_out
            ], check=True)

            het_df = pd.read_csv(het_out + ".het", delim_whitespace=True)
            pihat_df = pd.read_csv(genome_out + ".genome", delim_whitespace=True)

        return het_df[['#FID', 'IID', 'F']].rename(columns={"F": "HET_F"}), pihat_df

    # === Step 1: 清洗 info 表，标准化字段命名 ===
    info_cleaned = clean_info_file(info_file, bed_prefix)

    # === Step 2: 计算样本缺失率（SMISS） ===
    smiss_df = get_sample_missingness(plink2_path, bed_prefix)

    # === Step 3: 计算杂合性统计（HET_F）和亲缘关系（PI_HAT） ===
    het_df, pihat_df = compute_het_and_pi_hat_clean(plink2_path, plink1_path, bed_prefix, high_ld_file, threads=threads, mode="qc")

    # === Step 4: 合并 info、smiss、het 结果为一个样本质量控制表 ===
    sample_qc_summary = (
        info_cleaned
        .merge(smiss_df, on=["#FID", "IID"], how="inner")
        .merge(het_df, on=["#FID", "IID"], how="inner")
    )
    
    # 保存结果文件（可选），包括 sample_qc_summary 和 PI_HAT 表
    sample_qc_summary.to_csv(f"{bed_prefix}.sample_qc_summary.csv", index=False)
    pihat_df.to_csv(f"{bed_prefix}.pi_hat.csv", index=False)

    return sample_qc_summary, pihat_df

def compute_robust_z_by_group(
    df: pd.DataFrame,
    group_col: str,
    value_col: str,
    z_col_name: str = "ROBUST_Z",
    p_col_name: str = "P_ROBUST_Z",
    fdr_col_name: str = "FDR_ROBUST_Z",
    tail: Literal["two-side", "left", "right"] = "left"
) -> pd.DataFrame:
    """
    对指定分组内的数值列计算 Robust Z Score、p 值和 FDR。

    Parameters
    ----------
    df : pd.DataFrame
        输入数据框。
    group_col : str
        分组依据的列名。
    value_col : str
        要计算 z 分数的数值列。
    z_col_name : str, optional
        输出的 z 分数列名，默认 "ROBUST_Z"。
    p_col_name : str, optional
        输出的 p 值列名，默认 "P_ROBUST_Z"。
    fdr_col_name : str, optional
        输出的 FDR 列名，默认 "FDR_ROBUST_Z"。
    tail : {'two-side', 'left', 'right'}, optional
        指定 p 值计算的尾部类型：
        - 'two-side': 双尾检验
        - 'left': 检测偏低
        - 'right': 检测偏高
        默认 'left'。

    Returns
    -------
    pd.DataFrame
        输入数据框添加如下列：
            - z_col_name: Robust Z 分数
            - p_col_name: 单/双尾 p 值
            - fdr_col_name: 组内 FDR 多重比较校正结果

    Raises
    ------
    ValueError
        如果 tail 参数不在允许取值范围内，则抛出异常。
    """
    def calc_robust_z_with_p(group):
        # 计算当前分组中位数和中位绝对偏差（MAD），用于稳健 Z 分数
        median = group[value_col].median()
        mad = median_abs_deviation(group[value_col], scale='normal')

        # 如果 MAD 为 0，说明该组数据分布极为集中，中位数相同，避免除以零错误，赋值为 0
        if mad == 0:
            group[z_col_name] = 0
        else:
            group[z_col_name] = (group[value_col] - median) / mad

        z_scores = group[z_col_name]

        # 根据 tail 参数类型计算单尾或双尾 p 值
        if tail == "two-side":
            group[p_col_name] = 2 * (1 - norm.cdf(np.abs(z_scores)))
        elif tail == "left":
            group[p_col_name] = norm.cdf(z_scores)
        elif tail == "right":
            group[p_col_name] = 1 - norm.cdf(z_scores)
        else:
            raise ValueError("tail 参数必须是 'two-side', 'left', 或 'right'")

        # 组内执行 FDR 多重检验校正
        # 如果组内样本数少于2，则不进行校正以避免报错，直接赋值为 NaN
        if len(group) >= 2:
            _, qvals = fdrcorrection(group[p_col_name])
            group[fdr_col_name] = qvals
        else:
            group[fdr_col_name] = np.nan

        return group

    # 返回带有每个样本 z 分数、p 值与 FDR 的完整表格
    return df.groupby(group_col, group_keys=False).apply(calc_robust_z_with_p)