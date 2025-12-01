#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
基因型变体质量控制工具
用于计算变体的VMISS（总体和分组）以及MAF
"""

import pandas as pd
import numpy as np
import subprocess
import logging
import os
import tempfile
import multiprocessing as mp
from pathlib import Path
from typing import Optional, Literal

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def calculate_variant_metrics(
    bed_prefix: str,
    info_path: str = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/cteph_agp3k_jhrpv4.rev1.xlsx",
    output_path: Optional[str] = None,
    output_dir: Optional[str] = None,
    dp_column: str = "Target DP (JHRPv4)",
    id_column: str = "ID",
    maf_group: Literal["ctrl", "case", "all"] = "ctrl",
    pheno_column: str = "OUTCOME2",
    plink2_path: str = "/home/b/b37974/plink2",
    n_threads: int = 4
) -> str:
    """
    计算基因型文件中每个变体的VMISS和MAF
    
    参数:
        bed_prefix: plink2基因型文件的前缀路径（不含.bed/.bim/.fam后缀）
        info_path: 样本信息表路径，默认为cteph_agp3k项目的info文件
        output_path: 输出文件路径，如果指定则output_dir参数无效
        output_dir: 输出文件夹路径，默认为当前工作目录下的vmiss_results文件夹
        dp_column: 测序深度列名，默认为"Target DP (JHRPv4)"
        id_column: 样本ID列名，默认为"ID"
        maf_group: 计算MAF的样本组，可选"ctrl"/"case"/"all"，默认为"ctrl"
        pheno_column: 表型列名，默认为"OUTCOME2"（AGP3K=对照组，CTEPH=病例组）
        plink2_path: plink2可执行文件路径
        n_threads: 并行线程数
    
    返回:
        str: 输出的tab分割文件路径，包含以下列：
            - VARIANT_ID: 变体ID
            - VMISS: 所有样本的变体缺失率
            - VMISS_15X: 15X深度样本的变体缺失率
            - VMISS_30X: 30X深度样本的变体缺失率
            - MAF: 次要等位基因频率（Minor Allele Frequency），取min(ALT_FREQ, 1-ALT_FREQ)
    
    内存优化策略:
        1. 完全流式处理：从不将完整的plink2输出加载到内存
        2. 分批构建索引：使用chunksize=100000逐批读取构建字典
        3. 直接写入输出：逐行处理并直接写入最终文件
        4. 显式内存清理：使用gc.collect()主动释放内存
        
        这种设计可以处理数百万变体，内存占用主要取决于变体ID字典大小，
        而不是完整数据集大小。适合在内存受限的环境中运行。
    """
    
    import gc
    gc.collect()  # 开始前清理内存
    
    logger.info("="*60)
    logger.info("开始计算变体质量指标")
    logger.info(f"基因型文件前缀: {bed_prefix}")
    logger.info(f"样本信息文件: {info_path}")
    logger.info(f"MAF计算组别: {maf_group}")
    logger.info(f"线程数: {n_threads}")
    logger.info("="*60)
    
    # 检查文件是否存在
    _check_input_files(bed_prefix, info_path, plink2_path)
    
    # 读取样本信息
    logger.info("正在读取样本信息文件...")
    sample_info = _read_sample_info(info_path, id_column, dp_column, pheno_column)
    logger.info(f"成功读取 {len(sample_info)} 个样本的信息")
    
    # 读取fam文件获取样本列表
    logger.info("正在读取fam文件...")
    fam_samples = _read_fam_file(bed_prefix)
    logger.info(f"基因型文件包含 {len(fam_samples)} 个样本")
    
    # 合并样本信息和fam文件
    logger.info("正在匹配样本信息...")
    merged_samples = _merge_sample_info(fam_samples, sample_info, id_column, dp_column, pheno_column)
    logger.info(f"成功匹配 {len(merged_samples)} 个样本")
    
    # 按测序深度分组
    dp_groups = _group_samples_by_dp(merged_samples, dp_column)
    logger.info(f"测序深度分组: {dict(dp_groups['count'])}")
    
    # 设置输出路径
    if output_path is None:
        # 如果没有指定output_path，使用output_dir
        if output_dir is None:
            output_dir = os.path.join(os.getcwd(), "vmiss_results")
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        
        # 生成输出文件名（基于bed文件名）
        bed_name = Path(bed_prefix).name
        output_path = os.path.join(output_dir, f"{bed_name}_variant_metrics.tsv")
    else:
        # 如果指定了output_path，确保其目录存在
        output_dir = os.path.dirname(output_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
    
    logger.info(f"输出文件路径: {output_path}")
    
    # 创建临时目录
    with tempfile.TemporaryDirectory(prefix="vmiss_calc_") as tmpdir:
        logger.info(f"创建临时目录: {tmpdir}")
        
        # 生成样本列表文件
        sample_lists = _create_sample_list_files(
            tmpdir, merged_samples, dp_groups, maf_group, pheno_column
        )
        
        # 计算总体VMISS
        logger.info("正在计算总体VMISS...")
        vmiss_all_file = _calculate_vmiss_plink(
            bed_prefix, plink2_path, tmpdir, "all", 
            sample_list=None, n_threads=n_threads
        )
        logger.info(f"成功计算总体VMISS")
        
        # 计算15X组的VMISS
        if '15X' in sample_lists:
            logger.info("正在计算15X组VMISS...")
            vmiss_15x_file = _calculate_vmiss_plink(
                bed_prefix, plink2_path, tmpdir, "15x",
                sample_list=sample_lists['15X'], n_threads=n_threads
            )
            logger.info(f"成功计算15X组VMISS")
        else:
            logger.warning("未找到15X组样本")
            vmiss_15x_file = None
        
        # 计算30X组的VMISS
        if '30X' in sample_lists:
            logger.info("正在计算30X组VMISS...")
            vmiss_30x_file = _calculate_vmiss_plink(
                bed_prefix, plink2_path, tmpdir, "30x",
                sample_list=sample_lists['30X'], n_threads=n_threads
            )
            logger.info(f"成功计算30X组VMISS")
        else:
            logger.warning("未找到30X组样本")
            vmiss_30x_file = None
        
        # 计算MAF
        logger.info(f"正在计算{maf_group}组的MAF...")
        maf_file = _calculate_maf_plink(
            bed_prefix, plink2_path, tmpdir, maf_group,
            sample_list=sample_lists.get(maf_group), n_threads=n_threads
        )
        logger.info(f"成功计算MAF")
        
        # 合并结果并直接写入文件（流式处理，不加载到内存）
        logger.info("正在合并结果并写入文件（流式处理）...")
        n_variants = _merge_and_write_results_streaming(
            vmiss_all_file, vmiss_15x_file, vmiss_30x_file, maf_file, output_path
        )
    logger.info(f"成功写入 {n_variants} 个变体到文件")
    
    logger.info("="*60)
    logger.info(f"计算完成！结果已保存到: {output_path}")
    logger.info("="*60)
    
    return output_path


def _check_input_files(bed_prefix: str, info_path: str, plink2_path: str):
    """检查输入文件是否存在"""
    # 检查plink文件
    for ext in ['.bed', '.bim', '.fam']:
        file_path = bed_prefix + ext
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"找不到文件: {file_path}")
    
    # 检查info文件
    if not os.path.exists(info_path):
        raise FileNotFoundError(f"找不到样本信息文件: {info_path}")
    
    # 检查plink2
    if not os.path.exists(plink2_path):
        raise FileNotFoundError(f"找不到plink2: {plink2_path}")


def _read_sample_info(info_path: str, id_column: str, dp_column: str, pheno_column: str) -> pd.DataFrame:
    """读取样本信息文件"""
    # 根据文件扩展名选择读取方式
    if info_path.endswith('.xlsx'):
        df = pd.read_excel(info_path)
    elif info_path.endswith('.csv'):
        df = pd.read_csv(info_path)
    elif info_path.endswith('.tsv') or info_path.endswith('.txt'):
        df = pd.read_csv(info_path, sep='\t')
    else:
        raise ValueError(f"不支持的文件格式: {info_path}")
    
    # 检查必需的列是否存在
    required_cols = [id_column, dp_column]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"样本信息文件缺少必需的列: {missing_cols}\n可用的列: {df.columns.tolist()}")
    
    # 如果需要计算case/ctrl的MAF，检查表型列
    if pheno_column not in df.columns:
        logger.warning(f"未找到表型列'{pheno_column}'，将无法按病例对照分组")
    
    return df


def _read_fam_file(bed_prefix: str) -> pd.DataFrame:
    """读取fam文件"""
    fam_path = bed_prefix + '.fam'
    # fam文件格式: FID IID FATHER MOTHER SEX PHENO
    # 没有列名，使用列索引
    fam_df = pd.read_csv(
        fam_path, 
        sep=r'\s+', 
        header=None,
        names=['FID', 'IID', 'FATHER', 'MOTHER', 'SEX', 'PHENO']
    )
    return fam_df


def _merge_sample_info(fam_df: pd.DataFrame, info_df: pd.DataFrame, 
                       id_column: str, dp_column: str, pheno_column: str) -> pd.DataFrame:
    """合并fam文件和样本信息"""
    # 使用IID列匹配info文件的ID列
    merged = fam_df.merge(
        info_df[[id_column, dp_column] + ([pheno_column] if pheno_column in info_df.columns else [])],
        left_on='IID',
        right_on=id_column,
        how='left'
    )
    
    # 检查有多少样本成功匹配
    n_matched = merged[dp_column].notna().sum()
    if n_matched == 0:
        raise ValueError("没有样本成功匹配！请检查ID列名和样本ID格式")
    
    logger.info(f"成功匹配 {n_matched}/{len(fam_df)} 个样本的深度信息")
    
    return merged


def _group_samples_by_dp(merged_df: pd.DataFrame, dp_column: str) -> dict:
    """按测序深度分组样本"""
    # 提取深度值（假设格式为"30X"或"15X"）
    def extract_dp(dp_str):
        if pd.isna(dp_str):
            return None
        dp_str = str(dp_str).strip().upper()
        if '30' in dp_str or dp_str.startswith('30'):
            return '30X'
        elif '15' in dp_str or dp_str.startswith('15'):
            return '15X'
        else:
            return None
    
    merged_df['DP_GROUP'] = merged_df[dp_column].apply(extract_dp)
    
    # 统计各组样本数
    group_counts = merged_df['DP_GROUP'].value_counts()
    
    return {
        'data': merged_df,
        'count': group_counts
    }


def _create_sample_list_files(tmpdir: str, merged_df: pd.DataFrame, 
                               dp_groups: dict, maf_group: str, 
                               pheno_column: str) -> dict:
    """创建各组样本列表文件"""
    sample_lists = {}
    
    merged_df = dp_groups['data']
    
    # 创建15X组样本列表
    df_15x = merged_df[merged_df['DP_GROUP'] == '15X']
    if len(df_15x) > 0:
        list_15x = os.path.join(tmpdir, 'samples_15x.txt')
        df_15x[['FID', 'IID']].to_csv(list_15x, sep='\t', header=False, index=False)
        sample_lists['15X'] = list_15x
        logger.info(f"创建15X组样本列表: {len(df_15x)} 个样本")
    
    # 创建30X组样本列表
    df_30x = merged_df[merged_df['DP_GROUP'] == '30X']
    if len(df_30x) > 0:
        list_30x = os.path.join(tmpdir, 'samples_30x.txt')
        df_30x[['FID', 'IID']].to_csv(list_30x, sep='\t', header=False, index=False)
        sample_lists['30X'] = list_30x
        logger.info(f"创建30X组样本列表: {len(df_30x)} 个样本")
    
    # 创建MAF计算的样本列表
    if maf_group == "ctrl" and pheno_column in merged_df.columns:
        # 对照组：OUTCOME2列中值为"AGP3K"的样本
        df_ctrl = merged_df[
            merged_df[pheno_column].astype(str).str.upper() == 'AGP3K'
        ]
        if len(df_ctrl) > 0:
            list_ctrl = os.path.join(tmpdir, 'samples_ctrl.txt')
            df_ctrl[['FID', 'IID']].to_csv(list_ctrl, sep='\t', header=False, index=False)
            sample_lists['ctrl'] = list_ctrl
            logger.info(f"创建对照组样本列表 (AGP3K): {len(df_ctrl)} 个样本")
        else:
            logger.warning("未找到对照组样本 (AGP3K)，将使用所有样本计算MAF")
    
    elif maf_group == "case" and pheno_column in merged_df.columns:
        # 病例组：OUTCOME2列中值为"CTEPH"的样本
        df_case = merged_df[
            merged_df[pheno_column].astype(str).str.upper() == 'CTEPH'
        ]
        if len(df_case) > 0:
            list_case = os.path.join(tmpdir, 'samples_case.txt')
            df_case[['FID', 'IID']].to_csv(list_case, sep='\t', header=False, index=False)
            sample_lists['case'] = list_case
            logger.info(f"创建病例组样本列表 (CTEPH): {len(df_case)} 个样本")
        else:
            logger.warning("未找到病例组样本 (CTEPH)，将使用所有样本计算MAF")
    
    # maf_group为"all"时不需要创建额外列表
    
    return sample_lists


def _calculate_vmiss_plink(bed_prefix: str, plink2_path: str, tmpdir: str,
                            suffix: str, sample_list: Optional[str] = None,
                            n_threads: int = 4) -> str:
    """
    使用plink2计算VMISS
    返回输出文件路径而不是DataFrame，避免加载大数据到内存
    """
    output_prefix = os.path.join(tmpdir, f'vmiss_{suffix}')
    
    cmd = [
        plink2_path,
        '--bfile', bed_prefix,
        '--missing', 'variant-only',
        '--out', output_prefix,
        '--threads', str(n_threads)
    ]
    
    if sample_list is not None:
        cmd.extend(['--keep', sample_list])
    
    logger.info(f"执行命令: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        logger.debug(f"plink2标准输出:\n{result.stdout}")
    except subprocess.CalledProcessError as e:
        logger.error(f"plink2执行失败:\n{e.stderr}")
        raise
    
    # 返回结果文件路径而不是读取
    vmiss_file = output_prefix + '.vmiss'
    if not os.path.exists(vmiss_file):
        raise FileNotFoundError(f"找不到plink2输出文件: {vmiss_file}")
    
    # 只统计行数，不加载数据
    with open(vmiss_file, 'r') as f:
        n_lines = sum(1 for _ in f) - 1  # 减去表头
    logger.info(f"VMISS文件包含 {n_lines} 个变体")
    
    return vmiss_file


def _calculate_maf_plink(bed_prefix: str, plink2_path: str, tmpdir: str,
                         maf_group: str, sample_list: Optional[str] = None,
                         n_threads: int = 4) -> str:
    """
    使用plink2计算MAF
    返回输出文件路径而不是DataFrame，避免加载大数据到内存
    """
    output_prefix = os.path.join(tmpdir, f'freq_{maf_group}')
    
    cmd = [
        plink2_path,
        '--bfile', bed_prefix,
        '--freq',
        '--out', output_prefix,
        '--threads', str(n_threads)
    ]
    
    if sample_list is not None:
        cmd.extend(['--keep', sample_list])
    
    logger.info(f"执行命令: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        logger.debug(f"plink2标准输出:\n{result.stdout}")
    except subprocess.CalledProcessError as e:
        logger.error(f"plink2执行失败:\n{e.stderr}")
        raise
    
    # 返回结果文件路径而不是读取
    freq_file = output_prefix + '.afreq'
    if not os.path.exists(freq_file):
        raise FileNotFoundError(f"找不到plink2输出文件: {freq_file}")
    
    # 只统计行数，不加载数据
    with open(freq_file, 'r') as f:
        n_lines = sum(1 for _ in f) - 1  # 减去表头
    logger.info(f"频率文件包含 {n_lines} 个变体")
    
    return freq_file


def _merge_and_write_results_streaming(vmiss_all_file: str, vmiss_15x_file: Optional[str],
                                        vmiss_30x_file: Optional[str], maf_file: str,
                                        output_path: str) -> int:
    """
    完全流式处理：直接从plink2输出文件读取并合并写入
    极低内存占用，适合处理数百万变体
    
    参数:
        vmiss_all_file: 总体VMISS文件路径
        vmiss_15x_file: 15X组VMISS文件路径（可选）
        vmiss_30x_file: 30X组VMISS文件路径（可选）
        maf_file: MAF文件路径
        output_path: 输出文件路径
    
    返回:
        int: 写入的变体数量
    """
    import gc
    
    logger.info("开始流式构建索引字典...")
    
    # 分批构建字典的通用函数
    def build_dict_streaming(filepath, key_col, val_col, chunksize=50000):
        """流式构建字典，每次只处理一部分数据"""
        result_dict = {}
        if filepath is None or not os.path.exists(filepath):
            return result_dict
        
        chunk_count = 0
        for chunk in pd.read_csv(filepath, sep=r'\s+', chunksize=chunksize, 
                                  usecols=[key_col, val_col], dtype=str):  # 只读取需要的列，减少内存
            if key_col in chunk.columns and val_col in chunk.columns:
                chunk_dict = dict(zip(chunk[key_col].values, chunk[val_col].values))
                result_dict.update(chunk_dict)
                del chunk, chunk_dict
                chunk_count += 1
                if chunk_count % 5 == 0:  # 更频繁的垃圾回收
                    gc.collect()
        
        return result_dict
    
    # 构建各个字典
    logger.info("构建15X组索引...")
    vmiss_15x_dict = build_dict_streaming(vmiss_15x_file, 'ID', 'F_MISS')
    logger.info(f"15X组索引完成: {len(vmiss_15x_dict)} 个变体")
    
    logger.info("构建30X组索引...")
    vmiss_30x_dict = build_dict_streaming(vmiss_30x_file, 'ID', 'F_MISS')
    logger.info(f"30X组索引完成: {len(vmiss_30x_dict)} 个变体")
    
    logger.info("构建MAF索引...")
    maf_dict = build_dict_streaming(maf_file, 'ID', 'ALT_FREQS')
    logger.info(f"MAF索引完成: {len(maf_dict)} 个变体")
    
    # 流式读取主文件并写入结果
    logger.info("开始流式合并并写入...")
    n_variants = 0
    
    with open(output_path, 'w', buffering=8192*16) as outf:  # 增大写入缓冲区
        # 写入表头
        outf.write('VARIANT_ID\tVMISS\tVMISS_15X\tVMISS_30X\tMAF\n')
        
        # 逐块处理主VMISS文件，减小chunksize减少内存占用
        chunksize = 20000  # 从50000减小到20000
        for chunk in pd.read_csv(vmiss_all_file, sep=r'\s+', chunksize=chunksize):
            # 使用向量化操作而不是iterrows()，更快且内存效率更高
            for variant_id, vmiss in zip(chunk['ID'].values, chunk['F_MISS'].values):
                vmiss_15x = vmiss_15x_dict.get(variant_id, 'NA')
                vmiss_30x = vmiss_30x_dict.get(variant_id, 'NA')
                
                # 获取ALT频率并转换为真正的MAF（次要等位基因频率）
                alt_freq = maf_dict.get(variant_id, 'NA')
                if alt_freq != 'NA':
                    try:
                        alt_freq_val = float(alt_freq)
                        # MAF = min(ALT_FREQ, 1 - ALT_FREQ)
                        maf = min(alt_freq_val, 1 - alt_freq_val)
                    except (ValueError, TypeError):
                        maf = 'NA'
                else:
                    maf = 'NA'
                
                outf.write(f"{variant_id}\t{vmiss}\t{vmiss_15x}\t{vmiss_30x}\t{maf}\n")
                n_variants += 1
            
            del chunk
            if n_variants % 50000 == 0:
                logger.info(f"已处理 {n_variants} 个变体...")
                gc.collect()
    
    # 清理字典
    del vmiss_15x_dict, vmiss_30x_dict, maf_dict
    gc.collect()
    
    logger.info(f"流式处理完成，共写入 {n_variants} 个变体")
    return n_variants


def _merge_and_write_results(vmiss_all: pd.DataFrame, vmiss_15x: pd.DataFrame,
                             vmiss_30x: pd.DataFrame, maf_data: pd.DataFrame,
                             output_path: str) -> int:
    """
    使用完全流式处理合并结果，极低内存占用
    不在内存中保存完整数据，而是逐行读取和写入
    
    返回:
        int: 写入的变体数量
    """
    logger.info("正在准备数据字典（流式处理）...")
    
    # 分批创建字典以减少内存峰值
    # 使用chunksize分批读取（如果数据已经在内存中，先保存到临时文件）
    import gc
    
    # 创建临时文件保存各个数据
    tmpdir = os.path.dirname(output_path)
    tmp_vmiss_all = os.path.join(tmpdir, '.tmp_vmiss_all.txt')
    tmp_vmiss_15x = os.path.join(tmpdir, '.tmp_vmiss_15x.txt')
    tmp_vmiss_30x = os.path.join(tmpdir, '.tmp_vmiss_30x.txt')
    tmp_maf = os.path.join(tmpdir, '.tmp_maf.txt')
    
    # 保存数据到临时文件并立即释放内存
    logger.info("保存中间结果到临时文件...")
    vmiss_all[['ID', 'F_MISS']].to_csv(tmp_vmiss_all, sep='\t', index=False)
    del vmiss_all
    gc.collect()
    
    if len(vmiss_15x) > 0:
        vmiss_15x[['ID', 'F_MISS']].to_csv(tmp_vmiss_15x, sep='\t', index=False)
    del vmiss_15x
    gc.collect()
    
    if len(vmiss_30x) > 0:
        vmiss_30x[['ID', 'F_MISS']].to_csv(tmp_vmiss_30x, sep='\t', index=False)
    del vmiss_30x
    gc.collect()
    
    if len(maf_data) > 0:
        maf_data[['ID', 'ALT_FREQS']].to_csv(tmp_maf, sep='\t', index=False)
    del maf_data
    gc.collect()
    
    logger.info("构建索引字典（分批处理）...")
    
    # 分批读取构建字典，每次只处理一部分数据
    def build_dict_in_chunks(filepath, key_col, val_col, chunksize=100000):
        """分批构建字典以避免内存峰值"""
        result_dict = {}
        if not os.path.exists(filepath):
            return result_dict
        
        for chunk in pd.read_csv(filepath, sep='\t', chunksize=chunksize):
            chunk_dict = dict(zip(chunk[key_col], chunk[val_col]))
            result_dict.update(chunk_dict)
            del chunk, chunk_dict
            gc.collect()
        
        return result_dict
    
    vmiss_15x_dict = build_dict_in_chunks(tmp_vmiss_15x, 'ID', 'F_MISS')
    vmiss_30x_dict = build_dict_in_chunks(tmp_vmiss_30x, 'ID', 'F_MISS')
    maf_dict = build_dict_in_chunks(tmp_maf, 'ID', 'ALT_FREQS')
    
    logger.info(f"索引字典构建完成: 15X={len(vmiss_15x_dict)}, 30X={len(vmiss_30x_dict)}, MAF={len(maf_dict)}")
    
    # 流式处理：逐行读取、合并、写入
    logger.info("开始流式合并并写入结果...")
    n_variants = 0
    
    with open(output_path, 'w') as outf:
        # 写入表头
        outf.write('VARIANT_ID\tVMISS\tVMISS_15X\tVMISS_30X\tMAF\n')
        
        # 逐块读取并处理
        chunksize = 50000  # 每次处理5万行
        for chunk in pd.read_csv(tmp_vmiss_all, sep='\t', chunksize=chunksize):
            # 处理当前块
            for _, row in chunk.iterrows():
                variant_id = row['ID']
                vmiss = row['F_MISS']
                vmiss_15x = vmiss_15x_dict.get(variant_id, '')
                vmiss_30x = vmiss_30x_dict.get(variant_id, '')
                maf = maf_dict.get(variant_id, '')
                
                # 写入一行
                outf.write(f"{variant_id}\t{vmiss}\t{vmiss_15x}\t{vmiss_30x}\t{maf}\n")
                n_variants += 1
            
            # 每处理一块就清理一次
            del chunk
            gc.collect()
            
            if n_variants % 100000 == 0:
                logger.info(f"已处理 {n_variants} 个变体...")
    
    # 清理字典和临时文件
    del vmiss_15x_dict, vmiss_30x_dict, maf_dict
    gc.collect()
    
    # 删除临时文件
    for tmp_file in [tmp_vmiss_all, tmp_vmiss_15x, tmp_vmiss_30x, tmp_maf]:
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
    
    logger.info(f"流式处理完成，共写入 {n_variants} 个变体")
    
    return n_variants


# new function

def analyze_vmiss_thresholds(
    variant_metrics_file: str,
    output_dir: str,
    analysis_mode: Literal[1, 2, 3, 4] = 4,
    maf_thresholds: tuple = (0.01, 0.05),
    hist_step: float = 0.01,
    cdf_step: Optional[float] = None,
    cdf_step_15x: Optional[float] = None,
    cdf_step_30x: Optional[float] = None,
    bin_position: Literal['left', 'right'] = 'left',
    knee_curve: str = 'concave',
    knee_direction: str = 'increasing',
    knee_S: float = 1.0,
    knee_weight_x: float = 1.0,
    knee_weight_y: float = 1.0,
    dpi: int = 600,
    chunksize: int = 50000,
    n_threads: int = 4
) -> dict:
    """
    分析变体缺失率(VMISS)阈值，寻找最佳过滤参数
    
    参数:
        variant_metrics_file: 变体指标文件路径（calculate_variant_metrics的输出）
        output_dir: 输出目录路径
        analysis_mode: 分析模式
            1: 针对所有变体，识别VMISS（1个参数建议）
            2: 针对所有变体，识别VMISS_15X和VMISS_30X（2个参数建议）
            3: 针对三组MAF变体，识别VMISS（3个参数建议）
            4: 针对三组MAF变体，识别VMISS_15X和VMISS_30X（6个参数建议）
        maf_thresholds: MAF分组阈值，默认(0.01, 0.05)
        hist_step: 直方图bin的步长，默认0.01（固定，用于绘制分布）
        cdf_step: CDF计算的步长，默认None（自动根据模式设置）
            - 模式1/3: 默认0.001
            - 模式2/4: 由cdf_step_15x和cdf_step_30x分别指定
        cdf_step_15x: 15X组CDF计算的步长，默认0.001（仅模式2/4有效）
        cdf_step_30x: 30X组CDF计算的步长，默认0.01（仅模式2/4有效）
        bin_position: CDF的x轴坐标定义方式，默认'left'
            'left': 使用bin的左边界（CDF表示≤x的累积比例）
            'right': 使用bin的右边界（CDF表示<x的累积比例，与plink2 --geno行为一致）
        knee_curve: kneed曲线类型，默认'concave'
        knee_direction: kneed方向，默认'increasing'
        knee_S: kneed敏感度参数
        knee_weight_x: kneed x轴权重
        knee_weight_y: kneed y轴权重
        dpi: 图片分辨率，默认600
        chunksize: 流式读取的块大小
        n_threads: 并行线程数
    
    返回:
        dict: 包含所有分析结果的字典，包括拐点信息、图片路径、统计数据等
    """
    import gc
    import json
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')  # 非交互式后端
    from kneed import KneeLocator 
    
    logger.info("="*60)
    logger.info("开始VMISS阈值分析")
    logger.info(f"输入文件: {variant_metrics_file}")
    logger.info(f"分析模式: {analysis_mode}")
    logger.info(f"MAF分组阈值: {maf_thresholds}")
    logger.info(f"直方图步长(hist_step): {hist_step}")
    
    # 根据模式自动设置CDF步长
    if analysis_mode in [1, 3]:
        # 模式1和3：不分15X和30X，使用统一的cdf_step
        if cdf_step is None:
            cdf_step = 0.001  # 默认0.001
        logger.info(f"CDF步长: {cdf_step}")
        cdf_step_dict = {'default': cdf_step}
    elif analysis_mode in [2, 4]:
        # 模式2和4：分15X和30X，分别设置步长
        if cdf_step_15x is None:
            cdf_step_15x = 0.001  # 15X默认0.001
        if cdf_step_30x is None:
            cdf_step_30x = 0.01   # 30X默认0.01
        logger.info(f"CDF步长(15X): {cdf_step_15x}")
        logger.info(f"CDF步长(30X): {cdf_step_30x}")
        cdf_step_dict = {'15X': cdf_step_15x, '30X': cdf_step_30x}
    
    logger.info("="*60)
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 第一步：流式读取数据并分组统计（内存优化）
    logger.info("步骤1: 流式读取并分组数据...")
    group_stats = _streaming_group_variants(
        variant_metrics_file, 
        analysis_mode, 
        maf_thresholds, 
        chunksize
    )
    
    # 第二步：计算累积分布
    logger.info("步骤2: 计算VMISS累积分布...")
    cumulative_results = _calculate_cumulative_distribution(
        group_stats, 
        hist_step,
        cdf_step_dict,
        analysis_mode,
        bin_position
    )
    
    # 第三步：寻找累积分布的拐点
    logger.info("步骤3: 寻找累积分布拐点...")
    knee_results = _find_cumulative_knee_points(
        cumulative_results,
        knee_curve,
        knee_direction,
        knee_S,
        knee_weight_x,
        knee_weight_y
    )
    
    # 第四步：绘制分布图和累积分布图
    logger.info("步骤4: 绘制分布和累积分布图...")
    plot_path = _plot_vmiss_distributions(
        cumulative_results,
        knee_results,
        analysis_mode,
        output_dir,
        dpi
    )
    
    # 第五步：保存结果到JSON（排除临时数组数据）
    logger.info("步骤5: 保存结果到JSON...")
    json_path = os.path.join(output_dir, 'vmiss_threshold_analysis.json')
    
    # 清理group_stats中的临时数组（以_开头的键）
    group_stats_for_json = {}
    for key, value in group_stats.items():
        group_stats_for_json[key] = {k: v for k, v in value.items() if not k.startswith('_')}
    
    results = {
        'analysis_mode': analysis_mode,
        'maf_thresholds': maf_thresholds,
        'hist_step': hist_step,
        'cdf_step_dict': cdf_step_dict,
        'bin_position': bin_position,
        'knee_parameters': {
            'curve': knee_curve,
            'direction': knee_direction,
            'S': knee_S,
            'weight_x': knee_weight_x,
            'weight_y': knee_weight_y
        },
        'group_statistics': group_stats_for_json,
        'knee_points': knee_results,
        'plot_path': plot_path,
        'json_path': json_path
    }
    
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    
    logger.info("="*60)
    logger.info(f"分析完成！")
    logger.info(f"结果JSON: {json_path}")
    logger.info(f"分析图片: {plot_path}")
    logger.info("="*60)
    
    # 打印拐点建议
    _print_knee_recommendations(knee_results, analysis_mode)
    
    return results


def _streaming_group_variants(
    variant_file: str,
    analysis_mode: int,
    maf_thresholds: tuple,
    chunksize: int
) -> dict:
    """
    流式读取变体数据并分组统计
    内存高效：不一次性加载所有数据
    """
    import gc
    
    # 初始化分组统计
    group_stats = {}
    
    # 根据分析模式确定需要的列
    required_cols = ['VARIANT_ID', 'MAF']
    vmiss_cols = []
    
    if analysis_mode in [1, 3]:
        vmiss_cols = ['VMISS']
    elif analysis_mode in [2, 4]:
        vmiss_cols = ['VMISS_15X', 'VMISS_30X']
    
    required_cols.extend(vmiss_cols)
    
    # 定义MAF分组（如果需要）
    if analysis_mode in [3, 4]:
        maf_low, maf_mid = maf_thresholds
        maf_groups = {
            f'MAF<{maf_low}': lambda x: x < maf_low,
            f'MAF>={maf_low}&<={maf_mid}': lambda x: (x >= maf_low) & (x <= maf_mid),
            f'MAF>{maf_mid}': lambda x: x > maf_mid
        }
    else:
        maf_groups = {'All': lambda x: pd.Series([True] * len(x))}
    
    # 初始化每组的数据收集器
    for maf_group in maf_groups:
        for vmiss_col in vmiss_cols:
            group_key = f"{maf_group}_{vmiss_col}" if len(maf_groups) > 1 else vmiss_col
            group_stats[group_key] = {
                'vmiss_values': [],  # 存储VMISS值（用于后续阈值分析）
                'total_variants': 0,
                'vmiss_col': vmiss_col,
                'maf_group': maf_group
            }
    
    # 流式读取并分组
    logger.info(f"开始流式读取变体文件（chunksize={chunksize}）...")
    chunk_count = 0
    
    for chunk in pd.read_csv(variant_file, sep='\t', chunksize=chunksize, 
                              usecols=required_cols):
        chunk_count += 1
        
        # 转换MAF为数值类型，处理'NA'
        chunk['MAF'] = pd.to_numeric(chunk['MAF'], errors='coerce')
        
        # 对每个MAF分组处理
        for maf_group_name, maf_filter in maf_groups.items():
            if maf_group_name == 'All':
                # 对于'All'组，直接使用整个chunk，不需要过滤
                group_data = chunk
            else:
                # 对于特定MAF分组，使用过滤函数
                maf_mask = maf_filter(chunk['MAF'])
                group_data = chunk[maf_mask]
            
            # 对每个VMISS列处理
            for vmiss_col in vmiss_cols:
                group_key = f"{maf_group_name}_{vmiss_col}" if len(maf_groups) > 1 else vmiss_col
                
                # 转换VMISS为数值，过滤NA
                vmiss_series = pd.to_numeric(group_data[vmiss_col], errors='coerce')
                valid_vmiss = vmiss_series.dropna().values
                
                # 累积统计
                group_stats[group_key]['vmiss_values'].extend(valid_vmiss.tolist())
                group_stats[group_key]['total_variants'] += len(valid_vmiss)
        
        # 定期清理内存
        del chunk
        if chunk_count % 10 == 0:
            logger.info(f"已处理 {chunk_count} 个数据块...")
            gc.collect()
    
    # 转换为numpy数组以加速后续计算
    logger.info("转换数据为numpy数组...")
    for group_key in group_stats:
        vmiss_list = group_stats[group_key]['vmiss_values']
        vmiss_arr = np.array(vmiss_list, dtype=np.float32)
        
        # 将数组存储用于后续计算，但不保存到最终结果中
        group_stats[group_key]['_vmiss_array'] = vmiss_arr  # 下划线表示临时数据
        # 删除列表释放内存
        del group_stats[group_key]['vmiss_values']
        
        # 计算统计信息
        if len(vmiss_arr) > 0:
            # 使用实际的最小值和最大值，不要过度取整，这样可以保留分布的细节
            group_stats[group_key]['vmiss_min'] = float(vmiss_arr.min())
            group_stats[group_key]['vmiss_max'] = float(vmiss_arr.max())
            group_stats[group_key]['vmiss_mean'] = float(vmiss_arr.mean())
            group_stats[group_key]['vmiss_median'] = float(np.median(vmiss_arr))
        else:
            group_stats[group_key]['vmiss_min'] = 0
            group_stats[group_key]['vmiss_max'] = 0
            group_stats[group_key]['vmiss_mean'] = 0
            group_stats[group_key]['vmiss_median'] = 0
        
        logger.info(f"组 [{group_key}]: {group_stats[group_key]['total_variants']} 变体, "
                   f"VMISS范围 [{group_stats[group_key]['vmiss_min']:.2f}, "
                   f"{group_stats[group_key]['vmiss_max']:.2f}]")
    
    gc.collect()
    return group_stats


def _calculate_cumulative_distribution(
    group_stats: dict,
    hist_step: float,
    cdf_step_dict: dict,
    analysis_mode: int,
    bin_position: str = 'left'
) -> dict:
    """
    计算每组VMISS的分布和累积分布
    
    参数:
        group_stats: 分组统计数据
        hist_step: 直方图bin的步长（固定，用于绘制分布）
        cdf_step_dict: CDF计算步长字典
            - 模式1/3: {'default': 0.001}
            - 模式2/4: {'15X': 0.001, '30X': 0.01}
        analysis_mode: 分析模式（用于确定如何选择cdf_step）
        bin_position: CDF的x轴坐标定义方式 ('left', 'right')
    """
    cumulative_results = {}
    
    for group_key, stats in group_stats.items():
        logger.info(f"计算组 [{group_key}] 的累积分布...")
        
        vmiss_arr = stats['_vmiss_array']
        if len(vmiss_arr) == 0:
            logger.warning(f"组 [{group_key}] 没有有效数据，跳过")
            continue
        
        vmiss_min = stats['vmiss_min']
        vmiss_max = stats['vmiss_max']
        
        # 根据组名确定使用哪个CDF步长
        if analysis_mode in [1, 3]:
            # 模式1/3：所有组使用相同步长
            cdf_step = cdf_step_dict['default']
        elif analysis_mode in [2, 4]:
            # 模式2/4：根据组名判断是15X还是30X
            vmiss_col = stats.get('vmiss_col', '')
            if '15X' in vmiss_col.upper():
                cdf_step = cdf_step_dict['15X']
            elif '30X' in vmiss_col.upper():
                cdf_step = cdf_step_dict['30X']
            else:
                # 备用：如果无法判断，使用默认值
                cdf_step = cdf_step_dict.get('15X', 0.001)
                logger.warning(f"无法确定组 [{group_key}] 的步长，使用默认值 {cdf_step}")
        
        logger.info(f"  直方图步长: {hist_step}, CDF步长: {cdf_step}")
        
        # 将min和max对齐到cdf_step的倍数（CDF计算使用cdf_step）
        # 向下取整到最近的cdf_step倍数
        vmiss_min_aligned = np.floor(vmiss_min / cdf_step) * cdf_step
        # 向上取整到最近的cdf_step倍数
        vmiss_max_aligned = np.ceil(vmiss_max / cdf_step) * cdf_step
        
        # 对于right边界，需要额外扩展一个步长以确保最后一个bin的右边界包含所有数据
        if bin_position == 'right':
            vmiss_max_aligned += cdf_step
        
        # 确保范围至少包含3个bins
        if vmiss_max_aligned - vmiss_min_aligned < cdf_step * 3:
            # 向两边各扩展以保证至少3个bins
            center = (vmiss_min_aligned + vmiss_max_aligned) / 2
            vmiss_min_aligned = center - cdf_step * 1.5
            vmiss_max_aligned = center + cdf_step * 1.5
        
        # 确保在合理范围内（VMISS应该在[0, 1]之间）
        vmiss_min_aligned = max(0, vmiss_min_aligned)
        vmiss_max_aligned = min(1, vmiss_max_aligned)
        
        # 生成VMISS值序列（用于CDF计算，使用cdf_step）
        # 使用arange确保每个bin宽度都是cdf_step
        vmiss_bins = np.arange(vmiss_min_aligned, vmiss_max_aligned + cdf_step/2, cdf_step)
        vmiss_bins = np.round(vmiss_bins, 6)  # 提高精度以避免浮点误差
        
        # 计算直方图（VMISS分布）
        hist_counts, _ = np.histogram(vmiss_arr, bins=vmiss_bins)
        
        # 计算累积分布（从小到大累积）
        cumulative_counts = np.cumsum(hist_counts)
        cumulative_percentage = (cumulative_counts / len(vmiss_arr)) * 100
        
        # 根据bin_position参数确定x轴坐标
        if bin_position == 'left':
            # 左边界：使用bin的左边界（下界）
            # CDF含义：≤ x 的变体占比
            vmiss_x_values = vmiss_bins[:-1]
        elif bin_position == 'right':
            # 右边界：使用bin的右边界（上界）
            # CDF含义：< x 的变体占比（与plink2的--geno一致）
            vmiss_x_values = vmiss_bins[1:]
        else:
            raise ValueError(f"Invalid bin_position: {bin_position}. Must be 'left' or 'right'.")
        
        # 不强制添加起始点，直接使用实际数据的CDF
        hist_counts_padded = hist_counts
        
        cumulative_results[group_key] = {
            'vmiss_values': vmiss_x_values.tolist(),  # x轴：VMISS值（根据bin_position确定）
            'histogram': hist_counts_padded.tolist(),     # 直方图计数
            'cumulative_counts': cumulative_counts.tolist(),  # 累积计数
            'cumulative_percentage': cumulative_percentage.tolist(),  # 累积百分比
            'total_variants': stats['total_variants'],
            'vmiss_col': stats['vmiss_col'],
            'maf_group': stats['maf_group'],
            'vmiss_min': vmiss_min,
            'vmiss_max': vmiss_max,
            'vmiss_mean': stats['vmiss_mean'],
            'vmiss_median': stats['vmiss_median'],
            'hist_step': hist_step,  # 记录直方图步长
            'cdf_step': cdf_step,    # 记录CDF步长
            'bin_position': bin_position,  # 记录bin位置定义方式
            'n_bins': len(hist_counts)  # 记录bin数量
        }
        
        logger.info(f"  完成！原始范围: [{vmiss_min:.4f}, {vmiss_max:.4f}], "
                   f"对齐范围: [{vmiss_min_aligned:.4f}, {vmiss_max_aligned:.4f}], "
                   f"{len(hist_counts)} bins (CDF步长={cdf_step})")
    
    return cumulative_results


def _find_cumulative_knee_points(
    cumulative_results: dict,
    knee_curve: str,
    knee_direction: str,
    knee_S: float,
    knee_weight_x: float,
    knee_weight_y: float
) -> dict:
    """
    在累积分布曲线上寻找拐点
    
    优化策略：
    1. 根据CDF步长自适应调整knee_S（步长越小，knee_S越大）
    2. 过滤起始的低累积区域（避免误判陡峭起始部分）
    3. 使用更强的平滑参数减少噪声影响
    """
    from kneed import KneeLocator
    
    knee_results = {}
    
    for group_key, result in cumulative_results.items():
        logger.info(f"寻找组 [{group_key}] 累积分布的拐点...")
        
        vmiss_values = np.array(result['vmiss_values'])
        cumulative_pct = np.array(result['cumulative_percentage'])
        total = result['total_variants']
        cdf_step = result.get('cdf_step', 0.01)
        
        logger.info(f"  CDF数据范围: VMISS [{vmiss_values[0]:.3f}, {vmiss_values[-1]:.3f}], "
                   f"累积% [{cumulative_pct[0]:.1f}%, {cumulative_pct[-1]:.1f}%], "
                   f"{len(vmiss_values)} 个数据点, CDF步长={cdf_step}")
        
        # 检查数据点是否足够
        if len(vmiss_values) < 3:
            logger.warning(f"  组 [{group_key}] 数据点不足（<3），跳过拐点检测")
            knee_results[group_key] = {
                'knee_found': False,
                'message': '数据点不足'
            }
            continue
        
        # 自适应调整knee_S：CDF步长越小，需要更高的knee_S来避免噪声
        adaptive_knee_S = knee_S
        if cdf_step <= 0.001:
            # 超高精度：大幅增加knee_S
            adaptive_knee_S = knee_S * 3.0
            logger.info(f"  检测到超高精度CDF（步长={cdf_step}），自适应增加knee_S: {knee_S} → {adaptive_knee_S}")
        elif cdf_step <= 0.005:
            # 高精度：适度增加knee_S
            adaptive_knee_S = knee_S * 2.0
            logger.info(f"  检测到高精度CDF（步长={cdf_step}），自适应增加knee_S: {knee_S} → {adaptive_knee_S}")
        
        # 策略1：过滤起始的低累积区域，避免误判陡峭起始
        # 只保留累积百分比>=5%的数据点进行拐点检测
        min_cumulative_threshold = 5.0  # 最小累积百分比阈值
        valid_mask = cumulative_pct >= min_cumulative_threshold
        
        if valid_mask.sum() < 3:
            # 如果过滤后数据点不足，降低阈值重试
            min_cumulative_threshold = 1.0
            valid_mask = cumulative_pct >= min_cumulative_threshold
            logger.info(f"  降低累积阈值至{min_cumulative_threshold}%以保留足够数据点")
        
        if valid_mask.sum() < 3:
            logger.warning(f"  组 [{group_key}] 过滤后数据点仍不足，跳过拐点检测")
            knee_results[group_key] = {
                'knee_found': False,
                'message': f'过滤后数据点不足（累积<{min_cumulative_threshold}%）'
            }
            continue
        
        # 使用过滤后的数据
        vmiss_filtered = vmiss_values[valid_mask]
        cumulative_pct_filtered = cumulative_pct[valid_mask]
        
        logger.info(f"  过滤后数据: {len(vmiss_filtered)} 个点 "
                   f"(VMISS范围 [{vmiss_filtered[0]:.3f}, {vmiss_filtered[-1]:.3f}], "
                   f"累积范围 [{cumulative_pct_filtered[0]:.1f}%, {cumulative_pct_filtered[-1]:.1f}%])")
        
        try:
            # 在过滤后的累积分布曲线上寻找拐点
            # 使用自适应的knee_S参数
            kneedle = KneeLocator(
                vmiss_filtered,
                cumulative_pct_filtered,
                S=adaptive_knee_S,
                curve=knee_curve,
                direction=knee_direction,
                weight_x=knee_weight_x,
                weight_y=knee_weight_y
            )
            
            if kneedle.knee is not None:
                knee_vmiss = kneedle.knee
                # 在原始（未过滤）数据中找到对应的索引
                knee_index = np.argmin(np.abs(vmiss_values - knee_vmiss))
                knee_cumulative_pct = float(cumulative_pct[knee_index])
                knee_cumulative_count = int(result['cumulative_counts'][knee_index])
                
                knee_results[group_key] = {
                    'knee_found': True,
                    'knee_vmiss': float(knee_vmiss),  # 拐点对应的VMISS值
                    'knee_cumulative_percentage': knee_cumulative_pct,  # 累积百分比
                    'knee_cumulative_count': knee_cumulative_count,  # 累积变体数
                    'total_variants': total,
                    'vmiss_col': result['vmiss_col'],
                    'maf_group': result['maf_group'],
                    'adaptive_knee_S': adaptive_knee_S,  # 记录使用的knee_S
                    'min_cumulative_filter': min_cumulative_threshold  # 记录过滤阈值
                }
                
                logger.info(f"  ✓ 找到拐点: VMISS={knee_vmiss:.4f}, "
                           f"累积{knee_cumulative_pct:.1f}% ({knee_cumulative_count:,}变体), "
                           f"使用knee_S={adaptive_knee_S:.1f}")
            else:
                logger.warning(f"  未找到明显拐点")
                knee_results[group_key] = {
                    'knee_found': False,
                    'message': '未检测到明显拐点'
                }
        
        except Exception as e:
            logger.error(f"  拐点检测失败: {str(e)}")
            knee_results[group_key] = {
                'knee_found': False,
                'message': f'检测失败: {str(e)}'
            }
    
    return knee_results


def _plot_vmiss_distributions(
    cumulative_results: dict,
    knee_results: dict,
    analysis_mode: int,
    output_dir: str,
    dpi: int
) -> str:
    """
    绘制VMISS分布图和累积分布曲线
    每个子图包含:
    - 左Y轴: 直方图（VMISS分布）- 始终使用hist_step=0.01宽度显示
    - 右Y轴: 累积分布曲线 - 使用cdf_step精度计算
    - 标记拐点位置
    """
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    
    # 清空当前所有图形，避免内存泄漏
    plt.close('all')
    
    # 设置字体
    rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    rcParams['axes.unicode_minus'] = False
    
    # 确定子图布局
    n_groups = len(cumulative_results)
    if n_groups == 0:
        logger.warning("没有数据可绘制")
        return ""
    
    # 根据分析模式确定布局
    if analysis_mode == 1:
        nrows, ncols = 1, 1
    elif analysis_mode == 2:
        nrows, ncols = 1, 2
    elif analysis_mode == 3:
        nrows, ncols = 3, 1
    elif analysis_mode == 4:
        nrows, ncols = 3, 2
    
    # 创建图形
    fig, axes = plt.subplots(nrows, ncols, figsize=(8*ncols, 6*nrows), dpi=dpi)
    if n_groups == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    # 绘制每个组
    for idx, (group_key, result) in enumerate(cumulative_results.items()):
        ax1 = axes[idx]
        
        # 转换为numpy数组以确保正确绘图
        vmiss_values = np.array(result['vmiss_values'])  # CDF的x轴坐标（cdf_step精度）
        histogram = np.array(result['histogram'])  # 直方图计数（cdf_step精度）
        cumulative_pct = np.array(result['cumulative_percentage'])
        total = result['total_variants']
        
        # 使用存储的步长参数
        hist_step = result.get('hist_step', 0.01)  # 直方图步长（固定0.01，用于显示）
        cdf_step = result.get('cdf_step', 0.001)   # CDF步长（用于计算精度）
        
        # 关键修正：如果cdf_step != hist_step，需要重新聚合histogram数据
        # 这样才能保证所有模式的直方图视觉一致
        if abs(cdf_step - hist_step) > 1e-6:
            # cdf_step更细，需要按hist_step重新聚合
            vmiss_min = vmiss_values[0]
            vmiss_max = vmiss_values[-1]
            
            # 生成hist_step的bin边界
            hist_bins = np.arange(
                np.floor(vmiss_min / hist_step) * hist_step,
                np.ceil(vmiss_max / hist_step) * hist_step + hist_step/2,
                hist_step
            )
            
            # 将细粒度的histogram聚合到粗粒度
            hist_aggregated = []
            hist_x_values = []
            
            for i in range(len(hist_bins) - 1):
                bin_left = hist_bins[i]
                bin_right = hist_bins[i + 1]
                
                # 找到落在这个bin范围内的所有cdf数据点
                mask = (vmiss_values >= bin_left) & (vmiss_values < bin_right)
                bin_count = histogram[mask].sum()
                
                hist_aggregated.append(bin_count)
                hist_x_values.append(bin_left)
            
            # 使用聚合后的数据绘制
            hist_x_values = np.array(hist_x_values)
            hist_aggregated = np.array(hist_aggregated)
            
            logger.info(f"  [{group_key}] 重新聚合直方图: "
                       f"{len(histogram)}个点(步长{cdf_step}) → {len(hist_aggregated)}个点(步长{hist_step})")
        else:
            # cdf_step == hist_step，直接使用原始数据
            hist_x_values = vmiss_values
            hist_aggregated = histogram
        
        # 左Y轴: 绘制直方图（分布）
        # 现在histogram数据和bar宽度都基于hist_step，视觉完全一致
        color_hist = 'tab:blue'
        ax1.bar(hist_x_values, hist_aggregated, width=hist_step * 0.8,
               color=color_hist, alpha=0.6, label='Distribution')
        ax1.set_xlabel('VMISS', fontsize=12)
        ax1.set_ylabel('Variant Count', fontsize=12, color=color_hist)
        ax1.tick_params(axis='y', labelcolor=color_hist)
        
        # 格式化y轴为千位分隔
        from matplotlib.ticker import FuncFormatter
        ax1.yaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{int(x):,}'))
        
        # 右Y轴: 绘制累积分布曲线
        ax2 = ax1.twinx()
        color_cum = 'tab:orange'
        ax2.plot(vmiss_values, cumulative_pct, color=color_cum, linewidth=2.5, 
                label='Cumulative %', marker='o', markersize=4)
        ax2.set_ylabel('Cumulative Percentage (%)', fontsize=12, color=color_cum)
        ax2.tick_params(axis='y', labelcolor=color_cum)
        ax2.set_ylim([0, 105])
        
        # 标记拐点
        if group_key in knee_results and knee_results[group_key].get('knee_found'):
            knee_info = knee_results[group_key]
            knee_vmiss = knee_info['knee_vmiss']
            knee_cum_pct = knee_info['knee_cumulative_percentage']
            knee_cum_count = knee_info['knee_cumulative_count']
            
            # 在累积分布曲线上标记拐点（红色星星）
            ax2.scatter([knee_vmiss], [knee_cum_pct], c='red', marker='*', 
                       s=500, zorder=10, label='Knee Point')
            ax2.axvline(knee_vmiss, color='r', linestyle='--', linewidth=2, alpha=0.7)
            
            # 添加注释
            annotation_text = (f"VMISS threshold: {knee_vmiss:.3f}\n"
                             f"Cumulative: {knee_cum_pct:.1f}%\n"
                             f"Variant count: {knee_cum_count:,}\n"
                             f"Total variants: {total:,}")
            
            ax2.annotate(annotation_text,
                        xy=(knee_vmiss, knee_cum_pct),
                        xytext=(20, -40),
                        textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.7),
                        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0',
                                      color='red', lw=2),
                        fontsize=10,
                        ha='left')
        
        # 设置标题
        title = f"{result['maf_group']} - {result['vmiss_col']}"
        ax1.set_title(title, fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3, axis='x')
        
        # 合并图例，放在右下角
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='lower right', fontsize=10)
    
    # 隐藏多余的子图
    for idx in range(n_groups, len(axes)):
        axes[idx].set_visible(False)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    plot_filename = f'vmiss_distribution_analysis_mode{analysis_mode}.png'
    plot_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_path, dpi=dpi, bbox_inches='tight')
    
    # 关闭图形对象，释放内存
    plt.close(fig)
    plt.close('all')
    
    logger.info(f"分布分析图已保存: {plot_path}")
    return plot_path


def _print_knee_recommendations(knee_results: dict, analysis_mode: int):
    """
    打印基于累积分布拐点的VMISS阈值建议摘要
    """
    logger.info("\n" + "="*60)
    logger.info("VMISS阈值建议摘要 (基于累积分布拐点)")
    logger.info("="*60)
    
    for group_key, result in knee_results.items():
        if result.get('knee_found'):
            logger.info(f"\n组: {group_key}")
            logger.info(f"  推荐VMISS阈值: {result['knee_vmiss']:.3f}")
            logger.info(f"  累积百分比: {result['knee_cumulative_percentage']:.1f}%")
            logger.info(f"  累积变体数量: {result['knee_cumulative_count']:,}")
            logger.info(f"  总变体数量: {result['total_variants']:,}")
            logger.info(f"  解释: 在此VMISS值处，累积分布曲线出现明显拐点")
        else:
            logger.info(f"\n组: {group_key}")
            logger.info(f"  状态: {result.get('message', '未找到拐点')}")
    
    logger.info("\n" + "="*60)
