#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SNP基因型一致性分析工具
用于比较两个PLINK格式基因型文件中相同SNP的基因型一致性
"""

import subprocess
import tempfile
import os
import json
import logging
from pathlib import Path
from typing import Dict, Tuple, Optional
import numpy as np
import pandas as pd

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def compare_genotypes(
    prefix1: str,
    prefix2: str,
    snp_id: str,
    plink2_path: str = "/home/b/b37974/plink2",
    output_dir: str = "tmp",
    output_prefix: Optional[str] = None,
    missing_mode: str = "exclude_any"
) -> Dict:
    """
    比较两个PLINK格式基因型文件中指定SNP的基因型一致性
    
    参数:
        prefix1 (str): 第一个PLINK文件的前缀路径（不含.bed/.bim/.fam后缀）
        prefix2 (str): 第二个PLINK文件的前缀路径（不含.bed/.bim/.fam后缀）
        snp_id (str): 要比较的SNP ID
        plink2_path (str): plink2可执行文件的路径，默认为'/home/b/b37974/plink2'
        output_dir (str): 输出目录，默认为'tmp'
        output_prefix (str, optional): 输出文件的前缀，默认使用snp_id（冒号替换为下划线）
        missing_mode (str): 缺失值处理模式，可选值:
            - "exclude_any": 任意一个文件中有缺失就不纳入concordance计算（默认）
            - "exclude_file1": 文件1中有缺失就不纳入计算（包括两者都缺失）
            - "exclude_file2": 文件2中有缺失就不纳入计算（包括两者都缺失）
            - "include_all": 所有样本都纳入计算，包括缺失值
    
    返回:
        Dict: 包含混淆矩阵、AAF、缺失率、一致性率和共同样本数的字典
    
    示例:
        >>> result = compare_genotypes(
        ...     prefix1="/path/to/file1",
        ...     prefix2="/path/to/file2",
        ...     snp_id="chr17:13508135:G:A",
        ...     output_dir="results",
        ...     missing_mode="exclude_any"
        ... )
    """
    
    logger.info(f"=" * 80)
    logger.info(f"开始基因型一致性分析")
    logger.info(f"文件1: {prefix1}")
    logger.info(f"文件2: {prefix2}")
    logger.info(f"SNP ID: {snp_id}")
    logger.info(f"缺失值处理模式: {missing_mode}")
    logger.info(f"=" * 80)
    
    # 验证missing_mode参数
    valid_modes = ["exclude_any", "exclude_file1", "exclude_file2", "include_all"]
    if missing_mode not in valid_modes:
        raise ValueError(f"无效的missing_mode: {missing_mode}. 有效值: {valid_modes}")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"输出目录: {output_dir}")
    
    # 确定输出文件前缀
    if output_prefix is None:
        # 将冒号替换为下划线，使其成为合法的文件名
        output_prefix = snp_id.replace(':', '_')
    logger.info(f"输出文件前缀: {output_prefix}")
    
    # 检查PLINK文件是否存在
    logger.info("检查输入文件...")
    _check_plink_files(prefix1, "文件1")
    _check_plink_files(prefix2, "文件2")
    logger.info("输入文件检查通过")
    
    # 创建临时目录
    with tempfile.TemporaryDirectory() as tmpdir:
        logger.info(f"创建临时工作目录: {tmpdir}")
        
        # 提取指定SNP的基因型数据
        logger.info("提取文件1中的SNP基因型数据...")
        geno1_file = _extract_snp_genotypes(
            prefix1, snp_id, tmpdir, "file1", plink2_path
        )
        
        logger.info("提取文件2中的SNP基因型数据...")
        geno2_file = _extract_snp_genotypes(
            prefix2, snp_id, tmpdir, "file2", plink2_path
        )
        
        # 读取基因型数据
        logger.info("读取基因型数据...")
        df1 = pd.read_csv(geno1_file, sep='\t')
        df2 = pd.read_csv(geno2_file, sep='\t')
        
        # 合并数据找到共同样本
        logger.info("合并数据并找到共同样本...")
        merged_df = pd.merge(
            df1, df2,
            on=['IID'],
            how='inner',
            suffixes=('_1', '_2')
        )
        
        n_common_samples = len(merged_df)
        logger.info(f"共同样本数: {n_common_samples}")
        
        if n_common_samples == 0:
            logger.warning("警告: 两个文件中没有共同样本!")
            raise ValueError("两个基因型文件中没有共同样本")
        
        # 提取基因型列
        # PLINK2输出的列名格式: IID, 然后是SNP_ID (或者可能是SNP_ID_x, SNP_ID_y等)
        geno_cols_1 = [col for col in merged_df.columns if col.endswith('_1') and col != 'IID_1']
        geno_cols_2 = [col for col in merged_df.columns if col.endswith('_2') and col != 'IID_2']
        
        if len(geno_cols_1) == 0 or len(geno_cols_2) == 0:
            # 尝试直接使用SNP ID
            geno_col_1 = f"{snp_id}_1"
            geno_col_2 = f"{snp_id}_2"
            if geno_col_1 not in merged_df.columns or geno_col_2 not in merged_df.columns:
                # 使用第一个非IID列
                geno_col_1 = [col for col in df1.columns if col != 'IID'][0] + '_1'
                geno_col_2 = [col for col in df2.columns if col != 'IID'][0] + '_2'
        else:
            geno_col_1 = geno_cols_1[0]
            geno_col_2 = geno_cols_2[0]
        
        logger.info(f"使用基因型列: {geno_col_1} 和 {geno_col_2}")
        
        geno1 = merged_df[geno_col_1].values
        geno2 = merged_df[geno_col_2].values
        
        # 计算混淆矩阵
        logger.info("计算混淆矩阵...")
        confusion_matrix = _calculate_confusion_matrix(geno1, geno2)
        logger.info("混淆矩阵计算完成")
        
        # 计算AAF和缺失率
        logger.info("计算等位基因频率(AAF)和缺失率...")
        aaf1, miss_rate1 = _calculate_aaf_and_miss_rate(geno1)
        aaf2, miss_rate2 = _calculate_aaf_and_miss_rate(geno2)
        logger.info(f"文件1 - AAF: {aaf1:.4f}, 缺失率: {miss_rate1:.4f}")
        logger.info(f"文件2 - AAF: {aaf2:.4f}, 缺失率: {miss_rate2:.4f}")
        
        # 计算一致性率
        logger.info("计算基因型一致性率...")
        concordance_rate, n_concordant, n_compared = _calculate_concordance_rate(geno1, geno2, missing_mode)
        logger.info(f"一致性率: {concordance_rate:.4f} ({concordance_rate*100:.2f}%)")
        
        # 组装结果
        result = {
            "snp_id": snp_id,
            "file1_prefix": prefix1,
            "file2_prefix": prefix2,
            "n_common_samples": int(n_common_samples),
            "missing_mode": missing_mode,
            "confusion_matrix": confusion_matrix,
            "file1_stats": {
                "aaf": float(aaf1),
                "miss_rate": float(miss_rate1)
            },
            "file2_stats": {
                "aaf": float(aaf2),
                "miss_rate": float(miss_rate2)
            },
            "concordance_rate": float(concordance_rate),
            "n_concordant": int(n_concordant),
            "n_compared": int(n_compared)
        }
        
        # 保存JSON结果
        json_file = os.path.join(output_dir, f"{output_prefix}.json")
        logger.info(f"保存结果到JSON文件: {json_file}")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(result, f, indent=2, ensure_ascii=False)
        logger.info("JSON文件保存成功")
        
        # 保存混淆矩阵为tab分隔的CSV文件
        confusion_matrix_file = os.path.join(output_dir, f"{output_prefix}_confusion_matrix.csv")
        logger.info(f"保存混淆矩阵到CSV文件: {confusion_matrix_file}")
        _save_confusion_matrix(confusion_matrix, confusion_matrix_file)
        logger.info("混淆矩阵CSV文件保存成功")
        
        logger.info("=" * 80)
        logger.info("分析完成!")
        logger.info(f"输出文件:")
        logger.info(f"  - JSON: {json_file}")
        logger.info(f"  - 混淆矩阵CSV: {confusion_matrix_file}")
        logger.info("=" * 80)
        
        return result


def _check_plink_files(prefix: str, label: str) -> None:
    """检查PLINK文件是否存在"""
    required_extensions = ['.bed', '.bim', '.fam']
    for ext in required_extensions:
        filepath = prefix + ext
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"{label}缺少文件: {filepath}")


def _extract_snp_genotypes(
    prefix: str,
    snp_id: str,
    tmpdir: str,
    label: str,
    plink2_path: str
) -> str:
    """
    使用PLINK2提取指定SNP的基因型数据
    
    从SNP ID (格式: CHROM:POS:REF:ALT) 中解析ALT等位基因，
    并使用--alt1-allele强制指定ALT作为counted allele
    
    返回: 输出文件的路径
    """
    output_prefix = os.path.join(tmpdir, label)
    
    # 从SNP ID中解析ALT等位基因
    # 格式: CHROM:POS:REF:ALT
    alt_allele = None
    ref_allele = None
    if ':' in snp_id:
        parts = snp_id.split(':')
        if len(parts) == 4:
            ref_allele = parts[2]
            alt_allele = parts[3]
            logger.info(f"{label}: 从SNP ID解析 - REF={ref_allele}, ALT={alt_allele}")
    
    # 构建PLINK2命令
    cmd = [
        plink2_path,
        '--bfile', prefix,
        '--snp', snp_id,
        '--export', 'A-transpose',  # 导出转置的加性编码格式
        '--out', output_prefix
    ]
    
    # 如果成功解析了ALT等位基因，创建一个临时文件来指定ALT1-allele
    alt1_allele_file = None
    if alt_allele:
        # 创建临时文件指定ALT等位基因
        alt1_allele_file = os.path.join(tmpdir, f'{label}_alt1.txt')
        with open(alt1_allele_file, 'w') as f:
            # 格式: SNP_ID ALT_ALLELE
            f.write(f"{snp_id}\t{alt_allele}\n")
        
        # 添加--alt1-allele参数强制使用ALT作为counted allele
        cmd.extend(['--alt1-allele', alt1_allele_file])
        logger.info(f"{label}: 使用--alt1-allele强制ALT等位基因({alt_allele})作为counted allele")
    else:
        logger.warning(f"{label}: 无法从SNP ID解析ALT等位基因，将使用PLINK默认行为")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        logger.debug(f"PLINK2命令执行成功: {' '.join(cmd)}")
    except subprocess.CalledProcessError as e:
        logger.error(f"PLINK2命令执行失败!")
        logger.error(f"命令: {' '.join(cmd)}")
        logger.error(f"返回码: {e.returncode}")
        logger.error(f"标准输出: {e.stdout}")
        logger.error(f"标准错误: {e.stderr}")
        raise RuntimeError(f"PLINK2提取SNP失败: {snp_id}")
    
    output_file = output_prefix + '.traw'
    if not os.path.exists(output_file):
        raise FileNotFoundError(f"PLINK2未生成预期的输出文件: {output_file}")
    
    # 转换traw格式为更易处理的格式
    df = pd.read_csv(output_file, sep='\t')
    
    # traw格式: CHR SNP (cM) POS COUNTED ALT Sample1 Sample2 ...
    # 第5列(COUNTED)是被计数的等位基因，第6列(ALT)是替代等位基因
    sample_cols = df.columns[6:]  # 前6列是变异信息
    genotypes = df.iloc[0, 6:].values  # 第一行（也是唯一一行）的基因型
    
    # 检查COUNTED列和ALT列
    counted_allele = df.iloc[0, 4]  # COUNTED列
    alt_allele_from_file = df.iloc[0, 5]  # ALT列
    
    logger.info(f"{label}: PLINK输出 - COUNTED等位基因={counted_allele}, ALT等位基因={alt_allele_from_file}")
    
    # 创建新的DataFrame
    result_df = pd.DataFrame({
        'IID': sample_cols,
        snp_id: genotypes
    })
    
    # 将NA标记为NaN
    result_df[snp_id] = result_df[snp_id].replace('NA', np.nan)
    result_df[snp_id] = pd.to_numeric(result_df[snp_id], errors='coerce')
    
    # 验证：如果我们指定了--alt1-allele，COUNTED应该等于ALT
    if alt_allele and counted_allele != alt_allele_from_file:
        logger.warning(f"{label}: 警告！即使使用了--alt1-allele，COUNTED({counted_allele})仍不等于ALT({alt_allele_from_file})")
        logger.warning(f"{label}: 将反转基因型编码以确保计算ALT频率")
        # 反转基因型编码（0->2, 2->0, 1保持不变）
        result_df[snp_id] = result_df[snp_id].apply(
            lambda x: 2 - x if pd.notna(x) else x
        )
    elif counted_allele == alt_allele_from_file:
        logger.info(f"{label}: 验证通过 - COUNTED等位基因已正确设置为ALT")
    else:
        # 没有指定alt_allele的情况，检查是否需要反转
        if counted_allele != alt_allele_from_file:
            logger.info(f"{label}: COUNTED等位基因与ALT不同，反转基因型编码以确保计算ALT频率")
            result_df[snp_id] = result_df[snp_id].apply(
                lambda x: 2 - x if pd.notna(x) else x
            )
    
    # 保存处理后的文件
    processed_file = output_prefix + '_processed.txt'
    result_df.to_csv(processed_file, sep='\t', index=False)
    
    # 同时保存counted和alt等位基因信息
    allele_info_file = output_prefix + '_allele_info.txt'
    with open(allele_info_file, 'w') as f:
        f.write(f"COUNTED\t{counted_allele}\n")
        f.write(f"ALT\t{alt_allele_from_file}\n")
        if ref_allele:
            f.write(f"REF_from_ID\t{ref_allele}\n")
        if alt_allele:
            f.write(f"ALT_from_ID\t{alt_allele}\n")
    
    return processed_file


def _calculate_confusion_matrix(geno1: np.ndarray, geno2: np.ndarray) -> Dict:
    """
    计算混淆矩阵
    
    返回格式:
    {
        "GT_0": {"GT_0": count, "GT_1": count, "GT_2": count, "GT_missing": count},
        "GT_1": {"GT_0": count, "GT_1": count, "GT_2": count, "GT_missing": count},
        "GT_2": {"GT_0": count, "GT_1": count, "GT_2": count, "GT_missing": count},
        "GT_missing": {"GT_0": count, "GT_1": count, "GT_2": count, "GT_missing": count}
    }
    """
    # 定义基因型类别
    categories = [0, 1, 2, 'missing']
    
    # 初始化混淆矩阵
    confusion_matrix = {
        f"file1_GT_{cat}": {
            f"file2_GT_{cat2}": 0 
            for cat2 in categories
        }
        for cat in categories
    }
    
    # 填充混淆矩阵
    for g1, g2 in zip(geno1, geno2):
        # 确定基因型类别
        cat1 = 'missing' if pd.isna(g1) else int(g1)
        cat2 = 'missing' if pd.isna(g2) else int(g2)
        
        confusion_matrix[f"file1_GT_{cat1}"][f"file2_GT_{cat2}"] += 1
    
    return confusion_matrix


def _calculate_aaf_and_miss_rate(genotypes: np.ndarray) -> Tuple[float, float]:
    """
    计算等位基因频率(AAF)和缺失率
    
    注意: 基因型编码必须已经被标准化为以ALT等位基因为计数对象
    - GT=0: REF/REF (0个ALT等位基因)
    - GT=1: REF/ALT (1个ALT等位基因)
    - GT=2: ALT/ALT (2个ALT等位基因)
    
    返回: (aaf, miss_rate)
    """
    # 计算缺失率
    n_total = len(genotypes)
    n_missing = np.sum(pd.isna(genotypes))
    miss_rate = n_missing / n_total if n_total > 0 else 0.0
    
    # 计算AAF (Alternative Allele Frequency) - 针对ALT等位基因
    # GT=0 (REF/REF): 贡献 0 个 ALT 等位基因
    # GT=1 (REF/ALT): 贡献 1 个 ALT 等位基因
    # GT=2 (ALT/ALT): 贡献 2 个 ALT 等位基因
    valid_genotypes = genotypes[~pd.isna(genotypes)]
    n_valid = len(valid_genotypes)
    
    if n_valid == 0:
        aaf = 0.0
    else:
        # 等位基因总数 = 有效样本数 * 2
        n_alleles = n_valid * 2
        # ALT等位基因数 = GT=1的数量*1 + GT=2的数量*2
        n_alt_alleles = np.sum(valid_genotypes == 1) + np.sum(valid_genotypes == 2) * 2
        aaf = n_alt_alleles / n_alleles if n_alleles > 0 else 0.0
    
    return aaf, miss_rate


def _calculate_concordance_rate(geno1: np.ndarray, geno2: np.ndarray, missing_mode: str = "exclude_any") -> Tuple[float, int, int]:
    """
    计算基因型一致性率
    
    参数:
        geno1: 文件1的基因型数组
        geno2: 文件2的基因型数组
        missing_mode: 缺失值处理模式
            - "exclude_any": 任意一个文件中有缺失就不纳入concordance计算
            - "exclude_file1": 文件1中有缺失就不纳入计算
            - "exclude_file2": 文件2中有缺失就不纳入计算
            - "include_all": 所有样本都纳入计算，缺失值也算作一种基因型
    
    返回:
        Tuple[float, int, int]: (一致性率, 匹配数, 用于计算的样本数)
    """
    n_total = len(geno1)
    
    if missing_mode == "exclude_any":
        # 任意一个文件中有缺失就不纳入计算
        # 只在两个文件都有有效基因型的样本中计算
        both_valid = ~pd.isna(geno1) & ~pd.isna(geno2)
        valid_geno1 = geno1[both_valid]
        valid_geno2 = geno2[both_valid]
        n_valid = len(valid_geno1)
        
        if n_valid == 0:
            logger.warning("警告: 没有样本在两个文件中都有有效基因型")
            return 0.0, 0, 0
        
        n_match = np.sum(valid_geno1 == valid_geno2)
        concordance_rate = n_match / n_valid  # 分母应该是有效样本数
        
        logger.info(f"  排除模式: 任意缺失")
        logger.info(f"  有效样本数: {n_valid}/{n_total}")
        logger.info(f"  匹配样本数: {n_match}")
        logger.info(f"  一致性率: {n_match}/{n_valid} = {concordance_rate:.4f}")
        
        return concordance_rate, n_match, n_valid
        
    elif missing_mode == "exclude_file1":
        # 文件1中有缺失就不纳入计算（包括两者都缺失）
        # 只在文件1有有效基因型的样本中计算
        # 注意：file2可以是缺失，如果file2缺失则算作不匹配
        file1_valid = ~pd.isna(geno1)
        valid_geno1 = geno1[file1_valid]
        valid_geno2 = geno2[file1_valid]
        n_valid = len(valid_geno1)
        
        if n_valid == 0:
            logger.warning("警告: 文件1中没有有效基因型")
            return 0.0, 0, 0
        
        # 在file1有效的样本中计算匹配数
        # 只有当两者都有效且相等时才算匹配
        # 如果file2是缺失，算作不匹配（纳入分母，不纳入分子）
        n_match = 0
        for g1, g2 in zip(valid_geno1, valid_geno2):
            if pd.notna(g2) and g1 == g2:
                n_match += 1
        
        concordance_rate = n_match / n_valid  # 分母是file1有效的样本数
        
        logger.info(f"  排除模式: 文件1缺失")
        logger.info(f"  文件1有效样本数: {n_valid}/{n_total}")
        logger.info(f"  匹配样本数: {n_match} (file2缺失的情况算作不匹配)")
        logger.info(f"  一致性率: {n_match}/{n_valid} = {concordance_rate:.4f}")
        
        return concordance_rate, n_match, n_valid
        
    elif missing_mode == "exclude_file2":
        # 文件2中有缺失就不纳入计算（包括两者都缺失）
        # 只在文件2有有效基因型的样本中计算
        # 注意：file1可以是缺失，如果file1缺失则算作不匹配
        file2_valid = ~pd.isna(geno2)
        valid_geno1 = geno1[file2_valid]
        valid_geno2 = geno2[file2_valid]
        n_valid = len(valid_geno2)
        
        if n_valid == 0:
            logger.warning("警告: 文件2中没有有效基因型")
            return 0.0, 0, 0
        
        # 在file2有效的样本中计算匹配数
        # 只有当两者都有效且相等时才算匹配
        # 如果file1是缺失，算作不匹配（纳入分母，不纳入分子）
        n_match = 0
        for g1, g2 in zip(valid_geno1, valid_geno2):
            if pd.notna(g1) and g1 == g2:
                n_match += 1
        
        concordance_rate = n_match / n_valid  # 分母是file2有效的样本数
        
        logger.info(f"  排除模式: 文件2缺失")
        logger.info(f"  文件2有效样本数: {n_valid}/{n_total}")
        logger.info(f"  匹配样本数: {n_match} (file1缺失的情况算作不匹配)")
        logger.info(f"  一致性率: {n_match}/{n_valid} = {concordance_rate:.4f}")
        
        return concordance_rate, n_match, n_valid
        
    elif missing_mode == "include_all":
        # 所有样本都纳入计算，缺失值也算作一种基因型
        # 注意：两个都缺失不算作匹配
        # 只有当两个都有值且相等时才算匹配
        n_match = 0
        for g1, g2 in zip(geno1, geno2):
            if pd.notna(g1) and pd.notna(g2) and g1 == g2:
                # 只有两个都有值且相等才算匹配
                n_match += 1
            # 其他所有情况都是不匹配（包括两个都缺失、一个缺失一个不缺失、都有值但不相等）
        
        concordance_rate = n_match / n_total  # 分母是总样本数
        
        logger.info(f"  排除模式: 包含所有（缺失也算）")
        logger.info(f"  总样本数: {n_total}")
        logger.info(f"  匹配样本数: {n_match} (两个都缺失不算匹配)")
        logger.info(f"  一致性率: {n_match}/{n_total} = {concordance_rate:.4f}")
        
        return concordance_rate, n_match, n_total
    
    else:
        raise ValueError(f"无效的missing_mode: {missing_mode}")


def _save_confusion_matrix(confusion_matrix: Dict, output_file: str) -> None:
    """
    将混淆矩阵保存为tab分隔的CSV格式文件
    
    参数:
        confusion_matrix: 混淆矩阵字典
        output_file: 输出文件路径（.csv格式）
    """
    categories = [0, 1, 2, 'missing']
    
    # 准备数据用于DataFrame
    data = []
    for cat1 in categories:
        row_key = f"file1_GT_{cat1}"
        row = [f"GT_{cat1}"]  # 第一列是行标签
        for cat2 in categories:
            col_key = f"file2_GT_{cat2}"
            count = confusion_matrix[row_key][col_key]
            row.append(count)
        data.append(row)
    
    # 创建DataFrame
    columns = ['File1_File2'] + [f"GT_{cat}" for cat in categories]
    df = pd.DataFrame(data, columns=columns)
    
    # 保存为tab分隔的CSV文件
    df.to_csv(output_file, sep='\t', index=False, encoding='utf-8')



def print_result_summary(result: Dict) -> None:
    """
    打印结果摘要（格式化输出）
    
    参数:
        result: compare_genotypes函数返回的结果字典
    """
    print("\n" + "=" * 80)
    print("基因型一致性分析结果摘要")
    print("=" * 80)
    print(f"SNP ID: {result['snp_id']}")
    print(f"共同样本数: {result['n_common_samples']}")
    print(f"\n文件1统计:")
    print(f"  - 等位基因频率(AAF): {result['file1_stats']['aaf']:.4f}")
    print(f"  - 缺失率: {result['file1_stats']['miss_rate']:.4f}")
    print(f"\n文件2统计:")
    print(f"  - 等位基因频率(AAF): {result['file2_stats']['aaf']:.4f}")
    print(f"  - 缺失率: {result['file2_stats']['miss_rate']:.4f}")
    print(f"\n一致性统计:")
    print(f"  - 一致性率: {result['concordance_rate']:.4f} ({result['concordance_rate']*100:.2f}%)")
    print(f"  - 匹配样本数(分子): {result['n_concordant']}")
    print(f"  - 比较样本数(分母): {result['n_compared']}")
    print(f"  - 缺失值处理模式: {result['missing_mode']}")
    print(f"\n混淆矩阵:")
    print("-" * 80)
    
    # 打印混淆矩阵表格
    cm = result['confusion_matrix']
    categories = [0, 1, 2, 'missing']
    
    # 表头
    print(f"{'文件1 \\ 文件2':<20}", end="")
    for cat in categories:
        print(f"GT_{cat}".rjust(12), end="")
    print()
    print("-" * 80)
    
    # 表内容
    for cat1 in categories:
        row_key = f"file1_GT_{cat1}"
        print(f"GT_{cat1}".ljust(20), end="")
        for cat2 in categories:
            col_key = f"file2_GT_{cat2}"
            count = cm[row_key][col_key]
            print(f"{count}".rjust(12), end="")
        print()
    
    print("=" * 80 + "\n")


if __name__ == "__main__":
    # 示例用法
    import sys
    
    if len(sys.argv) < 4:
        print("用法: python snp_concordance_tools.py <prefix1> <prefix2> <snp_id> [output_dir] [output_prefix] [missing_mode]")
        print("示例: python snp_concordance_tools.py /path/to/file1 /path/to/file2 chr17:13508135:G:A tmp")
        print("      python snp_concordance_tools.py /path/to/file1 /path/to/file2 chr17:13508135:G:A results my_snp exclude_any")
        print("\nmissing_mode选项:")
        print("  exclude_any (默认): 任意一个文件中有缺失就不纳入concordance计算")
        print("  exclude_file1: 文件1中有缺失就不纳入计算")
        print("  exclude_file2: 文件2中有缺失就不纳入计算")
        print("  include_all: 所有样本都纳入计算，缺失值也算作一种基因型")
        sys.exit(1)
    
    prefix1 = sys.argv[1]
    prefix2 = sys.argv[2]
    snp_id = sys.argv[3]
    output_dir = sys.argv[4] if len(sys.argv) > 4 else "tmp"
    output_prefix = sys.argv[5] if len(sys.argv) > 5 else None
    missing_mode = sys.argv[6] if len(sys.argv) > 6 else "exclude_any"
    
    try:
        result = compare_genotypes(
            prefix1, 
            prefix2, 
            snp_id, 
            output_dir=output_dir,
            output_prefix=output_prefix,
            missing_mode=missing_mode
        )
        print_result_summary(result)
    except Exception as e:
        logger.error(f"错误: {str(e)}")
        sys.exit(1)

