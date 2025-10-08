#!/usr/bin/env python3
"""
RVTest 数据预处理主程序

本脚本是一个专业的命令行工具，用于自动化完成 RVTest（稀有变异关联分析）所需的全部数据预处理步骤。
主要功能包括：

1. 基因型数据转换：将 PLINK 格式(.bed/.bim/.fam)转换为标准化的 VCF 格式
2. 表型协变量重格式化：处理表型和协变量文件以适配 RVTest 输入要求
3. 基因注释文件处理：去除 refFlat 文件中的染色体前缀，确保格式一致性

脚本设计原则：
- 提供完整的命令行参数支持
- 详细的中文日志记录和错误处理
- 模块化设计，便于维护和扩展
- 支持并行处理以提高运行效率

作者: ZHAO TIE
创建时间: 2025年10月8日
版本: 1.0.0

使用示例:
    python rvtest_prepare_main.py \\
        --bed-prefix /path/to/plink_data \\
        --pheno-path /path/to/phenotype.csv \\
        --covar-path /path/to/covariate.csv \\
        --refflat-path /path/to/refFlat.txt.gz \\
        --threads 8
"""

import sys
import os
import argparse
import logging
from datetime import datetime
from typing import Tuple

# 添加模块路径
MODULE_DIR = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest/scripts"
if MODULE_DIR not in sys.path:
    sys.path.insert(0, MODULE_DIR)

import rvtest_prepare_tools


def setup_logging(verbose: bool = False) -> None:
    """
    配置日志记录系统
    
    参数
    ----
    verbose : bool
        是否启用详细日志模式（DEBUG级别）
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    
    logging.basicConfig(
        level=log_level,
        format=log_format,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(f'rvtest_prepare_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
        ]
    )


def parse_arguments() -> argparse.Namespace:
    """
    解析命令行参数
    
    返回
    ----
    argparse.Namespace
        解析后的命令行参数对象
    """
    parser = argparse.ArgumentParser(
        description='RVTest 数据预处理工具 - 自动化处理基因型、表型和注释数据',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本用法 - 处理所有必需文件
  %(prog)s --bed-prefix /path/to/data \\
           --pheno-path /path/to/pheno.csv \\
           --covar-path /path/to/covar.csv \\
           --refflat-path /path/to/refFlat.txt.gz

  # 高级用法 - 自定义线程数和输出格式
  %(prog)s --bed-prefix /path/to/data \\
           --pheno-path /path/to/pheno.csv \\
           --covar-path /path/to/covar.csv \\
           --refflat-path /path/to/refFlat.txt.gz \\
           --threads 16 \\
           --keep-chr-prefix \\
           --snps-only \\
           --verbose

注意事项:
  - 所有输入文件路径必须是绝对路径
  - 确保有足够的磁盘空间用于存储中间文件和最终输出
  - 建议在处理大型数据集时使用更多线程数
        """
    )
    
    # 必需参数组
    required_group = parser.add_argument_group('必需参数')
    required_group.add_argument(
        '--bed-prefix',
        type=str,
        required=True,
        help='PLINK 格式数据文件前缀（不包含.bed/.bim/.fam后缀）'
    )
    required_group.add_argument(
        '--pheno-path',
        type=str,
        required=True,
        help='表型数据文件路径（CSV 格式）'
    )
    required_group.add_argument(
        '--covar-path',
        type=str,
        required=True,
        help='协变量数据文件路径（CSV 格式）'
    )
    required_group.add_argument(
        '--refflat-path',
        type=str,
        required=True,
        help='refFlat 基因注释文件路径（.txt.gz 格式）'
    )
    
    # 可选参数组
    optional_group = parser.add_argument_group('可选参数')
    optional_group.add_argument(
        '--threads',
        type=int,
        default=8,
        help='并行处理线程数（默认: 8）'
    )
    optional_group.add_argument(
        '--snps-only',
        action='store_true',
        help='仅处理 A/C/G/T 标准 SNP 变异（默认: 处理所有变异类型）'
    )
    optional_group.add_argument(
        '--keep-chr-prefix',
        action='store_true',
        help='在 VCF 文件中保留染色体前缀 "chr"（默认: 去除前缀）'
    )
    optional_group.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='输出目录路径（默认: 当前工作目录）'
    )
    optional_group.add_argument(
        '--verbose',
        action='store_true',
        help='启用详细日志输出模式'
    )
    
    return parser.parse_args()


def validate_input_files(args: argparse.Namespace) -> None:
    """
    验证输入文件的存在性和可读性
    
    参数
    ----
    args : argparse.Namespace
        命令行参数对象
        
    抛出
    ----
    FileNotFoundError
        当必需的输入文件不存在时
    PermissionError
        当文件权限不足时
    """
    logging.info("正在验证输入文件...")
    
    # 检查 PLINK 文件
    plink_extensions = ['.bed', '.bim', '.fam']
    for ext in plink_extensions:
        file_path = args.bed_prefix + ext
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"PLINK 文件不存在: {file_path}")
        if not os.access(file_path, os.R_OK):
            raise PermissionError(f"无法读取文件: {file_path}")
    
    # 检查表型和协变量文件
    for file_path in [args.pheno_path, args.covar_path, args.refflat_path]:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"输入文件不存在: {file_path}")
        if not os.access(file_path, os.R_OK):
            raise PermissionError(f"无法读取文件: {file_path}")
    
    logging.info("✓ 所有输入文件验证通过")


def process_plink_to_vcf(args: argparse.Namespace) -> str:
    """
    执行 PLINK 到 VCF 的转换处理
    
    参数
    ----
    args : argparse.Namespace
        命令行参数对象
        
    返回
    ----
    str
        生成的 VCF 文件路径
    """
    logging.info("="*60)
    logging.info("步骤 1/3: 开始 PLINK 到 VCF 格式转换")
    logging.info("="*60)
    logging.info(f"输入文件前缀: {args.bed_prefix}")
    logging.info(f"使用线程数: {args.threads}")
    logging.info(f"仅处理 SNP: {'是' if args.snps_only else '否'}")
    logging.info(f"保留 chr 前缀: {'是' if args.keep_chr_prefix else '否'}")
    
    try:
        vcf_path = rvtest_prepare_tools.plink_to_vcf_raw(
            bed_prefix=args.bed_prefix,
            threads=args.threads,
            snps_only_just_acgt=args.snps_only,
            keep_chr_prefix=args.keep_chr_prefix,
        )
        logging.info(f"✓ VCF 转换完成，输出文件: {vcf_path}")
        return vcf_path
    except Exception as e:
        logging.error(f"✗ VCF 转换失败: {str(e)}")
        raise


def process_pheno_covar(args: argparse.Namespace) -> Tuple[str, str]:
    """
    执行表型和协变量文件的重格式化处理
    
    参数
    ----
    args : argparse.Namespace
        命令行参数对象
        
    返回
    ----
    Tuple[str, str]
        (新表型文件路径, 新协变量文件路径)
    """
    logging.info("="*60)
    logging.info("步骤 2/3: 开始表型和协变量文件重格式化")
    logging.info("="*60)
    logging.info(f"表型文件: {args.pheno_path}")
    logging.info(f"协变量文件: {args.covar_path}")
    
    try:
        new_pheno, new_covar = rvtest_prepare_tools.reformat_pheno_covar(
            pheno_path=args.pheno_path,
            covar_path=args.covar_path,
            out_dir=args.output_dir,
            out_sep="\t",  # RVTest 推荐使用制表符分隔
            force=True,
        )
        logging.info(f"✓ 表型文件重格式化完成: {new_pheno}")
        logging.info(f"✓ 协变量文件重格式化完成: {new_covar}")
        return new_pheno, new_covar
    except Exception as e:
        logging.error(f"✗ 表型/协变量文件处理失败: {str(e)}")
        raise


def process_refflat(args: argparse.Namespace) -> str:
    """
    执行 refFlat 注释文件的处理
    
    参数
    ----
    args : argparse.Namespace
        命令行参数对象
        
    返回
    ----
    str
        处理后的 refFlat 文件路径
    """
    logging.info("="*60)
    logging.info("步骤 3/3: 开始 refFlat 注释文件处理")
    logging.info("="*60)
    logging.info(f"输入文件: {args.refflat_path}")
    
    try:
        processed_refflat = rvtest_prepare_tools.reformat_refflat_remove_chr(
            refflat_path=args.refflat_path,
            out_dir=args.output_dir
        )
        logging.info(f"✓ refFlat 文件处理完成: {processed_refflat}")
        return processed_refflat
    except Exception as e:
        logging.error(f"✗ refFlat 文件处理失败: {str(e)}")
        raise


def print_summary(vcf_path: str, pheno_path: str, covar_path: str, refflat_path: str) -> None:
    """
    打印处理结果摘要
    
    参数
    ----
    vcf_path : str
        VCF 文件路径
    pheno_path : str
        表型文件路径
    covar_path : str
        协变量文件路径
    refflat_path : str
        refFlat 文件路径
    """
    logging.info("="*60)
    logging.info("RVTest 数据预处理完成！")
    logging.info("="*60)
    logging.info("处理结果摘要:")
    logging.info(f"  基因型数据 (VCF):     {vcf_path}")
    logging.info(f"  表型数据:             {pheno_path}")
    logging.info(f"  协变量数据:           {covar_path}")
    logging.info(f"  基因注释 (refFlat):   {refflat_path}")
    logging.info("")
    logging.info("接下来可以使用这些文件进行 RVTest 稀有变异关联分析。")
    logging.info("建议查看各输出文件以确认格式正确性。")


def main() -> None:
    """
    主程序入口点
    
    执行完整的 RVTest 数据预处理流程：
    1. 解析命令行参数
    2. 验证输入文件
    3. 转换 PLINK 到 VCF
    4. 重格式化表型和协变量
    5. 处理 refFlat 注释文件
    6. 输出处理结果摘要
    """
    # 解析命令行参数
    args = parse_arguments()
    
    # 设置日志记录
    setup_logging(args.verbose)
    
    # 打印程序开始信息
    logging.info("="*60)
    logging.info("RVTest 数据预处理工具启动")
    logging.info(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info("="*60)
    
    try:
        # 验证输入文件
        validate_input_files(args)
        
        # 步骤 1: PLINK 到 VCF 转换
        vcf_path = process_plink_to_vcf(args)
        
        # 步骤 2: 表型和协变量重格式化
        pheno_path, covar_path = process_pheno_covar(args)
        
        # 步骤 3: refFlat 文件处理
        refflat_path = process_refflat(args)
        
        # 打印结果摘要
        print_summary(vcf_path, pheno_path, covar_path, refflat_path)
        
        logging.info("="*60)
        logging.info("程序执行成功完成！")
        logging.info(f"结束时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logging.info("="*60)
        
    except Exception as e:
        logging.error("="*60)
        logging.error("程序执行过程中发生错误!")
        logging.error(f"错误信息: {str(e)}")
        logging.error("="*60)
        sys.exit(1)


if __name__ == "__main__":
    main()


