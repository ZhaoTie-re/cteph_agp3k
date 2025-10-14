#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VCF INFO字段过滤命令行工具
========================

基于 info_filter_tools 模块的命令行接口，用于根据INFO字段过滤VCF文件。

典型用法:
    python info_filter_main.py -i input.vcf.gz -k impact -v HIGH MODERATE -o filtered
    python info_filter_main.py -i input.vcf.gz -k impact -v HIGH -t 32 --check-chr-prefix

作者: ZHAO TIE  
日期: 2025-10-14
"""

import sys
import os
import importlib
import argparse
import logging
from pathlib import Path

# 动态加载自定义模块 (保持reload模式)
MODULE_DIR = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest/scripts"
if MODULE_DIR not in sys.path:
    sys.path.insert(0, MODULE_DIR)

import info_filter_tools
importlib.reload(info_filter_tools)  # 保持reload模式

from info_filter_tools import filter_vcf_by_info


def setup_logging(verbose: bool = False) -> None:
    """设置日志配置"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def parse_arguments() -> argparse.Namespace:
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='VCF INFO字段过滤工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  %(prog)s -i input.vcf.gz -k impact -v HIGH MODERATE -o filtered
  %(prog)s -i input.vcf.gz -k impact -v HIGH -t 32 --check-chr-prefix
  %(prog)s -i input.vcf.gz -k consequence -v "missense_variant" --no-keep-chr-prefix

注意事项:
  - 输入文件必须是bgzip压缩的VCF格式 (.vcf.gz)
  - values参数支持多个值，用空格分隔
  - 自动创建tabix索引文件
        """
    )
    
    # 必需参数
    parser.add_argument(
        '-i', '--input', '--anno-vcf',
        required=True,
        type=str,
        help='输入的注释VCF文件路径 (必须是.vcf.gz格式)',
        metavar='FILE'
    )
    
    parser.add_argument(
        '-k', '--info-key',
        required=True,
        type=str,
        help='要过滤的INFO字段名 (例如: impact, consequence)',
        metavar='KEY'
    )
    
    parser.add_argument(
        '-v', '--values',
        required=True,
        nargs='+',
        type=str,
        help='INFO字段的目标值 (支持多个值，例如: HIGH MODERATE)',
        metavar='VALUE'
    )
    
    # 可选参数
    parser.add_argument(
        '-o', '--out-prefix',
        type=str,
        default='filtered',
        help='输出文件前缀 (默认: filtered)',
        metavar='PREFIX'
    )
    
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=16,
        help='并行线程数 (默认: 16)',
        metavar='N'
    )
    
    parser.add_argument(
        '--check-chr-prefix',
        action='store_true',
        help='检查染色体前缀 (默认: False)'
    )
    
    parser.add_argument(
        '--keep-chr-prefix',
        action='store_true',
        help='保持染色体前缀 (默认: False)'
    )
    
    parser.add_argument(
        '--no-keep-chr-prefix',
        action='store_true',
        help='不保持染色体前缀 (与--keep-chr-prefix相反)'
    )
    
    # 程序选项
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='详细输出 (调试模式)'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='仅显示将要执行的命令，不实际运行'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 1.0.0'
    )
    
    return parser.parse_args()


def validate_arguments(args: argparse.Namespace) -> None:
    """验证命令行参数"""
    # 检查输入文件
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"输入文件不存在: {args.input}")
    
    if not args.input.endswith('.vcf.gz'):
        raise ValueError("输入文件必须是bgzip压缩的VCF格式 (.vcf.gz)")
    
    # 检查线程数
    if args.threads <= 0:
        raise ValueError(f"线程数必须大于0: {args.threads}")
    
    # 处理染色体前缀冲突
    if args.keep_chr_prefix and args.no_keep_chr_prefix:
        raise ValueError("--keep-chr-prefix 和 --no-keep-chr-prefix 不能同时使用")


def main() -> None:
    """主函数"""
    args = None
    try:
        # 解析参数
        args = parse_arguments()
        
        # 设置日志
        setup_logging(args.verbose)
        
        # 验证参数
        validate_arguments(args)
        
        # 确定keep_chr_prefix的值
        if args.no_keep_chr_prefix:
            keep_chr_prefix = False
        else:
            keep_chr_prefix = args.keep_chr_prefix
        
        logging.info(f"开始处理VCF文件: {args.input}")
        logging.info(f"INFO字段: {args.info_key}")
        logging.info(f"过滤值: {args.values}")
        logging.info(f"线程数: {args.threads}")
        logging.info(f"输出前缀: {args.out_prefix}")
        
        if args.dry_run:
            logging.info("DRY RUN模式 - 仅显示参数，不执行实际过滤")
            return
        
        # 执行过滤
        filtered_vcf = filter_vcf_by_info(
            vcf_path=args.input,
            info_key=args.info_key,
            values=args.values,
            threads=args.threads,
            check_chr_prefix=args.check_chr_prefix,
            keep_chr_prefix=keep_chr_prefix,
            out_prefix=args.out_prefix,
        )
        
        logging.info(f"过滤完成！输出文件: {filtered_vcf}")
        
        # 检查输出文件
        if os.path.exists(filtered_vcf):
            file_size = os.path.getsize(filtered_vcf)
            logging.info(f"输出文件大小: {file_size / (1024**2):.2f} MB")
        
    except KeyboardInterrupt:
        logging.error("用户中断操作")
        sys.exit(130)
    except Exception as e:
        logging.error(f"程序执行出错: {e}")
        if args and hasattr(args, 'verbose') and args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()


