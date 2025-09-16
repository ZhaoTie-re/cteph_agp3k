#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
基因型缺失偏差校正主程序 - FDR校正工具

本程序用于对基因型缺失偏差分析结果进行假阳性发现率(FDR)校正。
主要功能是读取合并的JSON结果文件，对其中的统计摘要数据进行
Benjamini-Hochberg多重检验校正，输出带有FDR校正p值的新文件。

作者: ZHAO TIE
创建日期: 2025年
更新日期: 2025年

主要功能:
- 读取基因型缺失偏差分析的合并JSON结果
- 对三个关键p值列进行独立的BH/FDR校正
- 生成包含校正后q值的新TSV文件
- 提供详细的运行日志和错误处理

使用方法:
    python geno_miss_bias_cor_main.py -i <input_json> [-o <output_tsv>] [-v]
"""

import argparse
import sys
import os
import logging
from pathlib import Path
from typing import Optional

# 添加脚本目录到Python路径，确保能正确导入自定义模块
SCRIPT_DIR = Path(__file__).parent.absolute()
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

try:
    from geno_miss_bias_tools import adjust_stat_summary_fdr
except ImportError as e:
    print(f"错误：无法导入必需的模块 'geno_miss_bias_tools': {e}", file=sys.stderr)
    print("请确保 geno_miss_bias_tools.py 文件位于同一目录下", file=sys.stderr)
    sys.exit(1)


def setup_logging(verbose: bool = False) -> None:
    """
    配置日志系统
    
    参数:
        verbose (bool): 是否启用详细输出模式
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    
    logging.basicConfig(
        level=log_level,
        format=log_format,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler('geno_miss_bias_cor.log', mode='a', encoding='utf-8')
        ]
    )


def validate_input_file(file_path: str) -> None:
    """
    验证输入JSON文件的有效性
    
    参数:
        file_path (str): JSON文件路径
        
    异常:
        FileNotFoundError: 文件不存在
        ValueError: 文件格式不正确
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"输入文件不存在: {file_path}")
    
    if not file_path.lower().endswith('.json'):
        logging.warning(f"输入文件可能不是JSON格式: {file_path}")
    
    # 检查文件是否可读
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            import json
            json.load(f)
    except json.JSONDecodeError as e:
        raise ValueError(f"JSON文件格式错误: {e}")
    except Exception as e:
        raise ValueError(f"无法读取文件: {e}")


def create_output_path(input_path: str, output_path: Optional[str] = None) -> str:
    """
    创建输出文件路径
    
    参数:
        input_path (str): 输入JSON文件路径
        output_path (Optional[str]): 用户指定的输出路径
        
    返回:
        str: 最终的输出文件路径
    """
    if output_path:
        # 确保输出目录存在
        output_dir = os.path.dirname(os.path.abspath(output_path))
        os.makedirs(output_dir, exist_ok=True)
        return output_path
    
    # 使用默认命名规则：在当前工作目录创建带时间戳的文件
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base_name = os.path.splitext(os.path.basename(input_path))[0]
    default_output = f"{base_name}_fdr_corrected_{timestamp}.tsv"
    
    return os.path.join(os.getcwd(), default_output)


def parse_arguments() -> argparse.Namespace:
    """
    解析命令行参数
    
    返回:
        argparse.Namespace: 解析后的参数对象
    """
    parser = argparse.ArgumentParser(
        description='基因型缺失偏差校正 - FDR多重检验校正工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本使用（输出到当前目录）
  %(prog)s -i /path/to/bias_results.merged.json
  
  # 指定输出路径
  %(prog)s -i /path/to/bias_results.merged.json -o /path/to/output.fdr.tsv
  
  # 启用详细输出
  %(prog)s -i /path/to/bias_results.merged.json -v
  
  # 查看帮助信息
  %(prog)s -h

注意事项:
  - 输入文件必须是由 merge_bias_results_json() 函数生成的JSON格式
  - 程序会对以下三个p值列进行FDR校正：
    * p_missing_to_missing_vs_other
    * p_not_missing_vs_drop_to_missing  
    * p_genotype_drop_composition
  - 校正后的q值会作为新列添加到输出文件中
  - 运行日志会同时输出到控制台和日志文件
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        metavar='JSON_FILE',
        help='输入的合并JSON结果文件路径（必需）'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        default=None,
        metavar='TSV_FILE',
        help='输出的FDR校正结果TSV文件路径（可选，默认在当前目录生成）'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='启用详细输出模式，显示更多调试信息'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='基因型缺失偏差校正工具 v1.0.0'
    )
    
    return parser.parse_args()


def main() -> int:
    """
    主函数 - 程序入口点
    
    返回:
        int: 程序退出码（0表示成功，非0表示失败）
    """
    try:
        # 解析命令行参数
        args = parse_arguments()
        
        # 设置日志系统
        setup_logging(args.verbose)
        
        # 记录程序开始运行
        logging.info("="*80)
        logging.info("基因型缺失偏差校正程序开始运行")
        logging.info(f"输入文件: {args.input}")
        logging.info(f"输出文件: {args.output if args.output else '自动生成'}")
        logging.info(f"详细模式: {'启用' if args.verbose else '禁用'}")
        
        # 验证输入文件
        logging.info("正在验证输入文件...")
        validate_input_file(args.input)
        logging.info("输入文件验证通过")
        
        # 确定输出路径
        output_path = create_output_path(args.input, args.output)
        logging.info(f"输出文件路径: {output_path}")
        
        # 执行FDR校正
        logging.info("开始执行FDR校正...")
        try:
            result_path = adjust_stat_summary_fdr(
                merged_json_path=args.input,
                out_path=output_path
            )
            
            # 验证输出文件是否成功生成
            if os.path.exists(result_path):
                file_size = os.path.getsize(result_path)
                logging.info(f"FDR校正完成！输出文件: {result_path}")
                logging.info(f"输出文件大小: {file_size:,} 字节")
                
                # 简单统计输出文件行数
                try:
                    with open(result_path, 'r') as f:
                        line_count = sum(1 for _ in f)
                    logging.info(f"输出文件行数: {line_count:,} 行（包含表头）")
                except Exception as e:
                    logging.warning(f"无法统计输出文件行数: {e}")
                
            else:
                logging.error("输出文件未成功生成")
                return 1
                
        except Exception as e:
            logging.error(f"FDR校正过程中发生错误: {e}")
            if args.verbose:
                import traceback
                logging.error(f"详细错误信息:\n{traceback.format_exc()}")
            return 1
        
        logging.info("程序执行完成")
        logging.info("="*80)
        return 0
        
    except KeyboardInterrupt:
        logging.info("程序被用户中断")
        return 130
    except Exception as e:
        logging.error(f"程序执行过程中发生未预期的错误: {e}")
        return 1


if __name__ == "__main__":
    """
    程序入口点
    
    当脚本被直接执行时调用main函数，并以其返回值作为程序退出码
    """
    exit_code = main()
    sys.exit(exit_code)


