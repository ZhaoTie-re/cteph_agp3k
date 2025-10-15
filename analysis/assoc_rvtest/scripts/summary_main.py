#!/usr/bin/env python3
"""
摘要主程序CLI工具
用于更新JSON清单文件的命令行界面工具
"""

import sys
import importlib
import argparse
import logging
from pathlib import Path

# 设置自定义库的模块路径
MODULE_DIR = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest/scripts"
if MODULE_DIR not in sys.path:
    sys.path.insert(0, MODULE_DIR)

# 导入并重载自定义库（开发模式）
import summary_tools
importlib.reload(summary_tools)

from summary_tools import update_json_manifest


def setup_logging(verbose=False):
    """设置日志配置"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="使用指定参数更新JSON清单文件",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--json-path',
        type=str,
        default='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest/results/05.json_manifest/cteph_agp3k_file_manifest.json',
        help='JSON清单文件的路径'
    )
    
    parser.add_argument(
        '--num-var-thr',
        type=int,
        default=3,
        help='变异数量阈值'
    )
    
    parser.add_argument(
        '--threads',
        type=int,
        default=8,
        help='使用的线程数'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='启用详细日志输出'
    )
    
    return parser.parse_args()


def validate_arguments(args):
    """验证命令行参数"""
    # 检查JSON文件是否存在
    json_path = Path(args.json_path)
    if not json_path.exists():
        raise FileNotFoundError(f"未找到JSON文件: {args.json_path}")
    
    # 验证正整数
    if args.num_var_thr < 0:
        raise ValueError("num_var_thr必须是非负整数")
    
    if args.threads < 1:
        raise ValueError("threads必须是正整数")
    
    logging.info(f"JSON路径: {args.json_path}")
    logging.info(f"变异数量阈值: {args.num_var_thr}")
    logging.info(f"线程数: {args.threads}")


def main():
    """主函数"""
    try:
        # 解析参数
        args = parse_arguments()
        
        # 设置日志
        setup_logging(args.verbose)
        
        # 验证参数
        validate_arguments(args)
        
        logging.info("开始更新JSON清单...")
        
        # 在开发模式下重载自定义库
        importlib.reload(summary_tools)
        logging.debug("已重载summary_tools模块")
        
        # 执行主函数
        update_json = update_json_manifest(
            json_path=args.json_path,
            num_var_thr=args.num_var_thr,
            threads=args.threads
        )
        
        logging.info("JSON清单更新成功完成")
        
    except Exception as e:
        logging.error(f"发生错误: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()


