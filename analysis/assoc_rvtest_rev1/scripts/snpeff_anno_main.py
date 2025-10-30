#!/usr/bin/env python3
"""
snpEff VCF注释工具 - 命令行接口

这个工具用于对VCF文件进行snpEff功能注释，支持并行处理以提高效率。
主要功能：
1. 自动为VCF添加chr前缀（如果需要）
2. 使用snpEff预计算的TSV注释文件为VCF添加功能注释
3. 支持并行和串行两种模式
4. 提供详细的日志记录和进度追踪

作者: ZHAO TIE
日期: 2025-10-12
版本: 1.0.0
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime

# 添加模块路径
MODULE_DIR = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest_rev1/scripts"
if MODULE_DIR not in sys.path:
    sys.path.insert(0, MODULE_DIR)

# 导入自定义模块
# import snpeff_anno_tools
from snpeff_anno_tools import add_chr_prefix_to_vcf, annotate_vcf_with_snpeff_tsv

# 开发模式：每次调用都重新加载自定义模块
import importlib

add_chr_prefix_to_vcf_module = sys.modules.get("snpeff_anno_tools")
if add_chr_prefix_to_vcf_module:
    importlib.reload(add_chr_prefix_to_vcf_module)
    from snpeff_anno_tools import add_chr_prefix_to_vcf, annotate_vcf_with_snpeff_tsv


def setup_argument_parser():
    """
    设置命令行参数解析器
    
    Returns:
        argparse.ArgumentParser: 配置好的参数解析器
    """
    parser = argparse.ArgumentParser(
        description="snpEff VCF注释工具 - 为VCF文件添加功能注释",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本用法（使用默认参数）
  python snpeff_anno_main.py
  
  # 自定义输入文件
  python snpeff_anno_main.py --vcf-path /path/to/your/file.vcf.gz
  
  # 自定义注释目录
  python snpeff_anno_main.py --snpeff-dir /path/to/annotations
  
  # 调整性能参数
  python snpeff_anno_main.py --max-workers 8 --threads 4
  
  # 使用串行模式
  python snpeff_anno_main.py --no-parallel
  
  # 保留缓存文件用于调试
  python snpeff_anno_main.py --keep-cache

注意事项:
  - 默认使用预设的文件路径，适合当前项目环境
  - 支持并行处理以提高大文件注释效率
  - 自动生成详细日志文件便于调试和监控
        """
    )
    
    # 输入文件选项
    parser.add_argument(
        "--vcf-path", 
        type=str,
        default="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest/results/01.rvtest_prepare/cteph_agp3k.rare.rand300x10000.all.nochr.norm.vcf.gz",
        help="输入VCF文件路径（默认：项目预设路径）"
    )
    
    parser.add_argument(
        "--snpeff-dir", 
        type=str,
        default="/LARGE1/gr10478/platform/JHRPv4/workspace/pipeline/output/snpEff.v4.index",
        help="snpEff注释文件目录路径（默认：项目预设路径）"
    )
    
    # 输出选项
    output_group = parser.add_argument_group("输出选项")
    output_group.add_argument(
        "-o", "--output", 
        type=str, 
        default=None,
        help="输出VCF文件路径（默认：工作目录下的*.snpeff.vcf.gz）"
    )
    output_group.add_argument(
        "--force", 
        action="store_true",
        default=True,
        help="强制覆盖已存在的输出文件（默认启用）"
    )
    
    # 性能选项
    performance_group = parser.add_argument_group("性能选项")
    performance_group.add_argument(
        "--parallel", 
        action="store_true", 
        default=True,
        help="启用并行模式（默认启用）"
    )
    performance_group.add_argument(
        "--no-parallel", 
        action="store_true",
        help="禁用并行模式，使用串行处理"
    )
    performance_group.add_argument(
        "--max-workers", 
        type=int, 
        default=16,
        help="最大并行worker数量（默认：16）"
    )
    performance_group.add_argument(
        "--threads", 
        type=int, 
        default=8,
        help="bcftools使用的线程数（默认：8）"
    )
    
    # 日志和调试选项
    debug_group = parser.add_argument_group("日志和调试选项")
    debug_group.add_argument(
        "--log-file", 
        type=str, 
        default=None,
        help="日志文件路径（默认：snpeff_annotation.log）"
    )
    debug_group.add_argument(
        "--log-level", 
        choices=["DEBUG", "INFO", "WARNING", "ERROR"], 
        default="INFO",
        help="日志级别（默认：INFO）"
    )
    debug_group.add_argument(
        "--keep-cache", 
        action="store_true",
        help="保留中间缓存文件用于调试（默认删除）"
    )
    
    # 版本信息
    parser.add_argument(
        "--version", 
        action="version", 
        version="snpEff VCF注释工具 v1.0.0"
    )
    
    return parser


def validate_inputs(args):
    """
    验证输入文件和目录的有效性
    
    Args:
        args: 命令行参数
        
    Raises:
        SystemExit: 当输入无效时退出程序
    """
    # 验证VCF文件
    vcf_path = Path(args.vcf_path)
    if not vcf_path.exists():
        print(f"❌ 错误: VCF文件不存在 - {vcf_path}")
        sys.exit(1)
    
    if not str(vcf_path).endswith('.vcf.gz'):
        print(f"❌ 错误: 输入文件必须是.vcf.gz格式 - {vcf_path}")
        sys.exit(1)
    
    # 验证注释目录
    snpeff_dir = Path(args.snpeff_dir)
    if not snpeff_dir.exists():
        print(f"❌ 错误: snpEff注释目录不存在 - {snpeff_dir}")
        sys.exit(1)
    
    # 检查注释文件
    tsv_files = list(snpeff_dir.glob("*.tsv.2.gz"))
    if not tsv_files:
        print(f"❌ 错误: 注释目录中未找到*.tsv.2.gz文件 - {snpeff_dir}")
        sys.exit(1)
    
    print(f"✓ 找到 {len(tsv_files)} 个注释文件")


def print_configuration(args):
    """
    打印任务配置信息
    
    Args:
        args: 命令行参数
    """
    print("=" * 70)
    print("snpEff VCF注释工具 - 任务配置")
    print("=" * 70)
    print(f"输入VCF文件    : {args.vcf_path}")
    print(f"注释文件目录    : {args.snpeff_dir}")
    print(f"输出文件        : {args.output if args.output else '自动生成'}")
    print(f"运行模式        : {'并行' if args.parallel else '串行'}")
    if args.parallel:
        print(f"最大worker数   : {args.max_workers}")
    print(f"bcftools线程数 : {args.threads}")
    print(f"强制覆盖        : {'是' if args.force else '否'}")
    print(f"保留缓存        : {'是' if args.keep_cache else '否'}")
    print(f"日志级别        : {args.log_level}")
    print("=" * 70)


def main():
    """
    主函数 - 执行完整的VCF注释流程
    """
    # 解析命令行参数
    parser = setup_argument_parser()
    args = parser.parse_args()
    
    # 处理并行选项冲突
    if args.no_parallel:
        args.parallel = False
    
    # 验证输入
    validate_inputs(args)
    
    # 打印配置
    print_configuration(args)
    
    try:
        total_start_time = datetime.now()
        
        # 步骤1: 添加chr前缀（如果需要）
        print("\n🔄 步骤1: 检查并添加chr前缀...")
        step1_start = datetime.now()
        
        vcf_with_chr = add_chr_prefix_to_vcf(
            vcf_path=args.vcf_path,
            force=False  # 如果已有chr前缀则跳过
        )
        
        step1_duration = datetime.now() - step1_start
        print(f"✓ chr前缀处理完成，耗时: {step1_duration}")
        print(f"  处理后VCF: {vcf_with_chr}")
        
        # 步骤2: snpEff功能注释
        print("\n🔄 步骤2: 执行snpEff功能注释...")
        step2_start = datetime.now()
        
        annotated_vcf = annotate_vcf_with_snpeff_tsv(
            vcf_path=vcf_with_chr,
            snpeff_tsv_dir=args.snpeff_dir,
            out_path=args.output,
            threads=args.threads,
            parallel=args.parallel,
            max_workers=args.max_workers,
            force=args.force,
            remove_cache=not args.keep_cache,
            log_file=args.log_file,
            log_level=args.log_level,
        )
        
        step2_duration = datetime.now() - step2_start
        total_duration = datetime.now() - total_start_time
        
        # 显示完成信息
        print("\n" + "=" * 70)
        print("🎉 任务完成摘要")
        print("=" * 70)
        print(f"✓ 总耗时: {total_duration}")
        print(f"  - chr前缀处理: {step1_duration}")
        print(f"  - snpEff注释 : {step2_duration}")
        print(f"✓ 最终输出文件: {annotated_vcf}")
        
        # 文件信息
        if annotated_vcf.exists():
            file_size = annotated_vcf.stat().st_size / (1024 * 1024)  # MB
            print(f"✓ 输出文件大小: {file_size:.2f} MB")
            
            # 检查索引文件
            tbi_file = Path(str(annotated_vcf) + ".tbi")
            if tbi_file.exists():
                print(f"✓ 索引文件已生成: {tbi_file}")
        
        # 日志文件信息
        log_path = args.log_file if args.log_file else "snpeff_annotation.log"
        print(f"✓ 详细日志: {log_path}")
        
        print("=" * 70)
        print("🎉 所有任务成功完成！")
        
    except KeyboardInterrupt:
        print("\n\n⚠️  用户中断操作")
        sys.exit(1)
        
    except Exception as e:
        print(f"\n\n❌ 任务执行失败: {e}")
        print(f"错误类型: {type(e).__name__}")
        log_path = args.log_file if args.log_file else "snpeff_annotation.log"
        print(f"请查看日志文件获取详细信息: {log_path}")
        sys.exit(1)


if __name__ == "__main__":
    main()


