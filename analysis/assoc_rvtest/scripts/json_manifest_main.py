#!/usr/bin/env python3
"""
JSON Manifest Generator for CTEPH AGP3K RVtest Analysis
Generate structured JSON manifest file for association analysis results
"""

import argparse
import sys
from pathlib import Path

# 导入自定义工具模块
from json_manifest_tools import ManifestGenerator, validate_input_files


def main():
    parser = argparse.ArgumentParser(description="Generate JSON manifest for RVtest association results")
    parser.add_argument("--vcf", required=True, help="Input VCF file path")
    parser.add_argument("--assoc-no-macmin", required=True, help="Association file without MAC minimum")
    parser.add_argument("--assoc-macmin-5", required=True, help="Association file with MAC minimum 5")
    parser.add_argument("--assoc-macmin-10", required=True, help="Association file with MAC minimum 10")
    parser.add_argument("--results-dir", required=True, help="Results directory path")
    parser.add_argument("--scripts-dir", required=True, help="Scripts directory path")
    parser.add_argument("--output", default="cteph_agp3k_rvtest_manifest.json", 
                       help="Output JSON manifest file name")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    
    args = parser.parse_args()
    
    # 创建manifest生成器
    generator = ManifestGenerator(verbose=args.verbose)
    
    if args.verbose:
        print("Creating JSON manifest for CTEPH AGP3K RVtest analysis...")
        print(f"VCF file: {args.vcf}")
        print(f"Output file: {args.output}")
    
    # 准备关联分析文件字典
    assoc_files = {
        "no_macmin": args.assoc_no_macmin,
        "macmin_5": args.assoc_macmin_5,
        "macmin_10": args.assoc_macmin_10
    }
    
    # 验证输入文件
    all_files = {"vcf": args.vcf, **assoc_files}
    if not validate_input_files(all_files):
        generator.log("Some input files are missing, but continuing with available files...")
    
    # 创建manifest数据
    manifest_data = generator.create_manifest_data(
        vcf_file=args.vcf,
        assoc_files=assoc_files,
        results_dir=args.results_dir,
        scripts_dir=args.scripts_dir
    )
    
    # 保存manifest文件
    success = generator.save_manifest(manifest_data, args.output)
    
    if success:
        # 打印摘要信息
        generator.print_summary(manifest_data)
        print(f"✓ Manifest file created successfully: {args.output}")
        sys.exit(0)
    else:
        print("✗ Failed to create manifest file")
        sys.exit(1)


if __name__ == "__main__":
    main()
