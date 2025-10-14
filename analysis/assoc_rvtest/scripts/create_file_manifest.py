#!/usr/bin/env python3
"""
Simple file path manifest generator
Store structured file paths in JSON format for easy access
"""

import json
import os
import argparse
from datetime import datetime


def create_file_manifest(vcf_file, assoc_files, results_dir, scripts_dir, output_file):
    """创建简单的文件路径manifest"""
    
    # 创建文件路径结构
    manifest = {
        "created_time": datetime.now().isoformat(),
        "input_files": {
            "vcf_file": os.path.abspath(vcf_file),
            "vcf_index": os.path.abspath(vcf_file + ".tbi") if os.path.exists(vcf_file + ".tbi") else None
        },
        "association_results": {
            "no_macmin": os.path.abspath(assoc_files["no_macmin"]),
            "macmin_5": os.path.abspath(assoc_files["macmin_5"]),
            "macmin_10": os.path.abspath(assoc_files["macmin_10"])
        },
        "directories": {
            "results_dir": os.path.abspath(results_dir),
            "scripts_dir": os.path.abspath(scripts_dir)
        }
    }
    
    # 保存JSON文件
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    
    print(f"File manifest saved to: {os.path.abspath(output_file)}")


def main():
    parser = argparse.ArgumentParser(description="Create file path manifest")
    parser.add_argument("--vcf", required=True, help="VCF file path")
    parser.add_argument("--assoc-no-macmin", required=True, help="Association file without MAC minimum")
    parser.add_argument("--assoc-macmin-5", required=True, help="Association file with MAC minimum 5")
    parser.add_argument("--assoc-macmin-10", required=True, help="Association file with MAC minimum 10")
    parser.add_argument("--results-dir", required=True, help="Results directory")
    parser.add_argument("--scripts-dir", required=True, help="Scripts directory")
    parser.add_argument("--output", default="file_manifest.json", help="Output JSON file")
    
    args = parser.parse_args()
    
    assoc_files = {
        "no_macmin": args.assoc_no_macmin,
        "macmin_5": args.assoc_macmin_5,
        "macmin_10": args.assoc_macmin_10
    }
    
    create_file_manifest(
        vcf_file=args.vcf,
        assoc_files=assoc_files,
        results_dir=args.results_dir,
        scripts_dir=args.scripts_dir,
        output_file=args.output
    )


if __name__ == "__main__":
    main()
