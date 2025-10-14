#!/usr/bin/env python3
"""
JSON Manifest Tools for CTEPH AGP3K RVtest Analysis
Utility functions for creating structured JSON manifest files
"""

import json
import os
import subprocess
from datetime import datetime
from pathlib import Path


class ManifestGenerator:
    """JSON Manifest生成器类"""
    
    def __init__(self, verbose=False):
        self.verbose = verbose
    
    def log(self, message):
        """打印日志信息"""
        if self.verbose:
            print(f"[INFO] {message}")
    
    def get_file_info(self, filepath):
        """获取文件基本信息"""
        if not os.path.exists(filepath):
            self.log(f"Warning: File not found - {filepath}")
            return None
        
        try:
            stat = os.stat(filepath)
            return {
                "path": os.path.abspath(filepath),
                "filename": os.path.basename(filepath),
                "size_bytes": stat.st_size,
                "size_mb": round(stat.st_size / (1024 * 1024), 2),
                "modified_time": datetime.fromtimestamp(stat.st_mtime).isoformat(),
                "exists": True
            }
        except Exception as e:
            self.log(f"Error getting file info for {filepath}: {e}")
            return {"path": filepath, "exists": False, "error": str(e)}
    
    def run_command(self, cmd, description="command"):
        """安全地运行shell命令"""
        try:
            result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
            return result.decode().strip()
        except subprocess.CalledProcessError as e:
            self.log(f"Error running {description}: {e}")
            return None
        except Exception as e:
            self.log(f"Unexpected error running {description}: {e}")
            return None
    
    def get_vcf_stats(self, vcf_file):
        """获取VCF文件统计信息"""
        self.log(f"Analyzing VCF file: {vcf_file}")
        
        stats = {
            "sample_count": "unknown",
            "variant_count": "unknown",
            "chromosomes": "unknown",
            "chromosome_count": "unknown"
        }
        
        # 检查bcftools是否可用
        if self.run_command("which bcftools") is None:
            self.log("Warning: bcftools not found in PATH")
            return stats
        
        # 获取样本数量
        sample_result = self.run_command(f"bcftools query -l {vcf_file} | wc -l", "sample count")
        if sample_result:
            try:
                stats["sample_count"] = int(sample_result)
            except ValueError:
                pass
        
        # 获取变异数量
        variant_result = self.run_command(f"bcftools view -H {vcf_file} | wc -l", "variant count")
        if variant_result:
            try:
                stats["variant_count"] = int(variant_result)
            except ValueError:
                pass
        
        # 获取染色体信息
        chr_result = self.run_command(f"bcftools view -H {vcf_file} | cut -f1 | sort | uniq", "chromosome list")
        if chr_result:
            chromosomes = [c for c in chr_result.split('\n') if c.strip()]
            stats["chromosomes"] = chromosomes
            stats["chromosome_count"] = len(chromosomes)
        
        return stats
    
    def get_assoc_stats(self, assoc_file):
        """获取关联分析结果统计信息"""
        self.log(f"Analyzing association file: {assoc_file}")
        
        stats = {
            "total_genes": "unknown",
            "significant_genes_p005": "unknown",
            "significant_genes_p001": "unknown",
            "minimum_pvalue": "unknown",
            "significance_rate": "unknown"
        }
        
        if not os.path.exists(assoc_file):
            return stats
        
        try:
            # 统计基因数量（除去表头）
            gene_result = self.run_command(f"tail -n +2 {assoc_file} | wc -l", "gene count")
            if gene_result:
                gene_count = int(gene_result)
                stats["total_genes"] = gene_count
                
                # 获取显著性结果数量 (p < 0.05)
                sig_result = self.run_command(f"tail -n +2 {assoc_file} | awk -F'\\t' '$NF < 0.05' | wc -l", "significant genes p<0.05")
                if sig_result:
                    sig_count = int(sig_result)
                    stats["significant_genes_p005"] = sig_count
                    stats["significance_rate"] = round(sig_count / gene_count * 100, 2) if gene_count > 0 else 0
                
                # 获取高度显著结果数量 (p < 0.01)
                highly_sig_result = self.run_command(f"tail -n +2 {assoc_file} | awk -F'\\t' '$NF < 0.01' | wc -l", "significant genes p<0.01")
                if highly_sig_result:
                    stats["significant_genes_p001"] = int(highly_sig_result)
                
                # 获取最小p值
                min_p_result = self.run_command(f"tail -n +2 {assoc_file} | awk -F'\\t' '{{print $NF}}' | sort -g | head -1", "minimum p-value")
                if min_p_result:
                    try:
                        stats["minimum_pvalue"] = float(min_p_result)
                    except ValueError:
                        pass
        
        except Exception as e:
            self.log(f"Error analyzing association file {assoc_file}: {e}")
        
        return stats
    
    def create_manifest_data(self, vcf_file, assoc_files, results_dir, scripts_dir):
        """创建manifest数据结构"""
        self.log("Creating manifest data structure...")
        
        # 解析关联分析文件
        assoc_data = {}
        for tag, filepath in assoc_files.items():
            self.log(f"Processing association file: {tag}")
            
            mac_value = "none" if tag == "no_macmin" else int(tag.split("_")[1])
            description = f"No MAC minimum threshold applied" if tag == "no_macmin" else f"MAC minimum threshold = {mac_value}"
            
            assoc_data[tag] = {
                "file_info": self.get_file_info(filepath),
                "parameters": {
                    "mac_minimum": mac_value,
                    "description": description
                },
                "statistics": self.get_assoc_stats(filepath)
            }
        
        # 创建完整的manifest数据结构
        manifest = {
            "analysis_info": {
                "analysis_type": "rare_variant_association",
                "method": "SKAT-O",
                "created_time": datetime.now().isoformat(),
                "pipeline_version": "cteph_agp3k_rvtest",
                "description": "Rare variant association analysis using SKAT-O method for CTEPH patients",
                "generator_version": "1.0.0"
            },
            "input_data": {
                "vcf_file": self.get_file_info(vcf_file),
                "vcf_stats": self.get_vcf_stats(vcf_file)
            },
            "association_results": assoc_data,
            "file_paths": {
                "results_directory": os.path.abspath(results_dir),
                "scripts_directory": os.path.abspath(scripts_dir)
            },
            "summary": {
                "total_analyses": len(assoc_files),
                "mac_thresholds_tested": [assoc_data[tag]["parameters"]["mac_minimum"] for tag in assoc_data.keys()],
                "analysis_tags": list(assoc_files.keys())
            }
        }
        
        return manifest
    
    def save_manifest(self, manifest_data, output_file):
        """保存manifest到JSON文件"""
        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(manifest_data, f, indent=2, ensure_ascii=False)
            
            self.log(f"Manifest saved to: {os.path.abspath(output_file)}")
            return True
        except Exception as e:
            self.log(f"Error saving manifest: {e}")
            return False
    
    def print_summary(self, manifest_data):
        """打印分析摘要"""
        print("\n=== Analysis Summary ===")
        print(f"Analysis Type: {manifest_data['analysis_info']['analysis_type']}")
        print(f"Method: {manifest_data['analysis_info']['method']}")
        print(f"Created: {manifest_data['analysis_info']['created_time']}")
        
        print(f"\nInput Data:")
        vcf_stats = manifest_data['input_data']['vcf_stats']
        print(f"  VCF Samples: {vcf_stats['sample_count']}")
        print(f"  VCF Variants: {vcf_stats['variant_count']}")
        
        print(f"\nAssociation Results:")
        for tag, data in manifest_data['association_results'].items():
            stats = data['statistics']
            if isinstance(stats['total_genes'], int):
                print(f"  {tag}: {stats['total_genes']} genes, "
                      f"{stats['significant_genes_p005']} significant (p<0.05), "
                      f"{stats['significance_rate']}% rate")
        print("========================\n")


def validate_input_files(files_dict):
    """验证输入文件是否存在"""
    missing_files = []
    for tag, filepath in files_dict.items():
        if not os.path.exists(filepath):
            missing_files.append(f"{tag}: {filepath}")
    
    if missing_files:
        print("Warning: The following files are missing:")
        for f in missing_files:
            print(f"  - {f}")
        return False
    return True
