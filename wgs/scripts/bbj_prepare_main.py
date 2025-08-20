"""
模块名称：bbj_prepare_main

概述：
    本脚本用于准备 BBJ（BioBank Japan）队列的基因型数据。通过调用 prepare_bbj_genotype 函数，对原始基因型文件进行预处理，包括染色体重命名、参考基因组对齐、MAF 过滤等操作，以便下游分析使用。

适用场景：
    适用于 BBJ 队列的基因型数据预处理，特别是在需要标准化染色体编号、筛选变异位点、以及保证数据格式统一时使用。

输入与输出（核心函数）：
    main() 函数解析命令行参数，并调用 prepare_bbj_genotype。参数说明如下：
        --bbj_bed_prefix   ：必需，BBJ bed 文件前缀路径。
        --rename_chr       ：可选，染色体重命名文件路径（默认：/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/rename_chr.txt）。
        --fasta_ref        ：可选，参考基因组 fasta 文件路径（默认：/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa）。
        --threads          ：可选，使用的线程数（默认：32）。
        --maf_threshold    ：可选，MAF 阈值（默认：0.05）。
    prepare_bbj_genotype 返回处理后数据的前缀路径。

依赖/环境：
    - Python 3.x
    - 依赖模块：bbj_projection_tools
    - 外部可执行：PLINK

约定与注意事项：
    - rename_chr、fasta_ref 默认值已内置，如需更改请显式指定参数。
    - MAF 默认阈值为 0.05，可根据需求调整。
    - 线程数默认 32，根据服务器资源酌情设置。

快速示例：
    python bbj_prepare_main.py --bbj_bed_prefix /path/to/bbj_data/bbj_genotype

作者: ZHAO TIE
"""
import importlib
import bbj_projection_tools
import argparse

# 强制重新加载模块（适用于开发调试阶段）
importlib.reload(bbj_projection_tools)

from bbj_projection_tools import (
    prepare_bbj_genotype,
)

def main():
    parser = argparse.ArgumentParser(description="Prepare BBJ genotype data")
    parser.add_argument('--bbj_bed_prefix', type=str, required=True, help='Path prefix for BBJ bed files')
    parser.add_argument('--rename_chr', type=str, default="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/rename_chr.txt", help='Path to rename_chr.txt file')
    parser.add_argument('--fasta_ref', type=str, default="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa", help='Path to fasta reference file')
    parser.add_argument('--threads', type=int, default=32, help='Number of threads to use')
    parser.add_argument('--maf_threshold', type=float, default=0.05, help='MAF threshold')

    args = parser.parse_args()

    bbj_prepared_prefix = prepare_bbj_genotype(
        bbj_bed_prefix=args.bbj_bed_prefix,
        threads=args.threads,
        rename_chr=args.rename_chr,
        fasta_ref=args.fasta_ref,
        maf=args.maf_threshold,
    )

if __name__ == "__main__":
    main()