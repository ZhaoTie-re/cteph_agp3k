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