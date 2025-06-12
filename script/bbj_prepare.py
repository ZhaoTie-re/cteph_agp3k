# %%
import pandas as pd
import numpy as np
import subprocess
import os
import argparse

# Set up environment
os.environ["PATH"] = "/home/b/b37974/:" + os.environ["PATH"]

# Argument parser for command line execution
parser = argparse.ArgumentParser(description="Prepare BBJ data")
parser.add_argument("--bbj_bed_prefix", type=str, help="Path to the BBJ bed prefix for plink2 operations (after sample QC & variant QC)")
parser.add_argument("--threads", type=int, default=8, help="Number of threads for plink2 operations")
args = parser.parse_args()

bbj_bed_prefix = args.bbj_bed_prefix

def run_command(command):
    """Run a shell command (tuple/list or string) and return the output."""
    if isinstance(command, (tuple, list)):
        result = subprocess.run(command, capture_output=True, text=True)
    else:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"Command failed:\n{command}\n{result.stderr}")
    return result.stdout.strip()

command = [
    "/home/b/b37974/plink2", 
    "--bfile", bbj_bed_prefix,
    "--export", "vcf", "bgz",
    "--maf", "0.01",
    "--out", "bbj.b38.auto.qc.maf0.01", 
    "--threads", str(args.threads)  # Convert threads to string for command
]

command_tbi = [
    "bcftools", "index", 
    "-t", "bbj.b38.auto.qc.maf0.01.vcf.gz", 
    "--threads", str(args.threads)
]

vcf_output = run_command(command)
vcf_tbi_output = run_command(command_tbi)
print(vcf_output)
print(vcf_tbi_output)

# %%
rename_chr = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/rename_chr.txt"

command_rechr = [
    "bcftools", "annotate",
    "bbj.b38.auto.qc.maf0.01.vcf.gz", 
    "--rename-chrs", rename_chr,
    "-Oz", 
    "-o", "bbj.b38.auto.sqc.vqc.rechr.vcf.gz", 
    "--threads", str(args.threads)
]

command_rechr_tbi = [
    "bcftools", "index", 
    "-t", "bbj.b38.auto.sqc.vqc.rechr.vcf.gz", 
    "--threads", str(args.threads)
]

rename_output = run_command(command_rechr)
rename_tbi_output = run_command(command_rechr_tbi)
print(rename_output)
print(rename_tbi_output)

# %%
fasta_ref = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa"

command_norm = [
    "bcftools", "norm",
    "--fasta-ref", fasta_ref,
    "-m-",
    "--check-ref", "s", "bbj.b38.auto.sqc.vqc.rechr.vcf.gz",
    "-Oz",
    "-o", "bbj.b38.auto.sqc.vqc.rechr.norm.vcf.gz",
    "--threads", str(args.threads)
]

command_norm_tbi = [
    "bcftools", "index", 
    "-t", "bbj.b38.auto.sqc.vqc.rechr.norm.vcf.gz", 
    "--threads", str(args.threads)
]

norm_output = run_command(command_norm)
norm_tbi_output = run_command(command_norm_tbi)
print(norm_output)
print(norm_tbi_output)

# %%
command_setids = [
    "bcftools", "annotate",
    "--set-id", "%CHROM:%POS:%REF:%ALT",
    "bbj.b38.auto.sqc.vqc.rechr.norm.vcf.gz",
    "-Oz",
    "-o", "bbj.b38.auto.sqc.vqc.rechr.norm.setid.vcf.gz",
    "--threads", str(args.threads)
]

command_setids_tbi = [
    "bcftools", "index", 
    "-t", "bbj.b38.auto.sqc.vqc.rechr.norm.setid.vcf.gz", 
    "--threads", str(args.threads)
]

setids_output = run_command(command_setids)
setids_tbi_output = run_command(command_setids_tbi)
print(setids_output)
print(setids_tbi_output)

# %%
command_vcf2bed = [
    "/home/b/b37974/plink",
    "--vcf", "bbj.b38.auto.sqc.vqc.rechr.norm.setid.vcf.gz", 
    "--make-bed",
    "--keep-allele-order", 
    "--double-id",
    "--out", "bbj.b38.auto.sqc.vqc.rechr.norm.setid"
]

vcf2bed_output = run_command(command_vcf2bed)
print(vcf2bed_output)

# fam_file = "bbj.b38.auto.sqc.vqc.rechr.norm.setid.fam"
# output_file = "update_fam.txt"

# fam_df = pd.read_csv(fam_file, delim_whitespace=True, header=None)

# update_df = pd.DataFrame({
#     "#OLD_FID": fam_df[0],
#     "OLD_IID": fam_df[1],
#     "NEW_FID": "bbj_" + fam_df[0].astype(str),
#     "NEW_IID": "bbj_" + fam_df[1].astype(str),
# })

# update_df.to_csv(output_file, sep="\t", header=True, index=False)

# command_update_fam = [
#     "/home/b/b37974/plink2",
#     "--bfile", "bbj.b38.auto.sqc.vqc.rechr.norm.setid",
#     "--update-ids", output_file,
#     "--make-bed",
#     "--out", "bbj.b38.auto.sqc.vqc.rechr.norm.setid.update_fam",
#     "--threads", str(args.threads)
# ]
