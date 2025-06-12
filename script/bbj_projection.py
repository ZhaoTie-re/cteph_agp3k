# %%
import pandas as pd
import numpy as np
import subprocess
import os
import argparse

# Set up environment
os.environ["PATH"] = "/home/b/b37974/:" + os.environ["PATH"]

# Argument parser for command line execution
parser = argparse.ArgumentParser(description="BBJ Projection")
parser.add_argument("--my_bed_prefix", type=str, help="Path to the my bed prefix for plink2 operations (after all QC steps)")
parser.add_argument("--bbj_bed_prefix", type=str, help="Path to the BBJ bed prefix for plink2 operations (after bbj_prepare.py & bbj_pca.py)")
parser.add_argument("--bbj_prune_in", type=str, default="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/bbj_projection/02.bbj_pca/bbj.b38.auto.sqc.vqc.norm.no_high_ld.prune.in", help="Path to the BBJ prune input file")
parser.add_argument("--my_prefix_out", type=str, default="cteph_agp3k", help="Output prefix for my bed file after pruning")
parser.add_argument("--bbj_prefix_out", type=str, default="bbj", help="Output prefix for BBJ bed file after pruning")
parser.add_argument("--bbj_pca_acount", type=str, default="bbj.b38.auto.sqc.vqc.norm.no_high_ld.prune.pca.acount.txt", help="Path to the BBJ PCA account file")
parser.add_argument("--bbj_pca_eigenvec_allele", type=str, default="bbj.b38.auto.sqc.vqc.norm.no_high_ld.prune.pca.eigenvec.allele", help="Path to the BBJ PCA eigenvec allele file")
parser.add_argument("--threads", type=int, default=8, help="Number of threads for plink2 operations")
args = parser.parse_args()


def run_command(command):
    """Run a shell command (tuple/list or string) and return the output."""
    if isinstance(command, (tuple, list)):
        result = subprocess.run(command, capture_output=True, text=True)
    else:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"Command failed:\n{command}\n{result.stderr}")
    return result.stdout.strip()

bed_my_prefix = args.my_bed_prefix
bed_bbj_prefix = args.bbj_bed_prefix
bbj_prune_in = args.bbj_prune_in

my_prefix = args.my_prefix_out
bbj_prefix = args.bbj_prefix_out

cmd_bed_prune = [
    "/home/b/b37974/plink2",
    "--bfile", bed_my_prefix,
    "--extract", bbj_prune_in,
    "--make-bed",
    "--out", f"{my_prefix}.bbj_pruned",
    "--threads", str(args.threads)  # Convert threads to string for command
]

run_command(cmd_bed_prune)

# %%
my_bbj_prune = pd.read_csv(f"{my_prefix}.bbj_pruned.bim", sep="\t", header=None)
my_bbj_prune[1].to_csv(f"{my_prefix}.bbj_pruned.snps.txt", index=False, header=False)

# %%
cmd_share_snps = [
    "/home/b/b37974/plink2",
    "--bfile", bed_bbj_prefix,
    "--extract", f"{my_prefix}.bbj_pruned.snps.txt",
    "--make-bed",
    "--out", f"{bbj_prefix}.bbj_pruned",
    "--threads", str(args.threads)  # Convert threads to string for command
]

run_command(cmd_share_snps)

# %%
cmd_merge = [
    "/home/b/b37974/plink",
    "--bfile", f"{my_prefix}.bbj_pruned",
    "--bmerge", f"{bbj_prefix}.bbj_pruned",
    "--keep-allele-order",
    "--make-bed",
    "--out", f"{my_prefix}.{bbj_prefix}.bbj_pruned.merged"
]

run_command(cmd_merge)

# %%

cmd_projection = [
    "/home/b/b37974/plink2",
    "--bfile", f"{my_prefix}.{bbj_prefix}.bbj_pruned.merged",
    "--read-freq", args.bbj_pca_acount,
    "--extract", f"{my_prefix}.bbj_pruned.snps.txt",
    "--score", args.bbj_pca_eigenvec_allele, "2", "6", "header-read", "no-mean-imputation", "variance-standardize", "list-variants",
    "--score-col-nums", "7-16",
    "--out", f"{my_prefix}.{bbj_prefix}.projection",
    "--threads", str(args.threads)  # Convert threads to string for command
]

run_command(cmd_projection)
