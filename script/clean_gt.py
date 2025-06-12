# %%
import pandas as pd
import numpy as np
import subprocess
import os
import argparse

# set up the environment
os.environ["PATH"] = "/home/b/b37974/:" + os.environ["PATH"]

# argument parser
parser = argparse.ArgumentParser(description="Prepare final genotype data, phnotype, and covariance files for CTEPH AGP3K project.")
parser.add_argument("--bed_prefix", type=str, required=True, help="Path to the bed prefix for plink2 operations (after QC)")
parser.add_argument("--projection_keep", type=str, required=True, help="Path to the projection keep file (after bbj_projection_keep.py)")
parser.add_argument("--kinship_remove", type=str, required=True, help="Path to the kinship remove file (after kinship.pi_hat.rev1.py)")
parser.add_argument("--sscore", type=str, required=True, help="Path to the sscore file (after bbj_projection.py)")
parser.add_argument("--out_prefix", type=str, default="cteph_agp3k.s_qc.gt_qc.v_qc.hwe.sample_keep", help="Prefix for the output files")
parser.add_argument("--threads", type=int, default=8, help="Number of threads for plink2 operations")
args = parser.parse_args()

projection_keep = pd.read_csv(args.projection_keep, header=0, sep='\t', dtype=str)
kinship_remove = pd.read_csv(args.kinship_remove, header=None, dtype=str)

projection_keep_filtered = projection_keep[~projection_keep['IID'].isin(kinship_remove[0])]
projection_keep_filtered.to_csv('sample_keep.txt',sep='\t', index=False, header=True)

# %%

def run_command(command):
    """Run a shell command (tuple/list or string) and return the output."""
    if isinstance(command, (tuple, list)):
        result = subprocess.run(command, capture_output=True, text=True)
    else:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"Command failed:\n{command}\n{result.stderr}")
    return result.stdout.strip()

bed_prefix = args.bed_prefix
out_prefix = args.out_prefix

cmd_sample_keep = [
    "/home/b/b37974/plink2",
    "--bfile", bed_prefix,
    "--keep", "sample_keep.txt",
    "--make-bed",
    "--out", out_prefix, 
    "--threads", str(args.threads)  # Convert threads to string for command
]

run_command(cmd_sample_keep)

# %%
fam_df = pd.read_csv(f'{out_prefix}.fam', sep='\t', header=None, dtype=str)
sex_df = fam_df[[0, 1, 4]].copy()
sex_df.columns = ['#FID', 'IID', 'SEX']

pheno_df = fam_df[[0, 1, 5]].copy()
pheno_df.columns = ['#FID', 'IID', 'PHENO1']
pheno_df.to_csv('pheno.txt', sep='\t', index=False, header=True)


# %%
sscore = pd.read_csv(args.sscore, header=0, sep='\t')

sscore_filtered = sscore[sscore['IID'].isin(projection_keep_filtered['IID'])].drop(sscore.columns[[2, 3, 4]], axis=1)
covariance_df = pd.merge(sex_df, sscore_filtered, on=['#FID', 'IID'], how='left')
covariance_df.to_csv('covariance.txt', sep='\t', index=False, header=True)



