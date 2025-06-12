# %%
import pandas as pd
import subprocess
import argparse
import os

parser = argparse.ArgumentParser(description='Update fam and pheno files for plink2')
parser.add_argument('--info_file', type=str, required=True, help='Path to the info file (excel)')
parser.add_argument('--bed_prefix', type=str, required=True, help='Path to the bed prefix')
parser.add_argument('--case_prefix', type=str, required=True, help='Prefix for case IDs')
parser.add_argument('--out', type=str, required=True, help='Output prefix for the updated files')
args = parser.parse_args()


os.environ["PATH"] = "/home/b/b37974/:" + os.environ["PATH"]

info_df = pd.read_excel(args.info_file, header=0)

plink_prefix = args.bed_prefix

# %%
def generate_sex_files(info_df):
    sex_files = pd.DataFrame({
        '#FID': info_df['ID'],
        'IID': info_df['ID'],
        'SEX': info_df['Sex']
    })
    return sex_files

sex_files = generate_sex_files(info_df)

# %%
def generate_pheno_files(info_df, prefix):
    pheno_files = pd.DataFrame({
        '#FID': info_df['ID'],
        'IID': info_df['ID'],
        'PHENO1': [2 if id.startswith(prefix) else 1 for id in info_df['ID']]
    })
    return pheno_files

prefix = args.case_prefix
pheno_files = generate_pheno_files(info_df, prefix)

# %%
sex_files.to_csv('update_sex.txt', sep='\t', index=False)
pheno_files.to_csv('update_pheno.txt', sep='\t', index=False)

# %%
# Define the command to run plink2
plink_command = [
    "/home/b/b37974/plink2",
    "--bfile", plink_prefix,
    "--pheno", "update_pheno.txt",
    "--update-sex", "update_sex.txt",
    "--make-bed",
    "--out", args.out
]

# Run the command using subprocess
try:
    subprocess.run(plink_command, check=True)
    print("Plink2 command executed successfully.")
except subprocess.CalledProcessError as e:
    print(f"Error occurred while running plink2: {e}")


