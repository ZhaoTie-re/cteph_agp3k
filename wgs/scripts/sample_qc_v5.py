#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
sample_qc_v5.py

Unified Sample QC Script for WGS Pipeline
Performs PLINK2-based QC (Missingness, Heterozygosity, Sex Check, Relatedness, PCA)
Generates publication-ready metrics and plots.

Author: ZHAO TIE
Date: 2026-01-06
"""

import argparse
import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging

# Configure Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger()

def run_command(cmd, shell=False):
    logger.info(f"Running command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    try:
        subprocess.run(cmd, shell=shell, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {e.stderr}")
        raise e

def parse_arguments():
    parser = argparse.ArgumentParser(description="WGS Sample QC Pipeline")
    parser.add_argument("--bed", required=True, help="Input BED file prefix")
    parser.add_argument("--sex-bed", help="Optional: Separate BED file prefix for Sex Check (e.g. combined X chr)")
    parser.add_argument("--info", required=True, help="Sample Info Excel file")
    parser.add_argument("--high-ld", required=True, help="High LD regions file (for PCA filtering)")
    parser.add_argument("--id-col", required=True, help="Column name for Sample ID in Info file")
    parser.add_argument("--platform-col", default="Platform", help="Column name for Platform (optional filtering)")
    parser.add_argument("--target-dp-col", default="Target_DP", help="Column for Target Depth")
    parser.add_argument("--mean-dp-col", default="DP", help="Column for Mean Depth")
    parser.add_argument("--sex-col", default="Sex", help="Column for Sex")
    parser.add_argument("--group-col", default="Outcome", help="Column for Phenotype/Group")
    parser.add_argument("--out-prefix", required=True, help="Output prefix")
    parser.add_argument("--threads", type=str, default="4", help="Number of threads")
    return parser.parse_args()

def load_fam(bed_prefix):
    fam_file = f"{bed_prefix}.fam"
    try:
        fam = pd.read_csv(fam_file, sep=r'\s+', header=None, 
                        names=['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO'])
        # Phenotype: 1=Control, 2=Case, 0/-9=Missing
        # Adjust mapping if needed. Typically 1=Control, 2=Case
        fam['Group_FAM'] = fam['PHENO'].map({1: 'Control', 2: 'Case', 0: 'Missing', -9: 'Missing'}).fillna('Unknown')
    except Exception as e:
        logger.error(f"Failed to read fam file: {fam_file}. Error: {e}")
        raise e
    return fam

class SampleQC:
    def __init__(self, args):
        self.args = args
        self.bed_prefix = args.bed
        self.out_prefix = args.out_prefix
    
    def run_qc_steps(self):
        # 1. Missingness (with MAF 0.05 filter as per old logic)
        run_command(['plink2', '--bfile', self.bed_prefix, '--missing', '--maf', '0.05', '--out', self.out_prefix, '--threads', self.args.threads])
        
        # 2. Heterozygosity
        run_command(['plink2', '--bfile', self.bed_prefix, '--het', '--out', self.out_prefix, '--threads', self.args.threads])
        
        # 3. Sex Check (Using plink2 alpha6 as requested)
        if self.args.sex_bed:
             logger.info(f"Using separate Sex BED for check-sex: {self.args.sex_bed}")
             run_command(['plink2', '--bfile', self.args.sex_bed, '--check-sex', '--out', self.out_prefix, '--threads', self.args.threads])
        else:
             run_command(['plink2', '--bfile', self.bed_prefix, '--check-sex', '--out', self.out_prefix, '--threads', self.args.threads])
        
        # 4. Relatedness (PLINK 1.9 Genome / Pi-hat)
        logger.info("Running Relatedness check using PLINK 1.9 --genome (Pi-hat)")
        
        # Pruning
        prune_prefix = f"{self.out_prefix}_prune"
        cmd_prune = ['plink2', '--bfile', self.bed_prefix, '--snps-only', 'just-acgt', 
                     '--indep-pairwise', '50', '5', '0.2', '--maf', '0.05', 
                     '--out', prune_prefix, '--threads', self.args.threads]
        
        if os.path.exists(self.args.high_ld):
            cmd_prune.extend(['--exclude', 'range', self.args.high_ld])
        else:
            logger.warning(f"High LD file {self.args.high_ld} not found. Skipping LD exclusion.")
            
        run_command(cmd_prune)
        
        # Extract Pruned Variants to temp bed
        pruned_bed = f"{self.out_prefix}.pruned"
        run_command(['plink2', '--bfile', self.bed_prefix, '--extract', f"{prune_prefix}.prune.in", 
                     '--make-bed', '--out', pruned_bed, '--threads', self.args.threads])
        
        # Calc Genome (PLINK 1.9)
        # Note: plink 1.9 must be in path or aliased. In NF script we export path.
        run_command(['plink', '--bfile', pruned_bed, '--genome', '--out', self.out_prefix])

        # 5. PCA (reuse pruned dataset)
        run_command(['plink2', '--bfile', pruned_bed, '--pca', '10', '--out', self.out_prefix, '--threads', self.args.threads])

        # Cleanup temp pruned bed
        for ext in ['.bed', '.bim', '.fam']:
            if os.path.exists(pruned_bed + ext):
                os.remove(pruned_bed + ext)

    def aggregate_data(self, fam_df):
        def read_plink_out(suffix, cols):
            fname = f"{self.out_prefix}.{suffix}"
            if os.path.exists(fname):
                return pd.read_csv(fname, sep=r'\s+', usecols=cols)
            return pd.DataFrame()

        # Load Metrics
        smiss = read_plink_out("smiss", ['#IID', 'F_MISS'])
        het = read_plink_out("het", ['#IID', 'F'])
        sexcheck = read_plink_out("sexcheck", ['#IID', 'F'])
        
        # Merge (Use #IID as key)
        merged = fam_df.copy()
        # Ensure IID is str
        merged['IID'] = merged['IID'].astype(str)
        
        if not smiss.empty:
            smiss['#IID'] = smiss['#IID'].astype(str)
            merged = merged.merge(smiss, left_on='IID', right_on='#IID', how='left').drop(columns=['#IID'])
        
        if not het.empty:
            het['#IID'] = het['#IID'].astype(str)
            merged = merged.merge(het, left_on='IID', right_on='#IID', how='left', suffixes=('', '_HET'))
            # Rename F to F_HET if collision happened or if suffix didn't apply
            if 'F' in merged.columns and 'F_HET' not in merged.columns:
                merged.rename(columns={'F': 'F_HET'}, inplace=True)
            elif 'F' in merged.columns and 'F_HET' in merged.columns:
                 # If original F exists (unlikely from fam), drop it or rename
                 pass

        if not sexcheck.empty:
            sexcheck['#IID'] = sexcheck['#IID'].astype(str)
            merged = merged.merge(sexcheck, left_on='IID', right_on='#IID', how='left', suffixes=('', '_SEX'))
            if 'F' in merged.columns: # This would be F_SEX coming in as F
                merged.rename(columns={'F': 'F_SEX'}, inplace=True)
            if 'F_SEX' not in merged.columns and 'F' in sexcheck.columns:
                # If merge resulted in F_SEX suffix, good. If not (unlikely due to F_HET), handle it.
                # Actually earlier merge renamed F to F_HET. So check-sex F will come in as F unless suffixed
                pass 
                
        # Fix Column names just in case
        # PLINK2 headers: .het -> F, .sexcheck -> F.
        # We need to distinguish Heterozygosity F and Sex Check F
        # Re-doing logic cleanly:
        # 1. Smiss -> F_MISS
        # 2. Het -> F_HET
        # 3. Sex -> F_SEX
        # The aggregation above might be sloppy with suffixes. Let's force rename.
        pass

        # Robust Merge
        final_df = fam_df.copy()
        final_df['IID'] = final_df['IID'].astype(str)
        
        if not smiss.empty: 
            smiss['#IID'] = smiss['#IID'].astype(str)
            final_df = final_df.merge(smiss.rename(columns={'F_MISS': 'MISS'}), left_on='IID', right_on='#IID', how='left')
        
        if not het.empty:
            het['#IID'] = het['#IID'].astype(str)
            final_df = final_df.merge(het.rename(columns={'F': 'F_HET'}), left_on='IID', right_on='#IID', how='left')
            
        if not sexcheck.empty:
            sexcheck['#IID'] = sexcheck['#IID'].astype(str)
            final_df = final_df.merge(sexcheck.rename(columns={'F': 'F_SEX'}), left_on='IID', right_on='#IID', how='left')

        # Load Info
        try:
            info_df = pd.read_excel(self.args.info)
            info_df[self.args.id_col] = info_df[self.args.id_col].astype(str)
            
            # Merge
            final_df = final_df.merge(info_df, left_on='IID', right_on=self.args.id_col, how='left')
        except Exception as e:
            logger.warning(f"Could not load or merge Excel file: {e}")

        return final_df

    def plot_metrics(self, df):
        sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
        fig, axes = plt.subplots(2, 3, figsize=(24, 16))
        
        # Prepare Plot Data
        plot_df = df.copy()
        sex_col = self.args.sex_col if self.args.sex_col in plot_df.columns else 'SEX'
        
        # 1. Missingness vs Heterozygosity
        ax = axes[0,0]
        if 'MISS' in plot_df.columns and 'F_HET' in plot_df.columns:
            sns.scatterplot(data=plot_df, x='MISS', y='F_HET', hue='Group_FAM', style=sex_col, ax=ax, s=60, alpha=0.7)
            ax.set_title('Heterozygosity vs Missingness')
            ax.set_xlabel('Missingness Rate')
            ax.set_ylabel('Heterozygosity F')
        
        # 2. Sex Check
        ax = axes[0,1]
        if 'F_SEX' in plot_df.columns:
            sns.boxplot(data=plot_df, x=sex_col, y='F_SEX', hue='Group_FAM', ax=ax)
            ax.axhline(0.2, color='red', linestyle='--', label='Female < 0.2')
            ax.axhline(0.8, color='blue', linestyle='--', label='Male > 0.8')
            ax.set_title('Sex Check (X-Chr F-stat)')
            ax.legend(loc='upper right')
        
        # 3. Depth (Check simple distribution)
        ax = axes[0,2]
        dp_col = self.args.mean_dp_col
        if dp_col in plot_df.columns:
            sns.histplot(data=plot_df, x=dp_col, hue='Group_FAM', element="step", ax=ax)
            ax.set_title('Mean Depth Distribution')
        
        # 4. PCA
        ax = axes[1,0]
        pca_file = f"{self.out_prefix}.eigenvec"
        if os.path.exists(pca_file):
            pca = pd.read_csv(pca_file, sep=r'\s+')
            pca['#IID'] = pca['#IID'].astype(str)
            # Merge PCA back to plot_df for consistent Hue
            pca_plot = plot_df.merge(pca, left_on='IID', right_on='#IID', how='inner')
            if 'PC1' in pca_plot.columns and 'PC2' in pca_plot.columns:
                sns.scatterplot(data=pca_plot, x='PC1', y='PC2', hue='Group_FAM', ax=ax, s=60, alpha=0.7)
                ax.set_title('PCA (PC1 vs PC2)')
        
        # 5. Relatedness (PI_HAT from .genome)
        ax = axes[1,1]
        genome_file = f"{self.out_prefix}.genome"
        if os.path.exists(genome_file):
            rel_df = pd.read_csv(genome_file, sep=r'\s+')
            if not rel_df.empty and 'PI_HAT' in rel_df.columns:
                 # Filter to show only related pairs to avoid crowding with 0s? 
                 # Or just show all? Usually mostly 0s. 
                 # Let's show > 0.1 to be useful or just the distribution.
                 # The user logic was "refer to old method".
                 # Old calculator just calculated it. 
                 # Let's plot histogram of PI_HAT > 0.1
                 rel_high = rel_df[rel_df['PI_HAT'] > 0.1]
                 if not rel_high.empty:
                     sns.histplot(data=rel_high, x='PI_HAT', bins=20, ax=ax)
                     ax.set_title('Relatedness (PI_HAT > 0.1)')
                 else:
                     ax.text(0.5, 0.5, "No pairs with PI_HAT > 0.1", ha='center')
                     ax.set_title('Relatedness (Low IBD)')
        elif os.path.exists(f"{self.out_prefix}.kin0"):
             # Fallback if old file exists for some reason
            king = pd.read_csv(f"{self.out_prefix}.kin0", sep=r'\s+')
            if not king.empty:
                sns.histplot(data=king, x='KINSHIP', bins=50, ax=ax, log_scale=(False, True))
                ax.set_title('Kinship Coefficients (KING)')
        
        # 6. DP vs Missingness
        ax = axes[1,2]
        if dp_col in plot_df.columns and 'MISS' in plot_df.columns:
            sns.scatterplot(data=plot_df, x=dp_col, y='MISS', hue='Group_FAM', ax=ax)
            ax.set_title('Depth vs Missingness')

        plt.tight_layout()
        plt.savefig(f"{self.out_prefix}.sample_qc_summary.png", dpi=300)
        plt.savefig(f"{self.out_prefix}.sample_qc_summary.pdf")
        plt.close()

    def generate_flags(self, df):
        flags = []
        
        # Thresholds
        miss_thresh = 0.05
        
        # Calc Het Mean/Std from Controls if possible, else all
        controls = df[df['Group_FAM'] == 'Control']
        if len(controls) > 10:
            het_mean = controls['F_HET'].mean()
            het_std = controls['F_HET'].std()
        else:
            het_mean = df['F_HET'].mean()
            het_std = df['F_HET'].std()
            
        for idx, row in df.iterrows():
            issues = []
            
            # Missingness
            if pd.notna(row.get('MISS')) and row['MISS'] > miss_thresh:
                issues.append(f"HighMissingness({row['MISS']:.3f})")
            
            # Het
            if pd.notna(row.get('F_HET')):
                z_het = (row['F_HET'] - het_mean) / (het_std if het_std > 0 else 1)
                if abs(z_het) > 4: # 4 SD
                    issues.append(f"HetOutlier(Z={z_het:.1f})")
            
            # Sex
            # Check reported sex col
            sex_val = str(row.get(self.args.sex_col, '')).lower().strip()
            f_sex = row.get('F_SEX')
            if pd.notna(f_sex):
                if f_sex > 0.8: # Genotypic Male
                    if sex_val in ['f', 'female', '2']:
                        issues.append(f"SexMismatch(Gen=M,Rep={sex_val})")
                elif f_sex < 0.2: # Genotypic Female
                    if sex_val in ['m', 'male', '1']:
                        issues.append(f"SexMismatch(Gen=F,Rep={sex_val})")
            
            if issues:
                flags.append({
                    'IID': row['IID'],
                    'Issues': "; ".join(issues)
                })
        
        flag_df = pd.DataFrame(flags)
        flag_df.to_csv(f"{self.out_prefix}.sample_qc_flags.csv", index=False)
        
        # Save Full Summary
        df.to_csv(f"{self.out_prefix}.sample_qc_summary.csv", index=False)

def main():
    args = parse_arguments()
    
    logger.info(f"Starting QC for {args.bed}")
    
    # Check High LD file
    if not os.path.exists(args.high_ld):
        logger.warning(f"High LD file not found at {args.high_ld}")
    
    pipeline = SampleQC(args)
    pipeline.run_qc_steps()
    
    fam_df = load_fam(args.bed)
    final_df = pipeline.aggregate_data(fam_df)
    
    pipeline.plot_metrics(final_df)
    pipeline.generate_flags(final_df)
    
    logger.info("QC Completed.")

if __name__ == "__main__":
    main()
