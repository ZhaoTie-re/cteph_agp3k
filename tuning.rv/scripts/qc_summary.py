import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
import argparse
import gzip
import sys
import subprocess

def count_lines_fast(fs_path):
    """
    Count lines in a gzipped file using zcat | wc -l which is 
    much faster than Python iteration for large files.
    """
    try:
        # Use zcat -f (force) to handle both compressed and uncompressed
        # On some systems zcat expects .Z, so use gzip -dc or zcat
        # wc -l counts newlines. A file with 1 line and no newline at EOF might count 0.
        # But Nextflow/Plink outputs are standard.
        cmd = f"gzip -dc {fs_path} | wc -l"
        result = subprocess.check_output(cmd, shell=True)
        return int(result.strip())
    except Exception as e:
        print(f"Warning: Fast line count failed ({e}), falling back to slow method.", file=sys.stderr)
        count = 0
        with gzip.open(fs_path, 'rt') as f:
            for i, _ in enumerate(f):
                count = i + 1 # count variants = lines
        return count

def analyze_qc(args):
    try:
        # 1. Variant Count
        # Count lines of compressed file, subtract 1 for header
        total_lines = count_lines_fast(args.variant_metrics)
        variant_count = max(0, total_lines - 1)
        
        # 2. Load Sample Metrics
        # pandas reads gzip automatically if filename ends in .gz
        df = pd.read_csv(args.sample_metrics, sep='\t')

        # Data Preprocessing for Regression
        # Filter for valid Group and Regression cols
        # Group Mapping: Case -> 1, Control -> 0, others -> drop/ignore for regression
        
        reg_df = df.copy()
        reg_df = reg_df[reg_df['Group'].isin(['Case', 'Control'])]
        reg_df['Group_Bin'] = reg_df['Group'].map({'Case': 1, 'Control': 0})
        
        # MeanDP Z-scale
        # Calculating Z-score on the subset or whole? user said "MeanDP要进行Z scale", usually implies the dataset used for regression.
        # But to be robust, let's use the distribution of the regression dataset.
        if len(reg_df) > 0:
            mu = reg_df['MeanDP'].mean()
            sigma = reg_df['MeanDP'].std()
            if sigma == 0:
                reg_df['MeanDP_Z'] = 0
            else:
                reg_df['MeanDP_Z'] = (reg_df['MeanDP'] - mu) / sigma
        else:
            reg_df['MeanDP_Z'] = np.nan

        # 3. Spearman Correlation (SMinAC vs MeanDP)
        # Using raw values from the full dataframe (excluding missing)
        corr_df = df.dropna(subset=['SMinAC', 'MeanDP'])
        if len(corr_df) > 1:
            spearman_rho, spearman_p = stats.spearmanr(corr_df['SMinAC'], corr_df['MeanDP'])
        else:
            spearman_rho, spearman_p = np.nan, np.nan

        # 4. Poisson Regression (SMinAC ~ Group + MeanDP_Z)
        beta, ci_lower, ci_upper, p_val = np.nan, np.nan, np.nan, np.nan
        
        reg_valid = reg_df.dropna(subset=['SMinAC', 'Group_Bin', 'MeanDP_Z'])
        if len(reg_valid) > 1:
            try:
                # Formula: SMinAC ~ Group_Bin + MeanDP_Z
                model = smf.glm(formula="SMinAC ~ Group_Bin + MeanDP_Z", 
                                data=reg_valid, 
                                family=sm.families.Poisson())
                result = model.fit()
                
                # Extract MeanDP_Z stats
                params = result.params
                conf = result.conf_int()
                pvalues = result.pvalues
                
                if 'MeanDP_Z' in params:
                    beta = params['MeanDP_Z']
                    ci_lower = conf.loc['MeanDP_Z', 0]
                    ci_upper = conf.loc['MeanDP_Z', 1]
                    p_val = pvalues['MeanDP_Z']
            except Exception as e:
                print(f"Warning: Poisson regression failed: {e}", file=sys.stderr)

        # 5. TargetDP Analysis (MWU & Wasserstein)
        # Identify groups. Expected 15x and 30x but handle generically.
        # We take the top 2 sorted unique values.
        # Convert to string to ensure consistent grouping
        df['TargetDP_Str'] = df['TargetDP'].astype(str)
        groups = df['TargetDP_Str'].unique()
        
        mwu_p = np.nan
        wass_dist = np.nan
        
        if len(groups) >= 2:
            # Sort to be deterministic (e.g. 15 < 30)
            groups = sorted(groups)
            # Pick first two
            g1_label = groups[0]
            g2_label = groups[1]
            
            vals1 = df[df['TargetDP_Str'] == g1_label]['MeanDP'].dropna().values
            vals2 = df[df['TargetDP_Str'] == g2_label]['MeanDP'].dropna().values
            
            if len(vals1) > 0 and len(vals2) > 0:
                try:
                    _, mwu_p = stats.mannwhitneyu(vals1, vals2, alternative='two-sided')
                except:
                    pass
                
                try:
                    wass_dist = stats.wasserstein_distance(vals1, vals2)
                except:
                    pass

        # 6. Write TSV
        cols = [
            "MinAC_Threshold", 
            "Variant_Count", 
            "Sample_Count",
            "Spearman_Rho_SMinAC_MeanDP", 
            "Spearman_P_SMinAC_MeanDP",
            "Poisson_MeanDP_Beta", 
            "Poisson_MeanDP_95CI_Lower", 
            "Poisson_MeanDP_95CI_Upper", 
            "Poisson_MeanDP_P",
            "MWU_P_TargetDP", 
            "Wasserstein_Dist_TargetDP"
        ]
        
        with open(args.out_tsv, 'w') as f:
            f.write("\t".join(cols) + "\n")
            vals = [
                str(args.min_ac),
                str(variant_count),
                str(len(df)),
                f"{spearman_rho:.6g}", f"{spearman_p:.6g}",
                f"{beta:.6g}", f"{ci_lower:.6g}", f"{ci_upper:.6g}", f"{p_val:.6g}",
                f"{mwu_p:.6g}", f"{wass_dist:.6g}"
            ]
            f.write("\t".join(vals) + "\n")
            
        print(f"Generated stats: {args.out_tsv}")

        # 7. Write MD
        with open(args.out_md, 'w') as f:
            f.write("# Quality Control Metrics Methodology and Interpretation\n\n")
            
            f.write("## 1. Basic Statistics\n")
            f.write(f"- **MinAC Threshold ($T_{{min\\_ac}}$)**: {args.min_ac}. The minimum allele count threshold applied for variant filtering in the current workflow.\n")
            f.write("- **Variant Count ($N_{var}$)**: The total number of variants retained in the dataset (calculated as total lines minus the header).\n\n")
            
            f.write("## 2. Sequencing Depth Bias Assessment\n")
            f.write("To evaluate whether sequencing depth systematically influences the detection of minor alleles (genotyping bias), we employ two statistical approaches.\n\n")
            
            f.write("### 2.1 Spearman's Rank Correlation\n")
            f.write("We assess the monotonic relationship between the **Total Burden of Minor Alleles ($S_{MinAC}$)** and **Mean Sequencing Depth ($D_{mean}$)** across all samples.\n\n")
            f.write("- **Metrics**: `Spearman_Rho_SMinAC_MeanDP` ($\\rho$), `Spearman_P_SMinAC_MeanDP` ($p$)\n")
            f.write("- **Hypothesis**: $H_0: \\rho = 0$ (No monotonic association) vs $H_1: \\rho \\neq 0$.\n")
            f.write("- **Interpretation**: A statistically significant positive $\\rho$ suggests that samples with higher sequencing depth tend to have higher called minor allele burdens, indicating potential depth-dependent sensitivity bias.\n\n")
            
            f.write("### 2.2 Poisson Regression Model\n")
            f.write("We model the count of minor alleles using a Generalized Linear Model (GLM) with a Poisson link function, adjusting for biological group covariance (Case/Control status).\n\n")
            f.write("- **Model Specification**:\n")
            f.write("  $$ \\ln(E[S_{MinAC}]) = \\beta_0 + \\beta_{group} \\cdot X_{group} + \\beta_{depth} \\cdot Z(D_{mean}) $$\n")
            f.write("  Where:\n")
            f.write("  - $X_{group}$: Indicator variable (1 for Case, 0 for Control).\n")
            f.write("  - $Z(D_{mean})$: Z-score standardized mean depth ($ \\frac{D - \\mu_D}{\\sigma_D} $).\n")
            f.write("- **Metrics**: `Poisson_MeanDP_Beta` ($\\beta_{depth}$), `95CI` (Confidence Interval), `P`.\n")
            f.write("- **Interpretation**: $\\beta_{depth}$ represents the log-count change in minor alleles per standard deviation increase in sequencing depth, holding phenotype constant. A value significantly different from 0 indicates an independent depth effect.\n\n")
            
            f.write("## 3. Depth Distribution Confounding Analysis\n")
            f.write("To ensure that sequencing depth is not confounded with the target depth design (e.g., 15x vs 30x batches), we compare the empirical distributions of observed $D_{mean}$ stratified by `TargetDP`.\n\n")
            
            f.write("### 3.1 Mann-Whitney U Test\n")
            f.write("- **Metric**: `MWU_P_TargetDP` ($p$)\n")
            f.write("- **Interpretation**: Non-parametric test for location shift. A significant p-value indicates the median depths differ between target groups.\n\n")
            
            f.write("### 3.2 Wasserstein Distance (Earth Mover's Distance)\n")
            f.write("- **Metric**: `Wasserstein_Dist_TargetDP` ($W_1$)\n")
            f.write("- **Intuitive Concept**: Often visualized as the \"Earth Mover's Distance\". If we treat the two depth distributions as piles of dirt, $W_1$ is the minimum average \"work\" (mass $\\times$ distance) required to shovel one pile into the shape of the other.\n")
            f.write("- **Algorithmic Detail**: For 1D data like sequencing depth, we avoid complex optimization by exploiting the closed-form dual representation involving Cumulative Distribution Functions (CDFs):\n")
            f.write("  $$ W_1(U, V) = \\int_{-\\infty}^{\\infty} |CDF_U(x) - CDF_V(x)| dx $$\n")
            f.write("  This measures the total area between the two cumulative frequency curves.\n")
            f.write("- **Interpretation**: It provides a dissimilarity score in the original units (coverage $x$). Unlike the P-value (which assesses significance), $W_1$ assesses the *magnitude* of the shift.\n\n")
            
            f.write("---\n\n")
            
            f.write("## 4. Example Interpretation\n\n")
            f.write("The following example illustrates the interpretation of these metrics using representative data:\n\n")
            
            f.write("| Category | Metric | Value |\n")
            f.write("| :--- | :--- | :--- |\n")
            f.write(f"| **Basic** | MinAC_Threshold | 0 |\n")
            f.write(f"| | Variant_Count | 20,838,333 |\n")
            f.write(f"| **Bias** | Spearman $\\rho$ | 0.339 |\n")
            f.write(f"| | Spearman P-value | $5.61 \\times 10^{-82}$ |\n")
            f.write(f"| | Poisson $\\beta_{{depth}}$ | 0.014 (95% CI: 0.0137 - 0.0142) |\n")
            f.write(f"| | Poisson P-value | 0.00 |\n")
            f.write(f"| **Distribution** | MWU P-value | $2.77 \\times 10^{-189}$ |\n")
            f.write(f"| | Wasserstein Dist | 8.50 |\n\n")
            
            f.write("**Detailed Interpretation:**\n\n")
            f.write("1.  **Sequencing Bias**: \n")
            f.write("    - The **Spearman correlation** ($\\rho=0.34, p < 10^{-80}$) indicates a moderate, statistically significant positive monotonic relationship: samples with deeper sequencing generally have higher minor allele burdens.\n")
            f.write("    - The **Poisson Regression** confirms this effect is independent of Case/Control status. The $\\beta_{depth} = 0.014$ suggests that for every **1 Standard Deviation increase** in sequencing depth, the expected minor allele count increases by approximately $e^{0.014} - 1 \\approx 1.4\\%$. While the effect size appears small, the exceedingly low p-values indicate a systematic bias across the dataset.\n\n")
            f.write("2.  **Batch/Design Confounding**: \n")
            f.write("    - The **Mann-Whitney U test** p-value ($< 10^{-180}$) confirms that the depth distributions between the TargetDP groups are statistically distinct.\n")
            f.write("    - The **Wasserstein Distance (8.50)**, also known as **Earth Mover's Distance**, quantifies the effort to align the distributions.\n")
            f.write("      **Intuitive Interpretation**: Imagine the depth distributions of the two groups as **piles of earth**. This metric measures the minimum \"work\" required to **push** the earth from one pile to reshape it into the other.\n")
            f.write("      **In Our Example**: A value of **8.50** means that the two sequencing batches are separated by an average of **8.5x coverage**. To make the two groups identical in depth profile, you would effectively need to \"shovel\" **8.5x worth of coverage depth** from the higher-depth group to fill in the lower-depth group. This confirms the batches are visibly and substantially different.\n")
            
        print(f"Generated documentation: {args.out_md}")

    except Exception as e:
        print(f"Error in analyze_qc: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample-metrics', required=True)
    parser.add_argument('--variant-metrics', required=True)
    parser.add_argument('--min-ac', required=True)
    parser.add_argument('--out-tsv', required=True)
    parser.add_argument('--out-md', required=True)
    args = parser.parse_args()
    
    analyze_qc(args)

if __name__ == "__main__":
    main()
