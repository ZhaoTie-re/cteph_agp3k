import pandas as pd
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Merge phenotype and covariate files.')
    parser.add_argument('--pheno', required=True, help='Phenotype file path')
    parser.add_argument('--cov', required=True, help='Covariate file path')
    parser.add_argument('--pheno_col', required=True, help='Phenotype column name')
    parser.add_argument('--cov_list', required=True, help='Comma-separated list of covariates')
    parser.add_argument('--sex_col', required=False, help='Sex column name (optional)')
    parser.add_argument('--out', required=True, help='Output file path')
    
    args = parser.parse_args()

    pheno_file = args.pheno
    cov_file = args.cov
    pheno_col = args.pheno_col
    cov_list_str = args.cov_list
    sex_col = args.sex_col
    output_file = args.out

    # Load data
    try:
        # Try reading with tab/whitespace separator since files seem to be PLINK-style despite .csv extension
        p = pd.read_csv(pheno_file, sep=r'\s+', engine='python')
        c = pd.read_csv(cov_file, sep=r'\s+', engine='python')
        
        # Clean up header if it starts with #
        p.rename(columns={'#FID': 'FID'}, inplace=True)
        c.rename(columns={'#FID': 'FID'}, inplace=True)
        
        # Check and recode phenotype if necessary (PLINK 1/2 -> SAIGE 0/1)
        if pheno_col in p.columns:
            unique_vals = set(p[pheno_col].dropna().unique())
            if unique_vals == {1, 2} or unique_vals == {1.0, 2.0}:
                print(f"Recoding phenotype '{pheno_col}' from 1/2 to 0/1 for SAIGE.")
                p[pheno_col] = p[pheno_col] - 1
        else:
             print(f"Warning: Phenotype column '{pheno_col}' not found in phenotype file.")

        # Check and recode sex if necessary (1/2 -> 0/1 for SAIGE: 1->0 Male, 2->1 Female)
        if sex_col and sex_col in c.columns:
             unique_sex = set(c[sex_col].dropna().unique())
             if unique_sex <= {1, 2, 1.0, 2.0}:
                  print(f"Recoding sex '{sex_col}' from 1/2 to 0/1 for SAIGE (Male=0, Female=1).")
                  c[sex_col] = c[sex_col] - 1
    except Exception as e:
        print(f"Error reading input files: {e}")
        sys.exit(1)

    # Parse covariates
    cov_cols = [x.strip() for x in cov_list_str.split(',')]
    if 'IID' not in cov_cols:
        cov_cols = ['IID'] + cov_cols

    # Check if columns exist
    missing_cols = [col for col in cov_cols if col not in c.columns]
    if missing_cols:
        print(f"Error: Missing columns in covariate file: {missing_cols}")
        sys.exit(1)

    # Select columns and merge
    c_sub = c[cov_cols]
    merged = pd.merge(p, c_sub, on='IID', how='inner')

    # Write output
    merged.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()
