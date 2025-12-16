#!/usr/bin/env python3
"""
Update PSAM file with sex and phenotype information from sample metadata.

This script updates a PLINK2 .psam file by:
1. Reading sample information from an Excel file
2. Updating SEX column (1=male, 2=female, 0=unknown)
3. Updating or creating phenotype column (2=case, 1=control, 0=missing)
4. Ensuring proper PLINK2 .psam format compliance
"""

import argparse
import pandas as pd
import sys


def convert_sex(iid, sex_map):
    """Convert sex value to PLINK2 encoding."""
    sex = sex_map.get(iid, 0)
    if pd.isna(sex):
        return 0
    sex_str = str(sex).upper()
    if sex_str in ['M', 'MALE', '1']:
        return 1
    elif sex_str in ['F', 'FEMALE', '2']:
        return 2
    else:
        return 0


def convert_pheno(iid, outcome_map, case_value):
    """Convert phenotype/outcome to PLINK2 encoding."""
    outcome = outcome_map.get(iid, 0)
    if pd.isna(outcome):
        return 0
    if str(outcome) == case_value:
        return 2
    else:
        return 1


def update_psam(psam_file, sample_info_file, id_col, sex_col, outcome_col, case_value, output_file):
    """Main function to update PSAM file."""
    
    # Read sample info from Excel
    sample_info = pd.read_excel(sample_info_file)
    
    # Read PSAM file (comment=None to preserve # in column names)
    psam = pd.read_csv(psam_file, sep='\t', comment=None)
    
    # Create mapping dictionaries
    sex_map = dict(zip(sample_info[id_col], sample_info[sex_col]))
    outcome_map = dict(zip(sample_info[id_col], sample_info[outcome_col]))
    
    # Detect FID and IID columns
    fid_col = None
    iid_col = None
    for col in psam.columns:
        if col in ['#FID', 'FID']:
            fid_col = col
        if col in ['#IID', 'IID']:
            iid_col = col
    
    if iid_col is None:
        raise ValueError("Cannot find IID column in PSAM file")
    
    # Extract IID values for mapping
    iid_values = psam[iid_col]
    
    # Ensure PAT and MAT columns exist
    if 'PAT' not in psam.columns:
        psam['PAT'] = '0'
    if 'MAT' not in psam.columns:
        psam['MAT'] = '0'
    
    # Update SEX column
    psam['SEX'] = iid_values.apply(lambda x: convert_sex(x, sex_map))
    
    # Handle phenotype column
    pheno_cols = [col for col in psam.columns if col.startswith('PHENO')]
    if pheno_cols:
        psam[pheno_cols[0]] = iid_values.apply(lambda x: convert_pheno(x, outcome_map, case_value))
    else:
        psam['PHENO1'] = iid_values.apply(lambda x: convert_pheno(x, outcome_map, case_value))
    
    # Reorder columns according to PLINK2 standard format
    if fid_col:
        standard_cols = [fid_col, iid_col, 'PAT', 'MAT', 'SEX']
    else:
        standard_cols = [iid_col, 'PAT', 'MAT', 'SEX']
    
    other_cols = [col for col in psam.columns if col not in standard_cols]
    final_cols = standard_cols + other_cols
    psam = psam[final_cols]
    
    # Write updated PSAM file
    psam.to_csv(output_file, sep='\t', index=False)
    print(f"Updated PSAM file written to: {output_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Update PLINK2 .psam file with sex and phenotype information',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--psam', required=True,
                        help='Input PSAM file path')
    parser.add_argument('--sample-info', required=True,
                        help='Sample information Excel file')
    parser.add_argument('--id-col', required=True,
                        help='Column name for sample IDs in Excel file')
    parser.add_argument('--sex-col', required=True,
                        help='Column name for sex information in Excel file')
    parser.add_argument('--outcome-col', required=True,
                        help='Column name for outcome/phenotype in Excel file')
    parser.add_argument('--case-value', required=True,
                        help='Value indicating case status in outcome column')
    parser.add_argument('--output', required=True,
                        help='Output PSAM file path')
    
    args = parser.parse_args()
    
    update_psam(
        args.psam,
        args.sample_info,
        args.id_col,
        args.sex_col,
        args.outcome_col,
        args.case_value,
        args.output
    )


if __name__ == '__main__':
    main()
