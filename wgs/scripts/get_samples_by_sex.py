import pandas as pd
import argparse
import os
import sys

def main():
    parser = argparse.ArgumentParser(description="Extract sample IDs by sex from Excel.")
    parser.add_argument("--excel", required=True, help="Path to Excel file")
    parser.add_argument("--id-col", required=True, help="Column name for Sample ID")
    parser.add_argument("--sex-col", required=True, help="Column name for Sex")
    parser.add_argument("--male-out", required=True, help="Output file for Male IDs")
    parser.add_argument("--female-out", required=True, help="Output file for Female IDs")
    
    args = parser.parse_args()
    
    try:
        df = pd.read_excel(args.excel)
    except Exception as e:
        print(f"Error reading Excel file: {e}", file=sys.stderr)
        sys.exit(1)
    
    if args.id_col not in df.columns:
        print(f"Error: ID column '{args.id_col}' not found in Excel file.", file=sys.stderr)
        sys.exit(1)
        
    if args.sex_col not in df.columns:
        print(f"Error: Sex column '{args.sex_col}' not found in Excel file.", file=sys.stderr)
        sys.exit(1)

    # Filter and select
    # Assuming 'M' and 'F' are the values. Strip whitespace just in case.
    df[args.sex_col] = df[args.sex_col].astype(str).str.strip()
    
    # Drop rows where ID is NaN or empty
    df = df.dropna(subset=[args.id_col])
    df[args.id_col] = df[args.id_col].astype(str).str.strip()
    df = df[df[args.id_col] != '']
    
    males = df[df[args.sex_col] == 'M'][args.id_col]
    females = df[df[args.sex_col] == 'F'][args.id_col]
    
    # Write to files
    males.to_csv(args.male_out, index=False, header=False)
    females.to_csv(args.female_out, index=False, header=False)
    
    print(f"Extracted {len(males)} males and {len(females)} females.", file=sys.stderr)

if __name__ == "__main__":
    main()
