import pandas as pd
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Generate PLINK sample info file from Excel.')
    parser.add_argument('--excel', required=True, help='Input Excel file path')
    parser.add_argument('--id-col', required=True, help='Column name for Sample ID')
    parser.add_argument('--sex-col', required=True, help='Column name for Sex')
    parser.add_argument('--pheno-col', required=True, help='Column name for Phenotype/Group')
    parser.add_argument('--out', required=True, help='Output file path')
    
    args = parser.parse_args()

    try:
        # Read Excel file
        df = pd.read_excel(args.excel)
        
        # Check if columns exist
        for col in [args.id_col, args.sex_col, args.pheno_col]:
            if col not in df.columns:
                print(f"Error: Column '{col}' not found in Excel file.", file=sys.stderr)
                sys.exit(1)

        # Select and rename columns
        df = df[[args.id_col, args.sex_col, args.pheno_col]].copy()
        df.columns = ['IID', 'Sex', 'Pheno']

        # Drop rows with missing IDs
        df = df.dropna(subset=['IID'])

        # Map Sex: M -> 1 (Male), F -> 2 (Female)
        # Handle case sensitivity and whitespace
        def map_sex(x):
            s = str(x).strip().upper()
            if s == 'M' or s == 'MALE': return 1
            if s == 'F' or s == 'FEMALE': return 2
            return 0 # Unknown
            
        df['SEX'] = df['Sex'].apply(map_sex)

        # Map Phenotype: PH -> 2 (Case), AGP3K -> 1 (Control)
        def map_pheno(x):
            s = str(x).strip()
            if s == 'PH': return 2
            if s == 'AGP3K': return 1
            return -9 # Missing
            
        df['PHENO'] = df['Pheno'].apply(map_pheno)

        # Set FID = IID (for --double-id)
        df['FID'] = df['IID']

        # Write to tab-separated file: FID IID SEX PHENO
        df[['FID', 'IID', 'SEX', 'PHENO']].to_csv(args.out, sep='\t', index=False, header=False)
        print(f"Successfully wrote {len(df)} samples to {args.out}")

    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
