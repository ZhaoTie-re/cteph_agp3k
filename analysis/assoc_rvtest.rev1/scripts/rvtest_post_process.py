import pandas as pd
import argparse
import sys
import re
from collections import Counter
from statsmodels.stats.multitest import multipletests

def parse_genomic_range(range_str):
    """
    Parses RANGE string to extract CHR, START, END.
    Handles multiple ranges and prioritizes regular chromosomes.
    """
    if pd.isna(range_str) or str(range_str).strip() == "":
        return None, None, None

    # Split by comma
    ranges = [r.strip() for r in str(range_str).split(',')]
    parsed = []

    # Regex for CHR:START-END
    pattern = re.compile(r'(.+):(\d+)-(\d+)')

    for r in ranges:
        m = pattern.match(r)
        if m:
            chrom = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3))
            parsed.append((chrom, start, end))

    if not parsed:
        return None, None, None

    # Regular chromosome pattern: 1-22, X, Y, M, MT (with or without chr)
    reg_pattern = re.compile(r'^(chr)?(\d{1,2}|[XYM]|MT)$', re.IGNORECASE)

    regular_parsed = [p for p in parsed if reg_pattern.match(p[0])]

    # Strategy: Use regular chromosomes if available
    if regular_parsed:
        target_list = regular_parsed
    else:
        target_list = parsed

    # Determine primary chromosome (majority vote if mixed, though unlikely for one gene)
    chrom_counts = Counter(p[0] for p in target_list)
    primary_chrom = chrom_counts.most_common(1)[0][0]

    # Filter for segments on primary chromosome
    final_segments = [p for p in target_list if p[0] == primary_chrom]

    min_start = min(p[1] for p in final_segments)
    max_end = max(p[2] for p in final_segments)

    return primary_chrom, min_start, max_end

def main():
    parser = argparse.ArgumentParser(description="Post-process rvtest results: filter by NumVar and calculate FDR.")
    parser.add_argument("--input", required=True, help="Input .assoc or .out file")
    parser.add_argument("--output", required=True, help="Output file")
    parser.add_argument("--num-var-threshold", type=int, required=True, help="Minimum NumVar threshold")
    
    args = parser.parse_args()
    
    input_path = args.input
    output_path = args.output
    threshold = args.num_var_threshold
    
    print(f"Processing {input_path} with NumVar threshold >= {threshold}")

    # Try reading as tab-separated or whitespace-separated
    try:
        # First attempt: assume tab-separated
        df = pd.read_csv(input_path, sep='\t')
        if df.shape[1] < 2:
            # Fallback: whitespace-separated
            df = pd.read_csv(input_path, delim_whitespace=True)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    print(f"Original columns: {df.columns.tolist()}")

    # Identify columns flexibly
    col_map = {c.lower(): c for c in df.columns}
    
    # Find P-value column
    p_col_candidates = ['pvalue', 'p_value', 'p.value', 'p']
    p_col = next((col_map[c] for c in p_col_candidates if c in col_map), None)
    
    # Find NumVar column
    # Common names: NumVar, N_VAR, N_MARKER, etc. rvtests usually uses 'NumVar'
    num_var_candidates = ['numvar', 'n_var', 'num_variant', 'n_marker']
    num_var_col = next((col_map[c] for c in num_var_candidates if c in col_map), None)

    if not p_col:
        print("Error: Could not find P-value column.")
        sys.exit(1)
        
    if not num_var_col:
        print(f"Error: Could not find NumVar column. Available: {df.columns.tolist()}")
        print("Cannot filter by NumVar. Exiting.")
        sys.exit(1)
    
    print(f"Using P-value column: '{p_col}' and NumVar column: '{num_var_col}'")

    # Parse RANGE column if exists
    range_candidates = ['range']
    range_col = next((col_map[c] for c in range_candidates if c in col_map), None)

    if range_col:
        print(f"Parsing genomic range from column: '{range_col}'")
        # Apply parsing
        parsed_data = df[range_col].apply(parse_genomic_range)
        
        # Create separate columns from the tuples
        # Using list comprehension to unzip
        parsed_list = parsed_data.tolist()
        chr_list = [x[0] for x in parsed_list]
        start_list = [x[1] for x in parsed_list]
        end_list = [x[2] for x in parsed_list]

        # Insert columns after RANGE
        loc_index = df.columns.get_loc(range_col) + 1
        df.insert(loc_index, 'CHR', chr_list)
        df.insert(loc_index + 1, 'START', start_list)
        df.insert(loc_index + 2, 'END', end_list)
    else:
        print("Warning: Range column not found. Skipping CHR/START/END extraction.")

    # Filter
    n_original = len(df)
    df_filtered = df[df[num_var_col] >= threshold].copy()
    n_filtered = len(df_filtered)
    print(f"Filtered rows: {n_original} -> {n_filtered} (removed {n_original - n_filtered} rows)")

    # Calculate FDR (Benjamini-Hochberg)
    # Remove NaNs from P-value for calculation
    df_valid_p = df_filtered.dropna(subset=[p_col])
    
    if not df_valid_p.empty:
        pvals = df_valid_p[p_col].values
        reject, qvals, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
        
        # Assign back aligning with index
        df_filtered.loc[df_valid_p.index, 'FDR'] = qvals
    else:
        df_filtered['FDR'] = pd.Series(dtype=float)

    # Sort by P-value for better readability
    df_filtered.sort_values(by=p_col, inplace=True)

    # Write output
    df_filtered.to_csv(output_path, sep='\t', index=False, float_format='%.6g')
    print(f"Written result to {output_path}")

if __name__ == "__main__":
    main()
