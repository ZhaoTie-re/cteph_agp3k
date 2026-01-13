#!/bin/bash

# ==============================================================================
# Script Name: plot_minac_dist.sh
# Description: Efficiently calculates and plots the distribution of MinAC 
#              from a large variant_metrics.txt.gz file.
#              Uses awk for memory-efficient counting and Python for plotting.
#
# Usage: ./plot_minac_dist.sh <input_file.txt.gz> <output_prefix>
#
# Example:
#   ./plot_minac_dist.sh ../results/00.qc_metrics/minac0/variant_metrics.txt.gz minac_dist_plot
# ==============================================================================

set -e

INPUT_FILE="$1"
OUTPUT_PREFIX="$2"

if [[ -z "$INPUT_FILE" || -z "$OUTPUT_PREFIX" ]]; then
    echo "Usage: $0 <input_file.txt.gz> <output_prefix>"
    exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' not found."
    exit 1
fi

echo "Processing $INPUT_FILE..."
echo "Step 1: Counting MinAC frequencies using awk..."

# Create a temporary file for the frequency data
FREQ_FILE="${OUTPUT_PREFIX}.counts.tsv"

# 1. zcat streams the compressed file
# 2. awk skips header (NR>1), increments count for the MinAC column ($8), then prints the map
#    MinAC is column 8 based on the header provided: #CHROM POS VariantID REF ALT RefAC AltAC MinAC ...
zcat "$INPUT_FILE" | awk '
BEGIN { OFS="\t" }
NR==1 {
    # Verify column index just in case, or default to 8
    for(i=1;i<=NF;i++) if($i=="MinAC") col=i;
    if(col=="") col=8; 
}
NR>1 { 
    count[$col]++; 
} 
END { 
    print "MinAC", "Count"
    for (val in count) print val, count[val] 
}' > "$FREQ_FILE"

echo "Frequency data saved to $FREQ_FILE"

echo "Step 2: Plotting using Python..."

# Create a temporary python script
PY_SCRIPT="${OUTPUT_PREFIX}.plot.py"

cat <<EOF > "$PY_SCRIPT"
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# Set style
sns.set_theme(style="ticks", context="paper", font_scale=1.2)

input_file = "${FREQ_FILE}"
output_prefix = "${OUTPUT_PREFIX}"

try:
    df = pd.read_csv(input_file, sep='\t')
    df = df.sort_values('MinAC')
    
    total_variants = df['Count'].sum()
    print(f"Total Variants: {total_variants:,}")

    # Plot 1: Linear Scale (MinAC vs Count)
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(df['MinAC'], df['Count'], color='#4c72b0', alpha=0.9, width=0.8)
    ax.set_xlabel("Minor Allele Count (MinAC)")
    ax.set_ylabel("Count")
    ax.set_title(f"MinAC Count Distribution (N={total_variants:,})")
    ax.grid(True, alpha=0.3)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}.dist_linear.pdf")
    plt.close()

    # Plot 2: Log Scale Y (MinAC vs Log Count)
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(df['MinAC'], df['Count'], color='#c44e52', alpha=0.9, width=0.8)
    ax.set_yscale('log')
    ax.set_xlabel("Minor Allele Count (MinAC)")
    ax.set_ylabel("Count (Log Scale)")
    ax.set_title("MinAC Count Distribution (Log Scale)")
    ax.grid(True, which="both", alpha=0.2)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}.dist_log.pdf")
    plt.close()
    
    print(f"Plots saved to {output_prefix}.dist_*.pdf")

except Exception as e:
    print(f"Error in plotting: {e}", file=sys.stderr)
    sys.exit(1)
EOF

python3 "$PY_SCRIPT"

# Cleanup
rm "$PY_SCRIPT"
# We keep the counts file as it might be useful
echo "Done."
