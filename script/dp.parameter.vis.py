# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Visualize the number of SNPs that meet the specified depth of coverage (DP) threshold for given sample codes")
parser.add_argument("--chrom", type=str, help="Chromosome number")
parser.add_argument("--count_file", type=str, help="Path to the count file")
parser.add_argument("--dp_summary", type=str, help="Path to the DP summary file")
parser.add_argument("--output_pdf", type=str, help="Output file name")
args = parser.parse_args()

chrom = args.chrom
chrom_count_file = args.count_file

with open(chrom_count_file, 'r') as file:
    for line in file:
        if line.startswith('Number of sites:'):
            num_sites = int(line.split(':')[1].strip())
            break

dp_summary = pd.read_csv(args.dp_summary, index_col=0)

dp_summary_transposed = dp_summary.T
dp_summary_transposed.plot(kind='line', figsize=(10, 6))
plt.axhline(y=num_sites, color='red', linestyle='--')
plt.style.use('default')
plt.title(f'Remaining Variants Across DP Thresholds ({chrom})')
plt.xlabel('DP Threshold')
plt.ylabel('Number of Variants')
plt.legend()
plt.grid(axis='x', linestyle='--', color='gray', alpha=0.7)
plt.xticks(ticks=range(len(dp_summary_transposed.index)), labels=dp_summary_transposed.index, rotation=45)
# plt.show()
output_pdf = args.output_pdf
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
print(f"Plot saved: {output_pdf}")
