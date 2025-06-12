import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

# Ensure bcftools path is available
os.environ["PATH"] = "/home/b/b37974/:" + os.environ["PATH"]

# Input VCF file
vcf_file = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/08.variant_filter/chr21.pass.mac1.het.smiss.v_filtered.vcf.gz"

# Step 1: Extract GQ matrix using bcftools query
gq_outfile = "gq_matrix.txt"
# query_cmd = [
#     "bcftools", "query",
#     "-f", "[%GQ\t]\n",  # tab-delimited matrix
#     vcf_file
# ]

# with open(gq_outfile, "w") as f:
#     subprocess.run(query_cmd, stdout=f)

# Step 2: Process GQ matrix in chunks to reduce memory usage
from collections import Counter

gq_counter = Counter()

with open(gq_outfile, "r") as f:
    for line in f:
        values = [v for v in line.strip().split("\t") if v != "."]
        chunk_values = pd.to_numeric(pd.Series(values), errors="coerce").dropna().astype(int)
        gq_counter.update(chunk_values)

# Step 3: Plot histogram
plt.figure(figsize=(10, 6))
plt.bar(gq_counter.keys(), gq_counter.values(), width=1.0, edgecolor='black')
plt.title("Distribution of FORMAT/GQ")
plt.xlabel("GQ")
plt.ylabel("Frequency")
plt.grid(True)
plt.tight_layout()

# Save plot to PDF
plt.savefig("gq_histogram.pdf")