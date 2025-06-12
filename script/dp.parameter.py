# %%
import pandas as pd
import numpy as np
import pysam
from concurrent.futures import ProcessPoolExecutor
import argparse

parser = argparse.ArgumentParser(description="Check the number of SNPs that meet the specified depth of coverage (DP) threshold for given sample codes")
parser.add_argument("--chrom", type=str, help="Chromosome number")
parser.add_argument("--vcf_path", type=str, help="Path to the VCF file")
parser.add_argument("--metadata", type=str, help="Path to the metadata file")
parser.add_argument("--output_csv", type=str, help="Output CSV file name")
args = parser.parse_args()

chrom = args.chrom
vcf_path = args.vcf_path

# %%
# read sample information from the metadata file
cteph_jhrp4 = pd.read_csv(args.metadata, sep=",")
cteph_30x = set(cteph_jhrp4[cteph_jhrp4["target_depth"] == "30x"]["sample_code"].values)
cteph_15x = set(cteph_jhrp4[cteph_jhrp4["target_depth"] == "15x"]["sample_code"].values)

def process_vcf(vcf_path, dp_thresholds, sample_codes):
    """
    Processes a VCF file to count the number of records that meet specified 
    depth of coverage (DP) thresholds for given sample codes.

    Args:
        vcf_path (str): Path to the input VCF file.
        dp_thresholds (list of int): A list of DP threshold values to evaluate.
        sample_codes (list of str): A list of sample codes to consider when 
            evaluating DP values.

    Returns:
        dict: A dictionary where keys are DP thresholds and values are the 
        counts of records that meet or exceed the corresponding DP threshold 
        for at least one of the specified sample codes.
    """
    vcf_in = pysam.VariantFile(vcf_path, 'r')  # open the VCF file for reading (r)
    counts = {dp: 0 for dp in dp_thresholds}  # initialize a dictionary to store the counts for each DP threshold

    for record in vcf_in.fetch():
        for dp in dp_thresholds:
            valid_sample = False
            for sample in record.samples.keys():
                if sample in sample_codes:
                    dp_value = record.samples[sample].get('DP', None)
                    if dp_value is not None and dp_value >= dp:
                        valid_sample = True
                        break
            if valid_sample:
                counts[dp] += 1

    vcf_in.close()
    return counts

def run_parallel_processing(vcf_path, dp_thresholds, sample_sets):
    """
    Executes parallel processing of VCF files for multiple sample sets using a process pool.

    Args:
        vcf_path (str): The file path to the VCF file to be processed.
        dp_thresholds (list): A list of depth thresholds to be applied during processing.
        sample_sets (list): A list of sample sets, where each sample set is processed independently.

    Returns:
        list: A list of results from processing each sample set. Each result corresponds to the output
        of the `process_vcf` function for a given sample set.
    """
    with ProcessPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(process_vcf, vcf_path, dp_thresholds, sample_set) for sample_set in sample_sets]
        results = [future.result() for future in futures]
    return results

# Define the DP thresholds to evaluate
dp_thresholds = list(range(1, 31))
# Run parallel processing
results = run_parallel_processing(vcf_path, dp_thresholds, [cteph_30x, cteph_15x])

# Display the results
print("Results for 30x samples:", results[0])
print("Results for 15x samples:", results[1])

# Save the results to a CSV file
results_df = pd.DataFrame(results, index=["30x", "15x"])
results_df.to_csv(args.output_csv)
print("Results saved to CSV file:", args.output_csv)
