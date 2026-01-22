import argparse
import pandas as pd
import subprocess
import json
import os
import sys
import datetime

def run_cmd(cmd, verbose=True):
    if verbose:
        print(f"[CMD] {cmd}")
    try:
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {cmd}")
        sys.exit(1)

def get_vcf_count(vcf_path):
    """Get variant count using bcftools index -n"""
    try:
        # Check if index exists, if not, return -1 or try to count
        # Assuming index exists as per pipeline flow
        output = subprocess.check_output(f"bcftools index -n {vcf_path}", shell=True).decode().strip()
        return int(output)
    except:
        return 0

def save_json_log(log_data, output_path):
    with open(output_path, 'w') as f:
        json.dump(log_data, f, indent=4)

def main():
    parser = argparse.ArgumentParser(description="Sensitivity Check Data Preparation")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--tbi", required=False) # Not explicitly used if bcftools auto-detects
    parser.add_argument("--summary", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    start_time = datetime.datetime.now().isoformat()
    log_data = {
        "timestamp": start_time,
        "input_file": args.vcf,
        "summary_file": args.summary,
        "stats": {}
    }

    # 1. Get Input Count
    print("Counting input variants...")
    input_count = get_vcf_count(args.vcf)
    log_data["stats"]["input_variants"] = input_count
    log_data["stats"]["input_variants_fmt"] = f"{input_count:,}"

    # 2. Load Summary and Prepare IDs
    print(f"Loading summary: {args.summary}")
    
    file_ids_stat1 = "temp.stat1.ids"
    file_ids_stat12 = "temp.stat12.ids"

    count_stat1 = 0
    count_stat12 = 0
    total_rare = 0

    try:
        chunk_size = 100000
        with open(file_ids_stat1, 'w') as f1, open(file_ids_stat12, 'w') as f12:
            # Iterate over chunks to save memory
            for chunk in pd.read_csv(args.summary, sep='\t', usecols=['VARIANT_ID', 'GROUP', 'FILTER_STAT'], chunksize=chunk_size):
                rare_chunk = chunk[chunk['GROUP'] == 'rare']
                if rare_chunk.empty:
                    continue
                
                total_rare += len(rare_chunk)

                # Set 1: Stat_1
                ids_stat1 = rare_chunk[rare_chunk['FILTER_STAT'] == 'Stat_1']['VARIANT_ID']
                if not ids_stat1.empty:
                    f1.write('\n'.join(ids_stat1) + '\n')
                    count_stat1 += len(ids_stat1)

                # Set 2: Stat_1 + Stat_2
                ids_stat12 = rare_chunk[rare_chunk['FILTER_STAT'].isin(['Stat_1', 'Stat_2'])]['VARIANT_ID']
                if not ids_stat12.empty:
                    f12.write('\n'.join(ids_stat12) + '\n')
                    count_stat12 += len(ids_stat12)
                    
    except Exception as e:
        print(f"Error loading summary or writing IDs: {e}")
        sys.exit(1)

    log_data["stats"]["group_rare_total"] = total_rare
    log_data["stats"]["ids_stat1_count"] = count_stat1
    log_data["stats"]["ids_stat12_count"] = count_stat12

    print(f"ID extraction complete. Stat1: {count_stat1}, Stat1+2: {count_stat12}")

    # 3. Process Stat 1
    out_stat1 = f"{args.out_prefix}.stat1.vcf.gz"
    print(f"Processing Stat1 -> {out_stat1}")
    cmd_stat1 = (f"bcftools view --threads {args.threads} "
                 f"--include 'ID=@{file_ids_stat1}' "
                 f"{args.vcf} -O z -o {out_stat1}")
    run_cmd(cmd_stat1)
    run_cmd(f"bcftools index -t {out_stat1}")
    
    count_stat1 = get_vcf_count(out_stat1)
    log_data["stats"]["output_stat1"] = {
        "file": out_stat1,
        "description": "GROUP=rare & FILTER_STAT=Stat_1",
        "variants_count": count_stat1,
        "variants_count_fmt": f"{count_stat1:,}"
    }

    # 4. Process Stat 1+2
    out_stat12 = f"{args.out_prefix}.stat1_stat2.vcf.gz"
    print(f"Processing Stat1+2 -> {out_stat12}")
    cmd_stat12 = (f"bcftools view --threads {args.threads} "
                 f"--include 'ID=@{file_ids_stat12}' "
                 f"{args.vcf} -O z -o {out_stat12}")
    run_cmd(cmd_stat12)
    run_cmd(f"bcftools index -t {out_stat12}")

    count_stat12 = get_vcf_count(out_stat12)
    log_data["stats"]["output_stat12"] = {
        "file": out_stat12,
        "description": "GROUP=rare & FILTER_STAT in [Stat_1, Stat_2]",
        "variants_count": count_stat12,
        "variants_count_fmt": f"{count_stat12:,}"
    }

    # Cleanup
    if os.path.exists(file_ids_stat1): os.remove(file_ids_stat1)
    if os.path.exists(file_ids_stat12): os.remove(file_ids_stat12)

    # Save Log
    log_file = f"{args.out_prefix}.sensitivity_extract.json"
    save_json_log(log_data, log_file)
    print(f"Log saved to {log_file}")

if __name__ == "__main__":
    main()
