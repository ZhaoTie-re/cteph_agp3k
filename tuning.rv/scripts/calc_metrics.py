import pandas as pd
import numpy as np
import argparse
import sys
import os
import multiprocessing
import logging
import time
from functools import partial

def setup_logging(log_file):
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    # Also log to stdout
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

METRIC_DESCRIPTIONS = {
    # Sample Metrics
    'SampleID': 'Unique identifier for the sample (from input VCF/Info file).',
    'Group': 'Sample Group (Case/Control) determined by info file outcome column.',
    'TargetDP': 'Target Depth (e.g., sequencing depth target from Info file).',
    'MeanDP': 'Mean Depth (Actual average coverage from Info file).',
    'SNumHomRef': 'Count of sites where the sample is Homozygous Reference.',
    'SNumHet': 'Count of sites where the sample is Heterozygous.',
    'SNumHomAlt': 'Count of sites where the sample is Homozygous Alternative (relative to REF genome).',
    'SMissCount': 'Count of sites where the sample has missing genotype calls.',
    'SMinAC': 'Sum of Minor Allele Counts (calculated as N_het + 2 * N_hom_minor using force-accepted major allele).',

    # Variant Metrics
    '#CHROM': 'Chromosome identifier.',
    'POS': 'Genomic position.',
    'VariantID': 'Variant identifier (format: CHROM:POS:REF:ALT).',
    'REF': 'Reference allele sequence.',
    'ALT': 'Alternative allele sequence.',
    'RefAC': 'Reference Allele Count. Total number of REF alleles observed.',
    'AltAC': 'Alternative Allele Count. Total number of ALT alleles observed.',
    'MinAC': 'Minor Allele Count. The count of the less frequent allele (min(RefAC, AltAC)).',
    'VNumHomRef': 'Number of samples that are Homozygous Reference (0/0).',
    'VNumHet': 'Number of samples that are Heterozygous (0/1).',
    'VNumHomAlt': 'Number of samples that are Homozygous Alternative (1/1).',
    'VMissCount': 'Number of samples where the genotype call is missing.'
}

def write_documentation(out_dir, sample_info, variant_info):
    """
    Generates a Markdown file describing the output metrics.
    """
    try:
        doc_path = os.path.join(out_dir, 'METRICS_DESCRIPTION.md')
        sample_rows, sample_cols = sample_info
        variant_rows, variant_cols = variant_info
        
        with open(doc_path, 'w') as f:
            f.write("# Metrics Output Documentation\n\n")
            f.write(f"Generated on {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write(f"## 1. Sample Metrics (sample_metrics.txt)\n")
            f.write(f"**Dimensions:** {sample_rows:,} samples x {len(sample_cols)} columns\n\n")
            f.write("| Column Name | Description |\n")
            f.write("| :--- | :--- |\n")
            for col in sample_cols:
                desc = METRIC_DESCRIPTIONS.get(col, 'Custom metric.')
                f.write(f"| **{col}** | {desc} |\n")
            f.write("\n")
            
            f.write(f"## 2. Variant Metrics (variant_metrics.txt)\n")
            f.write(f"**Dimensions:** {variant_rows:,} variants x {len(variant_cols)} columns\n\n")
            f.write("| Column Name | Description |\n")
            f.write("| :--- | :--- |\n")
            for col in variant_cols:
                desc = METRIC_DESCRIPTIONS.get(col, 'Custom metric.')
                f.write(f"| **{col}** | {desc} |\n")
                
        logging.info(f"Written documentation to {doc_path}")
    except Exception as e:
        logging.warning(f"Failed to write documentation: {e}")

def process_sample_metrics(args):
    logging.info("Processing Sample Metrics...")
    try:
        scount = pd.read_csv(args.sample_counts, sep='\t')
        
        # Handle ID column renaming
        if '#IID' in scount.columns:
            scount.rename(columns={'#IID': 'SampleID'}, inplace=True)
        elif 'IID' in scount.columns:
            scount.rename(columns={'IID': 'SampleID'}, inplace=True)
        else:
             logging.error(f"Error: Could not find IID or #IID column in sample counts. Columns: {scount.columns}")
             sys.exit(1)
        
        # Read Sample Missingness
        smiss = pd.read_csv(args.sample_missing, sep='\t')
        if '#IID' in smiss.columns:
            smiss.rename(columns={'#IID': 'SampleID'}, inplace=True)
        elif 'IID' in smiss.columns:
            smiss.rename(columns={'IID': 'SampleID'}, inplace=True)
        else:
             logging.error(f"Error: Could not find IID or #IID column in sample missingness. Columns: {smiss.columns}")
             sys.exit(1)
            
        # Merge missingness
        scount = pd.merge(scount, smiss[['SampleID', 'MISSING_CT']], on='SampleID', how='left')
        
        info = pd.read_excel(args.info)
        
        if args.id_col not in info.columns:
            logging.error(f"Error: ID column '{args.id_col}' not found in info file.")
            sys.exit(1)
            
        # Determine Group if columns provided
        group_col_present = False
        if args.group_col and args.group_col in info.columns:
            group_col_present = True
            
        cols_to_use = [args.id_col, args.tdp_col, args.mdp_col]
        if group_col_present:
            cols_to_use.append(args.group_col)
            
        info_subset = info[cols_to_use].copy()
        
        rename_dict = {
            args.id_col: 'SampleID',
            args.tdp_col: 'TargetDP',
            args.mdp_col: 'MeanDP'
        }
        
        info_subset.rename(columns=rename_dict, inplace=True)
        
        if group_col_present and args.case_value:
             # Process Group: map to Case/Control
             # Ensure string comparison
             info_subset['Group'] = info_subset[args.group_col].astype(str).apply(
                 lambda x: 'Case' if x == str(args.case_value) else 'Control'
             )
        
        sample_merged = pd.merge(scount, info_subset, on='SampleID', how='left')
        
        sample_final = pd.DataFrame()
        sample_final['SampleID'] = sample_merged['SampleID']
        
        if 'Group' in sample_merged.columns:
            sample_final['Group'] = sample_merged['Group']
        else:
            sample_final['Group'] = 'Unknown'
            
        sample_final['TargetDP'] = sample_merged['TargetDP']
        sample_final['MeanDP'] = sample_merged['MeanDP']
        
        if 'HOM_REF_CT' in sample_merged.columns:
            sample_final['SNumHomRef'] = sample_merged['HOM_REF_CT']
        
        if 'HET_CT' in sample_merged.columns:
            sample_final['SNumHet'] = sample_merged['HET_CT']
        elif 'HET_SNP_CT' in sample_merged.columns:
            sample_final['SNumHet'] = sample_merged['HET_SNP_CT']
            
        if 'HOM_ALT_CT' in sample_merged.columns:
            sample_final['SNumHomAlt'] = sample_merged['HOM_ALT_CT']
        elif 'TWO_ALT_CT' in sample_merged.columns:
             sample_final['SNumHomAlt'] = sample_merged['TWO_ALT_CT']
        elif 'HOM_ALT_SNP_CT' in sample_merged.columns:
             sample_final['SNumHomAlt'] = sample_merged['HOM_ALT_SNP_CT']
        
        # Handle Missing Count (Prioritize MISSING_CT from .smiss)
        if 'MISSING_CT' in sample_merged.columns:
            sample_final['SMissCount'] = sample_merged['MISSING_CT']
        elif 'MISSING_INCL_FEMALE_Y_CT' in sample_merged.columns:
            sample_final['SMissCount'] = sample_merged['MISSING_INCL_FEMALE_Y_CT']
        else:
            sample_final['SMissCount'] = 0 # Default or warning?
            logging.warning("Missing count column not found in sample metrics. Setting to 0.")

        # Process SMinAC if file is provided
        if args.sample_minac:
            try:
                logging.info(f"Processing SMinAC from {args.sample_minac}...")
                minac_df = pd.read_csv(args.sample_minac, sep='\t')
                
                # Align ID column
                if '#IID' in minac_df.columns:
                    minac_df.rename(columns={'#IID': 'SampleID'}, inplace=True)
                elif 'IID' in minac_df.columns:
                    minac_df.rename(columns={'IID': 'SampleID'}, inplace=True)
                
                # Ensure we have needed columns
                # We expect HET_CT and HOM_ALT_CT relative to Major-Ref alignment
                if 'HET_CT' in minac_df.columns and 'HOM_ALT_CT' in minac_df.columns:
                    # SMinAC = Heterozygous + 2 * Homozygous Minor (HOM_ALT in Maj-Ref aligned file)
                    minac_df['SMinAC'] = minac_df['HET_CT'] + 2 * minac_df['HOM_ALT_CT']
                    
                    # Merge SMinAC into sample_final
                    # sample_final has SampleID, so we can merge
                    sample_final = pd.merge(sample_final, minac_df[['SampleID', 'SMinAC']], on='SampleID', how='left')
                    
                    # Fill NaNs with 0 (implies 0 minor alleles if sample missing from file?)
                    if 'SMinAC' in sample_final.columns:
                        sample_final['SMinAC'] = sample_final['SMinAC'].fillna(0).astype(int)
                    
                else:
                    logging.warning(f"HET_CT or HOM_ALT_CT not found in sample-minac file ({minac_df.columns}). SMinAC will not be calculated.")

            except Exception as e:
                logging.warning(f"Failed to process SMinAC from {args.sample_minac}: {e}")

        # Reorder columns
        desired_cols = ['SampleID', 'Group', 'TargetDP', 'MeanDP', 'SNumHomRef', 'SNumHet', 'SNumHomAlt', 'SMissCount', 'SMinAC']
        final_cols = [c for c in desired_cols if c in sample_final.columns]
        sample_final = sample_final[final_cols]
        
        sample_final.to_csv(args.out_sample, sep='\t', index=False)
        logging.info(f"Written {args.out_sample}")
        
        return len(sample_final), list(sample_final.columns)
        
    except Exception as e:
        logging.error(f"Error processing sample metrics: {e}")
        sys.exit(1)

def process_variant_chunk(chunk_data):
    """
    Process a tuple of (gcount_chunk, acount_chunk, vmiss_chunk)
    """
    g_chunk, a_chunk, m_chunk = chunk_data
    
    # Optimization: Check if IDs are aligned (Fast Path)
    # This avoids the expensive merge operation
    # PLINK2 output is typically line-aligned
    ids_aligned = False
    if len(g_chunk) == len(a_chunk) == len(m_chunk):
        if g_chunk['ID'].equals(a_chunk['ID']) and g_chunk['ID'].equals(m_chunk['ID']):
            ids_aligned = True
            
    if ids_aligned:
        # FAST PATH: Direct assignment
        df = g_chunk
        
        # Vectorized calculations using numpy
        obs_ct = a_chunk['OBS_CT'].values
        alt_cts = a_chunk['ALT_CTS'].values
        
        df['RefAC'] = obs_ct - alt_cts
        df['AltAC'] = alt_cts
        df['MinAC'] = np.minimum(df['RefAC'], df['AltAC'])
        
        # Add missing count from m_chunk
        if 'MISSING_CT' in m_chunk.columns:
            df['VMissCount'] = m_chunk['MISSING_CT'].values
        else:
            df['VMissCount'] = 0
        
        # Rename
        rename_map = {
            'ID': 'VariantID',
            '#ID': 'VariantID',
            'HOM_REF_CT': 'VNumHomRef',
            'HET_CT': 'VNumHet',
            'HET_REF_ALT1_CT': 'VNumHet',
            'HET_REF_ALT_CTS': 'VNumHet',
            'HOM_ALT_CT': 'VNumHomAlt',
            'TWO_ALT_CT': 'VNumHomAlt',
            'HOM_ALT1_CT': 'VNumHomAlt',
            'CHROM': '#CHROM'
        }
        df.rename(columns=rename_map, inplace=True)
        
    else:
        # SLOW PATH: Merge
        # Normalize ID columns first
        for chunk in [g_chunk, a_chunk, m_chunk]:
            if '#ID' in chunk.columns:
                chunk.rename(columns={'#ID': 'ID'}, inplace=True)

        g_chunk.rename(columns={'ID': 'VariantID'}, inplace=True)
        a_chunk.rename(columns={'ID': 'VariantID'}, inplace=True)
        m_chunk.rename(columns={'ID': 'VariantID'}, inplace=True)
        
        df = pd.merge(g_chunk, a_chunk[['VariantID', 'ALT_CTS', 'OBS_CT']], on='VariantID', how='left')
        
        if 'MISSING_CT' in m_chunk.columns:
            df = pd.merge(df, m_chunk[['VariantID', 'MISSING_CT']], on='VariantID', how='left')
            df['VMissCount'] = df['MISSING_CT']
        else:
            df['VMissCount'] = 0
        
        df['AltAC'] = df['ALT_CTS']
        df['RefAC'] = df['OBS_CT'] - df['ALT_CTS']
        df['MinAC'] = np.minimum(df['RefAC'], df['AltAC'])
        
        rename_map = {
            'HOM_REF_CT': 'VNumHomRef',
            'HET_CT': 'VNumHet',
            'HET_REF_ALT1_CT': 'VNumHet',
            'HET_REF_ALT_CTS': 'VNumHet',
            'HOM_ALT_CT': 'VNumHomAlt',
            'TWO_ALT_CT': 'VNumHomAlt',
            'HOM_ALT1_CT': 'VNumHomAlt',
            'CHROM': '#CHROM'
        }
        df.rename(columns=rename_map, inplace=True)

    # Extract POS from VariantID (assuming format CHROM:POS:REF:ALT)
    # Use a vectorized string operation for speed if possible, or apply
    try:
        # Try vectorized split if VariantID is consistent
        df['POS'] = df['VariantID'].astype(str).str.split(':').str[1]
    except Exception:
        # Fallback
        df['POS'] = df['VariantID'].apply(lambda x: x.split(':')[1] if ':' in str(x) else 0)

    # Select output columns
    desired_cols = ['#CHROM', 'POS', 'VariantID', 'REF', 'ALT', 'RefAC', 'AltAC', 'MinAC', 'VNumHomRef', 'VNumHet', 'VNumHomAlt', 'VMissCount']
    final_cols = [c for c in desired_cols if c in df.columns]
    
    return df[final_cols]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample-counts', required=True)
    parser.add_argument('--sample-missing', required=True)
    parser.add_argument('--variant-counts', required=True)
    parser.add_argument('--variant-missing', required=True)
    parser.add_argument('--freq-counts', required=True)
    parser.add_argument('--info', required=True)
    parser.add_argument('--id-col', required=True)
    parser.add_argument('--group-col', required=False, help='Column name for Group/Outcome')
    parser.add_argument('--case-value', required=False, help='Value indicating Case')
    parser.add_argument('--tdp-col', required=True)
    parser.add_argument('--mdp-col', required=True)
    parser.add_argument('--out-sample', required=True)
    parser.add_argument('--out-variant', required=True)
    parser.add_argument('--sample-minac', required=False, help='Path to sample counts file with maj-ref alignment for correct SMinAC calculation')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads for processing')
    parser.add_argument('--log', required=True, help='Log file path')
    args = parser.parse_args()

    setup_logging(args.log)
    logging.info("Starting calculation metrics...")
    logging.info(f"Arguments: {args}")

    start_time = time.time()

    # 1. Process Sample Metrics (usually small enough for memory)
    sample_rows, sample_cols = process_sample_metrics(args)

    # 2. Process Variant Metrics (Chunked & Parallel)
    logging.info("Processing Variant Metrics (Chunked)...")
    
    chunk_size = 10000 # Adjust based on memory availability
    
    # Create iterators for reading files in chunks
    # Read first chunk to check columns
    try:
        g_preview = pd.read_csv(args.variant_counts, sep='\t', nrows=0)
        a_preview = pd.read_csv(args.freq_counts, sep='\t', nrows=0)
        m_preview = pd.read_csv(args.variant_missing, sep='\t', nrows=0)
        logging.info(f"Variant counts columns: {list(g_preview.columns)}")
        logging.info(f"Freq counts columns: {list(a_preview.columns)}")
        logging.info(f"Variant missing columns: {list(m_preview.columns)}")
        
        if 'ID' not in g_preview.columns:
             # Try to find ID column (e.g. #ID)
             possible_ids = [c for c in g_preview.columns if 'ID' in c]
             logging.warning(f"'ID' column not found in variant counts. Possible candidates: {possible_ids}")
             
        if 'ID' not in a_preview.columns:
             possible_ids = [c for c in a_preview.columns if 'ID' in c]
             logging.warning(f"'ID' column not found in freq counts. Possible candidates: {possible_ids}")

    except Exception as e:
        logging.warning(f"Could not preview columns: {e}")

    g_iter = pd.read_csv(args.variant_counts, sep='\t', chunksize=chunk_size)
    a_iter = pd.read_csv(args.freq_counts, sep='\t', chunksize=chunk_size)
    m_iter = pd.read_csv(args.variant_missing, sep='\t', chunksize=chunk_size)
    
    # Use multiprocessing pool
    # We need to zip the iterators to process corresponding chunks together
    # Note: This assumes PLINK outputs are line-aligned (which they are for same run)
    
    pool_size = args.threads if args.threads > 0 else multiprocessing.cpu_count()
    logging.info(f"Using {pool_size} threads for variant processing.")
    
    first_chunk = True
    chunk_count = 0
    total_variants = 0
    
    with multiprocessing.Pool(processes=pool_size) as pool:
        # Create a generator that yields tuples of (g_chunk, a_chunk, m_chunk)
        chunk_generator = zip(g_iter, a_iter, m_iter)
        
        # imap returns results in order, which is crucial for writing to file
        for result_df in pool.imap(process_variant_chunk, chunk_generator):
            mode = 'w' if first_chunk else 'a'
            header = first_chunk
            
            result_df.to_csv(args.out_variant, sep='\t', index=False, mode=mode, header=header)
            
            chunk_count += 1
            total_variants += len(result_df)
            if chunk_count % 10 == 0:
                logging.info(f"Processed {chunk_count} chunks ({total_variants:,} variants so far)...")
            
            first_chunk = False
            
    logging.info(f"Written {args.out_variant}")
    logging.info(f"Total variants processed: {total_variants:,}")
    
    # 3. Generate Documentation
    variant_cols = ['#CHROM', 'POS', 'VariantID', 'REF', 'ALT', 'RefAC', 'AltAC', 'MinAC', 'VNumHomRef', 'VNumHet', 'VNumHomAlt', 'VMissCount']
    out_dir = os.path.dirname(os.path.abspath(args.out_sample))
    write_documentation(out_dir, (sample_rows, sample_cols), (total_variants, variant_cols))
    
    end_time = time.time()
    duration = end_time - start_time
    logging.info(f"Completed successfully in {duration:.2f} seconds.")

if __name__ == "__main__":
    main()
