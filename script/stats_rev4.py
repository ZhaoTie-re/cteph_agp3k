# %%
# usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import argparse
import pickle
import os

parser = argparse.ArgumentParser(description="Calculate genotype confusion metrics & allele confusion matrix & summary statistics")
parser.add_argument("--detail_mx", type=str, help="Path to picard.vcf.GenotypeConcordanceDetailMetrics from GATK picard GenotypeConcordance tool")
parser.add_argument("--group_info", type=str, help="Path to the group info file (e.g. wgs_array_dp.csv)")
parser.add_argument("--DP", type=int, help="DP threshold")
parser.add_argument("--GQ", type=int, help="GQ threshold")
parser.add_argument("--common_ids", type=str, help="Path to the common IDs file (e.g. common_ids.txt)") 
args = parser.parse_args()

detail_mx = pd.read_csv(args.detail_mx, sep="\s+", header=1)


# %%
def add_hom_ref_rows(df, total_num):
    # Group by TRUTH_SAMPLE and calculate the sum of COUNT
    count_sums = df.groupby('TRUTH_SAMPLE')['COUNT'].sum()
    
    # Extract unique DP, GQ, and CALL_SAMPLE values for each TRUTH_SAMPLE
    dp_gq_call_values = df.groupby('TRUTH_SAMPLE')[['DP', 'GQ', 'CALL_SAMPLE']].first()
    
    # Create new rows for each TRUTH_SAMPLE
    new_rows = []
    for truth_sample, count_sum in count_sums.items():
        dp, gq, call_sample = dp_gq_call_values.loc[truth_sample]
        new_row = {
            'VARIANT_TYPE': 'SNP',
            'TRUTH_SAMPLE': truth_sample,
            'CALL_SAMPLE': call_sample,
            'TRUTH_STATE': 'HOM_REF',
            'CALL_STATE': 'HOM_REF',
            'COUNT': total_num - count_sum,
            'CONTINGENCY_VALUES': 'TN',
            'DP': dp,
            'GQ': gq
        }
        new_rows.append(new_row)
    
    # Append the new rows to the dataframe
    return pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)

common_ids = pd.read_csv(args.common_ids, header=None).values
# Define total_num
total_num = common_ids.shape[0] # This should be the total number of common SNPs in the VCF file

# Apply the function to detail_mx
detail_mx = add_hom_ref_rows(detail_mx, total_num)
detail_mx = detail_mx.sort_values(by='TRUTH_SAMPLE', ascending=True).reset_index(drop=True)

# %%
def reformat_detail_mx(df):
    # Rule 1: Update CONTINGENCY_VALUES to 'FN'
    condition_fn = (df['TRUTH_STATE'] == 'HOM_VAR1') & (df['CALL_STATE'] == 'HOM_REF')
    df.loc[condition_fn, 'CONTINGENCY_VALUES'] = 'FN'
    
    # Rule 2: Update CONTINGENCY_VALUES to 'FP'
    condition_fp = (df['TRUTH_STATE'] == 'HOM_REF') & (df['CALL_STATE'] == 'HOM_VAR1')
    df.loc[condition_fp, 'CONTINGENCY_VALUES'] = 'FP'
    
    # Rule 3: Update CONTINGENCY_VALUES to 'EMPTY'
    condition_empty = df['CALL_STATE'].str.contains('LOW_DP|LOW_GQ', regex=True)
    df.loc[condition_empty, 'CONTINGENCY_VALUES'] = 'EMPTY'
    
    return df

# Apply the function to detail_mx
detail_mx = reformat_detail_mx(detail_mx)

# %%
group_info = pd.read_csv(args.group_info)
group_info_30x = group_info[group_info['Target DP (JHRPv4)'] == '30x']
group_info_15x = group_info[group_info['Target DP (JHRPv4)'] == '15x']

detail_mx_all = detail_mx[detail_mx['TRUTH_SAMPLE'].isin(group_info['ID'])]
detail_mx_30x = detail_mx[detail_mx['TRUTH_SAMPLE'].isin(group_info_30x['ID'])]
detail_mx_15x = detail_mx[detail_mx['TRUTH_SAMPLE'].isin(group_info_15x['ID'])]

# %%
def calculate_confusion_summary_matrices(detail_mx):
    # Genotype confusion matrix
    # Ensure the categories match the actual values in the columns
    detail_mx['TRUTH_STATE'] = pd.Categorical(
        detail_mx['TRUTH_STATE'], 
        categories=['HOM_REF', 'HET_REF_VAR1', 'HOM_VAR1', 'NO_CALL'], 
        ordered=True
    )
    
    detail_mx['CALL_STATE'] = pd.Categorical(
        detail_mx['CALL_STATE'], 
        categories=['HOM_REF', 'HET_REF_VAR1', 'HOM_VAR1', 'NO_CALL', 'LOW_DP', 'LOW_GQ'], 
        ordered=True
    )
    
    genotype_confusion_matrix = detail_mx.pivot_table(
        index='TRUTH_STATE',
        columns='CALL_STATE',
        values='COUNT',
        aggfunc='sum',
        fill_value=0
    )

    # Rename rows and columns
    genotype_confusion_matrix.index = ['TRUE_HOM_REF', 'TRUE_HET_REF_VAR1', 'TRUE_HOM_VAR1', 'TRUE_NO_CALL']
    genotype_confusion_matrix.columns = ['CALL_HOM_REF', 'CALL_HET_REF_VAR1', 'CALL_HOM_VAR1', 'CALL_NO_CALL', 'CALL_LOW_DP', 'CALL_LOW_GQ']
    
    # Combine 'CALL_NO_CALL', 'CALL_LOW_DP', and 'CALL_LOW_GQ' into 'CALL_NO_CALL'
    genotype_confusion_matrix['CALL_NO_CALL'] += genotype_confusion_matrix[['CALL_LOW_DP', 'CALL_LOW_GQ']].sum(axis=1)
    
    # Drop the 'CALL_LOW_DP' and 'CALL_LOW_GQ' columns
    genotype_confusion_matrix = genotype_confusion_matrix.drop(columns=['CALL_LOW_DP', 'CALL_LOW_GQ'])

    # Allele confusion matrix
    allele_confusion_matrix = pd.DataFrame(index=['TRUE_REF', 'TRUE_ALT'], columns=['CALL_REF', 'CALL_ALT'])

    TN = (
        detail_mx.loc[detail_mx['CONTINGENCY_VALUES'] == 'TN', 'COUNT'].sum() * 2 +
        detail_mx.loc[detail_mx['CONTINGENCY_VALUES'] == 'FP,TN', 'COUNT'].sum() +
        detail_mx.loc[detail_mx['CONTINGENCY_VALUES'] == 'TN,FN', 'COUNT'].sum() +
        detail_mx.loc[detail_mx['CONTINGENCY_VALUES'] == 'TP,TN', 'COUNT'].sum()
    )

    FN = (
        detail_mx.loc[detail_mx['CONTINGENCY_VALUES'] == 'FN', 'COUNT'].sum() * 2 +
        detail_mx.loc[detail_mx['CONTINGENCY_VALUES'] == 'TN,FN', 'COUNT'].sum() +
        detail_mx.loc[detail_mx['CONTINGENCY_VALUES'] == 'TP,FN', 'COUNT'].sum()
    )

    TP = (
        detail_mx.loc[detail_mx['CONTINGENCY_VALUES'] == 'TP', 'COUNT'].sum() * 2 +
        detail_mx.loc[detail_mx['CONTINGENCY_VALUES'] == 'TP,FN', 'COUNT'].sum() +
        detail_mx.loc[detail_mx['CONTINGENCY_VALUES'] == 'TP,TN', 'COUNT'].sum() +
        detail_mx.loc[detail_mx['CONTINGENCY_VALUES'] == 'TP,FP', 'COUNT'].sum()
    )

    FP = (
        detail_mx.loc[detail_mx['CONTINGENCY_VALUES'] == 'FP', 'COUNT'].sum() * 2 +
        detail_mx.loc[detail_mx['CONTINGENCY_VALUES'] == 'FP,TN', 'COUNT'].sum() +
        detail_mx.loc[detail_mx['CONTINGENCY_VALUES'] == 'TP,FP', 'COUNT'].sum()
    )

    # Fill the allele confusion matrix
    allele_confusion_matrix.loc['TRUE_REF', 'CALL_REF'] = TN
    allele_confusion_matrix.loc['TRUE_REF', 'CALL_ALT'] = FP
    allele_confusion_matrix.loc['TRUE_ALT', 'CALL_REF'] = FN
    allele_confusion_matrix.loc['TRUE_ALT', 'CALL_ALT'] = TP
    
    # Calculate summary statistics table
    dp = detail_mx['DP'].iloc[0]
    gq = detail_mx['GQ'].iloc[0]    
    
    total_genotypes = genotype_confusion_matrix.loc[['TRUE_HOM_REF', 'TRUE_HET_REF_VAR1', 'TRUE_HOM_VAR1'], 
                                                ['CALL_HOM_REF', 'CALL_HET_REF_VAR1', 'CALL_HOM_VAR1']].values.sum()
    
    total_genotypes_no_call = genotype_confusion_matrix.values.sum()
    # Calculate FP RATE & FN RATE
    geno_fp = genotype_confusion_matrix.loc['TRUE_HOM_REF', ['CALL_HET_REF_VAR1', 'CALL_HOM_VAR1']].values.sum() + \
        genotype_confusion_matrix.loc['TRUE_HET_REF_VAR1', 'CALL_HOM_VAR1']
    
    geno_fn = genotype_confusion_matrix.loc['TRUE_HOM_VAR1', ['CALL_HOM_REF', 'CALL_HET_REF_VAR1']].values.sum() + \
        genotype_confusion_matrix.loc['TRUE_HET_REF_VAR1', 'CALL_HOM_REF']
    
    false_positive_rate = geno_fp / total_genotypes
    false_negative_rate = geno_fn / total_genotypes
    
    # Calculate GENOTYPE_CONCORDANCE & GENOTYPE_CONCORDANCE(WITH_EMPTY)
    concordant_genotypes = genotype_confusion_matrix.loc['TRUE_HOM_REF', 'CALL_HOM_REF'] + \
                       genotype_confusion_matrix.loc['TRUE_HET_REF_VAR1', 'CALL_HET_REF_VAR1'] + \
                       genotype_confusion_matrix.loc['TRUE_HOM_VAR1', 'CALL_HOM_VAR1']
    
    geno_concordance = concordant_genotypes / total_genotypes
    geno_concordance_no_call = concordant_genotypes / (total_genotypes + genotype_confusion_matrix['CALL_NO_CALL'].sum() - genotype_confusion_matrix.loc['TRUE_NO_CALL', 'CALL_NO_CALL'])
    
    # Calculate MISS_RATE
    miss_rate = (genotype_confusion_matrix['CALL_NO_CALL'].sum() - genotype_confusion_matrix.loc['TRUE_NO_CALL', 'CALL_NO_CALL']) / (total_genotypes + genotype_confusion_matrix['CALL_NO_CALL'].sum() - genotype_confusion_matrix.loc['TRUE_NO_CALL', 'CALL_NO_CALL'])
    
    # Calculate Individual misclassification components
    het2homvar_count = genotype_confusion_matrix.loc['TRUE_HET_REF_VAR1', 'CALL_HOM_VAR1']
    homref2homvar_count = genotype_confusion_matrix.loc['TRUE_HOM_REF', 'CALL_HOM_VAR1']
    homref2het_count = genotype_confusion_matrix.loc['TRUE_HOM_REF', 'CALL_HET_REF_VAR1']

    het2homvar_rate = het2homvar_count / total_genotypes
    homref2homvar_rate = homref2homvar_count / total_genotypes
    homref2het_rate = homref2het_count / total_genotypes
    
    
    # Create summary statistics DataFrame and fill it with values
    summary_statistics = pd.DataFrame(
    [[dp, gq, total_genotypes, total_genotypes_no_call, 
      geno_concordance, geno_concordance_no_call, miss_rate,
      false_positive_rate, false_negative_rate,
      het2homvar_count, homref2homvar_count, homref2het_count,
      het2homvar_rate, homref2homvar_rate, homref2het_rate]],
    columns=[
        'DP', 'GQ', 'TOTAL_GENOTYPE', 'TOTAL_GENOTYPE(WITH_EMPTY)',
        'GENOTYPE_CONCORDANCE', 'GENOTYPE_CONCORDANCE(WITH_EMPTY)', 'MISS_RATE',
        'FALSE_POSITIVE_RATE', 'FALSE_NEGATIVE_RATE',
        'HET>HOMVAR_COUNT', 'HOMREF>HOMVAR_COUNT', 'HOMREF>HET_COUNT',
        'HET>HOMVAR_RATE', 'HOMREF>HOMVAR_RATE', 'HOMREF>HET_RATE'
        ]
    )

    return genotype_confusion_matrix, allele_confusion_matrix, summary_statistics

# %%
geno_cmx_30x, allele_cmx_30x, summary_statistics_30x = calculate_confusion_summary_matrices(detail_mx_30x)
geno_cmx_15x, allele_cmx_15x, summary_statistics_15x = calculate_confusion_summary_matrices(detail_mx_15x)
geno_cmx_all, allele_cmx_all, summary_statistics_all = calculate_confusion_summary_matrices(detail_mx_all)

# %%
global_dp = args.DP
global_gq = args.GQ

set1 = frozenset([f'DP{global_dp}', f'GQ{global_gq}', '15X'])
set2 = frozenset([f'DP{global_dp}', f'GQ{global_gq}', '30X'])
set3 = frozenset([f'DP{global_dp}', f'GQ{global_gq}', 'ALL'])

results_dict = {
    set1: (geno_cmx_15x, allele_cmx_15x, summary_statistics_15x),
    set2: (geno_cmx_30x, allele_cmx_30x, summary_statistics_30x),
    set3: (geno_cmx_all, allele_cmx_all, summary_statistics_all)
}

# %%
# save the results_dict to a file

import pickle

with open(f'DP{global_dp}_GQ{global_gq}_dict.pkl', 'wb') as f:
    pickle.dump(results_dict, f)


