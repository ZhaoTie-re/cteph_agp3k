params.WGSPrefix = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/18.tommo_panel_filter/cteph_agp3k.lowfreq_common'
params.ArrayPrefix = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/wgs_vs_imputed_array/results/03.variant_qc_missing/cteph_agp3k.imputed_array.vmiss_qc'
params.OutputDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/wgs_vs_imputed_array_assoc/results'
params.Plink2Path = '/home/b/b37974/plink2'

// Create channel with WGS and Array data
Channel
    .from(
        ['wgs', file("${params.WGSPrefix}.bed"), file("${params.WGSPrefix}.bim"), file("${params.WGSPrefix}.fam")],
        ['array', file("${params.ArrayPrefix}.bed"), file("${params.ArrayPrefix}.bim"), file("${params.ArrayPrefix}.fam")]
    )
    .toList()
    .set { input_data }

process PrepareWGSandArrayCommon {
    executor 'slurm'
    queue 'gr10478b'
    time '6h'
    tag "wgs_vs_array: extract common samples and variants"
    publishDir "${params.OutputDir}/01.common_samples_variants", mode: 'symlink'

    input:
    val(data_list) from input_data

    output:
    tuple val("wgs"), path("wgs_vs_array.wgs_common.bed"), path("wgs_vs_array.wgs_common.bim"), path("wgs_vs_array.wgs_common.fam") into wgs_ch
    tuple val("array"), path("wgs_vs_array.array_common.bed"), path("wgs_vs_array.array_common.bim"), path("wgs_vs_array.array_common.fam") into array_ch
    path "wgs_array_intersection.log"
    path "common_samples.txt"
    path "common_variants.txt"

    script:
    def wgs_data = data_list.find { it[0] == 'wgs' }
    def arr_data = data_list.find { it[0] == 'array' }
    """
    set -euo pipefail
    echo "=== WGS vs Array: common samples & variants ===" > wgs_array_intersection.log
    echo "Date: \$(date)" >> wgs_array_intersection.log
    echo "WGS files: ${wgs_data[1].name}, ${wgs_data[2].name}, ${wgs_data[3].name}" >> wgs_array_intersection.log
    echo "Array files: ${arr_data[1].name}, ${arr_data[2].name}, ${arr_data[3].name}" >> wgs_array_intersection.log
    echo "" >> wgs_array_intersection.log
    
    # Set file paths using absolute paths from staged files
    wgs_bed="${wgs_data[1]}"
    wgs_bim="${wgs_data[2]}"
    wgs_fam="${wgs_data[3]}"
    arr_bed="${arr_data[1]}"
    arr_bim="${arr_data[2]}"
    arr_fam="${arr_data[3]}"
    
    wgs_pref=\${wgs_bed%.bed}
    arr_pref=\${arr_bed%.bed}
    wgs_out="wgs_vs_array.wgs_common"
    arr_out="wgs_vs_array.array_common"
    
    echo "WGS prefix: \${wgs_pref}" >> wgs_array_intersection.log
    echo "Array prefix: \${arr_pref}" >> wgs_array_intersection.log
    echo "" >> wgs_array_intersection.log

    echo "All input files found." >> wgs_array_intersection.log
    echo "" >> wgs_array_intersection.log

    # Identify common samples and variants using command line tools
    echo "Identifying common samples and variants..." >> wgs_array_intersection.log
    
    # Extract samples where FID == IID from both FAM files
    awk '\$1 == \$2 {print \$1}' \${wgs_fam} | sort > wgs_fid_eq_iid.txt
    awk '\$1 == \$2 {print \$1}' \${arr_fam} | sort > arr_fid_eq_iid.txt
    
    wgs_fid_eq_iid_cnt=\$(wc -l < wgs_fid_eq_iid.txt)
    arr_fid_eq_iid_cnt=\$(wc -l < arr_fid_eq_iid.txt)
    echo "WGS fam records (FID==IID): \${wgs_fid_eq_iid_cnt}" | tee -a wgs_array_intersection.log
    echo "Array fam records (FID==IID): \${arr_fid_eq_iid_cnt}" | tee -a wgs_array_intersection.log
    
    # Find common sample IDs (intersection)
    comm -12 wgs_fid_eq_iid.txt arr_fid_eq_iid.txt > common_ids.txt
    common_sample_cnt=\$(wc -l < common_ids.txt)
    echo "Common sample count: \${common_sample_cnt}" | tee -a wgs_array_intersection.log
    
    # Create PLINK2 --keep format file (FID IID)
    awk '{print \$1, \$1}' common_ids.txt > common_samples.txt
    
    # Extract variant IDs (column 2) from BIM files
    awk '{print \$2}' \${wgs_bim} | sort > wgs_variants.txt
    awk '{print \$2}' \${arr_bim} | sort > arr_variants.txt
    
    wgs_var_cnt=\$(wc -l < wgs_variants.txt)
    arr_var_cnt=\$(wc -l < arr_variants.txt)
    echo "WGS variants: \${wgs_var_cnt}" | tee -a wgs_array_intersection.log
    echo "Array variants: \${arr_var_cnt}" | tee -a wgs_array_intersection.log
    
    # Find common variants (intersection)
    comm -12 wgs_variants.txt arr_variants.txt > common_variants.txt
    common_var_cnt=\$(wc -l < common_variants.txt)
    echo "Common variants: \${common_var_cnt}" | tee -a wgs_array_intersection.log
    
    # Phenotype breakdown for common samples
    if [ \${common_sample_cnt} -gt 0 ]; then
        # Extract phenotype info for common samples in WGS
        grep -Ff common_ids.txt \${wgs_fam} | awk '\$1 == \$2' > wgs_common_samples.fam
        wgs_cases=\$(awk '\$6 == 2' wgs_common_samples.fam | wc -l)
        wgs_ctrls=\$(awk '\$6 == 1' wgs_common_samples.fam | wc -l)
        wgs_missing=\$(awk '\$6 == -9 || \$6 == 0' wgs_common_samples.fam | wc -l)
        echo "WGS pheno (cases,ctrls,missing) in common samples: (\${wgs_cases},\${wgs_ctrls},\${wgs_missing})" | tee -a wgs_array_intersection.log
        
        # Extract phenotype info for common samples in Array
        grep -Ff common_ids.txt \${arr_fam} | awk '\$1 == \$2' > arr_common_samples.fam
        arr_cases=\$(awk '\$6 == 2' arr_common_samples.fam | wc -l)
        arr_ctrls=\$(awk '\$6 == 1' arr_common_samples.fam | wc -l)
        arr_missing=\$(awk '\$6 == -9 || \$6 == 0' arr_common_samples.fam | wc -l)
        echo "Array pheno (cases,ctrls,missing) in common samples: (\${arr_cases},\${arr_ctrls},\${arr_missing})" | tee -a wgs_array_intersection.log
    else
        echo "No common samples to report phenotype stats." | tee -a wgs_array_intersection.log
    fi
    
    echo "Sample and variant identification completed." >> wgs_array_intersection.log
    echo "" >> wgs_array_intersection.log

    # Extract common samples and variants from both datasets using PLINK2
    echo "Subsetting WGS by common samples + variants..." >> wgs_array_intersection.log
    ${params.Plink2Path} --bfile \${wgs_pref} --keep common_samples.txt --extract common_variants.txt --make-bed --out \${wgs_out} --threads 8 >> wgs_array_intersection.log 2>&1

    echo "Subsetting Array by common samples + variants..." >> wgs_array_intersection.log
    ${params.Plink2Path} --bfile \${arr_pref} --keep common_samples.txt --extract common_variants.txt --make-bed --out \${arr_out} --threads 8 >> wgs_array_intersection.log 2>&1

    # Summarize results
    wgs_variants=\$(wc -l < \${wgs_out}.bim)
    wgs_samples=\$(wc -l < \${wgs_out}.fam)
    arr_variants=\$(wc -l < \${arr_out}.bim)
    arr_samples=\$(wc -l < \${arr_out}.fam)

    echo "WGS output: \${wgs_out}.bed / .bim / .fam -> variants: \${wgs_variants}, samples: \${wgs_samples}" >> wgs_array_intersection.log
    echo "Array output: \${arr_out}.bed / .bim / .fam -> variants: \${arr_variants}, samples: \${arr_samples}" >> wgs_array_intersection.log

    echo "Done." >> wgs_array_intersection.log
    """
}

// Mix the two output channels into one
wgs_array_common = wgs_ch.mix(array_ch)

params.PhenoFile = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/19.cov_pheno_prepare/cteph_agp3k.bbj.projection.pheno_df.csv'
params.CovarFile = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/19.cov_pheno_prepare/cteph_agp3k.bbj.projection.cov_df.no_age.csv'

// Perform additive model association analysis
process AssocPlink2 {
    executor 'slurm'
    queue 'gr10478b'
    time '6d'
    tag "assoc_plink2: ${dataset}_additive"

    publishDir "${params.OutputDir}/02.assoc_result/${dataset}", mode: 'symlink'

    input:
    tuple val(dataset), path(bed), path(bim), path(fam) from wgs_array_common

    output:
    path "*.log"
    tuple val(dataset), path("*.glm.logistic") into assoc_plink2_out
    script:
    def bed_prefix = bed.baseName
    def pheno_name = "PHENO1"
    def covar_name = "SEX,PC1_AVG-PC10_AVG"
    def out_prefix = "wgs_vs_array.${dataset}.sex.10pc.additive"

    """
    ${params.Plink2Path} \\
        --bfile ${bed_prefix} \\
        --pheno ${params.PhenoFile} \\
        --pheno-name ${pheno_name} \\
        --covar ${params.CovarFile} \\
        --covar-name ${covar_name} \\
        --glm cols=+beta omit-ref no-firth hide-covar \\
        --out ${out_prefix} \\
        --ci 0.95 \\
        --threads 16
    """
}

assoc_plink2_out.view()



