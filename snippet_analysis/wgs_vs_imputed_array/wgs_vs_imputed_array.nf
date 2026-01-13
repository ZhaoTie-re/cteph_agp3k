params.WGSPrefix = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/18.tommo_panel_filter/cteph_agp3k.lowfreq_common'
params.ArrayGenotypePrefix = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array_imputed_rev1/results/09.gt_hardcall/cteph_agp3k.imputed_array.sqc.vqc.gt'
params.ArrayDosagePrefix = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array_imputed_rev1/results/08.variant_qc/cteph_agp3k.imputed_array.sqc.vqc'
params.OutDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/wgs_vs_imputed_array/results'
params.ScriptsDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/wgs_vs_imputed_array/scripts'
params.Plink2 = '/home/b/b37974/plink2'
params.Tabix = '/home/b/b37974/htslib-1.9/tabix'

// Create channels for input files
Channel
    .fromPath("${params.WGSPrefix}.{bed,bim,fam}")
    .collect()
    .map { files -> tuple('wgs', params.WGSPrefix, files) }
    .set { wgs_ch }

Channel
    .fromPath("${params.ArrayGenotypePrefix}.{pgen,pvar,psam}")
    .collect()
    .map { files -> tuple('array_gt', params.ArrayGenotypePrefix, files) }
    .set { array_gt_ch }

Channel
    .fromPath("${params.ArrayDosagePrefix}.{pgen,pvar,psam}")
    .collect()
    .map { files -> tuple('array_ds', params.ArrayDosagePrefix, files) }
    .set { array_ds_ch }

// Combine all three datasets
wgs_ch
    .concat(array_gt_ch, array_ds_ch)
    .collect()
    .set { all_datasets_ch }

process findCommonSamplesVariants {
    executor 'slurm'
    queue 'gr10478b'
    time '6h'
    tag "find_common"
    
    publishDir "${params.OutDir}/01.common_samples_variants", mode: 'symlink'
    
    input:
    val(datasets) from all_datasets_ch
    
    output:
    file("common_samples.txt") into common_samples_ch
    file("common_variants.txt") into common_variants_ch
    file("intersection_summary.txt")
    file("*.sample_counts.txt")
    file("*.variant_counts.txt")
    
    script:
    """
    # Find common samples and variants across all three datasets
    python3 ${params.ScriptsDir}/find_common_samples_variants.py \
        --wgs-prefix ${params.WGSPrefix} \
        --array-gt-prefix ${params.ArrayGenotypePrefix} \
        --array-ds-prefix ${params.ArrayDosagePrefix} \
        --output-dir . \
        --summary intersection_summary.txt
    """
}

process extractCommonSubsets {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "extract_${dataset_type}"
    
    publishDir "${params.OutDir}/02.extracted_genotypes/${dataset_type}", mode: 'symlink', pattern: "*.{pgen,pvar,psam}"
    
    input:
    each dataset_type from Channel.from('wgs', 'array_gt', 'array_ds')
    file(common_samples) from common_samples_ch
    file(common_variants) from common_variants_ch
    
    output:
    tuple val(dataset_type), file("*.common.pgen"), file("*.common.pvar"), file("*.common.psam") into common_subsets_out
    
    script:
    if (dataset_type == 'wgs') {
        input_prefix = params.WGSPrefix
        input_format = 'bfile'
    } else if (dataset_type == 'array_gt') {
        input_prefix = params.ArrayGenotypePrefix
        input_format = 'pfile'
    } else {
        input_prefix = params.ArrayDosagePrefix
        input_format = 'pfile'
    }
    
    out_prefix = "cteph_agp3k.${dataset_type}.common"
    
    """
    # Count original samples and variants
    if [ "${input_format}" == "bfile" ]; then
        ORIG_SAMPLES=\$(wc -l < ${input_prefix}.fam)
        ORIG_VARIANTS=\$(wc -l < ${input_prefix}.bim)
    else
        ORIG_SAMPLES=\$(tail -n +2 ${input_prefix}.psam | wc -l)
        ORIG_VARIANTS=\$(grep -v '^#' ${input_prefix}.pvar | wc -l)
    fi
    
    # Display dataset information
    echo "=========================================="
    if [ "${dataset_type}" == "wgs" ]; then
        echo "Processing: WGS (Whole Genome Sequencing)"
    elif [ "${dataset_type}" == "array_gt" ]; then
        echo "Processing: Imputed Array - Hard-called Genotypes (GT)"
    else
        echo "Processing: Imputed Array - Dosage Information (DS)"
    fi
    echo "Original samples: \${ORIG_SAMPLES}"
    echo "Original variants: \${ORIG_VARIANTS}"
    echo "=========================================="
    
    # Extract common samples and variants
    ${params.Plink2} \
        --${input_format} ${input_prefix} \
        --keep ${common_samples} \
        --extract ${common_variants} \
        --make-pgen \
        --out ${out_prefix} \
        --threads 8 2>&1 | tee ${out_prefix}.plink2.log
    
    # Count extracted samples and variants
    EXTRACTED_SAMPLES=\$(tail -n +2 ${out_prefix}.psam | wc -l)
    EXTRACTED_VARIANTS=\$(grep -v '^#' ${out_prefix}.pvar | wc -l)
    
    # Append to extraction summary report
    bash ${params.ScriptsDir}/append_extraction_stats.sh \
        "${params.OutDir}/extraction_summary.txt" \
        "${dataset_type}" \
        "\${ORIG_SAMPLES}" \
        "\${ORIG_VARIANTS}" \
        "\${EXTRACTED_SAMPLES}" \
        "\${EXTRACTED_VARIANTS}" \
        "${out_prefix}.plink2.log"
    """
}

params.PhenoFile = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/19.cov_pheno_prepare/cteph_agp3k.bbj.projection.pheno_df.csv'
params.CovarFile = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/19.cov_pheno_prepare/cteph_agp3k.bbj.projection.cov_df.no_age.csv'

// Perform additive model association analysis for all three datasets
process AssocPlink2 {
    executor 'slurm'
    queue 'gr10478b'
    time '6d'
    tag "assoc_${dataset_type}"

    publishDir "${params.OutDir}/03.assoc_results/${dataset_type}", mode: 'symlink'

    input:
    tuple val(dataset_type), file(pgen), file(pvar), file(psam) from common_subsets_out

    output:
    path "*.log"
    tuple val(dataset_type), path("*.glm.logistic") into assoc_results_out, assoc_results_out_2
    
    script:
    def pgen_prefix = pgen.baseName.replaceAll(/\.pgen$/, '')
    def pheno_name = "PHENO1"
    def covar_name = "SEX,PC1_AVG-PC10_AVG"
    
    // Set output prefix based on dataset type
    def dataset_label = ""
    if (dataset_type == 'wgs') {
        dataset_label = "WGS"
    } else if (dataset_type == 'array_gt') {
        dataset_label = "Imputed_Array_GT"
    } else {
        dataset_label = "Imputed_Array_DS"
    }
    def out_prefix = "cteph_agp3k.${dataset_label}.sex.10pc.additive"

    """
    echo "=========================================="
    echo "Association Analysis: ${dataset_label}"
    echo "Input: ${pgen_prefix}"
    echo "Phenotype: ${pheno_name}"
    echo "Covariates: ${covar_name}"
    echo "=========================================="
    
    ${params.Plink2} \
        --pfile ${pgen_prefix} \
        --pheno ${params.PhenoFile} \
        --pheno-name ${pheno_name} \
        --covar ${params.CovarFile} \
        --covar-name ${covar_name} \
        --glm cols=+beta omit-ref no-firth hide-covar \
        --out ${out_prefix} \
        --ci 0.95 \
        --threads 16
    
    echo "Association analysis completed for ${dataset_label}"
    """
}

// assoc_results_out.subscribe {
//     dataset_type, result_file ->
//     println "Association analysis completed for ${dataset_type}: ${result_file}"
// }

// Plot Manhattan and QQ plots for association results
process PlotMQQ {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag "plot_${dataset_type}"
    
    publishDir "${params.OutDir}/04.mqq_plots", mode: 'symlink'
    
    input:
    tuple val(dataset_type), file(glm_file) from assoc_results_out
    
    output:
    tuple val(dataset_type), path("*.mqq.png")
    
    script:
    def dataset_label = ""
    if (dataset_type == 'wgs') {
        dataset_label = "WGS"
    } else if (dataset_type == 'array_gt') {
        dataset_label = "Imputed_Array_GT"
    } else {
        dataset_label = "Imputed_Array_DS"
    }
    def out_plot = "cteph_agp3k.${dataset_label}.mqq.png"
    
    """
    source activate gwaslab
    
    echo "=========================================="
    echo "Plotting MQQ for: ${dataset_label}"
    echo "Input: ${glm_file}"
    echo "Output: ${out_plot}"
    echo "=========================================="
    
    python3 ${params.ScriptsDir}/plot_mqq.py \
        --input ${glm_file} \
        --output ${out_plot} \
        --build 38 \
        --sig-level 5e-8 \
        --dpi 400
    
    echo "MQQ plot completed for ${dataset_label}"
    """
}

params.VariantStat = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array_imputed_rev1/results/06.variant_stats_annotated/cteph_agp3k.imputed_array.sqc.variant_stats.annotated.tsv.gz'

// Pairwise comparison of association results
assoc_results_out_2
    .toList()
    .map { results ->
        // Create pairwise combinations from collected results
        def comparisons = []
        for (int i = 0; i < results.size(); i++) {
            for (int j = i + 1; j < results.size(); j++) {
                comparisons.add([results[i], results[j]])
            }
        }
        return comparisons
    }
    .flatMap()
    .set { pairwise_comparisons_ch }

process PairwiseComparison {
    executor 'slurm'
    queue 'gr10478b'
    time '6h'
    tag "compare_${pair1[0]}_vs_${pair2[0]}"
    
    publishDir "${params.OutDir}/05.pairwise_comparison/${pair1[0]}_vs_${pair2[0]}", mode: 'symlink'
    
    input:
    tuple val(pair1), val(pair2) from pairwise_comparisons_ch
    
    output:
    path "*.png"
    path "*.log"
    
    script:
    def dataset1_type = pair1[0]
    def dataset1_file = pair1[1]
    def dataset2_type = pair2[0]
    def dataset2_file = pair2[1]
    
    // Set display names
    def name1 = ""
    def name2 = ""
    
    if (dataset1_type == 'wgs') {
        name1 = "WGS"
    } else if (dataset1_type == 'array_gt') {
        name1 = "Imputed_Array_GT"
    } else {
        name1 = "Imputed_Array_DS"
    }
    
    if (dataset2_type == 'wgs') {
        name2 = "WGS"
    } else if (dataset2_type == 'array_gt') {
        name2 = "Imputed_Array_GT"
    } else {
        name2 = "Imputed_Array_DS"
    }
    
    def out_prefix = "cteph_agp3k.${name1}_vs_${name2}"
    def log_file = "${out_prefix}.comparison.log"
    
    """
    echo "=========================================="
    echo "Pairwise Comparison: ${name1} vs ${name2}"
    echo "Dataset 1: ${dataset1_file}"
    echo "Dataset 2: ${dataset2_file}"
    echo "=========================================="
    
    source activate cteph_geno_pro
    python3 ${params.ScriptsDir}/compare_sumstats.py \
        --sumstat1 ${dataset1_file} \
        --sumstat2 ${dataset2_file} \
        --name1 ${name1} \
        --name2 ${name2} \
        --variant-stats ${params.VariantStat} \
        --color-columns MAF_ALL IMPUTED_MARKER IMPUTED_R2 HWE_CTRL \
        --output-prefix ${out_prefix} \
        --sig-level 5e-8 \
        --log-file ${log_file} \
        --threads 16 \
        --tabix ${params.Tabix}
    
    echo "Pairwise comparison completed: ${name1} vs ${name2}"
    """
}


