params.VcfPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/wgs/results/01.pass_norm_setid'
params.OutDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/tuning.variants.v5/results'
params.ScriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/tuning.variants.v5/scripts'
params.Bcftools = '/home/b/b37974/bcftools/bcftools'

params.MQThreshold = 58.75
params.VQSLODThreshold = 10


// Create channel for all chromosomes including X, Y, and PAR
Channel
    .from((1..22).collect { "chr${it}" } + ["chrX", "chrY", "PAR"])
    .set { chr_ch }


// Process 0: Extract MQ and VQSLOD values from VCF to TSV files
process extractMetrics {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"
    publishDir "${params.OutDir}/00.extracted_metrics", mode: 'symlink'
    
    input:
    val chr from chr_ch
    
    output:
    tuple val(chr), path("${chr}.metrics.tsv") into metrics_ch
    
    script:
    """
    vcf_file="${params.VcfPath}/${chr}.selected.pass.norm_split.setid.vcf.gz"
    
    python3 ${params.ScriptDir}/extract_metrics.py \\
        \$vcf_file \\
        ${chr} \\
        ${params.Bcftools} \\
        4 \\
        ${chr}.metrics.tsv
    """
}

// Split channel for parallel processing
metrics_ch.into { metrics_for_plot; metrics_for_count }


// Process 1: Generate distribution plots for MQ and VQSLOD from TSV
process plotDistributions {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"
    publishDir "${params.OutDir}/01.distribution_plots", mode: 'symlink'
    
    input:
    tuple val(chr), path(metrics_tsv) from metrics_for_plot
    
    output:
    tuple val(chr), path("${chr}.MQ_distribution.png"), path("${chr}.VQSLOD_distribution.png") into plot_results
    
    script:
    """
    python3 ${params.ScriptDir}/plot_distributions.py \\
        ${metrics_tsv} \\
        ${chr} \\
        ${params.MQThreshold} \\
        ${params.VQSLODThreshold} \\
        ${chr}
    """
}


// Process 2: Count variants passing thresholds from TSV
process countVariants {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"
    publishDir "${params.OutDir}/02.variant_counts", mode: 'symlink'
    
    input:
    tuple val(chr), path(metrics_tsv) from metrics_for_count
    
    output:
    path "${chr}.stats.txt" into individual_stats
    
    script:
    """
    python3 ${params.ScriptDir}/count_variants.py \\
        ${metrics_tsv} \\
        ${chr} \\
        ${params.MQThreshold} \\
        ${params.VQSLODThreshold} \\
        ${chr}.stats.txt
    """
}


// Process 3: Merge all chromosome statistics into a single summary file
process mergeChrStats {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    publishDir "${params.OutDir}/03.merged_summary", mode: 'symlink'
    
    input:
    path stats_files from individual_stats.collect()
    
    output:
    path "all_chromosomes_summary.txt"
    
    script:
    """
    python3 ${params.ScriptDir}/merge_stats.py \\
        all_chromosomes_summary.txt \\
        ${stats_files}
    """
}

