params.GTPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs'
params.InfoPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/cteph_agp3k_jhrpv4.rev1.xlsx'
params.OutDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/tuning.rv/results'
params.ScriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/tuning.rv/scripts'
params.Plink2 = '/home/b/b37974/plink2_alpha6/plink2'
params.Bcftools = '/home/b/b37974/bcftools/bcftools'
params.Tabix = '/home/b/b37974/htslib-1.9/tabix'
params.IDCol = 'ID'
params.GroupCol = 'OUTCOME1'
params.CaseValue = 'PH'
params.TdpCol = 'Target DP (JHRPv4)'
params.MdpCol = 'DP (JHRPv4)'

minac_ch = Channel.from(0..20)

process CALC_METRICS {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "minac:${min_ac}"
    publishDir "${params.OutDir}/00.qc_metrics/minac${min_ac}", mode: 'symlink'

    input:
    val min_ac from minac_ch

    output:
    tuple val(min_ac), file("sample_metrics.txt.gz"), file("variant_metrics.txt.gz"), file("variant_metrics.txt.gz.tbi") into mertics_ch
    file "calc_metrics.log"
    file "*.md"

    script:
    bed_prefix = "${params.GTPath}/18.tommo_panel_filter/cteph_agp3k.rare"
    info_file = params.InfoPath

    """
    bash ${params.ScriptDir}/run_calc_metrics.sh \\
        --bed-prefix "${bed_prefix}" \\
        --info-file "${info_file}" \\
        --out-sample "sample_metrics.txt" \\
        --out-variant "variant_metrics.txt" \\
        --script-dir "${params.ScriptDir}" \\
        --plink2 "${params.Plink2}" \\
        --tabix "${params.Tabix}" \\
        --id-col "${params.IDCol}" \\
        --group-col "${params.GroupCol}" \\
        --case-value "${params.CaseValue}" \\
        --tdp-col "${params.TdpCol}" \\
        --mdp-col "${params.MdpCol}" \\
        --min-ac ${min_ac} \\
        --threads 4
    """
}

process QC_SUMMARY {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "minac:${min_ac}"
    publishDir "${params.OutDir}/01.qc_summary/minac${min_ac}", mode: 'symlink'

    input:
    tuple val(min_ac), file(sample_metrics), file(variant_metrics), file(variant_metrics_tbi) from mertics_ch

    output:
    tuple val(min_ac), file("qc_summary_stats.minac${min_ac}.tsv") into qc_summary_ch
    file "qc_summary_plots.pdf"
    file "*.md"

    script:
    """
    python3 ${params.ScriptDir}/qc_summary.py \
        --sample-metrics ${sample_metrics} \
        --variant-metrics ${variant_metrics} \
        --min-ac ${min_ac} \
        --out-tsv qc_summary_stats.minac${min_ac}.tsv \
        --out-md qc_summary_methods.md

    python3 ${params.ScriptDir}/plot_qc.py \
        --sample-metrics ${sample_metrics} \
        --variant-metrics ${variant_metrics} \
        --qc-stats qc_summary_stats.minac${min_ac}.tsv \
        --min-ac ${min_ac} \
        --out-pdf qc_summary_plots.pdf \
        --hist-stat density
    """
}

process MERGE_ALL_SUMMARIES {
    executor 'slurm'
    queue 'gr10478b'
    time '10m'
    publishDir "${params.OutDir}/02.qc_collect", mode: 'symlink'

    input:
    file tsvs from qc_summary_ch.collect { it[1] }

    output:
    file "all_qc_summary_stats.tsv" into all_qc_summary_stats_ch

    script:
    """
    # 1. Merge all TSVs, keeping the header from the first file only
    awk 'NR==1 || FNR > 1' ${tsvs} > temp_merged.tsv

    # 2. Sort by MinAC_Threshold (ascending) using Python
    python3 -c "import pandas as pd; df = pd.read_csv('temp_merged.tsv', sep='\\t'); df.sort_values('MinAC_Threshold').to_csv('all_qc_summary_stats.tsv', sep='\\t', index=False)"
    """
}

process PLOT_SUMMARY_TREND {
    executor 'slurm'
    queue 'gr10478b'
    time '10m'
    publishDir "${params.OutDir}/02.qc_collect", mode: 'symlink'

    input:
    file qc_stats from all_qc_summary_stats_ch

    output:
    file "qc_trend_plots.pdf"

    script:
    fam_file = "${params.GTPath}/18.tommo_panel_filter/cteph_agp3k.rare.fam"
    """
    N_SAMPLES=\$(wc -l < ${fam_file})
    python3 ${params.ScriptDir}/plot_summary_trend.py \
        --qc-stats ${qc_stats} \
        --out-pdf qc_trend_plots.pdf \
        --sample-n \${N_SAMPLES}
    """
}


