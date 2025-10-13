params.gtPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs'
params.snpeffDir = '/LARGE1/gr10478/platform/JHRPv4/workspace/pipeline/output/snpEff.v4.index'
params.scriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest/scripts'
params.outDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest/results'

process RVtestPrepare {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outDir}/01.rvtest_prepare", mode: 'symlink'

    input:
    val gtPath from params.gtPath

    output:
    file('*.log')
    tuple file("*.vcf.gz"), file("*.vcf.gz.tbi") into rvtest_vcf_ch
    tuple file("*.pheno_df.csv"), file("*.cov_df.no_age.csv") into rvtest_pheno_covar_ch
    file("refFlat.hg38.nochr.txt.gz") into rvtest_refflat_ch

    script:
    bed_prefix = "${gtPath}/19.tommo_panel_filter/cteph_agp3k.rare"
    pheno = "${gtPath}/20.cov_pheno_prepare/cteph_agp3k.bbj.projection.pheno_df.csv"
    covar = "${gtPath}/20.cov_pheno_prepare/cteph_agp3k.bbj.projection.cov_df.no_age.csv"
    refflat = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/anno_raw/refFlat.hg38.txt.gz"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/rvtest_prepare_main.py \
        --bed-prefix ${bed_prefix} \
        --pheno-path ${pheno} \
        --covar-path ${covar} \
        --refflat-path ${refflat} \
        --threads 32 \
        --verbose
    """
}

process SnpEffAnnotate {
    executor 'slurm'
    queue 'gr10478b'
    time '48h'

    publishDir "${params.outDir}/02.snpeff_annotate", mode: 'symlink'

    input:
    val snpeffDir from params.snpeffDir
    tuple file(vcf), file(vcf_tbi) from rvtest_vcf_ch

    output:
    file('*.log')
    tuple file("*.snpeff.vcf.gz"), file("*.snpeff.vcf.gz.tbi") into snpeff_annotated_vcf_ch

    script:
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/snpeff_anno_main.py \
        --vcf-path ${vcf} \
        --snpeff-dir ${snpeffDir} \
        --parallel \
        --max-workers 16 \
        --threads 16 \
        --keep-cache \
    """
}




// macmin_opts = Channel.from([
//     ["no_macmin", ""],               // 不指定--siteMACMin
//     ["macmin_5",  "--siteMACMin 5"], // 指定为5
//     ["macmin_10", "--siteMACMin 10"] // 指定为10
// ])

// process RVtestRun {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '48h'
//     tag "${macmin_tag}"

//     publishDir "${params.outDir}/02.rvtest_run/${macmin_tag}", mode: 'symlink'

//     input:
//     tuple val(macmin_tag), val(macmin_opt) from macmin_opts
//     tuple file(vcf), file(vcf_tbi) from rvtest_vcf_ch
//     tuple file(pheno), file(covar) from rvtest_pheno_covar_ch
//     file(refflat) from rvtest_refflat_ch

//     output:
//     file("*.log")
//     tuple val(macmin_tag), file("cteph_agp3k.rare.${macmin_tag}.SkatO.assoc") into skato_result_ch

//     script:
//     """
//     export PATH=/home/b/b37974/rvtests/executable/:$PATH
//     rvtest \
//         --inVcf ${vcf} \
//         --pheno ${pheno} --pheno-name pheno1 \
//         --covar ${covar} --covar-name sex,pc1_avg,pc2_avg,pc3_avg,pc4_avg,pc5_avg,pc6_avg,pc7_avg,pc8_avg,pc9_avg,pc10_avg \
//         --geneFile ${refflat} \
//         --kernel skato \
//         --out cteph_agp3k.rare.${macmin_tag} \
//         --noweb \
//         --numThread 16 \
//         ${macmin_opt}
//     """
// }