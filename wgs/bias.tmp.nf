params.vcfFilePath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/06.variant_filter'
params.plinkFilePath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/19.tommo_panel_filter'
params.infoPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info'
params.scriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/scripts'
params.outDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs'

Channel
    .from((1..22).collect { "chr${it}" })
    .set { chr_ch }

process BiasCalculation {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outDir}/21.geno_bias_calc", mode: 'symlink'

    input:
    val(chr) from chr_ch
    val(infoPath) from params.infoPath
    val(vcfFilePath) from params.vcfFilePath
    val(plinkFilePath) from params.plinkFilePath

    output:
    file("*.log")
    tuple val(chr), file("${chr}.coverage_transitions.tsv"), file("${chr}.stat_summary.tsv") into bias_result_ch

    script:
    bed_prefix = "${plinkFilePath}/cteph_agp3k.lowfreq_common"
    vcf = "${vcfFilePath}/${chr}.pass.mac1.vfilter.vcf.gz"
    info = "${infoPath}/cteph_agp3k_jhrpv4.xlsx"
    """
    source activate compute_env
    python ${params.scriptDir}/geno_miss_bias_chr_main.py \
        --plink-prefix ${bed_prefix} \
        --vcf ${vcf} \
        --info-xls ${info} \
        --chr ${chr} \
        --n-chunks 12 \
        --max-parallel 8 \
        --plink-threads 8 \
        --bcftools-threads 8 \
        --keep-temp
    """
}
