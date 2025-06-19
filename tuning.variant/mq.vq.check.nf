params.vcfDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/05.gt_norm'
params.scriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/script'
params.outdir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/tuning.variant'

Channel
    .from((1..22).collect { "chr${it}" } + ["PAR"])
    .set { chr_ch }

// Channel
//     .from((1..22).collect { "chr${it}" })
//     .set { chr_ch }

process mq_vq_check {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    conda '/home/b/b37974/anaconda3/envs/cteph_geno_pro'

    publishDir "${params.outdir}/tmp", mode: 'symlink'

    input:
    val(vcfDir) from params.vcfDir
    val(scriptDir) from params.scriptDir
    val(chr) from chr_ch

    output:
    tuple chr, file(pdf) into mq_vq_check_out

    script:
    vcf = "${vcfDir}/${chr}.pass.mac1.setid.gt_af.norm.vcf.gz"
    pdf = "${chr}.mq.vq.check.pdf"
    """
    python ${scriptDir}/variant.parameter.py --chrom ${chr} --vcf_path ${vcf} --MQ 58.75 --VQSLOD 10 --output_pdf ${pdf}
    """
}

