params.vcfDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/02.mac1'
params.scriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/script'
params.outdir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/tuning.genotype'

Channel
    .from((21..22).collect { "chr${it}" } + ["PAR"])
    .set { chr_ch_dp }

Channel
    .from((21..22).collect { "chr${it}" } + ["PAR"])
    .set { chr_ch_gq }

process dp_check {
    executor 'slurm'
    queue 'gr10478b'
    time '6d'
    tag "${chr}"

    conda '/home/b/b37974/anaconda3/envs/cteph_geno_pro'

    publishDir "${params.outdir}/tmp", mode: 'symlink'

    input:
    val(vcfDir) from params.vcfDir
    val(scriptDir) from params.scriptDir
    val(chr) from chr_ch_dp

    output:
    tuple chr, file(csv) into dp_check_out

    script:
    vcf = "${vcfDir}/${chr}.pass.mac1.vcf.gz"
    csv = "${chr}.dp.check.csv"
    cteph_mata = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/cteph_jhrp4_info.csv"
    """
    source activate cteph_geno_pro
    python ${scriptDir}/dp.parameter.py --chrom ${chr} --vcf_path ${vcf} --metadata ${cteph_mata} --output_csv ${csv}
    """
}

process gq_check {
    executor 'slurm'
    queue 'gr10478b'
    time '6d'
    tag "${chr}"

    conda '/home/b/b37974/anaconda3/envs/cteph_geno_pro'

    publishDir "${params.outdir}/tmp", mode: 'symlink'

    input:
    val(vcfDir) from params.vcfDir
    val(scriptDir) from params.scriptDir
    val(chr) from chr_ch_gq

    output:
    tuple chr, file(csv) into gq_check_out

    script:
    vcf = "${vcfDir}/${chr}.pass.mac1.vcf.gz"
    csv = "${chr}.gq.check.csv"
    cteph_mata = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/cteph_jhrp4_info.csv"
    """
    source activate cteph_geno_pro
    python ${scriptDir}/gq.parameter.py --chrom ${chr} --vcf_path ${vcf} --metadata ${cteph_mata} --output_csv ${csv}
    """
}

process dp_check_vis {
    // executor 'slurm'
    // queue 'gr10478b'
    // time '6d'
    tag "${chr}"

    conda '/home/b/b37974/anaconda3/envs/cteph_geno_pro'

    publishDir "${params.outdir}/tmp", mode: 'symlink'

    input:
    tuple val(chr), file(csv) from dp_check_out

    output:
    file(pdf) into dp_check_vis_out

    script:
    pdf = "${chr}.dp.check.pdf"
    count_file = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/my_tools/countsSNP/counts/${chr}.pass.mac1.vcf.gz.count.txt"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/dp.parameter.vis.py --chrom ${chr} --count_file ${count_file} --dp_summary ${csv} --output_pdf ${pdf}
    """
}

process gq_check_vis {
    // executor 'slurm'
    // queue 'gr10478b'
    // time '6d'
    tag "${chr}"

    conda '/home/b/b37974/anaconda3/envs/cteph_geno_pro'

    publishDir "${params.outdir}/tmp", mode: 'symlink'

    input:
    tuple val(chr), file(csv) from gq_check_out

    output:
    file(pdf) into gq_check_vis_out

    script:
    pdf = "${chr}.gq.check.pdf"
    count_file = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/my_tools/countsSNP/counts/${chr}.pass.mac1.vcf.gz.count.txt"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/gq.parameter.vis.rev3.py --chrom ${chr} --count_file ${count_file} --gq_summary ${csv} --output_pdf ${pdf}
    """
}