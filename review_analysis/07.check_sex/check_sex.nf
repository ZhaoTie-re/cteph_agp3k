params.filePath = '/LARGE1/gr10478/platform/JHRPv4/workspace/pipeline/output/VQSR.v4'
params.sampleList = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/00.info/cteph_agp3k.lowfreq_common.samples.txt'
params.outDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/07.check_sex'
params.refGenome = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa'
params.plink2 = '/home/b/b37974/plink2'

chrXY = Channel.of('chrX', 'chrY')

process selectPASS {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outDir}/01.pass", mode: 'symlink'

    input:
    val(filePath) from params.filePath
    val(chr) from chrXY
    path(sampleList) from params.sampleList

    output:
    tuple chr, file(pass_vcf), file(pass_vcf_tbi) into selectPASS_out

    script:
    vcf = "${filePath}/all.VQSR3.${chr}.vcf.gz"
    pass_vcf = "${chr}.pass.vcf.gz"
    pass_vcf_tbi = "${chr}.pass.vcf.gz.tbi"

    """
    bcftools view ${vcf} --threads 8 -S ${sampleList} --force-samples -Ou | bcftools view --threads 8 -f"PASS" -Oz -o ${pass_vcf}
    bcftools index --threads 8 -t ${pass_vcf}
    """

}

process normAndSplit {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outDir}/02.norm_split", mode: 'symlink'

    input:
    tuple val(chr), path(vcf), path(vcf_tbi) from selectPASS_out
    val(refGenome) from params.refGenome

    output:
    tuple chr, file(norm_vcf), file(norm_vcf_tbi) into normAndSplit_out

    script:
    norm_vcf = "${chr}.pass.norm.vcf.gz"
    norm_vcf_tbi = "${chr}.pass.norm.vcf.gz.tbi"

    """
    bcftools norm -f ${refGenome} -m -any --threads 8 ${vcf} -Ou | \
    bcftools annotate --threads 8 --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o ${norm_vcf}
    bcftools index --threads 8 -t ${norm_vcf}
    """

}



// process vcfToPlink {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "${chr}"

//     publishDir "${params.outDir}/03.plink", mode: 'symlink'

//     input:
//     tuple val(chr), path(vcf), path(vcf_tbi) from normAndSplit_out
//     val(plink2) from params.plink2

//     output:
//     tuple val(chr), file("${chr}.pass.norm.bed"), file("${chr}.pass.norm.bim"), file("${chr}.pass.norm.fam") into vcfToPlink_out

//     script:
//     prefix = "${chr}.pass.norm"

//     """
//     ${plink2} --vcf ${vcf} --make-bed --out ${prefix} --threads 8 --allow-extra-chr
//     """

// }

// process mergePlink {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'

//     publishDir "${params.outDir}/04.merged", mode: 'symlink'

//     input:
//     val(plink2) from params.plink2
//     val(plink_files) from vcfToPlink_out.toSortedList { a, b -> 
//         def order = ['chrX': 0, 'chrY': 1]
//         order[a[0]] <=> order[b[0]]
//     }

//     output:
//     tuple file("chrXY.merged.bed"), file("chrXY.merged.bim"), file("chrXY.merged.fam") into mergePlink_out

//     script:
//     chrX_prefix = plink_files.find { it[0] == 'chrX' }[1].baseName
//     chrY_prefix = plink_files.find { it[0] == 'chrY' }[1].baseName
    
//     """
//     echo "${chrY_prefix}" > merge_list.txt
//     ${plink2} --bfile ${chrX_prefix} --merge-list merge_list.txt --make-bed --out chrXY.merged --threads 8
//     """

// }



