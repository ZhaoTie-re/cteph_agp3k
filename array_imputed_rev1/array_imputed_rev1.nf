nextflow.enable.dsl=1

params.ImputedArrayPath = '/LARGE0/gr10478/project/Pulmonary_Hypertension/DATA_SOURCE/impute20251029/JSA'
params.OutDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array_imputed_rev1/results'
params.ScriptsDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array_imputed_rev1/scripts'
params.SampleSelectList = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array_rev1/info/cteph_agp3k.v4.array.ls'
params.NagasakiPipelinePath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline'
params.Bcftools = '/home/b/b37974/bcftools/bcftools'
params.Plink = '/home/b/b37974/plink'
params.Plink2 = '/home/b/b37974/plink2'
params.Tabix = '/home/b/b37974/htslib-1.9/tabix'
params.SampleInfo = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/cteph_agp3k_jhrpv4.rev1.xlsx'
params.SampleIDColumn = 'ID'
params.SampleSexColumn = 'Sex'
params.OutcomeColumn = 'OUTCOME2'
params.CaseValue = 'CTEPH'
params.PhenoColumn = 'PHENO1'

// params.SampleQCMode = 'TRUE'
// params.VariantQCMode = 'TRUE'
// params.GTHardCallMode = 'TRUE'

// Create channel for chr1-chr22 imputed array BCF files
imputed_array_ch = Channel
    .from(1..22)
    .map { chr -> 
        def chrName = "chr${chr}"
        def bcfFile = file("${params.ImputedArrayPath}/impute.${chrName}.imputed.bcf")
        return tuple(chrName, bcfFile)
    }

process indexBcf {
    tag "${chr}"
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    
    publishDir "${params.OutDir}/00.raw", mode: 'symlink'
    
    input:
    set val(chr), path(bcf) from imputed_array_ch
    
    output:
    set val(chr), path(bcf), path("${bcf}.csi") into indexed_array_ch
    
    script:
    """
    ${params.Bcftools} index -f -c --threads 4 ${bcf}
    """
}

process normalizeBcf {
    tag "${chr}"
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    
    publishDir "${params.OutDir}/01.normalized", mode: 'symlink'
    
    input:
    set val(chr), path(bcf), path(csi) from indexed_array_ch
    
    output:
    set val(chr), path("${chr}.normalized.vcf.gz"), path("${chr}.normalized.vcf.gz.tbi") into normalized_ch, normalized_ch_2
    
    script:
    """
    ${params.Bcftools} view \
        -S ${params.SampleSelectList} \
        --force-samples \
        --threads 8 \
        ${bcf} \
    | ${params.Bcftools} norm \
        --multiallelics -any \
        --fasta-ref ${params.NagasakiPipelinePath}/data/hs38DH.fa \
        --check-ref s \
        --threads 8 \
        -Ou \
    | ${params.Bcftools} annotate \
        --set-id '%CHROM:%POS:%REF:%ALT' \
        --threads 8 \
        -Oz \
        -o ${chr}.normalized.vcf.gz
    
    ${params.Bcftools} index -t --threads 8 ${chr}.normalized.vcf.gz
    """
}

process convertToPlink {
    tag "${chr}"
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    
    publishDir "${params.OutDir}/02.plink", mode: 'symlink'
    
    input:
    set val(chr), path(vcf), path(tbi) from normalized_ch
    
    output:
    set val(chr), path("${chr}.pgen"), path("${chr}.pvar"), path("${chr}.psam") into plink_ch
    
    script:
    """
    # Convert VCF to PLINK2 pgen format with both GT and DS
    # --double-id: Set FID = IID for all samples
    ${params.Plink2} \
        --vcf ${vcf} dosage=DS \
        --make-pgen \
        --double-id \
        --out ${chr} \
        --threads 8
    
    # Update sex and pheno information using external Python script
    python3 ${params.ScriptsDir}/update_psam.py \
        --psam ${chr}.psam \
        --sample-info ${params.SampleInfo} \
        --id-col ${params.SampleIDColumn} \
        --sex-col ${params.SampleSexColumn} \
        --outcome-col ${params.OutcomeColumn} \
        --case-value ${params.CaseValue} \
        --output ${chr}.psam
    """
}

plink_ch
    .toList()
    .map { entries ->
        def chr_order = (1..22).collect { "chr${it}" }
        def filtered = entries.findAll { it[0] in chr_order }
        def sorted = filtered.sort { a, b ->
            Integer.parseInt(a[0].replace("chr", "")) <=> Integer.parseInt(b[0].replace("chr", ""))
        }
        def paths = sorted.collect { it[1].toString().replaceAll(/\.pgen$/, '') }
        return paths.join('\n')
    }
    .set { pmerge_list_content }

process writePmergeList {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    publishDir "${params.OutDir}/03.plink_merge", mode: 'symlink'

    input:
    val(list_text) from pmerge_list_content

    output:
    file("pmerge_lst.txt") into pmerge_lst_ch

    script:
    """
    echo "${list_text}" > pmerge_lst.txt
    """
}

process pmerge {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "pmerge"

    publishDir "${params.OutDir}/03.plink_merge", mode: 'symlink'

    input:
    file(pmerge_list) from pmerge_lst_ch

    output:
    tuple file("${out_prefix}.pgen"), file("${out_prefix}.pvar"), file("${out_prefix}.psam") into pmerge_out_qc

    script:
    out_prefix = "cteph_agp3k.imputed_array"
    """
    ${params.Plink2} \
        --pmerge-list ${pmerge_list} pfile \
        --threads 8 \
        --make-pgen \
        --out ${out_prefix}
    """
}

process sampleQC {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "sample_qc"
    
    // when:
    // params.SampleQCMode == 'TRUE'
    
    publishDir "${params.OutDir}/04.sample_qc", mode: 'symlink'
    
    input:
    tuple file(pgen), file(pvar), file(psam) from pmerge_out_qc
    
    output:
    tuple file("*.sqc.pgen"), file("*.sqc.pvar"), file("*.sqc.psam") into sample_qc_out, sample_qc_out_2, sample_qc_out_3
    file("*.log")
    file("*.smiss")
    
    script:
    in_prefix = pgen.baseName
    out_prefix = "${in_prefix}.sqc"
    """
    # Calculate sample missingness
    ${params.Plink2} \
        --pfile ${in_prefix} \
        --missing sample-only \
        --threads 8 \
        --out ${in_prefix}
    
    # Filter samples based on missingness threshold
    ${params.Plink2} \
        --pfile ${in_prefix} \
        --mind 0.02 \
        --make-pgen \
        --threads 8 \
        --out ${out_prefix}
    """
}

process variantStats {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    tag "variant_stats"
    
    publishDir "${params.OutDir}/05.variant_stats", mode: 'symlink'
    
    input:
    tuple file(pgen), file(pvar), file(psam) from sample_qc_out
    
    output:
    tuple file("*.variant_stats.tsv.gz"), file("*.variant_stats.tsv.gz.tbi") into variant_stats_out
    
    script:
    in_prefix = pgen.baseName
    out_prefix = "${in_prefix}.variant_stats"
    """
    # Run variant statistics calculation script
    python3 ${params.ScriptsDir}/calc_variant_stats.py \
        --pfile ${in_prefix} \
        --plink2 ${params.Plink2} \
        --tabix ${params.Tabix} \
        --pheno-col ${params.PhenoColumn} \
        --threads 12 \
        --output ${out_prefix}.tsv.gz
    """
}

process addImputationInfo {
    executor 'slurm'
    queue 'gr10478b'
    time '48h'
    tag "add_imputation"
    
    publishDir "${params.OutDir}/06.variant_stats_annotated", mode: 'symlink'
    
    input:
    tuple file(stats_tsv), file(stats_tbi) from variant_stats_out
    tuple file(pgen), file(pvar), file(psam) from sample_qc_out_2
    
    output:
    tuple file("*.annotated.tsv.gz"), file("*.annotated.tsv.gz.tbi") into variant_stats_annotated
    
    script:
    in_prefix = stats_tsv.baseName.replaceAll(/\.variant_stats\.tsv$/, '')
    out_prefix = "${in_prefix}.variant_stats.annotated"
    """
    # Add imputation information from pvar file
    python3 ${params.ScriptsDir}/add_imputation_info.py \
        --variant-stats ${stats_tsv} \
        --pvar ${pvar} \
        --tabix ${params.Tabix} \
        --threads 12 \
        --output ${out_prefix}.tsv.gz
    """
}

process filterVariants {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    tag "filter_variants"
    
    publishDir "${params.OutDir}/07.variant_filter", mode: 'symlink'
    
    input:
    tuple file(annotated_tsv), file(annotated_tbi) from variant_stats_annotated
    
    output:
    file("*.exclude.txt") into variant_exclude_list
    file("*.upset.png") optional true
    file("*.upset.pdf") optional true
    file("*.summary.txt")
    
    script:
    in_prefix = annotated_tsv.baseName.replaceAll(/\.annotated\.tsv$/, '')
    out_prefix = "${in_prefix}.filtered"
    """
    source activate cteph_geno_pro 
    # Filter variants based on quality metrics
    python3 ${params.ScriptsDir}/filter_variants.py \
        --input ${annotated_tsv} \
        --config ${params.ScriptsDir}/variant_filter_config.json \
        --output-prefix ${out_prefix} \
        --chunk-size 500000
    """
}

process variantQC {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "exclude_variants"
    
    publishDir "${params.OutDir}/08.variant_qc", mode: 'symlink'
    
    input:
    tuple file(pgen), file(pvar), file(psam) from sample_qc_out_3
    file(exclude_list) from variant_exclude_list
    
    output:
    tuple file("*.vqc.pgen"), file("*.vqc.pvar"), file("*.vqc.psam") into variant_qc_out
    file("*.log")
    
    script:
    in_prefix = pgen.baseName.replaceAll(/\.pgen$/, '')
    out_prefix = "${in_prefix}.vqc"
    """
    # Exclude variants from plink fileset
    ${params.Plink2} \
        --pfile ${in_prefix} \
        --exclude ${exclude_list} \
        --make-pgen \
        --out ${out_prefix} \
        --threads 8
    """
}

process eraseGenotypeHardcall {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "gt_hardcall"
    
    publishDir "${params.OutDir}/09.gt_hardcall", mode: 'symlink'
    
    input:
    tuple file(pgen), file(pvar), file(psam) from variant_qc_out
    
    output:
    tuple file("*.gt.pgen"), file("*.gt.pvar"), file("*.gt.psam") into gt_hardcall_out
    file("*.log")
    
    script:
    in_prefix = pgen.baseName.replaceAll(/\.pgen$/, '')
    out_prefix = "${in_prefix}.gt"
    """
    # Erase dosage information, keep only hard-called genotypes (GT)
    ${params.Plink2} \
        --pfile ${in_prefix} \
        --make-pgen erase-dosage \
        --out ${out_prefix} \
        --threads 8
    """
}


