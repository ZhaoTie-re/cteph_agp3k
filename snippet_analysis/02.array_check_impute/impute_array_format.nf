params.arrayPath = '/LARGE0/gr10478/project/pulmonary_hypertension/DATA_SOURCE/impute20251029/JSA'
params.scriptPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/scripts'
params.outPath = '//LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/02.array_check_impute'
params.refGenome = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa'

Channel
    .from((1..22).collect { "chr${it}" })
    .set { chr_ch }

process normalizeAndSetId {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outPath}/01.norm_setid", mode: 'symlink'

    input:
    val(chr) from chr_ch
    val(arrayPath) from params.arrayPath
    val(refGenome) from params.refGenome

    output:
    tuple val(chr), file(output_vcf), file(output_vcf_tbi) into norm_setid_out

    script:
    input_bcf = "${arrayPath}/impute.${chr}.imputed.bcf"
    norm_bcf = "${chr}.norm.bcf"
    output_vcf = "${chr}.norm.setid.vcf.gz"
    output_vcf_tbi = "${chr}.norm.setid.vcf.gz.tbi"

    """
    # Step 1: Normalize BCF file
    bcftools norm -m- --fasta-ref ${refGenome} -c s --threads 8 ${input_bcf} -Ob -o ${norm_bcf}
    
    # Step 2: Set variant IDs and convert to VCF.gz
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' ${norm_bcf} -Oz -o ${output_vcf}
    
    # Step 3: Index the output VCF.gz file
    bcftools index --threads 4 -t ${output_vcf}
    """
}

process vcfToPlink {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outPath}/02.vcf_to_plink", mode: 'symlink'

    input:
    tuple val(chr), file(vcf), file(vcf_tbi) from norm_setid_out

    output:
    tuple val(chr), file("${chr}.bed"), file("${chr}.bim"), file("${chr}.fam") into plink_chr_out

    script:
    """
    # Convert VCF to PLINK format with double-id (FID = IID)
    export PATH=/home/b/b37974/:$PATH
    plink2 \
        --vcf ${vcf} \
        --make-bed \
        --double-id \
        --out ${chr} \
        --threads 8
    """
}

// Collect all chromosome files for merging
plink_chr_out
    .toList()
    .map { entries ->
        def chr_order = (1..22).collect { "chr${it}" }
        def filtered = entries.findAll { it[0] in chr_order }
        def sorted = filtered.sort { a, b ->
            Integer.parseInt(a[0].replace("chr", "")) <=> Integer.parseInt(b[0].replace("chr", ""))
        }
        def paths = sorted.collect { it[1].toString().replaceAll(/\.bed$/, '') }
        return paths.join('\n')
    }
    .set { pmerge_list_content }

process writePmergeList {
    executor 'local'
    time '1h'
    
    publishDir "${params.outPath}/02.vcf_to_plink", mode: 'copy'

    input:
    val(list_text) from pmerge_list_content

    output:
    file("pmerge_list.txt") into pmerge_list_ch

    script:
    """
    echo "${list_text}" > pmerge_list.txt
    """
}

process mergePlink {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "merge_all_chr"

    publishDir "${params.outPath}/02.vcf_to_plink", mode: 'symlink'

    input:
    file(pmerge_list) from pmerge_list_ch

    output:
    tuple file("cteph_agp3k.array.bed"), file("cteph_agp3k.array.bim"), file("cteph_agp3k.array.fam") into merged_plink_out

    script:
    """
    # Merge all chromosome PLINK files
    export PATH=/home/b/b37974/:$PATH
    plink2 \
        --pmerge-list ${pmerge_list} bfile \
        --make-bed \
        --out cteph_agp3k.array \
        --threads 16
    """
}

process extractCommonVariantsAndSamples {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "extract_common"

    publishDir "${params.outPath}/03.extract_common", mode: 'symlink'

    input:
    tuple file(array_bed), file(array_bim), file(array_fam) from merged_plink_out

    output:
    file('*.log') into extract_log_out
    tuple file('cteph_agp3k.array.common.bed'), file('cteph_agp3k.array.common.bim'), file('cteph_agp3k.array.common.fam') into array_common_out
    tuple file('cteph_agp3k.wgs.common.bed'), file('cteph_agp3k.wgs.common.bim'), file('cteph_agp3k.wgs.common.fam') into wgs_common_out

    script:
    wgs_prefix = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/19.tommo_panel_filter/cteph_agp3k.lowfreq_common"
    case_prefix = "PHOM"
    log_file = "extract_common_variants_samples.log"
    
    """
    # Set PATH for plink2
    export PATH=/home/b/b37974/:\$PATH
    
    # Create log file
    echo "=== Array vs WGS Common Variants and Samples Extraction ===" > ${log_file}
    echo "Date: \$(date)" >> ${log_file}
    echo "Array input prefix: ${array_bed.baseName}" >> ${log_file}
    echo "WGS input prefix: ${wgs_prefix}" >> ${log_file}
    echo "Case prefix: ${case_prefix}" >> ${log_file}
    echo "" >> ${log_file}
    
    # Check input files
    echo "--- Input File Check ---" >> ${log_file}
    echo "Array BIM file: \$(wc -l ${array_bim} | cut -d' ' -f1) variants" >> ${log_file}
    echo "Array FAM file: \$(wc -l ${array_fam} | cut -d' ' -f1) samples" >> ${log_file}
    echo "WGS BIM file: \$(wc -l ${wgs_prefix}.bim | cut -d' ' -f1) variants" >> ${log_file}
    echo "WGS FAM file: \$(wc -l ${wgs_prefix}.fam | cut -d' ' -f1) samples" >> ${log_file}
    echo "" >> ${log_file}
    
    # Extract variant IDs
    echo "--- Extracting Variant IDs ---" >> ${log_file}
    cut -f2 ${array_bim} | sort > array_variants.txt
    cut -f2 ${wgs_prefix}.bim | sort > wgs_variants.txt
    
    # Find common variants
    comm -12 array_variants.txt wgs_variants.txt > common_variants.txt
    comm -23 array_variants.txt wgs_variants.txt > array_only_variants.txt
    comm -13 array_variants.txt wgs_variants.txt > wgs_only_variants.txt
    
    echo "Array unique variants: \$(wc -l array_variants.txt | cut -d' ' -f1)" >> ${log_file}
    echo "WGS unique variants: \$(wc -l wgs_variants.txt | cut -d' ' -f1)" >> ${log_file}
    echo "Common variants: \$(wc -l common_variants.txt | cut -d' ' -f1)" >> ${log_file}
    echo "Array-only variants: \$(wc -l array_only_variants.txt | cut -d' ' -f1)" >> ${log_file}
    echo "WGS-only variants: \$(wc -l wgs_only_variants.txt | cut -d' ' -f1)" >> ${log_file}
    echo "" >> ${log_file}
    
    # Extract sample IDs (using FID)
    echo "--- Extracting Sample IDs ---" >> ${log_file}
    cut -f1 ${array_fam} | sort > array_samples.txt
    cut -f1 ${wgs_prefix}.fam | sort > wgs_samples.txt
    
    # Find common samples
    comm -12 array_samples.txt wgs_samples.txt > common_samples.txt
    comm -23 array_samples.txt wgs_samples.txt > array_only_samples.txt
    comm -13 array_samples.txt wgs_samples.txt > wgs_only_samples.txt
    
    echo "Array unique samples: \$(wc -l array_samples.txt | cut -d' ' -f1)" >> ${log_file}
    echo "WGS unique samples: \$(wc -l wgs_samples.txt | cut -d' ' -f1)" >> ${log_file}
    echo "Common samples: \$(wc -l common_samples.txt | cut -d' ' -f1)" >> ${log_file}
    echo "Array-only samples: \$(wc -l array_only_samples.txt | cut -d' ' -f1)" >> ${log_file}
    echo "WGS-only samples: \$(wc -l wgs_only_samples.txt | cut -d' ' -f1)" >> ${log_file}
    echo "" >> ${log_file}
    
    # Count cases and controls in common samples
    echo "--- Case/Control Analysis in Common Samples ---" >> ${log_file}
    grep "^${case_prefix}" common_samples.txt | wc -l > common_cases_count.txt
    grep -v "^${case_prefix}" common_samples.txt | wc -l > common_controls_count.txt
    
    echo "Common cases (${case_prefix}*): \$(cat common_cases_count.txt)" >> ${log_file}
    echo "Common controls: \$(cat common_controls_count.txt)" >> ${log_file}
    echo "" >> ${log_file}
    
    # Create keep files for both datasets
    echo "--- Creating Keep Files ---" >> ${log_file}
    # For PLINK, keep file needs FID and IID (both columns)
    awk '{print \$1, \$1}' common_samples.txt > common_samples_keep.txt
    echo "Keep file created with \$(wc -l common_samples_keep.txt | cut -d' ' -f1) samples" >> ${log_file}
    echo "" >> ${log_file}
    
    # Extract common variants and samples from Array data
    echo "--- Extracting Array Common Data ---" >> ${log_file}
    plink2 \
        --bfile ${array_bed.baseName} \
        --extract common_variants.txt \
        --keep common_samples_keep.txt \
        --make-bed \
        --out cteph_agp3k.array.common \
        --threads 8
    
    echo "Array extraction completed" >> ${log_file}
    echo "Final Array variants: \$(wc -l cteph_agp3k.array.common.bim | cut -d' ' -f1)" >> ${log_file}
    echo "Final Array samples: \$(wc -l cteph_agp3k.array.common.fam | cut -d' ' -f1)" >> ${log_file}
    echo "" >> ${log_file}
    
    # Extract common variants and samples from WGS data
    echo "--- Extracting WGS Common Data ---" >> ${log_file}
    plink2 \
        --bfile ${wgs_prefix} \
        --extract common_variants.txt \
        --keep common_samples_keep.txt \
        --make-bed \
        --out cteph_agp3k.wgs.common \
        --threads 8
    
    echo "WGS extraction completed" >> ${log_file}
    echo "Final WGS variants: \$(wc -l cteph_agp3k.wgs.common.bim | cut -d' ' -f1)" >> ${log_file}
    echo "Final WGS samples: \$(wc -l cteph_agp3k.wgs.common.fam | cut -d' ' -f1)" >> ${log_file}
    echo "" >> ${log_file}
    
    # Final summary
    echo "--- Final Summary ---" >> ${log_file}
    echo "Process completed successfully at: \$(date)" >> ${log_file}
    echo "Output files:" >> ${log_file}
    echo "  - cteph_agp3k.array.common.{bed,bim,fam}" >> ${log_file}
    echo "  - cteph_agp3k.wgs.common.{bed,bim,fam}" >> ${log_file}
    echo "Log file: ${log_file}" >> ${log_file}
    """
}

process AssocPlink2ArrayModel {
    executor 'slurm'
    queue 'gr10478b'
    time '6d'
    tag "assoc_plink2_array: additive"

    publishDir "${params.outPath}/04.array_assoc_result", mode: 'symlink'

    input:
    tuple file(array_bed), file(array_bim), file(array_fam) from array_common_out

    output:
    file("*.log")
    tuple val("additive"), file("*.glm.logistic") into array_assoc_out

    script:
    PHENO_FILE = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/20.cov_pheno_prepare/cteph_agp3k.bbj.projection.pheno_df.csv"
    PHENO_NAME = "PHENO1"
    COVAR_FILE = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/20.cov_pheno_prepare/cteph_agp3k.bbj.projection.cov_df.no_age.csv"
    COVAR_NAME = "SEX, PC1_AVG-PC10_AVG"
    OUT_PR = "cteph_agp3k.array.sex.10pc.additive"

    """
    export PATH=/home/b/b37974/:\$PATH
    plink2 \\
        --bfile ${array_bed.baseName} \\
        --pheno ${PHENO_FILE} \\
        --pheno-name ${PHENO_NAME} \\
        --covar ${COVAR_FILE} \\
        --covar-name ${COVAR_NAME} \\
        --glm omit-ref no-firth hide-covar \\
        --out ${OUT_PR} \\
        --ci 0.95 \\
        --threads 16
    """
}




