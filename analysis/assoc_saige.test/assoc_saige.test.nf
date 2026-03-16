nextflow.enable.dsl=2

// Define pipeline parameters
params.WGSPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs'
params.InfoPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info'
params.ContainerPath = '/home/b/b37974/simg'
params.ScriptsPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_saige.test/scripts'
params.OutPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_saige.test/results'
params.Covariates = 'SEX,PC1_AVG,PC2_AVG,PC3_AVG,PC4_AVG,PC5_AVG,PC6_AVG,PC7_AVG,PC8_AVG,PC9_AVG,PC10_AVG'
params.GRM = 'full' // Options: 'sparse', 'full'
params.FilterMAC = true // Options: true, false
params.MACThreshold = 20

/*
 * Process 0: Filter by MAC (Optional)
 * Inputs: PLINK files
 * Outputs: Filtered PLINK files
 */
process FILTER_BY_MAC {
    executor 'slurm'
    queue 'gr10478b'
    time '6h'
    publishDir "${params.OutPath}/saige_pre_filter", mode: 'symlink'

    input:
    tuple val(prefix), path(bed), path(bim), path(fam)

    output:
    tuple val(prefix), path("*.bed"), path("*.bim"), path("*.fam"), emit: filtered_plink

    script:
    def out_prefix = "${prefix}.mac${params.MACThreshold}"
    """
    # Filter by MAC using PLINK2
    /home/b/b37974/plink2_alpha6/plink2 \
        --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --mac ${params.MACThreshold} \
        --make-bed \
        --out ${out_prefix} \
        --threads 8
    """
}

/*
 * Process 0.5: LD Pruning
 * Inputs: PLINK files
 * Outputs: Pruned marker list (prune.in)
 */
process LD_PRUNING {
    executor 'slurm'
    queue 'gr10478b'
    time '6h'
    publishDir "${params.OutPath}/saige_pre_filter", mode: 'symlink'

    input:
    tuple val(prefix), path(bed), path(bim), path(fam)
    path(high_ld)

    output:
    path("*.prune.in"), emit: prune_in

    script:
    def out_prefix = "${prefix}.ld_pruning"
    """
    # Check chromosome format in BIM file
    # Get the detailed first chromosome name (e.g. '1' or 'chr1')
    BIM_CHR=\$(head -n 1 ${bim} | awk '{print \$1}')
    
    # Check chromosome format in High LD file
    # High LD file typically has columns: CHR START END ...
    # We assume standard headerless or specific format. Usually first column is CHR.
    LD_CHR=\$(head -n 1 ${high_ld} | awk '{print \$1}')
    
    # Determine if modification is needed
    # Case 1: BIM has 'chr', LD does not -> Add 'chr' to LD file
    # Case 2: BIM does not have 'chr', LD has 'chr' -> Remove 'chr' from LD file
    # Case 3: Both match -> Do nothing (copy LD file)
    
    cp ${high_ld} high_ld_regions.txt
    
    if [[ "\$BIM_CHR" == chr* ]] && [[ "\$LD_CHR" != chr* ]]; then
        echo "BIM has chr prefix, LD file does not. Adding chr prefix to LD file."
        awk '{print "chr"\$0}' ${high_ld} > high_ld_regions.txt
    elif [[ "\$BIM_CHR" != chr* ]] && [[ "\$LD_CHR" == chr* ]]; then
        echo "BIM does not have chr prefix, LD file has it. Removing chr prefix from LD file."
        sed 's/^chr//' ${high_ld} > high_ld_regions.txt
    else
        echo "Chromosome formats match or no adjustment needed."
    fi

    # LD Pruning using PLINK2
    
    /home/b/b37974/plink2_alpha6/plink2 \
        --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --exclude range high_ld_regions.txt \
        --maf 0.05 \
        --snps-only just-acgt \
        --indep-pairwise 50 5 0.2 \
        --out ${out_prefix} \
        --threads 8
    """
}

/*
 * Process 1: Calculate Sparse GRM
 * Inputs: PLINK file set (bed, bim, fam), Pruned marker list
 * Outputs: Sparse GRM matrix and sample IDs
 */
process CALC_SPARSE_GRM {
    executor 'slurm'
    queue 'gr10478b'
    time '48h'
    publishDir "${params.OutPath}/saige_sparse_grm/00.sparse_grm", mode: 'symlink'

    input:
    tuple val(prefix), path(bed), path(bim), path(fam)
    path(prune_in)

    output:
    tuple path("*.sparseGRM.mtx"), path("*.sparseGRM.mtx.sampleIDs.txt"), emit: sparse_grm

    script:
    def out_prefix = prefix
    """
    # Create pruned dataset
    /home/b/b37974/plink2_alpha6/plink2 \
        --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --extract ${prune_in} \
        --make-bed \
        --out ${out_prefix}.pruned \
        --threads 8

    source activate saige
    
    # Calculate sparse GRM using R script
    createSparseGRM.R \
        --plinkFile=${out_prefix}.pruned \
        --nThreads=8 \
        --outputPrefix=${out_prefix} \
        --numRandomMarkerforSparseKin=2000 \
        --relatednessCutoff=0.125 \
        --isDiagofKinSetAsOne=TRUE
    """
}

/*
 * Process 2: Fit Null GLMM
 * Inputs: Sparse GRM from previous step
 * Outputs: Null model file and variance ratio file
 */
process FIT_NULL_GLMM_SPARSE_GRM {
    executor 'slurm'
    queue 'gr10478b'
    time '48h'
    publishDir "${params.OutPath}/saige_sparse_grm/01.null_model", mode: 'symlink'

    input:
    tuple path(sparse_grm), path(sparse_grm_id)
    tuple val(plink_prefix), path(plink_bed), path(plink_bim), path(plink_fam)
    path(pheno_file)
    path(cov_file)

    output:
    tuple path("*.rda"), path("*.varianceRatio.txt"), emit: null_model

    script:
    def out_prefix = "cteph_agp3k.saige.null"
    def merge_script = "${params.ScriptsPath}/merge_pheno_cov.py"
    def plink_prefix_real = plink_bed.baseName
    """
    # Merge phenotype and covariates using python script
    python3 ${merge_script} \\
        --pheno "${pheno_file}" \\
        --cov "${cov_file}" \\
        --pheno_col "PHENO1" \\
        --cov_list "${params.Covariates}" \\
        --sex_col "SEX" \\
        --out "merged_pheno_cov.txt"

    # Run SAIGE Step 1: Fit Null GLMM
    source activate saige
    
    step1_fitNULLGLMM.R \
        --plinkFile=${plink_prefix_real} \
        --phenoFile=merged_pheno_cov.txt \
        --phenoCol=PHENO1 \
        --covarColList=${params.Covariates} \
        --sexCol=SEX \
        --sampleIDColinphenoFile=IID \
        --traitType=binary \
        --outputPrefix=${out_prefix} \
        --nThreads=32 \
        --numRandomMarkerforVarianceRatio=200 \
        --skipVarianceRatioEstimation=FALSE \
        --useSparseGRMtoFitNULL=TRUE  \
        --sparseGRMFile=${sparse_grm} \
        --sparseGRMSampleIDFile=${sparse_grm_id} \
        --IsOverwriteVarianceRatioFile=TRUE
    """
}

/*
 * Process 3: Association Test
 * Inputs: Null model, Sparse GRM, PLINK files
 * Outputs: Association results file
 */
process ASSOC_TEST_SPARSE_GRM {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    publishDir "${params.OutPath}/saige_sparse_grm/02.saige_assoc", mode: 'symlink'

    input:
    tuple path(model_file), path(variance_ratio)
    tuple path(sparse_grm_file), path(sparse_grm_id)
    tuple val(plink_prefix), path(plink_bed), path(plink_bim), path(plink_fam)

    output:
    path("*.txt"), emit: assoc_results

    script:
    def out_file = "cteph_agp3k.saige.assoc.txt"
    
    """
    source activate saige
    
    # Run SAIGE Step 2: Association Test
    step2_SPAtests.R \
        --bedFile=${plink_bed} \
        --bimFile=${plink_bim} \
        --famFile=${plink_fam} \
        --AlleleOrder=alt-first \
        --SAIGEOutputFile=${out_file} \
        --GMMATmodelFile=${model_file} \
        --varianceRatioFile=${variance_ratio} \
        --LOCO=FALSE \
        --is_output_moreDetails=TRUE \
        --sparseGRMFile=${sparse_grm_file} \
        --sparseGRMSampleIDFile=${sparse_grm_id} \
        --is_fastTest=TRUE
    """
}

/*
 * Process 4: Plot Results
 * Inputs: Association results
 * Outputs: Manhattan and QQ plots
 */
process PLOT_MANHATTAN_QQ_SPARSE_GRM {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    publishDir "${params.OutPath}/saige_sparse_grm/03.plots", mode: 'symlink'

    input:
    path(assoc_file)

    output:
    path("*.png"), emit: plots

    script:
    def script_path = "${params.ScriptsPath}/saige_manhattan_qq.py"
    def out_prefix = assoc_file.baseName
    """
    source activate cteph_geno_pro

    python3 ${script_path} \\
        --input ${assoc_file} \\
        --output-prefix ${out_prefix} \\
        --title "SAIGE Results (Sparse GRM)"
    """
}

/*
 * Process 5: Split PLINK by Chromosome and Convert to BGEN
 * Inputs: PLINK files
 * Outputs: BGEN files separated by chromosome with indexes and .sample files
 */
process SPLIT_PLINK_TO_BGEN {
    executor 'slurm'
    queue 'gr10478b'
    time '6h'
    publishDir "${params.OutPath}/saige_full_grm/00.prepared_bgen", mode: 'symlink'

    input:
    tuple val(prefix), path(bed), path(bim), path(fam)
    each chr

    output:
    tuple val(chr), path("*.bgen"), path("*.bgen.bgi"), path("*.sample"), emit: bgen_files

    script:
    def out_file = "${prefix}.chr${chr}"
    """
    # Set PATH for plink2
    export PATH=/home/b/b37974:\$PATH

    # Split by chromosome and convert to BGEN
    # Using --chr ${chr} to extract specific chromosome
    # Using --export bgen-1.2 bits=8 to get BGEN format
    # Using ref-first to ensure alignment with reference genome (good for SAIGE)
    
    plink2 \
        --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --chr ${chr} \
        --export bgen-1.2 bits=8 ref-first id-paste=iid\
        --out ${out_file} \
        --threads 4
        
    # Index the BGEN file using bgenix
    # BGEN file will be ${out_file}.bgen
    # Output index will be ${out_file}.bgen.bgi
    
    # Indexing is necessary for SAIGE to efficiently read the BGEN files during association testing
    source activate saige
    bgenix -g ${out_file}.bgen -index -clobber
    """
}

/*
 * Process 6: Fit Null GLMM (Full GRM)
 * Inputs: PLINK files, Phenotype, Covariates
 * Outputs: Null model file and variance ratio file
 */
process FIT_NULL_GLMM_FULL_GRM {
    executor 'slurm'
    queue 'gr10478b'
    time '72h'
    publishDir "${params.OutPath}/saige_full_grm/01.null_model", mode: 'symlink'

    input:
    tuple val(prefix), path(bed), path(bim), path(fam)
    path(prune_in)
    path(pheno_file)
    path(cov_file)

    output:
    tuple path("*.rda"), path("*.varianceRatio.txt"), emit: null_model

    script:
    def out_prefix = "cteph_agp3k.saige.full_grm.null"
    def merge_script = "${params.ScriptsPath}/merge_pheno_cov.py"
    """
    # Create pruned dataset
    /home/b/b37974/plink2_alpha6/plink2 \
        --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --extract ${prune_in} \
        --make-bed \
        --out ${prefix}.pruned \
        --threads 8

    # Merge phenotype and covariates using python script
    python3 ${merge_script} \\
        --pheno "${pheno_file}" \\
        --cov "${cov_file}" \\
        --pheno_col "PHENO1" \\
        --cov_list "${params.Covariates}" \\
        --sex_col "SEX" \\
        --out "merged_pheno_cov.txt"

    # Run SAIGE Step 1: Fit Null GLMM (Full GRM)
    source activate saige
    
    step1_fitNULLGLMM.R \
        --plinkFile=${prefix}.pruned \
        --phenoFile=merged_pheno_cov.txt \
        --phenoCol=PHENO1 \
        --covarColList=${params.Covariates} \
        --qCovarColList=SEX \
        --sexCol=SEX \
        --sampleIDColinphenoFile=IID \
        --traitType=binary \
        --outputPrefix=${out_prefix} \
        --nThreads=32 \
        --isDiagofKinSetAsOne=TRUE \
        --useSparseGRMtoFitNULL=FALSE \
        --IsOverwriteVarianceRatioFile=TRUE
    """
}

/*
 * Process 7: Association Test (Full GRM)
 * Inputs: Null model (Full GRM), BGEN files
 * Outputs: Association results file per chromosome
 */
process ASSOC_TEST_FULL_GRM {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    publishDir "${params.OutPath}/saige_full_grm/02.saige_assoc", mode: 'symlink'

    input:
    tuple path(model_file), path(variance_ratio), val(chr), path(bgen_file), path(bgen_index), path(sample_file)

    output:
    tuple val(chr), path("*.assoc.txt"), emit: assoc_results

    script:
    def out_file = "cteph_agp3k.saige.full_grm.chr${chr}.assoc.txt"
    
    """
    source activate saige
    
    # SAIGE requires a specific sample file format for BGEN
    # It often fails with standard Oxford format (5 cols). 
    # We create a simple headerless ID file from the .sample file provided by PLINK.
    # The PLINK sample file starts with 2 header lines.
    # Column 1 is ID_1.
    
    awk 'NR>2 {print \$1}' ${sample_file} > sample_ids.txt
    
    # Run SAIGE Step 2: Association Test using Full GRM
    # Note: For Full GRM, LOCO=TRUE is typically recommended and enabled by default in step 1, 
    # but here we must specify it again or ensure consistency.
    # Also using BGEN input instead of PLINK.
    
    step2_SPAtests.R \
        --bgenFile=${bgen_file} \
        --bgenFileIndex=${bgen_index} \
        --sampleFile=sample_ids.txt \
        --AlleleOrder=ref-first \
        --SAIGEOutputFile=${out_file} \
        --chrom=${chr} \
        --GMMATmodelFile=${model_file} \
        --varianceRatioFile=${variance_ratio} \
        --is_Firth_beta=FALSE \
        --LOCO=TRUE \
        --is_output_moreDetails=TRUE
    """
}

/*
 * Process 8: Merge Association Results
 * Inputs: Association results from all chromosomes
 * Outputs: Single merged association file
 */
process MERGE_ASSOC_RESULTS_FULL_GRM {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    publishDir "${params.OutPath}/saige_full_grm/03.saige_assoc_merged", mode: 'symlink'

    input:
    path(assoc_files)

    output:
    path("cteph_agp3k.saige.full_grm.assoc.txt"), emit: merged_results

    script:
    """
    # Create the header from the first file
    # We find the file for chr1 specifically to use its header, or just pick the first available one
    # Note: assoc_files is a list of files. 
    
    # Extract header from one file
    head -n 1 \$(ls *.assoc.txt | head -n 1) > cteph_agp3k.saige.full_grm.assoc.txt
    
    # Append content from all files, skipping headers, sorted by chromosome number if possible
    # A simple cat would work but order might be random.
    # We loop through 1 to 22.
    
    for i in {1..22}; do
        # Find the file corresponding to this chromosome
        # The filename pattern is cteph_agp3k.saige.full_grm.chr\${i}.assoc.txt
        file="cteph_agp3k.saige.full_grm.chr\${i}.assoc.txt"
        
        if [ -f "\$file" ]; then
            tail -n +2 "\$file" >> cteph_agp3k.saige.full_grm.assoc.txt
        fi
    done
    """
}

/*
 * Process 9: Plot Full GRM Results
 * Inputs: Merged association results
 * Outputs: Manhattan and QQ plots
 */
process PLOT_MANHATTAN_QQ_FULL_GRM {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    publishDir "${params.OutPath}/saige_full_grm/04.plots", mode: 'symlink'

    input:
    path(assoc_file)

    output:
    path("*.png"), emit: plots

    script:
    def script_path = "${params.ScriptsPath}/saige_manhattan_qq.py"
    def out_prefix = assoc_file.baseName
    """
    source activate cteph_geno_pro

    python3 ${script_path} \\
        --input ${assoc_file} \\
        --output-prefix ${out_prefix} \\
        --title "SAIGE Results (Full GRM)"
    """
}

/*
 * Main Workflow
 */
workflow {
    // Define input channels
    plink_file_prefix = "${params.WGSPath}/15.run_variant_qc/cteph_agp3k.sqc.vqc"
    
    // Create a channel for the PLINK fileset (bed, bim, fam)
    // Using fromFilePairs to group them together
    plink_ch = channel.fromFilePairs("${plink_file_prefix}.{bed,bim,fam}", size: 3)
        .map { id, files -> tuple(id, files[0], files[1], files[2]) } 
    
    // Channels for Phenotype and Covariate files
    pheno_ch = channel.fromPath("${params.WGSPath}/19.cov_pheno_prepare/cteph_agp3k.bbj.projection.pheno_df.csv")
    cov_ch = channel.fromPath("${params.WGSPath}/19.cov_pheno_prepare/cteph_agp3k.bbj.projection.cov_df.no_age.csv")

    // Chromosome channel (1-22)
    chr_ch = channel.of(1..22)
    
    // High LD regions file
    high_ld_ch = channel.fromPath("${params.InfoPath}/high-LD-regions-hg38-GRCh38.txt")

    // Optional MAC Filtering
    if (params.FilterMAC) {
        FILTER_BY_MAC(plink_ch)
        target_plink_ch = FILTER_BY_MAC.out.filtered_plink
    } else {
        target_plink_ch = plink_ch
    }

    // Perform LD Pruning on the target PLINK dataset (original or MAC-filtered)
    LD_PRUNING(target_plink_ch, high_ld_ch)
    prune_in_ch = LD_PRUNING.out.prune_in

    // Run processes based on GRM type
    if (params.GRM == 'sparse') {
        CALC_SPARSE_GRM(target_plink_ch, prune_in_ch)
        
        FIT_NULL_GLMM_SPARSE_GRM(
            CALC_SPARSE_GRM.out.sparse_grm,
            target_plink_ch,
            pheno_ch,
            cov_ch
        )
        
        ASSOC_TEST_SPARSE_GRM(
            FIT_NULL_GLMM_SPARSE_GRM.out.null_model,
            CALC_SPARSE_GRM.out.sparse_grm,
            target_plink_ch
        )
        
        PLOT_MANHATTAN_QQ_SPARSE_GRM(ASSOC_TEST_SPARSE_GRM.out.assoc_results)
    } else if (params.GRM == 'full') {
        SPLIT_PLINK_TO_BGEN(target_plink_ch, chr_ch)
        
        // Fit Null GLMM using Full GRM
        FIT_NULL_GLMM_FULL_GRM(target_plink_ch, prune_in_ch, pheno_ch, cov_ch)
        
        ASSOC_TEST_FULL_GRM(
            FIT_NULL_GLMM_FULL_GRM.out.null_model.combine(SPLIT_PLINK_TO_BGEN.out.bgen_files)
        )
        
        // Collect all association results and merge them
        // .collect() gathers all emitted files into a single list
        MERGE_ASSOC_RESULTS_FULL_GRM(
            ASSOC_TEST_FULL_GRM.out.assoc_results.map{ x -> x[1] }.collect()
        )
        
        // Plot Manhattan/QQ for Full GRM
        PLOT_MANHATTAN_QQ_FULL_GRM(MERGE_ASSOC_RESULTS_FULL_GRM.out.merged_results)
    } else {
        error "Invalid GRM parameter: ${params.GRM}. Options are 'sparse' or 'full'."
    }
}

