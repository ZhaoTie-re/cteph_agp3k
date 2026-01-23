#!/bin/bash
# Wrapper Script for Gene Detail Check
# Usage: ./check_gene.sh <Gene_Name>

GENE_NAME=$1

# ----------------- Configuration -----------------
# Determine script directory
SCRIPT_DIR=$(dirname "$0")/scripts
WORK_DIR=$(dirname "$0")

# Input File Paths (Modify these defaults if needed)
ASSOC_FILE="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest.rev1/results/07.sensitivity_check/01.rvtest_run/impact_moderate_high.stat1_stat2/skato/cteph_agp3k.rare.impact_moderate_high.stat1_stat2.skato.SkatO.assoc" 
# Assuming one assoc file. If parameterization needed, add argument.
# But request says "RVTest results file (parameterized)" so let's allow override via env var or arg?
# Ideally, hardcode the main one or detect? 
# "Provide a gene name, tool finds in rvtest result file (parameterized)" -> Input argument.
# Let's make ASSOC_FILE the second argument.

if [ -z "$2" ]; then
    echo "Usage: $0 <Gene_Name> <Assoc_File> [VCF_File] [Plink_Prefix] [Tommo_VCF]"
    echo "Using default paths for missing arguments."
else
    ASSOC_FILE=$2
fi

# Determine VCF File based on Assoc File Path if possible (Correct logic for pipeline)
# Path usually contains: .../impact_moderate_high.stat1_stat2/...
# VCFs are in .../03.info_filter/cteph_agp3k.rare.{IMPACT}.vcf.gz

# Default fallback
DEFAULT_VCF="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest.rev1/results/07.sensitivity_check/00.data_prepare/cteph_agp3k.rare.impact_low_moderate_high.stat1_stat2.vcf.gz"
DETECTED_VCF=""

if [[ "$ASSOC_FILE" == *"impact_high"* || "$ASSOC_FILE" == *"impact"* ]]; then
    # Dynamic VCF Detection Logic
    DATA_PREPARE_DIR="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest.rev1/results/07.sensitivity_check/00.data_prepare"
    
    # Get basename of assoc file
    ASSOC_BASENAME=$(basename "$ASSOC_FILE")
    
    # Iterate over VCFs in data prepare dir
    BEST_MATCH=""
    MAX_LEN=0
    
    for vcf in "${DATA_PREPARE_DIR}"/*.vcf.gz; do
        if [ ! -f "$vcf" ]; then continue; fi
        
        # Get base name of VCF (remove extension)
        # e.g. cteph_agp3k.rare.impact_moderate_high.stat1_stat2
        VCF_BASE=$(basename "$vcf" .vcf.gz)
        
        # Check if Assoc filename starts with VCF base name
        if [[ "$ASSOC_BASENAME" == "$VCF_BASE"* ]]; then
             LEN=${#VCF_BASE}
             if (( LEN > MAX_LEN )); then
                 MAX_LEN=$LEN
                 BEST_MATCH="$vcf"
             fi
        fi
    done
    
    if [ -n "$BEST_MATCH" ]; then
        DETECTED_VCF="$BEST_MATCH"
    else
        # Fallback to previous logic if no match found
        if [[ "$ASSOC_FILE" == *"impact_moderate_high"* ]]; then
            DETECTED_VCF="${DATA_PREPARE_DIR}/cteph_agp3k.rare.impact_moderate_high.stat1_stat2.vcf.gz"
        elif [[ "$ASSOC_FILE" == *"impact_low_moderate_high"* ]]; then
            DETECTED_VCF="${DATA_PREPARE_DIR}/cteph_agp3k.rare.impact_low_moderate_high.stat1_stat2.vcf.gz"
        elif [[ "$ASSOC_FILE" == *"impact_high"* ]]; then
            DETECTED_VCF="${DATA_PREPARE_DIR}/cteph_agp3k.rare.impact_high.stat1_stat2.vcf.gz"
        fi
    fi
fi

# Use detected if exists, else argument or default
if [ -z "$3" ] && [ -f "$DETECTED_VCF" ]; then
   VCF_FILE=$DETECTED_VCF
   echo "Auto-detected VCF: $VCF_FILE"
else
   VCF_FILE=${3:-"$DEFAULT_VCF"}
fi

# Determine Post Process Files (Burden & SKATO with FDR)
BURDEN_FILE=""
SKATO_FILE=""

if [[ "$ASSOC_FILE" == *"/01.rvtest_run/"* ]]; then
    # Extract Group Name (folder name after 01.rvtest_run)
    GROUP_NAME=$(echo "$ASSOC_FILE" | awk -F'/01.rvtest_run/' '{print $2}' | cut -d'/' -f1)
    
    # Extract Results Root
    RESULTS_ROOT=$(echo "$ASSOC_FILE" | awk -F'/01.rvtest_run/' '{print $1}')
    POST_PROCESS_DIR="${RESULTS_ROOT}/02.post_process/${GROUP_NAME}"
    
    if [ -d "$POST_PROCESS_DIR" ]; then
        # Format: cteph_agp3k.rare.{GROUP}.burden.CMC.filtered.fdr.assoc
        #         cteph_agp3k.rare.{GROUP}.skato.SkatO.filtered.fdr.assoc
        
        BURDEN_CANDIDATE="${POST_PROCESS_DIR}/burden/cteph_agp3k.rare.${GROUP_NAME}.burden.CMC.filtered.fdr.assoc"
        SKATO_CANDIDATE="${POST_PROCESS_DIR}/skato/cteph_agp3k.rare.${GROUP_NAME}.skato.SkatO.filtered.fdr.assoc"
        
        if [ -f "$BURDEN_CANDIDATE" ]; then BURDEN_FILE="$BURDEN_CANDIDATE"; fi
        if [ -f "$SKATO_CANDIDATE" ]; then SKATO_FILE="$SKATO_CANDIDATE"; fi
    fi
fi

if [ -n "$BURDEN_FILE" ]; then echo "Auto-detected Burden File: $BURDEN_FILE"; fi
if [ -n "$SKATO_FILE" ]; then echo "Auto-detected SKAT-O File: $SKATO_FILE"; fi

PLINK_PREFIX=${4:-"/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest.rev1/results/00.pre_step/cteph_agp3k.rare.mac2.rm_samples"}
TOMMO_VCF=${5:-"/LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz"}
PHENO_FILE=${6:-"/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest.rev1/results/01.rvtest_prepare/cteph_agp3k.bbj.projection.pheno_df.csv"}
# REFFLAT (Optional, hardcode path for now as per user request context or add arg)
REFFLAT_FILE="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest.rev1/results/01.rvtest_prepare/refFlat.hg38.nochr.txt.gz"

# Plink2 Path
PLINK2_PATH="/home/b/b37974/plink2_alpha6/plink2"

# Output
TMP_DIR="${WORK_DIR}/tmp/${GENE_NAME}_$(date +%s)"
mkdir -p ${TMP_DIR}
LOG_FILE="${WORK_DIR}/${GENE_NAME}.detail.log"

# Clean up previous log to ensure we don't see stale results if python fails
if [ -f "${LOG_FILE}" ]; then
    rm "${LOG_FILE}"
fi

echo "=========================================="
echo "Checking Gene: ${GENE_NAME}"
echo "Assoc File: ${ASSOC_FILE}"
echo "VCF File: ${VCF_FILE}"
echo "Plink Prefix: ${PLINK_PREFIX}"
echo "Tommo VCF: ${TOMMO_VCF}"
echo "Pheno File: ${PHENO_FILE}"
echo "RefFlat File: ${REFFLAT_FILE}"
echo "Temp Dir: ${TMP_DIR}"
echo "=========================================="

source activate cteph_geno_pro

python ${SCRIPT_DIR}/check_gene_detail.py \
    --gene ${GENE_NAME} \
    --assoc-file ${ASSOC_FILE} \
    --burden-file "${BURDEN_FILE}" \
    --skato-file "${SKATO_FILE}" \
    --vcf-file ${VCF_FILE} \
    --plink-prefix ${PLINK_PREFIX} \
    --tommo-vcf ${TOMMO_VCF} \
    --pheno-file ${PHENO_FILE} \
    --refflat-file ${REFFLAT_FILE} \
    --plink2-path ${PLINK2_PATH} \
    --out-dir ${TMP_DIR} \
    --out-log ${LOG_FILE}

# Clean up temp
# rm -rf ${TMP_DIR}
echo ""
echo "Done. Log saved to ${LOG_FILE}"
echo "------------------------------------------"
cat ${LOG_FILE}
echo "------------------------------------------"
