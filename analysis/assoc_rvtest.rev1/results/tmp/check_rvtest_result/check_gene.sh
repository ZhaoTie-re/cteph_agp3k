#!/bin/bash
# Wrapper Script for Gene Detail Check (Smart & Robust Version)
# Usage: ./check_gene.sh <Gene_Name> [Group_Name] [Mode]
# Mode options: 'sensitivity' (default) or 'original'

GENE_NAME=$1
GROUP_NAME=${2} # User provided group, or auto-set below
MODE=${3:-"sensitivity"}

# Default Group Logic depends on Mode
if [ -z "$GROUP_NAME" ]; then
    if [[ "$MODE" == "original" || "$MODE" == "main" ]]; then
        GROUP_NAME="impact_moderate_high"
    else
        GROUP_NAME="impact_moderate_high.stat1_stat2"
    fi
else
    # Automatic suffix handling for Sensitivity Mode
    if [[ "$MODE" != "original" && "$MODE" != "main" ]]; then
        if [[ "$GROUP_NAME" != *".stat1_stat2" && "$GROUP_NAME" != *"stat"* ]]; then
            GROUP_NAME="${GROUP_NAME}.stat1_stat2"
            echo "[INFO] Sensitivity Mode: Auto-appended .stat1_stat2 to group name."
        fi
    fi
fi

if [ -z "$GENE_NAME" ]; then
    echo "Usage: $0 <Gene_Name> [Group_Name] [Mode]"
    echo "  Mode: 'sensitivity' (default) or 'original'"
    echo "  Default Group (sensitivity): impact_moderate_high.stat1_stat2"
    echo "  Default Group (original):    impact_moderate_high"
    exit 1
fi

# ----------------- Path Auto-Detection -----------------
# Determine script location and project root
SCRIPT_DIR=$(dirname "$0")/scripts
WORK_DIR=$(dirname "$0")

# Anchor Path
ANCHOR_PATH="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest.rev1/results"

# Set Base Directories based on Mode
if [[ "$MODE" == "original" || "$MODE" == "main" ]]; then
    # Original Results
    # VCFs: 03.info_filter
    # Assocs: 04.rvtest_run (Inferred)
    # FDRs: 05.post_process
    
    VCF_BASE="${ANCHOR_PATH}/03.info_filter"
    ASSOC_BASE="${ANCHOR_PATH}/04.rvtest_run"
    POST_PROCESS_BASE="${ANCHOR_PATH}/05.post_process"
    
    echo "Mode: ORIGINAL RESULTS"
else
    # Sensitivity Results (Default)
    SENSITIVITY_DIR="${ANCHOR_PATH}/07.sensitivity_check"
    
    VCF_BASE="${SENSITIVITY_DIR}/00.data_prepare"
    ASSOC_BASE="${SENSITIVITY_DIR}/01.rvtest_run"
    POST_PROCESS_BASE="${SENSITIVITY_DIR}/02.post_process"
    
    echo "Mode: SENSITIVITY CHECK"
fi

# 1. Locate Raw Assoc File (Mainly for Fallback/Range)
# Pattern: .../{ASSOC_BASE}/{GROUP}/skato/cteph_agp3k.rare.{GROUP}.skato.SkatO.assoc
ASSOC_FILE="${ASSOC_BASE}/${GROUP_NAME}/skato/cteph_agp3k.rare.${GROUP_NAME}.skato.SkatO.assoc"

# 2. Locate Burden FDR File
# Pattern: .../{POST_PROCESS_BASE}/{GROUP}/burden/cteph_agp3k.rare.{GROUP}.burden.CMC.filtered.fdr.assoc
BURDEN_FILE="${POST_PROCESS_BASE}/${GROUP_NAME}/burden/cteph_agp3k.rare.${GROUP_NAME}.burden.CMC.filtered.fdr.assoc"

# 3. Locate SKAT-O FDR File
# Pattern: .../{POST_PROCESS_BASE}/{GROUP}/skato/cteph_agp3k.rare.${GROUP_NAME}.skato.SkatO.filtered.fdr.assoc
SKATO_FILE="${POST_PROCESS_BASE}/${GROUP_NAME}/skato/cteph_agp3k.rare.${GROUP_NAME}.skato.SkatO.filtered.fdr.assoc"

# 4. Locate VCF File
# Pattern: .../{VCF_BASE}/cteph_agp3k.rare.{GROUP}.vcf.gz
# Note: In 03.info_filter, naming might differ slightly if group name doesn't match perfectly.
# Standard assumption: cteph_agp3k.rare.{GROUP}.vcf.gz
# However, for 'original', user said 'original does not contain stat'.
# If Group is 'impact_moderate_high', VCF is usually 'cteph_agp3k.rare.impact_moderate_high.vcf.gz'.
# This fits the pattern.
VCF_FILE="${VCF_BASE}/cteph_agp3k.rare.${GROUP_NAME}.vcf.gz"

# ----------------- Validation -----------------
echo "=========================================="
echo "Smart Gene Check: ${GENE_NAME}"
echo "Analysis Group  : ${GROUP_NAME}"
echo "=========================================="

MISSING=0
if [ ! -f "$ASSOC_FILE" ]; then 
    echo "[WARN] Raw Assoc file not found: $ASSOC_FILE"; 
    # Try finding ANY assoc file in that folder if strict name fails?
    ASSOC_DIR=$(dirname "$ASSOC_FILE")
    FOUND=$(find "$ASSOC_DIR" -maxdepth 1 -name "*.assoc" | head -n 1)
    if [ -n "$FOUND" ]; then
        echo "[INFO] Using alternative: $FOUND"
        ASSOC_FILE="$FOUND"
    else
        MISSING=1
    fi
fi

if [ ! -f "$VCF_FILE" ]; then 
    echo "[WARN] VCF file not found: $VCF_FILE"
    # Try finding partial name match in VCF Dir
    echo "  Searching for partial match in ${VCF_BASE}..."
    FOUND=$(find "$VCF_BASE" -maxdepth 1 -name "*${GROUP_NAME}*.vcf.gz" | head -n 1)
    if [ -n "$FOUND" ]; then
         echo "[INFO] Using alternative VCF: $FOUND"
         VCF_FILE="$FOUND"
    else
         MISSING=1
    fi
fi

if [ ! -f "$BURDEN_FILE" ]; then echo "[WARN] Burden FDR file not found (Optional): $BURDEN_FILE"; fi
if [ ! -f "$SKATO_FILE" ]; then echo "[WARN] SkatO FDR file not found (Optional): $SKATO_FILE"; fi

if [ $MISSING -eq 1 ]; then
    echo "------------------------------------------"
    echo "[ERROR] Critical files missing. Check Group Name and Mode."
    exit 1
fi

# ----------------- Execution -----------------

# Static Paths (usually constant)
PLINK_PREFIX="${ANCHOR_PATH}/00.pre_step/cteph_agp3k.rare.mac2.rm_samples"
PHENO_FILE="${ANCHOR_PATH}/01.rvtest_prepare/cteph_agp3k.bbj.projection.pheno_df.csv"
TOMMO_VCF="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz"
REFFLAT_FILE="${ANCHOR_PATH}/01.rvtest_prepare/refFlat.hg38.nochr.txt.gz"
PLINK2_PATH="/home/b/b37974/plink2_alpha6/plink2"

# Output
TMP_DIR="${WORK_DIR}/tmp/${GENE_NAME}_${GROUP_NAME}_$(date +%s)"
mkdir -p ${TMP_DIR}
LOG_FILE="${WORK_DIR}/${GENE_NAME}.${GROUP_NAME}.detail.log"

if [ -f "${LOG_FILE}" ]; then rm "${LOG_FILE}"; fi

source activate cteph_geno_pro

python ${SCRIPT_DIR}/check_gene_detail.py \
    --gene ${GENE_NAME} \
    --assoc-file "${ASSOC_FILE}" \
    --burden-file "${BURDEN_FILE}" \
    --skato-file "${SKATO_FILE}" \
    --vcf-file "${VCF_FILE}" \
    --plink-prefix "${PLINK_PREFIX}" \
    --tommo-vcf "${TOMMO_VCF}" \
    --pheno-file "${PHENO_FILE}" \
    --refflat-file "${REFFLAT_FILE}" \
    --plink2-path "${PLINK2_PATH}" \
    --out-dir "${TMP_DIR}" \
    --out-log "${LOG_FILE}"

echo ""
echo "Done. Log saved to ${LOG_FILE}"
echo "------------------------------------------"
cat ${LOG_FILE}
echo "------------------------------------------"
