#!/bin/bash
# Wrapper Script for Gene Detail Check (Smart & Robust Version)
# Usage: ./check_gene.sh <Gene_Name> [Group_Name] [Mode] [Sample_Group]

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m' # No Color

GENE_NAME=$1
GROUP_NAME=${2} # User provided group, or auto-set below
MODE=${3:-"sensitivity"}
SAMPLE_GROUP=${4:-"both"}

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
            echo -e "${BLUE}[INFO] Sensitivity Mode: Auto-appended .stat1_stat2 to group name.${NC}"
        fi
    fi
fi

if [ -z "$GENE_NAME" ]; then
    echo -e "${RED}Usage: $0 <Gene_Name> [Group_Name] [Mode] [Sample_Group]${NC}"
    echo "  Mode: 'sensitivity' (default) or 'original'"
    echo "  Sample_Group: 'both' (default), 'case', or 'control'"
    exit 1
fi

# ----------------- Path Auto-Detection -----------------
SCRIPT_DIR=$(dirname "$0")/scripts
WORK_DIR=$(dirname "$0")
ANCHOR_PATH="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest.rev1/results"

# Set Base Directories based on Mode
if [[ "$MODE" == "original" || "$MODE" == "main" ]]; then
    # Original Results
    VCF_BASE="${ANCHOR_PATH}/03.info_filter"
    ASSOC_BASE="${ANCHOR_PATH}/04.rvtest_run"
    POST_PROCESS_BASE="${ANCHOR_PATH}/05.post_process"
    echo -e "Mode: ${BOLD}ORIGINAL RESULTS${NC}"
else
    # Sensitivity Results (Default)
    SENSITIVITY_DIR="${ANCHOR_PATH}/07.sensitivity_check"
    VCF_BASE="${SENSITIVITY_DIR}/00.data_prepare"
    ASSOC_BASE="${SENSITIVITY_DIR}/01.rvtest_run"
    POST_PROCESS_BASE="${SENSITIVITY_DIR}/02.post_process"
    echo -e "Mode: ${BOLD}SENSITIVITY CHECK${NC}"
fi

# 1. Locate Raw Assoc File
ASSOC_FILE="${ASSOC_BASE}/${GROUP_NAME}/skato/cteph_agp3k.rare.${GROUP_NAME}.skato.SkatO.assoc"

# 2. Locate Burden FDR File
BURDEN_FILE="${POST_PROCESS_BASE}/${GROUP_NAME}/burden/cteph_agp3k.rare.${GROUP_NAME}.burden.CMC.filtered.fdr.assoc"

# 3. Locate SKAT-O FDR File
SKATO_FILE="${POST_PROCESS_BASE}/${GROUP_NAME}/skato/cteph_agp3k.rare.${GROUP_NAME}.skato.SkatO.filtered.fdr.assoc"

# 4. Locate VCF File
VCF_FILE="${VCF_BASE}/cteph_agp3k.rare.${GROUP_NAME}.vcf.gz"

# ----------------- Validation -----------------
echo -e "${CYAN}==========================================${NC}"
echo -e "${BOLD}Target Gene     : ${GREEN}${GENE_NAME}${NC}"
echo -e "${BOLD}Analysis Group  : ${BLUE}${GROUP_NAME}${NC}"
echo -e "${CYAN}==========================================${NC}"

MISSING=0
if [ ! -f "$ASSOC_FILE" ]; then 
    echo -e "${YELLOW}[WARN] Raw Assoc file not found: $ASSOC_FILE${NC}"; 
    ASSOC_DIR=$(dirname "$ASSOC_FILE")
    FOUND=$(find "$ASSOC_DIR" -maxdepth 1 -name "*.assoc" 2>/dev/null | head -n 1)
    if [ -n "$FOUND" ]; then
        echo -e "${BLUE}[INFO] Using alternative: $FOUND${NC}"
        ASSOC_FILE="$FOUND"
    else
        MISSING=1
    fi
fi

if [ ! -f "$VCF_FILE" ]; then 
    echo -e "${YELLOW}[WARN] VCF file not found: $VCF_FILE${NC}"
    echo "  Searching for partial match in ${VCF_BASE}..."
    FOUND=$(find "$VCF_BASE" -maxdepth 1 -name "*${GROUP_NAME}*.vcf.gz" 2>/dev/null | head -n 1)
    if [ -n "$FOUND" ]; then
         echo -e "${BLUE}[INFO] Using alternative VCF: $FOUND${NC}"
         VCF_FILE="$FOUND"
    else
         MISSING=1
    fi
fi

if [ $MISSING -eq 1 ]; then
    echo -e "${RED}------------------------------------------${NC}"
    echo -e "${RED}[ERROR] Critical files missing. Check Group Name and Mode.${NC}"
    exit 1
fi

# ----------------- Execution -----------------
PLINK_PREFIX="${ANCHOR_PATH}/00.pre_step/cteph_agp3k.rare.mac2.rm_samples"
PHENO_FILE="${ANCHOR_PATH}/01.rvtest_prepare/cteph_agp3k.bbj.projection.pheno_df.csv"
TOMMO_VCF="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz"
REFFLAT_FILE="${ANCHOR_PATH}/01.rvtest_prepare/refFlat.hg38.nochr.txt.gz"
PLINK2_PATH="/home/b/b37974/plink2_alpha6/plink2"

OUT_BASE="${WORK_DIR}/output/${GENE_NAME}"
mkdir -p "${OUT_BASE}"

TMP_DIR="${OUT_BASE}/tmp/${GROUP_NAME}_$(date +%s)"
mkdir -p ${TMP_DIR}

LOG_FILE="${OUT_BASE}/${GENE_NAME}.${GROUP_NAME}.summary.txt"
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
    --out-log "${LOG_FILE}" \
    --sample-group "${SAMPLE_GROUP}"

# Clean up
rm -rf "${TMP_DIR}"

echo -e "${GREEN}Done. Full Log saved to: ${LOG_FILE}${NC}"
