#!/bin/bash
# Wrapper Script for Gene Detail Check (Structured & Robust)
# Usage: ./check_gene.sh <Gene> [Mode] [Impact] [Sample_Group]

# --- Colors ---
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# --- Arguments ---
GENE_NAME=$1
MODE=${2:-"sens1"}          # main, sens1, sens2
IMPACT=${3:-"moderate_high"} # high, moderate_high, low_moderate_high
SAMPLE_GROUP=${4:-"both"}    # both, case, control

if [ -z "$GENE_NAME" ]; then
    echo -e "${RED}Usage: $0 <Gene_Name> [Mode] [Impact] [Sample_Group]${NC}"
    echo "  Mode:"
    echo "    main  : Original Results (05.post_process)"
    echo "    sens1 : Sensitivity Type 1 (Pass ToMMo QC) -> .stat1_stat2 suffix"
    echo "    sens2 : Sensitivity Type 2 (Type 1 + AF Consistency) -> .stat1 suffix"
    echo "  Impact:"
    echo "    high"
    echo "    moderate_high (default)"
    echo "    low_moderate_high"
    exit 1
fi

# --- Base Directories ---
ROOT_DIR="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest.rev1/results"
MAIN_POST_DIR="${ROOT_DIR}/05.post_process"
SENS_POST_DIR="${ROOT_DIR}/07.sensitivity_check/02.post_process"

# --- Logic Parsing ---

# 1. Determine Suffix and Base Paths based on Mode
if [[ "$MODE" == "main" ]]; then
    # Original Analysis
    MODE_DESC="Original Results (Main)"
    echo -e "Mode: ${BOLD}${MODE_DESC}${NC}"
    
    # Path Config
    VCF_BASE="${ROOT_DIR}/03.info_filter"
    ASSOC_BASE="${ROOT_DIR}/04.rvtest_run"
    POST_PROCESS_BASE="${MAIN_POST_DIR}"
    
    # Group naming for Main usually strictly follows impact
    GROUP_SUFFIX="" 

elif [[ "$MODE" == "sens1" ]]; then
    # Sensitivity 1: "stat1_stat2"
    MODE_DESC="Sensitivity 1 (stat1_stat2)"
    echo -e "Mode: ${BOLD}${MODE_DESC}${NC}"
    
    VCF_BASE="${ROOT_DIR}/07.sensitivity_check/00.data_prepare"
    ASSOC_BASE="${ROOT_DIR}/07.sensitivity_check/01.rvtest_run"
    POST_PROCESS_BASE="${SENS_POST_DIR}"
    
    GROUP_SUFFIX=".stat1_stat2"

elif [[ "$MODE" == "sens2" ]]; then
    # Sensitivity 2: "stat1" (Based on directory structure on disk)
    MODE_DESC="Sensitivity 2 (stat1)"
    echo -e "Mode: ${BOLD}${MODE_DESC}${NC}"
    
    VCF_BASE="${ROOT_DIR}/07.sensitivity_check/00.data_prepare"
    ASSOC_BASE="${ROOT_DIR}/07.sensitivity_check/01.rvtest_run"
    POST_PROCESS_BASE="${SENS_POST_DIR}"
    
    GROUP_SUFFIX=".stat1"

else
    echo -e "${RED}[ERROR] Invalid Mode: $MODE${NC}"
    echo "Valid modes: main, sens1, sens2"
    exit 1
fi

# 2. Determine Full Group Name
# Normalize Impact Input
case "$IMPACT" in
    "high")
        IMPACT_STR="impact_high"
        ;;
    "moderate_high"|"mod_high")
        IMPACT_STR="impact_moderate_high"
        ;;
    "low_moderate_high"|"all")
        IMPACT_STR="impact_low_moderate_high"
        ;;
    "impact_"*)
        IMPACT_STR="$IMPACT" # User entered full string
        ;;
    *)
        echo -e "${YELLOW}[WARN] Unknown impact format '$IMPACT', assuming 'impact_$IMPACT'${NC}"
        IMPACT_STR="impact_$IMPACT"
        ;;
esac

GROUP_NAME="${IMPACT_STR}${GROUP_SUFFIX}"

# Define Output Tag for naming files (Include Mode explicitly)
OUTPUT_TAG="${IMPACT_STR}.${MODE}"

# ----------------- File Locations -----------------

# 1. Raw Assoc File (skato output)
ASSOC_FILE="${ASSOC_BASE}/${GROUP_NAME}/skato/cteph_agp3k.rare.${GROUP_NAME}.skato.SkatO.assoc"

# 2. Burden FDR File
BURDEN_FILE="${POST_PROCESS_BASE}/${GROUP_NAME}/burden/cteph_agp3k.rare.${GROUP_NAME}.burden.CMC.filtered.fdr.assoc"

# 3. SKAT-O FDR File
SKATO_FILE="${POST_PROCESS_BASE}/${GROUP_NAME}/skato/cteph_agp3k.rare.${GROUP_NAME}.skato.SkatO.filtered.fdr.assoc"

# 4. VCF File
VCF_FILE="${VCF_BASE}/cteph_agp3k.rare.${GROUP_NAME}.vcf.gz"

# ----------------- Validation & Execution -----------------

echo -e "${CYAN}==========================================${NC}"
echo -e "${BOLD}Target Gene     : ${GREEN}${GENE_NAME}${NC}"
echo -e "${BOLD}Analysis Group  : ${BLUE}${GROUP_NAME}${NC}"
echo -e "${BOLD}Sample Analysis : ${BLUE}${SAMPLE_GROUP}${NC}"
echo -e "${CYAN}==========================================${NC}"

# Check Critical Files
MISSING=0
if [ ! -f "$ASSOC_FILE" ]; then 
    echo -e "${YELLOW}[WARN] Raw Assoc file not found: $ASSOC_FILE${NC}"
    # Approximate search
    PARENT=$(dirname "$ASSOC_FILE")
    FOUND=$(find "$PARENT" -maxdepth 1 -name "*.assoc" 2>/dev/null | head -n 1)
    if [ -n "$FOUND" ]; then
        echo -e "${BLUE}[INFO] Using alternative: $FOUND${NC}"
        ASSOC_FILE="$FOUND"
    else
        MISSING=1
    fi
fi

if [ ! -f "$VCF_FILE" ]; then 
    echo -e "${YELLOW}[WARN] VCF file not found: $VCF_FILE${NC}"
     # Approximate search in VCF_BASE
    FOUND=$(find "$VCF_BASE" -maxdepth 1 -name "*${GROUP_NAME}*.vcf.gz" 2>/dev/null | head -n 1)
    if [ -n "$FOUND" ]; then
         echo -e "${BLUE}[INFO] Using alternative VCF: $FOUND${NC}"
         VCF_FILE="$FOUND"
    else
         MISSING=1
    fi
fi

if [ $MISSING -eq 1 ]; then
    echo -e "${RED}[ERROR] Critical files missing for group '$GROUP_NAME'.${NC}"
    echo "  - Check if the Mode/Impact combination exists."
    exit 1
fi

# Paths
SCRIPT_DIR=$(dirname "$0")/scripts
WORK_DIR=$(dirname "$0")

PLINK_PREFIX="${ROOT_DIR}/00.pre_step/cteph_agp3k.rare.mac2.rm_samples"
PHENO_FILE="${ROOT_DIR}/01.rvtest_prepare/cteph_agp3k.bbj.projection.pheno_df.csv"
TOMMO_VCF="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz"
REFFLAT_FILE="${ROOT_DIR}/01.rvtest_prepare/refFlat.hg38.nochr.txt.gz"
PLINK2_PATH="/home/b/b37974/plink2_alpha6/plink2"

# Output
OUT_BASE="${WORK_DIR}/output/${GENE_NAME}"
mkdir -p "${OUT_BASE}"
TMP_DIR="${OUT_BASE}/tmp/${GROUP_NAME}_$(date +%s)"
mkdir -p ${TMP_DIR}
LOG_FILE="${OUT_BASE}/${GENE_NAME}.${OUTPUT_TAG}.summary.txt"
if [ -f "${LOG_FILE}" ]; then rm "${LOG_FILE}"; fi

source activate cteph_geno_pro

# Execute Python
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
    --sample-group "${SAMPLE_GROUP}" \
    --group-name "${OUTPUT_TAG}"

# Clean
rm -rf "${TMP_DIR}"

# --- Final Report ---
echo -e ""
echo -e "${CYAN}============================================================${NC}"
echo -e "${BOLD}                   ANALYSIS COMPLETE                        ${NC}"
echo -e "${CYAN}============================================================${NC}"
echo -e "${BOLD}[ANALYSIS CONTEXT]${NC}"
echo -e "  > Gene Target   : ${GREEN}${GENE_NAME}${NC}"
echo -e "  > Variant Group : ${BLUE}${GROUP_NAME}${NC}"
echo -e "  > Analysis Mode : ${YELLOW}${MODE_DESC}${NC}"
echo -e "  > Sample Scope  : ${YELLOW}${SAMPLE_GROUP}${NC}"
echo -e ""
echo -e "${BOLD}[OUTPUT RESULTS]${NC}"
echo -e "  > Directory     : ${OUT_BASE}"
echo -e "  > Summary Log   : ${GREEN}${LOG_FILE}${NC}"

# Check for Sample Details
CASE_DETAILS="${OUT_BASE}/${GENE_NAME}.${OUTPUT_TAG}.case.sample_details.tsv"
CTRL_DETAILS="${OUT_BASE}/${GENE_NAME}.${OUTPUT_TAG}.control.sample_details.tsv"

if [ -f "$CASE_DETAILS" ]; then
    echo -e "  > Case Details  : ${GREEN}${CASE_DETAILS}${NC}"
fi
if [ -f "$CTRL_DETAILS" ]; then
    echo -e "  > Ctrl Details  : ${GREEN}${CTRL_DETAILS}${NC}"
fi

echo -e "${CYAN}============================================================${NC}"
echo -e ""
