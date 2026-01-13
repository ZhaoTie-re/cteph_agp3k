#!/bin/bash

# ==============================================================================
# Script Name: run_calc_metrics.sh
# Description: Wrapper script to run PLINK2 and Python metric calculation.
#              This script encapsulates the logic for generating sample and
#              variant metrics from PLINK binary files.
#
# Usage: ./run_calc_metrics.sh [options]
#
# Options:
#   --bed-prefix    Path prefix for PLINK binary files (.bed, .bim, .fam)
#   --info-file     Path to sample information Excel file
#   --out-sample    Output filename for sample metrics (uncompressed)
#   --out-variant   Output filename for variant metrics (uncompressed)
#   --script-dir    Directory containing the python script (calc_metrics.py)
#   --plink2        Path to plink2 executable
#   --tabix         Path to tabix executable
#   --id-col        Column name for Sample ID in info file
#   --group-col     Column name for Group/Outcome in info file
#   --case-value    Value in group column that indicates Case status (others are Control)
#   --tdp-col       Column name for Target Depth in info file
#   --mdp-col       Column name for Mean Depth in info file
#   --threads       Number of threads to use (default: 1)
# ==============================================================================

set -e  # Exit immediately if a command exits with a non-zero status

# Default values
THREADS=1
MIN_AC=0

# specialized parsing for named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --bed-prefix) BED_PREFIX="$2"; shift ;;
        --info-file) INFO_FILE="$2"; shift ;;
        --out-sample) OUT_SAMPLE="$2"; shift ;;
        --out-variant) OUT_VARIANT="$2"; shift ;;
        --script-dir) SCRIPT_DIR="$2"; shift ;;
        --plink2) PLINK2="$2"; shift ;;
        --tabix) TABIX="$2"; shift ;;
        --id-col) ID_COL="$2"; shift ;;
        --group-col) GROUP_COL="$2"; shift ;;
        --case-value) CASE_VALUE="$2"; shift ;;
        --tdp-col) TDP_COL="$2"; shift ;;
        --mdp-col) MDP_COL="$2"; shift ;;
        --min-ac) MIN_AC="$2"; shift ;;
        --threads) THREADS="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Check required arguments
if [[ -z "$BED_PREFIX" || -z "$INFO_FILE" || -z "$OUT_SAMPLE" || -z "$OUT_VARIANT" || -z "$SCRIPT_DIR" || -z "$PLINK2" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 --bed-prefix <path> --info-file <path> --out-sample <path> --out-variant <path> --script-dir <path> --plink2 <path> [options]"
    exit 1
fi

echo "=========================================================="
echo "Starting Metrics Calculation Pipeline"
echo "Date: $(date)"
echo "PLINK2: $PLINK2"
echo "Input BED: $BED_PREFIX"
echo "MinAC Cutoff: $MIN_AC"
echo "Threads: $THREADS"
echo "=========================================================="

# Step 0: Filter by MinAC (if > 0)
CURRENT_BED="$BED_PREFIX"

if [[ "$MIN_AC" -ne 0 ]]; then
    echo "[Step 0] Filtering variants with MinAC >= $MIN_AC..."
    
    FILTERED_PREFIX="cteph_agp3k.minac${MIN_AC}"
    
    $PLINK2 --threads $THREADS --bfile "$BED_PREFIX" \
        --mac $MIN_AC \
        --make-bed \
        --out "$FILTERED_PREFIX" \
        --silent
    
    # Check if filtering was successful (file exists)
    if [[ -f "${FILTERED_PREFIX}.bed" ]]; then
         # Update prefix to use the filtered files (relative path in current execution dir)
         CURRENT_BED="./$FILTERED_PREFIX"
         echo "         Filtering completed. Working with: $CURRENT_BED"
    else
         echo "Error: Filtering failed. Output file ${FILTERED_PREFIX}.bed not found."
         exit 1
    fi
else
    echo "[Step 0] MinAC is 0, skipping filtering (using original data)."
fi

# 1. Generate standard sample counts (HomRef, Het, HomAlt based on REF genome)
# cols=maybefid,homref,het,homalt -> Requests explicit counts including homozygous reference
echo "[Step 1] Running PLINK2 sample counts..."
$PLINK2 --threads $THREADS --bfile "$CURRENT_BED" \
    --sample-counts cols=maybefid,homref,het,homalt \
    --missing sample-only \
    --out temp_sample

# 2. Generate SMinAC (Sum of Minor Allele Counts)
# Requires aligning REF to Major Allele so that ALT becomes Minor Allele.
# 'maj-ref force' forces this assignment.
# We create a temporary bed file for this purpose because modifiers like --maj-ref 
# cannot be used directly with --sample-counts in some context/versions or purely for safety.
echo "[Step 2] Calculating SMinAC (Minor Allele Counts)..."
$PLINK2 --threads $THREADS --bfile "$CURRENT_BED" \
    --maj-ref force \
    --make-bed \
    --out temp_aligned

$PLINK2 --threads $THREADS --bfile temp_aligned \
    --sample-counts cols=maybefid,het,homalt \
    --out temp_sample_minac

# Clean up temp aligned files
rm temp_aligned.*

# 3. Generate Variant Counts and Missingness
echo "[Step 3] Running PLINK2 variant counts..."
$PLINK2 --threads $THREADS --bfile "$CURRENT_BED" \
    --geno-counts cols=chrom,pos,ref,alt,homref,refalt1,homalt1 \
    --missing variant-only \
    --out temp_variant

# 4. Generate Allele Frequencies
echo "[Step 4] Running PLINK2 allele frequencies..."
$PLINK2 --threads $THREADS --bfile "$CURRENT_BED" \
    --freq counts \
    --out temp_freq

# 5. Run Python Aggregation Script
echo "[Step 5] Running Python aggregation script..."
python3 "$SCRIPT_DIR/calc_metrics.py" \
    --sample-counts temp_sample.scount \
    --sample-missing temp_sample.smiss \
    --sample-minac temp_sample_minac.scount \
    --variant-counts temp_variant.gcount \
    --variant-missing temp_variant.vmiss \
    --freq-counts temp_freq.acount \
    --info "$INFO_FILE" \
    --id-col "$ID_COL" \
    --group-col "${GROUP_COL}" \
    --case-value "${CASE_VALUE}" \
    --tdp-col "$TDP_COL" \
    --mdp-col "$MDP_COL" \
    --out-sample "$OUT_SAMPLE" \
    --out-variant "$OUT_VARIANT" \
    --threads "$THREADS" \
    --log calc_metrics.log

# 6. Compress and Index Outputs
echo "[Step 6] Compressing and indexing outputs..."
# Use bgzip if available in path, or just assume it is (Nextflow env usually provides it)
# If tabix path is provided explicitly, use it, otherwise assume 'tabix'
TABIX_CMD=${TABIX:-tabix}

if command -v bgzip &> /dev/null; then
    bgzip -f "$OUT_SAMPLE"
    bgzip -f "$OUT_VARIANT"
    
    # Index variant metrics (Sample metrics are not strictly suitable for tabix indexing as they are not genomic range data)
    $TABIX_CMD -f -s 1 -b 2 -e 2 "${OUT_VARIANT}.gz"
else
    echo "Warning: bgzip not found. Outputs are left uncompressed."
fi

# 7. Cleanup Intermediate Files
echo "[Step 7] Cleaning up intermediate files..."

# Remove filtered PLINK files if they were created
if [[ "$MIN_AC" -ne 0 ]]; then
    rm -f "cteph_agp3k.minac${MIN_AC}.bed" "cteph_agp3k.minac${MIN_AC}.bim" "cteph_agp3k.minac${MIN_AC}.fam" "cteph_agp3k.minac${MIN_AC}.log"    
fi

# Remove temporary metric files
rm -f temp_sample.scount temp_sample.smiss temp_sample.log \
      temp_sample_minac.scount temp_sample_minac.log \
      temp_variant.gcount temp_variant.vmiss temp_variant.log \
      temp_freq.acount temp_freq.log

echo "=========================================================="
echo "Pipeline Completed Successfully"
echo "=========================================================="
