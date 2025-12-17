#!/bin/bash
# Append extraction statistics to summary file
# Usage: append_extraction_stats.sh <summary_file> <dataset> <orig_samples> <orig_variants> <ext_samples> <ext_variants> <plink2_log>

SUMMARY_FILE="$1"
DATASET="$2"
ORIG_SAMPLES="$3"
ORIG_VARIANTS="$4"
EXT_SAMPLES="$5"
EXT_VARIANTS="$6"
PLINK2_LOG="$7"

# Calculate percentages
SAMPLE_PCT=$(awk "BEGIN {printf \"%.2f\", ($EXT_SAMPLES/$ORIG_SAMPLES)*100}")
VARIANT_PCT=$(awk "BEGIN {printf \"%.2f\", ($EXT_VARIANTS/$ORIG_VARIANTS)*100}")

# Create directory if not exists
mkdir -p "$(dirname "$SUMMARY_FILE")"

# Use flock to prevent concurrent write conflicts
(
    flock -x 200

    # If file doesn't exist or is empty, create header
    if [ ! -s "$SUMMARY_FILE" ]; then
        REPORT_TIME=$(date '+%Y-%m-%d %H:%M:%S')
        cat > "$SUMMARY_FILE" << EOF
====================================================================================================
                         WGS vs IMPUTED ARRAY COMPARISON REPORT
                              EXTRACTION STATISTICS SUMMARY
====================================================================================================

Report generated: $REPORT_TIME

This report summarizes the extraction of common samples and variants across three datasets:
  - WGS:              Whole Genome Sequencing data
  - Imputed Array GT: Array-based imputed genotypes (hard-called, 0/1/2)
  - Imputed Array DS: Array-based imputed genotypes (dosage, 0.0-2.0)

All datasets have been subset to include only the samples and variants present in all three datasets.

----------------------------------------------------------------------------------------------------
Dataset              Original Samples    Original Variants    Extracted Samples   Extracted Variants   Sample %    Variant %
----------------------------------------------------------------------------------------------------
EOF
    fi
    
    # Append data line with proper dataset name
    DATASET_DISPLAY="$DATASET"
    if [ "$DATASET" == "wgs" ]; then
        DATASET_DISPLAY="WGS"
    elif [ "$DATASET" == "array_gt" ]; then
        DATASET_DISPLAY="Imputed Array GT"
    elif [ "$DATASET" == "array_ds" ]; then
        DATASET_DISPLAY="Imputed Array DS"
    fi
    
    printf "%-20s %15s %20s %19s %20s %10s%% %12s%%\n" \
        "$DATASET_DISPLAY" "$ORIG_SAMPLES" "$ORIG_VARIANTS" "$EXT_SAMPLES" "$EXT_VARIANTS" "$SAMPLE_PCT" "$VARIANT_PCT" >> "$SUMMARY_FILE"
    
    # Append detailed PLINK2 log information
    if [ -f "$PLINK2_LOG" ]; then
        echo "" >> "$SUMMARY_FILE"
        echo "  [$DATASET_DISPLAY] PLINK2 Execution Details:" >> "$SUMMARY_FILE"
        echo "  $(date '+%Y-%m-%d %H:%M:%S') - Processing $DATASET_DISPLAY dataset" >> "$SUMMARY_FILE"
        
        # Extract key information from PLINK2 log
        if grep -q "Error:" "$PLINK2_LOG"; then
            echo "  ⚠ ERROR detected in PLINK2 execution" >> "$SUMMARY_FILE"
            grep "Error:" "$PLINK2_LOG" | sed 's/^/    /' >> "$SUMMARY_FILE"
        else
            # Extract processing time
            START_TIME=$(grep "Start time:" "$PLINK2_LOG" | tail -1)
            END_TIME=$(grep "End time:" "$PLINK2_LOG" | tail -1)
            if [ -n "$START_TIME" ] && [ -n "$END_TIME" ]; then
                echo "  $START_TIME" >> "$SUMMARY_FILE"
                echo "  $END_TIME" >> "$SUMMARY_FILE"
            fi
            
            # Extract key statistics if available
            if grep -q "variants loaded" "$PLINK2_LOG"; then
                grep "variants loaded" "$PLINK2_LOG" | tail -1 | sed 's/^/  /' >> "$SUMMARY_FILE"
            fi
            if grep -q "samples" "$PLINK2_LOG"; then
                grep "samples" "$PLINK2_LOG" | head -1 | sed 's/^/  /' >> "$SUMMARY_FILE"
            fi
            
            echo "  ✓ Extraction completed successfully" >> "$SUMMARY_FILE"
        fi
        echo "" >> "$SUMMARY_FILE"
    fi
    
    # Check if all three datasets are completed (use display names)
    LINE_COUNT=$(grep -c "^WGS\|^Imputed Array GT\|^Imputed Array DS" "$SUMMARY_FILE" 2>/dev/null || echo 0)
    if [ "$LINE_COUNT" -eq 3 ]; then
        # Add footer with output locations
        cat >> "$SUMMARY_FILE" << 'EOF'
----------------------------------------------------------------------------------------------------

EXTRACTED GENOTYPE FILES LOCATION:
  
  WGS genotypes (from whole genome sequencing):
    02.extracted_genotypes/wgs/cteph_agp3k.wgs.common.{pgen,pvar,psam}
    
  Imputed Array GT genotypes (hard-called, 0/1/2):
    02.extracted_genotypes/array_gt/cteph_agp3k.array_gt.common.{pgen,pvar,psam}
    
  Imputed Array DS genotypes (with dosage, 0.0-2.0):
    02.extracted_genotypes/array_ds/cteph_agp3k.array_ds.common.{pgen,pvar,psam}

ADDITIONAL INFORMATION:
  - Common samples list:     01.common_samples_variants/common_samples.txt
  - Common variants list:    01.common_samples_variants/common_variants.txt
  - Detailed intersection:   01.common_samples_variants/intersection_summary.txt

PROCESSING LOG:
  All PLINK2 execution details are recorded above in the respective dataset sections.
  Check for any warnings or errors marked with ⚠ symbol.

NOTES:
  All extracted genotype files are in PLINK2 format (pgen/pvar/psam).
  Files contain only the common samples and variants identified across all three datasets.
  
  DATA TYPE CLARIFICATION:
    - WGS: Direct sequencing data from whole genome sequencing
    - Imputed Array GT: Array genotyping followed by imputation, hard-called genotypes
    - Imputed Array DS: Array genotyping followed by imputation, includes dosage information
  
  Use these files for downstream comparative analyses between WGS and imputed array genotyping.

====================================================================================================
                                    END OF REPORT
====================================================================================================
EOF
    fi

) 200>"$SUMMARY_FILE.lock"

# Clean up lock file
rm -f "$SUMMARY_FILE.lock"
