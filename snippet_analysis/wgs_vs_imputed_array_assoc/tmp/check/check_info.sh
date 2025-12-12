#!/bin/bash

set -e
set -u

IMPUTED_ARRAY_ASSOC='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/wgs_vs_imputed_array_assoc/results/02.assoc_result/array/wgs_vs_array.array.sex.10pc.additive.PHENO1.glm.logistic'
WGS_ASSOC='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/wgs_vs_imputed_array_assoc/results/02.assoc_result/wgs/wgs_vs_array.wgs.sex.10pc.additive.PHENO1.glm.logistic'
IMPUTED_VCF_PATH='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/wgs_vs_imputed_array/results/01.normalized'
TARGET_LOCI='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/wgs_vs_imputed_array_assoc/tmp/check/peak.tsv'

OUTPUT_FILE="$(dirname $TARGET_LOCI)/peak_info.tsv"

echo "=========================================="
echo "Extracting variant information"
echo "=========================================="
echo "Target loci file: $TARGET_LOCI"
echo "WGS association:  $WGS_ASSOC"
echo "Array association: $IMPUTED_ARRAY_ASSOC"
echo "VCF directory:    $IMPUTED_VCF_PATH"
echo "Output file:      $OUTPUT_FILE"
echo ""

# Check if input files exist
for file in "$TARGET_LOCI" "$WGS_ASSOC" "$IMPUTED_ARRAY_ASSOC"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: File not found: $file"
        exit 1
    fi
done

if [ ! -d "$IMPUTED_VCF_PATH" ]; then
    echo "ERROR: Directory not found: $IMPUTED_VCF_PATH"
    exit 1
fi

# Write header
echo -e "VARIANT_ID\tWGS_BETA\tARRAY_BETA\tWGS_SE\tARRAY_SE\tWGS_P\tARRAY_P\tWGS_AF\tARRAY_AF\tIMPUTED_MARKER_TYPE\tIMPUTED_R2" > "$OUTPUT_FILE"

echo "Processing variants..."

# Get column numbers from association files
wgs_header=$(head -n 1 "$WGS_ASSOC")
array_header=$(head -n 1 "$IMPUTED_ARRAY_ASSOC")

# Function to get column number by name
get_col_num() {
    local header="$1"
    local col_name="$2"
    echo "$header" | tr '\t' '\n' | nl -v 1 | grep -w "$col_name" | awk '{print $1}'
}

WGS_BETA_COL=$(get_col_num "$wgs_header" "BETA")
WGS_SE_COL=$(get_col_num "$wgs_header" "SE")
WGS_P_COL=$(get_col_num "$wgs_header" "P")
WGS_AF_COL=$(get_col_num "$wgs_header" "A1_FREQ")

ARRAY_BETA_COL=$(get_col_num "$array_header" "BETA")
ARRAY_SE_COL=$(get_col_num "$array_header" "SE")
ARRAY_P_COL=$(get_col_num "$array_header" "P")
ARRAY_AF_COL=$(get_col_num "$array_header" "A1_FREQ")

echo "Column positions detected:"
echo "  WGS: BETA=$WGS_BETA_COL, SE=$WGS_SE_COL, P=$WGS_P_COL, A1_FREQ=$WGS_AF_COL"
echo "  Array: BETA=$ARRAY_BETA_COL, SE=$ARRAY_SE_COL, P=$ARRAY_P_COL, A1_FREQ=$ARRAY_AF_COL"
echo ""

# Read TARGET_LOCI file (all lines, no header)
# Use || [ -n "$variant_id" ] to handle last line without newline
while IFS=$'\t' read -r variant_id rest || [ -n "$variant_id" ]; do
    # Skip empty lines
    [ -z "$variant_id" ] && continue
    
    echo "  Processing: $variant_id"
    
    # Extract chromosome, position, ref, alt from VARIANT_ID (format: chr1:146895349:A:G)
    chrom=$(echo "$variant_id" | cut -d':' -f1)
    pos=$(echo "$variant_id" | cut -d':' -f2)
    ref=$(echo "$variant_id" | cut -d':' -f3)
    alt=$(echo "$variant_id" | cut -d':' -f4)
    
    # Construct VCF file path directly (more reliable than find)
    vcf_file="${IMPUTED_VCF_PATH}/${chrom}.normalized.vcf.gz"
    
    if [ ! -f "$vcf_file" ]; then
        echo "    WARNING: VCF file not found: $vcf_file"
        marker_type="NA"
        r2="NA"
    else
        # Use tabix to query VCF by position (try with chr prefix first)
        vcf_line=$(tabix "$vcf_file" "${chrom}:${pos}-${pos}" 2>/dev/null)
        
        # If not found, try without chr prefix
        if [ -z "$vcf_line" ]; then
            vcf_line=$(tabix "$vcf_file" "${chrom#chr}:${pos}-${pos}" 2>/dev/null)
        fi
        
        if [ -z "$vcf_line" ]; then
            echo "    WARNING: No variants found at position ${chrom}:${pos}"
            marker_type="NA"
            r2="NA"
        else
            # Find the matching variant by comparing REF and ALT (columns 4 and 5)
            matching_line=$(echo "$vcf_line" | awk -v ref="$ref" -v alt="$alt" '$4==ref && $5==alt')
            
            if [ -z "$matching_line" ]; then
                echo "    WARNING: Variant ${ref}>${alt} not found at ${chrom}:${pos}"
                echo "    Available variants:"
                echo "$vcf_line" | awk '{print "      " $4 ">" $5}'
                marker_type="NA"
                r2="NA"
            else
                # Extract INFO field (8th column)
                info_field=$(echo "$matching_line" | cut -f8)
                
                # Check for IMPUTED or TYPED flags
                has_typed=$(echo "$info_field" | grep -o "TYPED" || echo "")
                has_imputed=$(echo "$info_field" | grep -o "IMPUTED" || echo "")
                
                if [ -n "$has_typed" ] && [ -n "$has_imputed" ]; then
                    marker_type="TYPED-IMPUTED"  # Both flags present
                elif [ -n "$has_typed" ]; then
                    marker_type="TYPED"
                elif [ -n "$has_imputed" ]; then
                    marker_type="IMPUTED"
                else
                    marker_type="UNKNOWN"
                fi
                
                # Extract R2 value (split by semicolon, match field starting with R2=)
                r2=$(echo "$info_field" | awk -F';' '{for(i=1;i<=NF;i++) if($i ~ /^R2=/) {sub(/^R2=/,"",$i); print $i; exit}}')
                [ -z "$r2" ] && r2="NA"
            fi
        fi
    fi
    
    # Extract WGS association data
    wgs_data=$(grep -w "$variant_id" "$WGS_ASSOC" | head -n 1)
    if [ -z "$wgs_data" ]; then
        echo "    WARNING: Variant not found in WGS association: $variant_id"
        wgs_beta="NA"
        wgs_se="NA"
        wgs_p="NA"
        wgs_af="NA"
    else
        wgs_beta=$(echo "$wgs_data" | awk -v col="$WGS_BETA_COL" '{print $col}')
        wgs_se=$(echo "$wgs_data" | awk -v col="$WGS_SE_COL" '{print $col}')
        wgs_p=$(echo "$wgs_data" | awk -v col="$WGS_P_COL" '{print $col}')
        wgs_af=$(echo "$wgs_data" | awk -v col="$WGS_AF_COL" '{print $col}')
    fi
    
    # Extract Array association data
    array_data=$(grep -w "$variant_id" "$IMPUTED_ARRAY_ASSOC" | head -n 1)
    if [ -z "$array_data" ]; then
        echo "    WARNING: Variant not found in Array association: $variant_id"
        array_beta="NA"
        array_se="NA"
        array_p="NA"
        array_af="NA"
    else
        array_beta=$(echo "$array_data" | awk -v col="$ARRAY_BETA_COL" '{print $col}')
        array_se=$(echo "$array_data" | awk -v col="$ARRAY_SE_COL" '{print $col}')
        array_p=$(echo "$array_data" | awk -v col="$ARRAY_P_COL" '{print $col}')
        array_af=$(echo "$array_data" | awk -v col="$ARRAY_AF_COL" '{print $col}')
    fi
    
    # Write to output file
    echo -e "${variant_id}\t${wgs_beta}\t${array_beta}\t${wgs_se}\t${array_se}\t${wgs_p}\t${array_p}\t${wgs_af}\t${array_af}\t${marker_type}\t${r2}" >> "$OUTPUT_FILE"
done < "$TARGET_LOCI"

echo ""
echo "=========================================="
echo "Analysis completed!"
echo "Output file: $OUTPUT_FILE"
echo "Total variants: $(tail -n +2 "$OUTPUT_FILE" | wc -l)"
echo "=========================================="
