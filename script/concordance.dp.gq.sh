#usr/bin/zsh

PICARD_SIG="/home/b/b37974/simg/picard.sif"
TRUTH_VCF=$1
CALL_VCF=$2

DP=$3
GQ=$4

# limit the number of parallel jobs
MAX_JOBS=8
CURRENT_JOBS=0

for SAMPLE in $(bcftools query -l ${TRUTH_VCF} | grep -Fxf - <(bcftools query -l ${CALL_VCF})); do
    # execute the command in the background and increase the job count by 1
    (
        singularity exec -B /LARGE0/gr10478/b37974:/LARGE0/gr10478/b37974 ${PICARD_SIG} java -XX:ParallelGCThreads=4 -jar /home/b/b37974/anaconda3/envs/picard_env/share/picard-3.3.0-0/picard.jar GenotypeConcordance \
            -CALL_VCF ${CALL_VCF} \
            -CALL_SAMPLE ${SAMPLE} \
            --MIN_DP ${DP} \
            --MIN_GQ ${GQ} \
            -TRUTH_VCF ${TRUTH_VCF} \
            -TRUTH_SAMPLE ${SAMPLE} \
            --OUTPUT test_concordance_${SAMPLE}_
    ) &

    # increase the job count
    CURRENT_JOBS=$((CURRENT_JOBS + 1))

    # if the maximum number of parallel jobs is reached, wait for all jobs to finish
    if [ ${CURRENT_JOBS} -ge ${MAX_JOBS} ]; then
        wait
        CURRENT_JOBS=0
    fi
done

# wait for the remaining jobs to finish
wait

# Combine the results into a single file for each metric in order
for SAMPLE in $(bcftools query -l ${TRUTH_VCF} | grep -Fxf - <(bcftools query -l ${CALL_VCF})); do
    for METRIC in genotype_concordance_contingency_metrics genotype_concordance_detail_metrics genotype_concordance_summary_metrics; do
        if [ -f test_concordance_${SAMPLE}_.${METRIC} ]; then
            if [ ! -f DP${DP}_GQ${GQ}.${METRIC} ]; then
                # extract the header and add new columns
                awk 'BEGIN {OFS="\t"} /^## METRICS CLASS/ {print $0} /^VARIANT_TYPE/ {print $0, "DP", "GQ"}' test_concordance_${SAMPLE}_.${METRIC} > DP${DP}_GQ${GQ}.${METRIC}
            fi
            # add the DP and GQ columns to the existing file
            awk -v dp=${DP} -v gq=${GQ} 'BEGIN {OFS="\t"} /^VARIANT_TYPE/ {next} /^#/ {next} NF {print $0, dp, gq}' test_concordance_${SAMPLE}_.${METRIC} >> DP${DP}_GQ${GQ}.${METRIC}
            # remove the temporary file
            rm -f test_concordance_${SAMPLE}_.${METRIC}
        fi
    done
done