#!/bin/bash
#SBATCH -p gr10478b
#SBATCH -t 7-00:00:00
#SBATCH --rsc p=1:t=1:c=4:m=8G
#SBATCH -o nextflow_head.log
#SBATCH -e nextflow_head.err

# ==============================================================================
# Nextflow Pipeline Launcher
# ==============================================================================
# Author:
#   ZHAO TIE
#
# Description:
#   Sets up the environment and submits the Nextflow pipeline.
#   Includes memory configuration and run artifact cleanup.
#
# Usage:
#   sbatch run_nextflow.sh [script.nf] [clean]
#   - script.nf: Optional Nextflow script (default: tuning.concordance.v5.nf)
#   - 'clean'  : Run without -resume (start fresh).
#   Order of arguments does not matter.
# ==============================================================================

# Activate environment
source activate dsl1

# Increase Nextflow Java Heap memory to prevent OOM when scheduling many tasks
# Physical memory 8G, leave 1G for system, max Java Heap set to 7G
export NXF_OPTS='-Xms2g -Xmx7g'

# Default configurations
NF_SCRIPT="tuning.concordance.v5.nf"
RESUME_FLAG="-resume"

# Parse arguments
for arg in "$@"; do
    if [[ "$arg" == *.nf ]]; then
        NF_SCRIPT="$arg"
    elif [[ "$arg" == "clean" ]]; then
        echo "Starting fresh run (cleaning cache)..."
        RESUME_FLAG=""
    fi
done

if [ -n "$RESUME_FLAG" ]; then
    echo "Resuming from last checkpoint..."
fi

# Detect and remove existing trace.txt and report.html to prevent conflicts
if [ -f trace.txt ]; then
    rm trace.txt
fi
if [ -f report.html ]; then
    rm report.html
fi

# Run Nextflow
# -resume: Resume interrupted tasks (default behavior unless 'clean' is specified)
# -with-trace: Generate real-time task tracing file 'trace.txt'
# -with-report: Generate execution report 'report.html'
nextflow run $NF_SCRIPT $RESUME_FLAG -with-trace trace.txt -with-report report.html
