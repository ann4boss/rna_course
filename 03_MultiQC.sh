#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=multiqc
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_multiqc_%A_%a.e


# Ensure the script stops on error
set -euo pipefail

# Check for correct usage
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

# Directory with soft links to RNA-seq paired-end reads
INPUT_DIR=$1
mkdir -p "${INPUT_DIR}"

# Directory for FastQC output
OUTPUT_DIR=$2
mkdir -p "${OUTPUT_DIR}"

# Samples
SAMPLES=("HER21" "HER22" "HER23" "NonTNBC1" "NonTNBC2" "NonTNBC3" "Normal1" "Normal2" "Normal3" "TNBC1" "TNBC2" "TNBC3")

# Apptainer paths
APPTAINER_MULTIQC="/containers/apptainer/multiqc-1.19.sif"

#-------------------------MultiQC-----------------------------------------------------
# Run MultiQC for all FastQC files to compare using Apptainer
apptainer exec --bind ${INPUT_DIR} ${APPTAINER_MULTIQC} multiqc ${OUTPUT_DIR} -o ${OUTPUT_DIR}
