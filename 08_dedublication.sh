#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=8000
#SBATCH --time=02:00:00
#SBATCH --array=0-11
#SBATCH --job-name=deduplication
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_dedup_%A_%a.e
#SBATCH --output=/data/users/aboss/rna_course/output_dedup_%A_%a.o
#SBATCH --partition=pcoursea

# Paths
INPUT_DIR=/data/users/aboss/rna_course/06_mapping_results
OUTPUT_DIR=/data/users/aboss/rna_course/08_deduplicated_bam
# Ensure output directory exists
mkdir -p ${OUTPUT_DIR}
# Samples
SAMPLES=("HER21" "HER22" "HER23" "NonTNBC1" "NonTNBC2" "NonTNBC3" "Normal1" "Normal2" "Normal3" "TNBC1" "TNBC2" "TNBC3")
# Path to apptainer
APPTAINER_IMAGE=/containers/apptainer/umi_tools-1.1.3.sif



# Get the current sample based on the SLURM array index
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
INPUT_BAM=${INPUT_DIR}/${SAMPLE}_aligned.sorted.bam
OUTPUT_BAM=${OUTPUT_DIR}/${SAMPLE}_deduplicated.bam
LOG_FILE=${OUTPUT_DIR}/${SAMPLE}_deduplication.log
STATS_FILE=${OUTPUT_DIR}/${SAMPLE}_deduplication_stats.txt

# Deduplication with UMI-tools
apptainer --bind ${INPUT_DIR} exec ${APPTAINER_IMAGE} umi_tools dedup \
    --stdin=${INPUT_BAM} \
    --stdout=${OUTPUT_BAM} \
    --log=${LOG_FILE} \
    --output-stats=${STATS_FILE} \
    --paired

echo "Deduplication completed for ${SAMPLE}. Output saved to ${OUTPUT_BAM}"
