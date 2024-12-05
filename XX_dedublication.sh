#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=8000
#SBATCH --time=02:00:00
#SBATCH --array=0-11
#SBATCH --job-name=dedup
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_dedup_%A_%a.e
#SBATCH --output=/data/users/aboss/rna_course/output_dedup_%A_%a.o
#SBATCH --partition=pibu_el8

# Paths
#CLEANED_DIR=./analysis/02_cleaning_results
CLEANED_DIR=/data/users/aboss/rna_course/03_cleaning_results
OUTPUT_DIR=./analysis/03_deduplication_results
mkdir -p ${OUTPUT_DIR}
# Samples
SAMPLES=("HER21" "HER22" "HER23" "NonTNBC1" "NonTNBC2" "NonTNBC3" "Normal1" "Normal2" "Normal3" "TNBC1" "TNBC2" "TNBC3")
# Path to apptainer
APPTAINER=/containers/apptainer/fastqc-0.12.1.sif
# samtools for deduplication -> RMDUP

# Create input file for FastUniq
FASTUNIQ_INPUT=${OUTPUT_DIR}/fastuniq_input.txt
echo "Preparing FastUniq input list..."
> ${FASTUNIQ_INPUT}
for SAMPLE in "${SAMPLES[@]}"; do
    echo "${CLEANED_DIR}/${SAMPLE}_cleaned_R1.fastq.gz ${CLEANED_DIR}/${SAMPLE}_cleaned_R2.fastq.gz" >> ${FASTUNIQ_INPUT}
done

# Run FastUniq for deduplication
echo "Running FastUniq for deduplication..."
apptainer exec --bind ${CLEANED_DIR},${OUTPUT_DIR} ${APPTAINER_FASTUNIQ} fastuniq \
    -i ${FASTUNIQ_INPUT} \
    -o ${OUTPUT_DIR}/deduplicated_R1.fastq.gz \
    -p ${OUTPUT_DIR}/deduplicated_R2.fastq.gz

# Loop to handle individual sample outputs
echo "Splitting deduplicated reads into individual samples..."
for SAMPLE in "${SAMPLES[@]}"; do
    zcat ${OUTPUT_DIR}/deduplicated_R1.fastq.gz | grep -A 3 "^@.*${SAMPLE}_" | gzip > ${OUTPUT_DIR}/${SAMPLE}_deduplicated_R1.fastq.gz
    zcat ${OUTPUT_DIR}/deduplicated_R2.fastq.gz | grep -A 3 "^@.*${SAMPLE}_" | gzip > ${OUTPUT_DIR}/${SAMPLE}_deduplicated_R2.fastq.gz
done
