#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --time=0:30:00
#SBATCH --job-name=copyfiles
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_file_copy_%A_%a.e
#SBATCH --partition=pibu_el8

# Reads source directory
SOURCE_DIR=/data/courses/rnaseq_course/breastcancer_de/reads
# Define the directory to copy RNA-seq reads
READS_DIR=./data/01_reads_copy
mkdir -p ${READS_DIR}

# Define the sample names
SAMPLES=("HER21" "HER22" "HER23" "NonTNBC1" "NonTNBC2" "NonTNBC3" "Normal1" "Normal2" "Normal3" "TNBC1" "TNBC2" "TNBC3")

# Loop through each sample to copy files
for SAMPLE in "${SAMPLES[@]}"; do
    # Define the R1 and R2 file paths
    R1_FILE="${SAMPLE}_R1.fastq.gz"
    R2_FILE="${SAMPLE}_R2.fastq.gz"
    
    # Copy the R1 file
    cp ${SOURCE_DIR}/${R1_FILE} ${READS_DIR}/${R1_FILE}
    # Copy the R2 file
    cp ${SOURCE_DIR}/${R2_FILE} ${READS_DIR}/${R2_FILE}
done

echo "File copying completed."
