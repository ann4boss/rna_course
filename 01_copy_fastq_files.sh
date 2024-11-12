#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=linkfiles
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_fastq_%A_%a.e

#Reads source directory
SOURCE_DIR="/data/courses/rnaseq_course/breastcancer_de/reads"
# Define the directory containing the soft link to RNA-seq reads
READS_DIR="/data/users/aboss/rna_course/reads"
mkdir -p ${READS_DIR}
# Directory for FastQC output    
OUTPUT_DIR="/data/users/aboss/rna_course/fastqc_results"
mkdir -p ${OUTPUT_DIR}
# strains
STRAINS=("HER21" "HER22" "HER23" "NonTNBC1" "NonTNBC2" "NonTNBC3" "Normal1" "Normal2" "Normal3" "TNBC1" "TNBC2" "TNBC3")

# link to reads
# Loop through each strain
for STRAIN in "${STRAINS[@]}"; do
    # Define the R1 and R2 file paths
    R1_FILE="${STRAIN}_R1.fastq.gz"
    R2_FILE="${STRAIN}_R2.fastq.gz"
    # Create the soft link for the R1 file
    ln -s ${SOURCE_DIR}/${R1_FILE} ${READS_DIR}/${R1_FILE}
    # Create the soft link for the R2 file
    ln -s ${SOURCE_DIR}/${R2_FILE} ${READS_DIR}/${R2_FILE}
done