#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=4000
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=count_read_per_gene
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_countread_%A_%a.e

# Working directory
WORK_DIR=/data/users/aboss/rna_course/
# Directory containing the processed BAM files
BAM_DIR=/data/users/aboss/rna_course/06_mapping_results
# Annotation file path (e.g., GTF file)
ANNOTATION_FILE=/data/users/aboss/rna_course/04_reference_genome/Homo_sapiens.GRCh38.113.gtf
# Output directory and file
OUTPUT_DIR=/data/users/aboss/rna_course/07_featureCounts_results
mkdir -p ${OUTPUT_DIR}
OUTPUT_FILE=${OUTPUT_DIR}/gene_counts.txt
# Apptainer paths
APPTAINER=/containers/apptainer/subread_2.0.1--hed695b0_0.sif


# List all BAM files in the BAM directory
BAM_FILES=(${BAM_DIR}/*.bam)

# Run featureCounts using Apptainer
apptainer exec --bind ${WORK_DIR} ${APPTAINER} \
    featureCounts -T 4 \
                  -a ${ANNOTATION_FILE} \
                  -o ${OUTPUT_FILE} \
                  -g gene_id \
                  -t exon \
                  -p \
                  "${BAM_FILES[@]}"

echo "featureCounts completed. Output saved to ${OUTPUT_FILE}"
