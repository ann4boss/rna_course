#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=4000
#SBATCH --time=13:00:00
#SBATCH --job-name=count_read_per_gene
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_countread_%A_%a.e
#SBATCH --partition=pcoursea

# Working directory
WORK_DIR=/data/users/aboss/rna_course/
# Directory containing the processed BAM files
BAM_DIR=/data/users/aboss/rna_course/06_mapping_results_uncleaned_stranded
# Annotation file path (e.g., GTF file)
ANNOTATION_FILE=/data/users/aboss/rna_course/04_reference_genome/Homo_sapiens.GRCh38.113.gtf
# Output directory and file
OUTPUT_DIR=/data/users/aboss/rna_course/09_featureCounts_results_stranded
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


# Analyze the output
SUMMARY_FILE=${OUTPUT_DIR}/gene_counts_summary.txt
echo "Processing featureCounts output to generate summary..."

# Extract summary information
awk 'NR > 2 && $1 == "Assigned" {assigned_reads+=$2}
     NR > 2 && $1 == "Ambiguity" {ambiguous_reads+=$2; ambiguous_count++}
     END {
         print "Total Assigned Reads:", assigned_reads
         print "Total Ambiguous Reads:", ambiguous_reads
         print "Average Ambiguous Reads Per Sample:", ambiguous_reads/ambiguous_count
     }' ${OUTPUT_FILE}.summary > ${SUMMARY_FILE}

echo "Summary saved to ${SUMMARY_FILE}"