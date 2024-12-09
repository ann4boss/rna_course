#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=10000
#SBATCH --time=13:00:00
#SBATCH --job-name=count_read_per_gene
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=./error_countread_%A_%a.e
#SBATCH --partition=pibu_el8

# Directory containing the processed BAM files
BAM_DIR=$(realpath ./analysis/04_mapping_results)
ANNOTATION_FILE=$(realpath ./data/02_reference_genome/Homo_sapiens.GRCh38.113.gtf)
# Output directory and file
OUTPUT_DIR=$(realpath ./data/04_featureCounts_results)
mkdir -p ${OUTPUT_DIR}
OUTPUT_FILE=${OUTPUT_DIR}/gene_counts.txt
SUMMARY_FILE=${OUTPUT_FILE}.summary
# Apptainer paths
APPTAINER=/containers/apptainer/subread_2.0.6.sif

# List all BAM files in the BAM directory
BAM_FILES=(${BAM_DIR}/*.bam)

# Run featureCounts
# -g: Specifies the attribute "gene_id" in the annotation file (GTF/GFF) to be used for counting
# -t: Reads are only counted if they overlap features labeled as exon in the annotation file
# -Q: Sets a minimum mapping quality score for a read to be included
# -p: Ensures that paired-end reads are treated as a single unit
# --countReadPairs: Counts paired-end reads as pairs instead of individual reads
# --minOverlap 1: Sets the minimum number of overlapping bases between a read and a feature to 1 for the read to be assigned -> For short reads 50 bp) a small value like 1–3 is appropriate
apptainer exec --bind ${BAM_DIR}:${BAM_DIR} \
               --bind ${OUTPUT_DIR}:${OUTPUT_DIR} \
               --bind $(dirname ${ANNOTATION_FILE}):$(dirname ${ANNOTATION_FILE}) \
               ${APPTAINER} \
    featureCounts \
        -T 4 \
        -a ${ANNOTATION_FILE} \
        -o ${OUTPUT_FILE} \
        -g gene_id \
        -t exon \
        -Q 10 \
        -p \
        --countReadPairs \
        --minOverlap 1 \
        "${BAM_FILES[@]}" 2> ${OUTPUT_DIR}/featureCounts_log.log

echo "featureCounts completed. Output saved to ${OUTPUT_FILE}"

# Adjust the header of the summary file (gene_counts.txt.summary) to include only sample names
awk 'BEGIN {OFS="\t"} 
    NR == 1 {
        for (i = 2; i <= NF; i++) {
            n = split($i, a, "/")
            m = split(a[n], b, "_")
            $i = b[1]
        }
    }
    {print}' ${SUMMARY_FILE} > ${OUTPUT_DIR}/adjusted_gene_counts.txt.summary

# Adjust the header of the counts file (gene_counts.txt) to include only sample names
awk 'BEGIN {OFS="\t"} 
    NR == 2 { 
        for (i = 7; i <= NF; i++) {
            n = split($i, a, "/")
            m = split(a[n], b, "_")
            $i = b[1]
        }
    }
    {print}' ${OUTPUT_FILE} > ${OUTPUT_DIR}/adjusted_gene_counts.txt