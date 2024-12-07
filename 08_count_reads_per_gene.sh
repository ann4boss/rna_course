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
OUTPUT_DIR=$(realpath ./analysis/06_featureCounts_results)
mkdir -p ${OUTPUT_DIR}
OUTPUT_FILE=${OUTPUT_DIR}/gene_counts.txt
# Apptainer paths
APPTAINER=/containers/apptainer/subread_2.0.6.sif

# List all BAM files in the BAM directory
BAM_FILES=(${BAM_DIR}/*.bam)


# Run featureCounts
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
        -M 20 \
        "${BAM_FILES[@]}" 2> ${OUTPUT_DIR}/featureCounts_log.log

echo "featureCounts completed. Output saved to ${OUTPUT_FILE}"

# Adjust the header of the output file to include only sample names
awk 'BEGIN {OFS="\t"} 
    NR==1 {
        for (i=2; i<=NF; i++) {
            split($i, a, "/")
            split(a[length(a)], b, "_")
            $i = b[1]
        }
    }
    {print}' gene_counts.txt.summary > gene_counts.txt.summary

# Analyze the output
SUMMARY_FILE=${OUTPUT_DIR}/gene_counts_summary.txt

# Extract summary information
awk 'NR > 2 && $1 == "Assigned" {assigned_reads+=$2}
     NR > 2 && $1 == "Ambiguity" {ambiguous_reads+=$2; ambiguous_count++}
     END {
         print "Total Assigned Reads:", assigned_reads
         print "Total Ambiguous Reads:", ambiguous_reads
         print "Average Ambiguous Reads Per Sample:", ambiguous_reads/ambiguous_count
     }' ${OUTPUT_FILE}.summary > ${SUMMARY_FILE}
