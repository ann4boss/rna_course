#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=8000
#SBATCH --time=13:00:00
#SBATCH --array=0-11
#SBATCH --job-name=mapping
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_mapping_%A_%a.e
#SBATCH --output=/data/users/aboss/rna_course/output_mapping_%A_%a.o
#SBATCH --partition=pcoursea


# Directory with input files fastq from cleaned files
#READS_DIR=/data/users/aboss/rna_course/03_cleaning_results
READS_DIR=/data/users/aboss/rna_course/03_cleaning_results
# Directory for mapping output    
OUTPUT_DIR=/data/users/aboss/rna_course/06_mapping_results_stranded
mkdir -p ${OUTPUT_DIR}
# Directory with reference genome and indexed files
REFERENCE_GENOME_DIR=/data/users/aboss/rna_course/04_reference_genome
BASENAME='Homo_sapiens.GRCh38_indexed'
HISAT2_INDEX=/data/users/aboss/rna_course/05_indexing/${BASENAME}
# Samples
SAMPLES=("HER21" "HER22" "HER23" "NonTNBC1" "NonTNBC2" "NonTNBC3" "Normal1" "Normal2" "Normal3" "TNBC1" "TNBC2" "TNBC3")
# Apptainer paths
APPTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif


# Get the current sample based on the job array index
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
READ1="${READS_DIR}/${SAMPLE}_cleaned_R1.fastq.gz"
READ2="${READS_DIR}/${SAMPLE}_cleaned_R2.fastq.gz"

#--------- Mapping with HISAT2 and Post-Processing (SAM to BAM, Sort, Index) -----------------
# Output SAM and BAM files directly in OUTPUT_DIR
BAM_OUTPUT="${OUTPUT_DIR}/${SAMPLE}_aligned.sorted.bam"
    
# Run HISAT2 (for aligning) and samtools (for converting to BAM, sorting, and indexing) using Apptainer
apptainer exec --bind ${REFERENCE_GENOME_DIR} --bind ${READS_DIR} --bind ${OUTPUT_DIR} ${APPTAINER} bash -c "
hisat2 -x ${HISAT2_INDEX} -1 ${READ1} -2 ${READ2} | \
samtools view -b - | \
samtools sort -o ${BAM_OUTPUT} - && \
samtools index ${BAM_OUTPUT} && \
samtools flagstat ${BAM_OUTPUT} > ${OUTPUT_DIR}/${SAMPLE}_flagstat.txt && \
samtools idxstats ${BAM_OUTPUT} > ${OUTPUT_DIR}/${SAMPLE}_idxstats.txt && \
samtools depth ${BAM_OUTPUT} > ${OUTPUT_DIR}/${SAMPLE}_coverage.txt
" 2> ${OUTPUT_DIR}/${SAMPLE}_mapping.log

# Calculate average coverage
awk '{sum+=$3} END {if (NR > 0) print "Average Coverage: " sum/NR; else print "Average Coverage: 0"}' \
    ${OUTPUT_DIR}/${SAMPLE}_coverage.txt > ${OUTPUT_DIR}/${SAMPLE}_coverage_summary.txt

echo "Mapping completed for Sample: ${SAMPLE}"