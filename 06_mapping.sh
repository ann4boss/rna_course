#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=8000
#SBATCH --time=13:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=mapping
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_mapping_%A_%a.e


# Directory with input files fastq from !!! copy and uncleaned!!!
READS_DIR=/data/users/aboss/rna_course/readscopy
# Directory for mapping output    
OUTPUT_DIR=/data/users/aboss/rna_course/06_mapping_results
mkdir -p ${OUTPUT_DIR}
# Directory for processed BAM files
PROCESSED_BAM_DIR=${OUTPUT_DIR}/processed_bam_files
mkdir -p ${PROCESSED_BAM_DIR}
# Directory with reference genome and indexed files
REFERENCE_GENOME_DIR=/data/users/aboss/rna_course/reference_genome
BASENAME='Homo_sapiens.GRCh38_indexed'
HISAT2_INDEX=${REFERENCE_GENOME_DIR}/indexing/${BASENAME}
# SAMPLES
#SAMPLES=("HER21" "HER22" "HER23" "NonTNBC1" "NonTNBC2" "NonTNBC3" "Normal1" "Normal2" "Normal3" "TNBC1" "TNBC2" "TNBC3")
SAMPLES=("NonTNBC1")
# Apptainer paths
APPTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif



#---------Mapping with HISAT2 and Post-Processing (SAM to BAM, Sort, Index)-----------------
# Loop through each sample
for SAMPLE in "${SAMPLES[@]}"; do
    # R1 and R2 read files
    READ1="${READS_DIR}/${SAMPLE}_R1.fastq.gz"
    READ2="${READS_DIR}/${SAMPLE}_R2.fastq.gz"
    
    # Check if both READ1 and READ2 files exist - could be deleted
    if [[ -e "${READ1}" && -e "${READ2}" ]]; then
        echo "Running mapping for Sample: ${SAMPLE}"
        
        # Output SAM file name
        SAM_OUTPUT="${OUTPUT_DIR}/${SAMPLE}_aligned.sam"
        
        # Run HISAT2 (for aligning) and samtools (for converting to BAM, sorting, and indexing) using Apptainer
        apptainer exec --bind ${REFERENCE_GENOME_DIR} --bind ${READS_DIR} --bind ${PROCESSED_BAM_DIR} ${APPTAINER} \
        bash -c "
        hisat2 -x ${HISAT2_INDEX} -1 ${READ1} -2 ${READ2} | \
        samtools view -b - | \
        samtools sort -o ${PROCESSED_BAM_DIR}/${SAMPLE}_aligned.sorted.bam - && \
        samtools index ${PROCESSED_BAM_DIR}/${SAMPLE}_aligned.sorted.bam
        " 2> ${OUTPUT_DIR}/${SAMPLE}_mapping.log
    else
        echo "Warning: One or both read files for ${SAMPLE} are missing. Skipping this sample."
    fi
done

echo "Mapping completed for all Samples."
