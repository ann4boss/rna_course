#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=8000
#SBATCH --time=13:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=mapping
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_mapping_%A_%a.e



# Directory for mapping output    
OUTPUT_DIR=/data/users/aboss/rna_course/mapping_results
mkdir -p ${OUTPUT_DIR}
# Directory for processed BAM files
PROCESSED_BAM_DIR=${OUTPUT_DIR}/processed_bam_files
mkdir -p ${PROCESSED_BAM_DIR}
# SAMPLES
#SAMPLES=("HER21" "HER22" "HER23" "NonTNBC1" "NonTNBC2" "NonTNBC3" "Normal1" "Normal2" "Normal3" "TNBC1" "TNBC2" "TNBC3")
SAMPLES=("HER21")
# Apptainer paths
APPTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif



#---------Mapping with HISAT2 and Post-Processing (SAM to BAM, Sort, Index)-----------------
# Loop through each sample
for SAMPLE in "${SAMPLES[@]}"; do
    # Convert SAM to BAM, sort, and index using Samtools
    # Convert SAM to BAM
    apptainer exec ${APPTAINER} samtools view -Sb "${SAM_OUTPUT}" > "${PROCESSED_BAM_DIR}/${SAMPLE}_aligned.bam"
    # Sort the BAM file by genomic coordinates
    apptainer exec ${APPTAINER} samtools sort "${PROCESSED_BAM_DIR}/${SAMPLE}_aligned.bam" -o "${PROCESSED_BAM_DIR}/${SAMPLE}_aligned.sorted.bam"
    # Index the sorted BAM file
    apptainer exec ${APPTAINER} samtools index "${PROCESSED_BAM_DIR}/${SAMPLE}_aligned.sorted.bam"
    # Remove the intermediate unsorted BAM file to save space
    rm "${PROCESSED_BAM_DIR}/${SAMPLE}_aligned.bam"

    echo "BAM file sorted and indexed for ${SAMPLE}."
else
    echo "Error running HISAT2 for ${SAMPLE}. Skipping BAM conversion and sorting."
fi
done

echo "BAM conversion, sorting, and indexing completed for all SAMPLES."
