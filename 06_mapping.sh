#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --time=03:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=mapping
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_mapping_%A_%a.e


# Directory with input files fastq from copy and uncleaned!!!
READS_DIR="/data/users/aboss/rna_course/readscopy"
# Directory for mapping output    
OUTPUT_DIR="/data/users/aboss/rna_course/mapping_results"
mkdir -p ${OUTPUT_DIR}
# Directory for processed BAM files
PROCESSED_BAM_DIR="${OUTPUT_DIR}/processed_bam_files"
mkdir -p ${PROCESSED_BAM_DIR}
# Directory with reference genome and indexed files
REFERENCE_GENOME_DIR="/data/users/aboss/rna_course/reference_genome"
HISAT2_INDEX="${REFERENCE_GENOME_DIR}/indexing/"
# Strains
#STRAINS=("HER21" "HER22" "HER23" "NonTNBC1" "NonTNBC2" "NonTNBC3" "Normal1" "Normal2" "Normal3" "TNBC1" "TNBC2" "TNBC3")
STRAINS=("HER21")
# Apptainer paths
APPTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"



#---------Mapping with HISAT2 and Post-Processing (SAM to BAM, Sort, Index)-----------------
# Loop through each strain
for STRAIN in "${STRAINS[@]}"; do
    # R1 and R2 read files based on the strain
    READ1="${READS_DIR}/${STRAIN}_R1.fastq.gz"
    READ2="${READS_DIR}/${STRAIN}_R2.fastq.gz"
    
    # Check if both READ1 and READ2 files exist
    if [[ -e "${READ1}" && -e "${READ2}" ]]; then
        echo "Running HISAT2 mapping for strain: ${STRAIN}"
        
        # Output SAM file name
        SAM_OUTPUT="${OUTPUT_DIR}/${STRAIN}_aligned.sam"
        
        # Run HISAT2 using Apptainer
        apptainer exec ${APPTAINER} hisat2 -x ${HISAT2_INDEX} -1 ${READ1} -2 ${READ2} -S ${SAM_OUTPUT}
        
        # Check if HISAT2 was successful
        if [[ $? -eq 0 ]]; then
            echo "HISAT2 mapping completed for ${STRAIN}. Now converting to BAM format, sorting, and indexing."

            # Convert SAM to BAM, sort, and index using Samtools
            # Convert SAM to BAM
            apptainer exec ${APPTAINER} samtools view -Sb "${SAM_OUTPUT}" > "${PROCESSED_BAM_DIR}/${STRAIN}_aligned.bam"
            # Sort the BAM file by genomic coordinates
            apptainer exec ${APPTAINER} samtools sort "${PROCESSED_BAM_DIR}/${STRAIN}_aligned.bam" -o "${PROCESSED_BAM_DIR}/${STRAIN}_aligned.sorted.bam"
            # Index the sorted BAM file
            apptainer exec ${APPTAINER} samtools index "${PROCESSED_BAM_DIR}/${STRAIN}_aligned.sorted.bam"
            # Remove the intermediate unsorted BAM file to save space
            rm "${PROCESSED_BAM_DIR}/${STRAIN}_aligned.bam"

            echo "BAM file sorted and indexed for ${STRAIN}."
        else
            echo "Error running HISAT2 for ${STRAIN}. Skipping BAM conversion and sorting."
        fi
    else
        # If either READ1 or READ2 file is missing, print a warning
        echo "Missing files for strain: ${STRAIN}. Skipping HISAT2 mapping."
    fi
done

echo "HISAT2 mapping, BAM conversion, sorting, and indexing completed for all strains."
