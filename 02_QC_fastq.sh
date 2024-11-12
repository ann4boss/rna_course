#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastq
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_fastq_%A_%a.e


# Directory with soft link to RNA-seq paired-end reads (try with readycopy)
READS_DIR="/data/users/aboss/rna_course/readscopy"
mkdir -p ${READS_DIR}
# Directory for FastQC output    
OUTPUT_DIR="/data/users/aboss/rna_course/fastqc_results"
mkdir -p ${OUTPUT_DIR}
# strains
STRAINS=("Normal3" "TNBC1" "TNBC2" "TNBC3")
#STRAINS=("HER21" "HER22" "HER23" "NonTNBC1" "NonTNBC2" "NonTNBC3" "Normal1" "Normal2" "Normal3" "TNBC1" "TNBC2" "TNBC3")
# apptainer paths
APPTAINER_fastqc="/containers/apptainer/fastqc-0.12.1.sif"
APPTAINER_multiqc="/containers/apptainer/multiqc-1.19.sif"

# Loop through all strains
for STRAIN in "${STRAINS[@]}"; do
    # Define the R1 and R2 read files based on the strain
    READ1="${READS_DIR}/${STRAIN}_R1.fastq.gz"
    READ2="${READS_DIR}/${STRAIN}_R2.fastq.gz"
    
    # Check if both READ1 and READ2 files exist
    if [[ -e "${READ1}" && -e "${READ2}" ]]; then
        echo "Running FastQC for strain: ${STRAIN}"
        
        # Run FastQC using Apptainer
        apptainer exec ${APPTAINER_fastqc} fastqc -o ${OUTPUT_DIR} ${READ1} ${READ2}
        
        # Check if FastQC was successful
        if [[ $? -eq 0 ]]; then
            echo "FastQC completed for ${STRAIN}."
        else
            echo "Error running FastQC for ${STRAIN}."
        fi
    else
        # If either READ1 or READ2 file is missing, print a warning
        echo "Missing files for strain: ${STRAIN}. Skipping FastQC."
    fi
done


#-------------------------MultiQC-----------------------------------------------------
#run multiQC for all fastqc files to compare ussing Apptainer
apptainer exec ${APPTAINER_multiqc} multiqc ${OUTPUT_DIR} -o ${OUTPUT_DIR}







                          