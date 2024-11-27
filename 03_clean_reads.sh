#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=8000
#SBATCH --time=1:00:00
#SBATCH --job-name=fastp_cleaning
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_cleaning_%A_%a.e
#SBATCH --output=/data/users/aboss/rna_course/output_cleaning_%A_%a.o
#SBATCH --array=0-11
#SBATCH --ntasks=1
#SBATCH --partition=pcoursea


# Directory with soft link to RNA-seq paired-end reads (-> link not working, tried with readycopy)
READS_DIR=/data/users/aboss/rna_course/readscopy
# Directory for FastQC output    
OUTPUT_DIR=/data/users/aboss/rna_course/03_cleaning_results
mkdir -p ${OUTPUT_DIR}
# Path to uncleaned fastqc files
UNCLEANED_FASTQC=/data/users/aboss/rna_course/02_fastqc_results
# Samples
SAMPLES=("HER21" "HER22" "HER23" "NonTNBC1" "NonTNBC2" "NonTNBC3" "Normal1" "Normal2" "Normal3" "TNBC1" "TNBC2" "TNBC3")
# apptainer paths
APPTAINER_fastp=/containers/apptainer/fastp_0.23.2--h5f740d0_3.sif
APPTAINER_multiqc=/containers/apptainer/multiqc-1.19.sif
APPTAINER_fastqc=/containers/apptainer/fastqc-0.12.1.sif

# Get the sample for the current job array index
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
READ1="${READS_DIR}/${SAMPLE}_R1.fastq.gz"
READ2="${READS_DIR}/${SAMPLE}_R2.fastq.gz"

# Run fastp for cleaning
echo "Running fastp for sample: ${SAMPLE}"
apptainer exec --bind ${READS_DIR},${OUTPUT_DIR} ${APPTAINER_fastp} \
    fastp \
    -i ${READ1} -I ${READ2} \
    -o ${OUTPUT_DIR}/${SAMPLE}_cleaned_R1.fastq.gz -O ${OUTPUT_DIR}/${SAMPLE}_cleaned_R2.fastq.gz \
    -j 4 \
    --detect_adapter_for_pe \
    --low_complexity_filter \
    --complexity_threshold 30 \
    --overrepresentation_analysis \
    --qualified_quality_phred 15 \
    --cut_by_quality3 -q 20 \
    --unqualified_percent_limit 40 \
    --length_required 50 \
    --n_base_limit 1 \
    --correction \
    --html ${OUTPUT_DIR}/${SAMPLE}_report.html


# Run FastQC on cleaned reads
echo "Running FastQC for cleaned reads of sample: ${SAMPLE}"
apptainer exec --bind ${OUTPUT_DIR} ${APPTAINER_fastqc} fastqc -o ${OUTPUT_DIR} ${OUTPUT_DIR}/${SAMPLE}_cleaned_R1.fastq.gz ${OUTPUT_DIR}/${SAMPLE}_cleaned_R2.fastq.gz


# this needs to change otherwise it runs for each array!! -> create separate script
# Run MultiQC for all FastQC results
apptainer exec --bind ${OUTPUT_DIR},${UNCLEANED_FASTQC} ${APPTAINER_multiqc} multiqc ${OUTPUT_DIR} -o ${OUTPUT_DIR}
