#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=8000
#SBATCH --time=1:00:00
#SBATCH --job-name=fastp_cleaning
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=./error_cleaning_%A_%a.e
#SBATCH --output=./output_cleaning_%A_%a.o
#SBATCH --array=0-11
#SBATCH --ntasks=1
#SBATCH --partition=pibu_el8


# Directory with raw RNA-seq paired-end reads
READS_DIR=$(realpath ./data/01_reads_copy)
# Directory for FastP output    
mkdir -p ./analysis/02_cleaning_results
OUTPUT_DIR=$(realpath ./analysis/02_cleaning_results)
# Samples
SAMPLES=("HER21" "HER22" "HER23" "NonTNBC1" "NonTNBC2" "NonTNBC3" "Normal1" "Normal2" "Normal3" "TNBC1" "TNBC2" "TNBC3")
# apptainer paths
APPTAINER_fastp=/containers/apptainer/fastp_0.23.2--h5f740d0_3.sif
APPTAINER_fastqc=/containers/apptainer/fastqc-0.12.1.sif

# Get the sample for the current job array index
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
READ1="${READS_DIR}/${SAMPLE}_R1.fastq.gz"
READ2="${READS_DIR}/${SAMPLE}_R2.fastq.gz"

# Run fastp for cleaning with following settings:
# detect_adapter_for_pe: Automatically detects and removes adapter sequences for paired-end reads
# low_complexity_filter: Filters out reads with low complexity, such as repetitive genomic regions, Sequencing artifacts or errors, or Adapter contamination or other technical issues
# complexity_threshold: Sets the threshold for low-complexity filtering, default 30% complexity
# overrepresentation_analysis: Analyzes and reports overrepresented sequences, such as PCR duplicates or contaminants
# qualified_quality_phred: Sets the Phred quality score threshold for base qualification (bases with quality < 15 are considered unqualified)
# cut_tail -q 20: Trims low-quality bases from the tail of reads (quality < 20).
# unqualified_percent_limit 40: Discards reads if more than 40% of their bases are unqualified (quality < 15).
# length_required 50: Discards reads shorter than 50 bases after trimming -> minimizes ambiguity in mapping
# n_base_limit 1: Discards reads with more than 1 ambiguous base (N)
# dedup: Removes duplicate reads
# correction: Enables base correction for overlapping regions in paired-end reads.
echo "Running fastp for sample: ${SAMPLE}"
apptainer exec --bind ${READS_DIR},${OUTPUT_DIR} ${APPTAINER_fastp} \
    fastp \
    -i ${READ1} -I ${READ2} \
    -o ${OUTPUT_DIR}/${SAMPLE}_cleaned_R1.fastq.gz -O ${OUTPUT_DIR}/${SAMPLE}_cleaned_R2.fastq.gz \
    -t 4 \
    --detect_adapter_for_pe \
    --low_complexity_filter \
    --complexity_threshold 30 \
    --overrepresentation_analysis \
    --qualified_quality_phred 15 \
    --cut_tail -q 20 \
    --unqualified_percent_limit 40 \
    --length_required 50 \
    --n_base_limit 1 \
    --dedup \
    --correction \
    --html ${OUTPUT_DIR}/${SAMPLE}_fastp_report.html


# Run FastQC on cleaned reads
echo "Running FastQC for cleaned reads of sample: ${SAMPLE}"
apptainer exec --bind ${OUTPUT_DIR} ${APPTAINER_fastqc} fastqc -o ${OUTPUT_DIR} ${OUTPUT_DIR}/${SAMPLE}_cleaned_R1.fastq.gz ${OUTPUT_DIR}/${SAMPLE}_cleaned_R2.fastq.gz
