#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastq
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_fastq_%A_%a.e


mkdir /data/users/aboss/rna_course/reads_point
cd /data/users/aboss/rna_course/reads_point
FASTQC_OUTPUT=/data/users/aboss/rna_course/fastq_output
SAMPLE='HER21'

#reads
ln -s /data/users/aboss/rna_course/reads${SAMPLE}_R1.fastq.gz ${SAMPLE}_R1.fastq.gz
ln -s /data/users/aboss/rna_course/reads${SAMPLE}_R2.fastq.gz ${SAMPLE}_R2.fastq.gz



apptainer exec --bind ${SAMPLE}_R1.fastq.gz /containers/apptainer/fastqc-0.12.1.sif fastqc -t -o ${FASTQC_OUTPUT}








                          