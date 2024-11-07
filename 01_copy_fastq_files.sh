#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastq
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_fastq_%A_%a.e


mkdir /data/users/aboss/rna_course/reads
OUTPUT_FOLDER=/data/users/aboss/rna_course/reads

#copy all reads to my repository
cp -r /data/courses/rnaseq_course/breastcancer_de/reads/ $OUTPUT_FOLDER