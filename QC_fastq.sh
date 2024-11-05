#!/bin bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastq
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/aboss/slurm_tutorial/output_alignment_%A_%a.o
#SBATCH --error=/data/users/aboss/slurm_tutorial/error_alignment_%A_%a.e

READS_DIR=/data/courses/rnaseq_course/breastcancer_de
FASTQ_OUTPUT=/data/users/aboss/slurm_tutorial/ecoli/alignment

mkdir -p $ALIGN_DIR

cd $READS_DIR

FILES=(*.fastq)

apptainer exec \
--bind /data/courses \
/data/courses/HPC_tutorial/containers/minimap2.sif \
minimap2 \
-a \
-x sr \
$REFERENCE_DIR/ecoli-strK12-MG1655.fasta \
${FILES[$SLURM_ARRAY_TASK_ID]} \
> $ALIGN_DIR/${FILES[$SLURM_ARRAY_TASK_ID]}.sam 