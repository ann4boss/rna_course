#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1000MG
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastq
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_fastq_%A_%a.e

#create variables and links 
READS_DIR=/data/courses/rnaseq_course/breastcancer_de
FASTQ_OUTPUT=/data/users/aboss/rna_course
group='HER2'

mkdir /data/users/aboss/fastq_output


#reads
ln -s /data/courses/rnaseq_course/breastcancer_de/reads/${group}1_R1.fastq.gz ${group}1_R1.fastq.gz
ln -s /data/courses/rnaseq_course/breastcancer_de/reads/${group}1_R2.fastq.gz ${group}1_R2.fastq.gz

for k in `ls -1 ${group}*.fastq.gz`;
do fastqc -t 2 ${k};
done


FILES=(*.fastq)

apptainer exec \
--bind /data/courses \
/containers/apptainer/fastqc-0.12.1.siff \
fastq \
-a \
-x sr \
$REFERENCE_DIR/ecoli-strK12-MG1655.fasta \
${FILES[$SLURM_ARRAY_TASK_ID]} \
> $ALIGN_DIR/${FILES[$SLURM_ARRAY_TASK_ID]}.sam 






                          