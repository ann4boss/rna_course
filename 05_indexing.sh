#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --time=03:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=indexing
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_indexing_%A_%a.e


# Directory with reference genome
REFERENCE_GENOME_DIR="/data/users/aboss/rna_course/reference_genome"
# output directory for indexed ref genome
HISAT2_INDEX="${REFERENCE_GENOME_DIR}/indexing/Homo_sapiens.GRCh38_indexed"
mkdir -p ${REFERENCE_GENOME_DIR}/indexing/
# apptainer paths
APPTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"


###------------file preparation for HISAT2 - indexing---------------------
# input fasta file
GENOME_FA="${REFERENCE_GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

# indexing, using --verbose to get updates
apptainer exec --bind ${REFERENCE_GENOME_DIR} ${APPTAINER} hisat2-build --verbose ${GENOME_FA} ${HISAT2_INDEX}
