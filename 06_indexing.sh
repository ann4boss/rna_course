#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --time=03:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=indexing
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=./error_indexing_%A_%a.e
#SBATCH --output=./output_indexing_%A_%a.o


# Directory with reference genome
REFERENCE_GENOME_DIR=$(realpath ./data/02_reference_genome)
# output directory for indexed ref genome
BASENAME='Homo_sapiens.GRCh38_indexed'
HISAT2_INDEX=./data/03_indexing/${BASENAME}
mkdir -p ./data/03_indexing
# Apptainer paths
APPTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif


###------------file preparation for HISAT2 - indexing---------------------
# Unzip files
# Paths for unzipped files
GENOME_FA=${REFERENCE_GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa
ANNOTATION_GTF=${REFERENCE_GENOME_DIR}/Homo_sapiens.GRCh38.113.gtf

# Unzip the reference genome and annotation file if not already unzipped
if [ ! -f "$GENOME_FA" ]; then
    gunzip -c ${REFERENCE_GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > "$GENOME_FA"
fi
if [ ! -f "$ANNOTATION_GTF" ]; then
    gunzip -c ${REFERENCE_GENOME_DIR}/Homo_sapiens.GRCh38.113.gtf.gz > "$ANNOTATION_GTF"
fi


# input unzipped fasta file
GENOME_FA=${REFERENCE_GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# indexing, using --verbose to get updates
apptainer exec --bind ${REFERENCE_GENOME_DIR}:/mnt ${APPTAINER} hisat2-build --verbose ${GENOME_FA} ${HISAT2_INDEX}
