#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=download_genome
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=./error_downloadgenome_%A_%a.e
#SBATCH --output=./output_downloadgenome_%A_%a.o

#-----------------Get Ref Genome from Ensembl ftp site------------------------------------
# Species of interest -> breast cancer in human specimens; homo sapiens, status 12Nov2024: latest build -> 113
SPECIES="homo_sapiens"
ASSEMBLY="GRCh38"
BUILD="113"
# URL for reference genome
URL_REFERENCE_GENOME="ftp://ftp.ensembl.org/pub/release-${BUILD}/fasta/${SPECIES}/dna/"
URL_ANNOTATION="ftp://ftp.ensembl.org/pub/release-${BUILD}/gtf/${SPECIES}/"
# Files
REFERENCE_GENOME_FILE="${SPECIES^}.${ASSEMBLY}.dna.primary_assembly.fa.gz"
ANNOTATION_FILE="${SPECIES^}.${ASSEMBLY}.${BUILD}.gtf.gz"
CHECKSUMS_FILE="CHECKSUMS"

# Directory for downloaded files
DOWNLOAD_DIR=./data/02_reference_genome
mkdir -p "${DOWNLOAD_DIR}"

# Save the download date
DATE=$(date '+%Y-%m-%d %H:%M:%S')
echo "Download date: ${DATE}" > "${DOWNLOAD_DIR}/reference_genome_log.txt"

# Move to the download directory
cd "${DOWNLOAD_DIR}"

# Download reference genome sequence and annotation
wget "${URL_REFERENCE_GENOME}${REFERENCE_GENOME_FILE}"
wget "${URL_ANNOTATION}${ANNOTATION_FILE}"

# Download the CHECKSUMS files to verify the downloaded files
wget -O CHECKSUMS_REF_GENOME "${URL_REFERENCE_GENOME}${CHECKSUMS_FILE}"
wget -O CHECKSUMS_ANNOTATION_GENOME "${URL_ANNOTATION}${CHECKSUMS_FILE}"


#---------------Verification of Dowload with CHECKSUM-----------------
# Calculate the sum of the reference genome file
CALCULATED_CHECKSUM=$(sum "${REFERENCE_GENOME_FILE}" | awk '{print $1}')
# Extract the expected checksum from the CHECKSUMS_REF_GENOME file
EXPECTED_CHECKSUM=$(grep "${REFERENCE_GENOME_FILE}" CHECKSUMS_REF_GENOME | awk '{print $1}')

# Compare the calculated checksum with the expected checksum
if [ "$CALCULATED_CHECKSUM" != "$EXPECTED_CHECKSUM" ]; then
    echo "Checksum verification failed for reference genome file."
    exit 1
fi
echo "Reference genome checksum verification passed." >> reference_genome_log.txt

# Verify the integrity of the downloaded annotation genome using the CHECKSUM file
# Calculate the checksum of the annotation file
CALCULATED_CHECKSUM=$(sum "${ANNOTATION_FILE}" | awk '{print $1}')
# Extract the expected checksum from the CHECKSUMS_ANNOTATION_FILE file
EXPECTED_CHECKSUM=$(grep "${ANNOTATION_FILE}" CHECKSUMS_ANNOTATION_GENOME | awk '{print $1}')
# Compare the calculated checksum with the expected checksum
if [ "$CALCULATED_CHECKSUM" != "$EXPECTED_CHECKSUM" ]; then
    echo "Checksum verification failed for annotation file."
    exit 1
fi
echo "Annotation genome checksum verification passed." >> reference_genome_log.txt