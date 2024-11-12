#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=download_genome
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_downloadgenome_%A_%a.e


#-----------------Get Ref Genome form Ensemble ftp site------------------------------------
# Specise of interest -> breast cancer in human specimens; homo sapiens, status 12Nov2024: latest build -> 113
SPECIES="homo_sapiens"
ASSEMBLY="GRCh38"
BUILD="113"
# URL for reference genome
URL_REFERENCE_GENOME="ftp://ftp.ensembl.org/pub/release-${BUILD}/fasta/${SPECIES}/dna/"
URL_ANNOTATION="ftp://ftp.ensembl.org/pub/release-${BUILD}/gtf/${SPECIES}/"
# files
REFERENCE_GENOME_FILE="${SPECIES^}.${ASSEMBLY}.dna.primary_assembly.fa.gz"
ANNOTATION_FILE="${SPECIES^}.${ASSEMBLY}.${BUILD}.gtf.gz"
CHECKSUMS_FILE="CHECKSUMS"

# Directory for downloaded files
DOWNLOAD_DIR="/data/users/aboss/rna_course/reference_genome" 
mkdir -p "${DOWNLOAD_DIR}"
# Move to the download directory
cd "${DOWNLOAD_DIR}"

# Download reference genome sequence and annotation
wget "${URL_REFERENCE_GENOME}${REFERENCE_GENOME_FILE}"
wget "${URL_ANNOTATION}${ANNOTATION_FILE}"

# Download the CHECKSUMS files to verify the downloaded files
echo "Downloading CHECKSUMS file..."
wget -O CHECKSUMS_REF_GENOME "${URL_REFERENCE_GENOME}${CHECKSUMS_FILE}"
wget -O CHECKSUMS_ANNOTATION_GENOME "${URL_ANNOTATION}${CHECKSUMS_FILE}"

###-------------Verification with CHECKSUM file------------------
##not working

# Verify the integrity of the downloaded genome sequence file
echo "Verifying genome file checksum..."
EXPECTED_GENOME_CHECKSUM=$(grep "$REFERENCE_GENOME_FILE" CHECKSUMS_REF_GENOME | awk '{print $1}')
ACTUAL_GENOME_CHECKSUM=$(md5sum "$REFERENCE_GENOME_FILE" | awk '{print $1}')

if [[ "$EXPECTED_GENOME_CHECKSUM" == "$ACTUAL_GENOME_CHECKSUM" ]]; then
    echo "Genome file checksum verified successfully."
else
    echo "Error: Genome file checksum does not match!"
    exit 1
fi

# Verify the integrity of the downloaded annotation file
echo "Verifying annotation file checksum..."
EXPECTED_ANNOTATION_CHECKSUM=$(grep "$ANNOTATION_FILE" CHECKSUMS_ANNOTATION_GENOME | awk '{print $1}')
ACTUAL_ANNOTATION_CHECKSUM=$(md5sum "$ANNOTATION_FILE" | awk '{print $1}')

if [[ "$EXPECTED_ANNOTATION_CHECKSUM" == "$ACTUAL_ANNOTATION_CHECKSUM" ]]; then
    echo "Annotation file checksum verified successfully."
else
    echo "Error: Annotation file checksum does not match!"
    exit 1
fi

echo "All files downloaded and verified successfully."