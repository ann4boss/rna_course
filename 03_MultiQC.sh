#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=multiqc
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=./error_multiqc_%A_%a.e
#SBATCH --output=./output_multiqc_%A_%a.o

# Apptainer paths
APPTAINER_multiqc=/containers/apptainer/multiqc-1.19.sif

# Check if at least 3 arguments are passed (FILE_EXTENSION, at least one input directory, and output directory)
if [ $# -lt 3 ]; then
    echo "Usage: $0 <FILE_EXTENSION> <input_directory_1> [input_directory_2 ...] <output_directory>"
    exit 1
fi

# Extract the file extension (first argument)
FILE_EXTENSION=$1
# Remove the file extension argument from the list
shift

# Extract input directories (all arguments except the last one and the file extension)
INPUT_DIRS=()
while [[ $# -gt 1 ]]; do
    INPUT_DIRS+=("$1")
    shift
done

# Extract the output directory
OUTPUT_DIR=$1
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")
mkdir -p ${OUTPUT_DIR}

# Prepare the input directories for MultiQC
INPUT_ARGS=""
for DIR in "${INPUT_DIRS[@]}"; do
    if [ ! -d "$DIR" ]; then
        echo "Error: Directory $DIR does not exist."
        exit 1
    fi
    INPUT_ARGS+=" $DIR/$FILE_EXTENSION"
done

#-------------------------MultiQC-----------------------------------------------------
# Run MultiQC
echo "Running MultiQC..."
apptainer exec ${APPTAINER_multiqc} multiqc ${INPUT_ARGS} -o ${OUTPUT_DIR}

echo "MultiQC analysis complete. Results saved to $OUTPUT_DIR"