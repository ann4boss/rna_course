#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=multiqc
#SBATCH --mail-user=anna.boss@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --error=/data/users/aboss/rna_course/error_multiqc_%A_%a.e


# Check if at least one input path and output path are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_directory_1> [input_directory_2 ...] <output_directory>"
    exit 1
fi

# Extract the output directory (the last argument)
OUTPUT_DIR=${@: -1}
# Ensure the output directory is an absolute path
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")
mkdir -p ${OUTPUT_DIR}

# Extract input directories (all arguments except the last one)
INPUT_DIRS=("${@:1:$#-1}")

# Ensure all input directories are absolute paths
for i in "${!INPUT_DIRS[@]}"; do
    INPUT_DIRS[$i]=$(realpath "${INPUT_DIRS[$i]}")
done

# Create the input directories if they don't exist
for INPUT_DIR in ${INPUT_DIRS[@]}; do
    if [ ! -d "$INPUT_DIR" ]; then
        echo "Warning: Input directory $INPUT_DIR does not exist. Skipping."
    else
        mkdir -p "$INPUT_DIR"
    fi
done

# Apptainer paths
APPTAINER_multiqc="/containers/apptainer/multiqc-1.19.sif"

#-------------------------MultiQC-----------------------------------------------------
# Debugging: Print input directories and files
echo "Input directories: ${INPUT_DIRS[@]}"
for DIR in "${INPUT_DIRS[@]}"; do
    echo "Listing files in ${DIR}:"
    ls -l ${DIR}
done

# Run MultiQC for all input directories together and generate a single report
INPUT_PATHS="${INPUT_DIRS[@]}"  # Concatenate all input paths into a single string
echo "${INPUT_PATHS}"

echo "Running MultiQC for all input directories..."
apptainer exec --bind ${INPUT_PATHS} ${APPTAINER_multiqc} multiqc ${INPUT_PATHS} -o ${OUTPUT_DIR}

echo "MultiQC analysis complete. Results are in ${OUTPUT_DIR}"