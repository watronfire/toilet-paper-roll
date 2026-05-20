#!/bin/bash
#SBATCH --job-name=seqtk_subsample
#SBATCH --output=logs/subsample_%j.log
#SBATCH --error=logs/subsample_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --partition=short  # Adjust to your cluster's partition name

# 1. Load the seqtk module (if available) or ensure it's in your PATH
source /PHShome/nm104/mambaforge/etc/profile.d/conda.sh
conda activate vibecheck
cd /data/wohllab/2025.07.10_fp/

# 2. Path to your text file containing the list of .fastq.gz files
FILE_LIST="subsample_files.txt"
OUT_DIR="./input/terra_upload/"

# 3. Target number of reads
READ_COUNT=500000

# 4. Loop through each file in the list
while read -r FASTQ; do

    FILENAME=$(basename "$FASTQ")

    # Define the output filename (e.g., sample_sub.fastq.gz)
    OUT_FILE="${OUT_DIR}/${FILENAME%.fastq.gz}_sub.fastq.gz"

    echo "Processing $FASTQ to $OUT_FILE..."

    # seqtk sample -s100 uses a seed of 100 for reproducibility
    # We pipe to gzip to keep the output compressed
    seqtk sample -s100 "$FASTQ" "$READ_COUNT" | gzip > "$OUT_FILE"

done < "$FILE_LIST"

echo "Subsampling complete."
