#!/bin/bash
#SBATCH --job-name=filter_snakemake
#SBATCH --output=logs/filter_%j.log
#SBATCH --error=logs/filter_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=72:00:00
#SBATCH --partition=normal  # Adjust to your cluster's partition name

# 1. Load the seqtk module (if available) or ensure it's in your PATH
source /PHShome/nm104/mambaforge/etc/profile.d/conda.sh
conda activate filter
cd /data/wohllab/2025.07.10_fp/


TARGET_FILE="filter.complete"
SNAKE_CMD="snakemake --profile profiles/slurm"

echo "Starting Snakemake orchestration loop at $(date)"

# Loop until the target file exists
until [ -f "$TARGET_FILE" ]; do
    echo "Running Snakemake..."
    $SNAKE_CMD

    # Check if Snakemake succeeded or failed
    if [ $? -eq 0 ] && [ -f "$TARGET_FILE" ]; then
        echo "Pipeline finished successfully."
        break
    else
        echo "Snakemake exited with an error or file not found. Restarting in 30 seconds..."
        sleep 30
    fi
done

echo "Workflow complete at $(date)"
