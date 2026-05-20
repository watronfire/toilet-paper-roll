#!/bin/bash
#SBATCH --job-name=upload_terrrra
#SBATCH --output=logs/upload_%j.log
#SBATCH --error=logs/upload_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --partition=normal  # Adjust to your cluster's partition name

# 1. Load the seqtk module (if available) or ensure it's in your PATH
module load gcloud

cd /data/wohllab/2025.07.10_fp/

# 2. Path to your text file containing the list of .fastq.gz files
FILE_LIST="to_upload.txt"

cat $FILE_List | gsutil -m cp -I gs://fc-c5000c68-c6db-4966-9d2f-5fcd712c96d8/test/
