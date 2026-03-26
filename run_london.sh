#!/bin/bash
#SBATCH --job-name=filter_paper
#SBATCH --time=24:0:0
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --partition=normal
#SBATCH -o out.reads
#SBATCH -e err.reads

## This is a comment
source /PHShome/nm104/mambaforge/etc/profile.d/conda.sh
conda activate bacpage

cd /data/wohllab/2025.07.10_fp

FILE=results/variants_subsample.csv
while [ ! -f $FILE ]
do
    snakemake --rerun-incomplete --keep-going --profile profiles/slurm/ $FILE
done

#snakemake --rerun-incomplete --keep-going --profile profiles/slurm/ intermediates/illumina/assembly/B1.assembly.fasta
