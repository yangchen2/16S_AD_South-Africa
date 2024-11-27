#!/bin/bash -l
# SLURM job script to run Snakemake
# Author: Yang Chen
# Date: 11-27-2024
# Description: Executes Xebec Snakemake workflow using SLURM on 16S AD South Africa dataset.

# SLURM directives
#SBATCH --job-name=smk-xebec              # Job name
#SBATCH --output=logs/smk-xebec_%j.out    # Standard output log
#SBATCH --error=logs/smk-xebec_%j.err     # Standard error log
#SBATCH --time=1:00:00                    # Maximum runtime (HH:MM:SS)
#SBATCH --nodes=1                         # Number of nodes
#SBATCH --ntasks=1                        # Number of tasks (keep as 1 for Snakemake)
#SBATCH --cpus-per-task=4                 # Number of cores per task
#SBATCH --mem=16G                         # Memory allocation
#SBATCH --partition=short                 # Partition/queue name
#SBATCH --mail-type=END,FAIL              # Send email on job completion or failure
#SBATCH --mail-user=yac027@ucsd.edu       # Email address for notifications

# Activate conda environment
echo 'Activating xebec_python39 environment'
source activate xebec_python39

# Run Snakemake
echo 'Running Snakemake'
snakemake --use-conda --cores $SLURM_CPUS_PER_TASK

# Job completion message
echo 'Xebec workflow completed'

