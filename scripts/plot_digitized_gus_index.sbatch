#!/bin/bash
#SBATCH --job-name=gus-index-figure
#SBATCH --account=publicgrp
#SBATCH --partition=high
#SBATCH --time=00:03:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/gus-index-figure-%j.out
#SBATCH --error=logs/gus-index-figure-%j.err
source /etc/profile.d/modules.sh
set -euo pipefail
module load R/4.4.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"
Rscript scripts/plot_digitized_gus_index.R
