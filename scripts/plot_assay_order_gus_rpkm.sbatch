#!/bin/bash
#SBATCH --job-name=gus-rpkm-plot
#SBATCH --account=publicgrp
#SBATCH --partition=high
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/gus-rpkm-plot-%j.out
#SBATCH --error=logs/gus-rpkm-plot-%j.err
source /etc/profile.d/modules.sh
set -euo pipefail
module load R/4.4.2
cd "$SLURM_SUBMIT_DIR"
Rscript scripts/plot_assay_order_gus_rpkm.R
