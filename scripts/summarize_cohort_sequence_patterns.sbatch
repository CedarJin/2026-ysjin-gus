#!/bin/bash
#SBATCH --job-name=gus-sequence-patterns
#SBATCH --account=publicgrp
#SBATCH --partition=high
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/gus-sequence-patterns-%j.out
#SBATCH --error=logs/gus-sequence-patterns-%j.err
source /etc/profile.d/modules.sh
set -euo pipefail
module load R/4.4.2
cd "$SLURM_SUBMIT_DIR"
Rscript scripts/summarize_cohort_sequence_patterns.R
