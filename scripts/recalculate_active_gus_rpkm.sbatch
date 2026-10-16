#!/bin/bash
#SBATCH --job-name=active_gus_rpkm
#SBATCH --account=publicgrp
#SBATCH --partition=high
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=results/summary/pident95_cov90/abundance_all_subject/rpkm-%j.out
#SBATCH --error=results/summary/pident95_cov90/abundance_all_subject/rpkm-%j.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
python3 scripts/recalculate_active_gus_rpkm.py
