#!/bin/bash
#SBATCH --job-name=gus-site-conservation
#SBATCH --account=publicgrp
#SBATCH --partition=high
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/gus-site-conservation-%j.out
#SBATCH --error=logs/gus-site-conservation-%j.err
source /etc/profile.d/modules.sh
set -euo pipefail
module load mafft/7.505
cd "$SLURM_SUBMIT_DIR"
/usr/bin/python3 scripts/analyze_reviewer_gus_sites.py
