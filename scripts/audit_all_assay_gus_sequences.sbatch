#!/bin/bash
#SBATCH --job-name=gus-cohort-sequences
#SBATCH --account=publicgrp
#SBATCH --partition=high
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/gus-cohort-sequences-%j.out
#SBATCH --error=logs/gus-cohort-sequences-%j.err
source /etc/profile.d/modules.sh
set -euo pipefail
module load mafft/7.505
cd "$SLURM_SUBMIT_DIR"
/usr/bin/python3 scripts/audit_all_assay_gus_sequences.py
