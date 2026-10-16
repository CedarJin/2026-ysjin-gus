#!/bin/bash
#SBATCH --job-name=gus-seq-compare
#SBATCH --account=publicgrp
#SBATCH --partition=high
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/gus-seq-compare-%j.out
#SBATCH --error=logs/gus-seq-compare-%j.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
/usr/bin/python3 scripts/compare_reviewer_gus_sequences.py
