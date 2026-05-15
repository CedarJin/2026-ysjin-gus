#!/bin/bash

#SBATCH --job-name=megahit_cleandata
#SBATCH --output=slurm_megahit_cleandata_%j.out
#SBATCH --error=slurm_megahit_cleandata_%j.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=3-00:00:00

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/home/jys0914/2026-ysjin-gus}"
SMK_FILE="scripts/assemble_megahit_cleandata.smk"

cd "${PROJECT_DIR}"

export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}"
mkdir -p assembly/megahit logs/megahit_cleandata

# Pre-create Megahit output dirs (fewer mkdir races on shared storage)
python - <<'PY'
import glob
import os

for r1 in glob.glob("cleandata/*/*_1.fastq"):
    sid = os.path.basename(os.path.dirname(r1))
    os.makedirs(os.path.join("assembly", "megahit", sid), exist_ok=True)
PY

# Initialize conda
module load conda/base/latest
source "$(conda info --base)/etc/profile.d/conda.sh"

conda activate assemble

if [[ -f ".venv/bin/activate" ]]; then
  source ".venv/bin/activate"
fi

echo "[$(date)] Starting Snakemake megahit cleandata job ${SLURM_JOB_ID} on ${SLURMD_NODENAME:-?}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}, MEM: ${SLURM_MEM_PER_NODE:-unknown}"

python -m snakemake \
  -s "${SMK_FILE}" \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  --config megahit_threads="${MEGAHIT_THREADS:-8}" min_contig_len="${MIN_CONTIG_LEN:-1000}" \
  --rerun-incomplete \
  --nolock \
  --printshellcmds

echo "[$(date)] Snakemake megahit cleandata job finished"
