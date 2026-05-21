#!/bin/bash
#SBATCH --job-name=gus-abund
#SBATCH --output=logs/gus_abundance-%j.out
#SBATCH --error=logs/gus_abundance-%j.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --partition=low
#SBATCH --account=publicgrp

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"
ROOT_DIR="$(pwd)"
mkdir -p logs

GUS_ENV="${ROOT_DIR}/envs/gus-map"
SMK_FILE="${ROOT_DIR}/scripts/gus_abundance.smk"
CONFIG="${ROOT_DIR}/config/gus_abundance.yaml"

# Initialize conda (same pattern as submit_assemble_megahit_cleandata.sbatch)
module load conda/base/latest
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

conda activate "${GUS_ENV}"

python -m snakemake --version >/dev/null || {
  echo "ERROR: snakemake not importable in gus-map env." >&2
  exit 1
}

if [[ "${SNAKEMAKE_SKIP_UNLOCK:-0}" != "1" ]]; then
  python -m snakemake --directory "${ROOT_DIR}" --unlock
fi

exec python -m snakemake \
  --directory "${ROOT_DIR}" \
  -s "${SMK_FILE}" \
  -c "${SLURM_CPUS_PER_TASK:-8}" \
  --configfile "${CONFIG}" \
  --rerun-incomplete \
  --printshellcmds \
  ${SNAKEMAKE_EXTRA:-}
