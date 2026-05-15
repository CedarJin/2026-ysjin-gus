#!/usr/bin/env bash
#SBATCH --job-name=cleandata
#SBATCH --output=logs/cleandata-%j.out
#SBATCH --error=logs/cleandata-%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=low
#SBATCH --account=publicgrp

#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL
#
# Run scripts/cleandata.smk on all 18 PE samples in rawdata/.
#
# Pipeline:
#   FastQC -> MultiQC -> fastp -> FastQC -> MultiQC -> Bowtie2 host removal (T2T)
#
# Submit from repo root:
#   sbatch scripts/submit_cleandata.sbatch
#
# Useful overrides:
#   # Dry-run / print DAG only
#   sbatch --export=ALL,SNAKEMAKE_EXTRA='-n -p' scripts/submit_cleandata.sbatch
#
#   # Different conda env name (defaults to `qc`)
#   sbatch --export=ALL,QC_ENV=my-qc-env scripts/submit_cleandata.sbatch
#
#   # Skip Snakemake --unlock at startup (only if you know no other run is active)
#   sbatch --export=ALL,SNAKEMAKE_SKIP_UNLOCK=1 scripts/submit_cleandata.sbatch

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"
ROOT_DIR="$(pwd)"
mkdir -p logs

QC_ENV="${QC_ENV:-qc}"
SMK_FILE="${SMK_FILE:-${ROOT_DIR}/scripts/cleandata.smk}"

# ---- Activate the `qc` conda environment ------------------------------------
# Non-interactive batch jobs do not run `conda init`; load the shell hook first.
module load conda/latest 2>/dev/null || true

if command -v conda >/dev/null 2>&1; then
  _conda_base="$(conda info --base 2>/dev/null || true)"
  if [[ -n "${_conda_base}" && -f "${_conda_base}/etc/profile.d/conda.sh" ]]; then
    # shellcheck disable=SC1091
    source "${_conda_base}/etc/profile.d/conda.sh"
  else
    eval "$(conda shell.bash hook 2>/dev/null)" || true
  fi
  if ! conda activate "${QC_ENV}"; then
    echo "ERROR: conda activate ${QC_ENV} failed." >&2
    echo "       Make sure the env exists: 'conda env list'" >&2
    exit 1
  fi
else
  echo "ERROR: conda not on PATH." >&2
  exit 1
fi

echo "[$(date -Is)] Using conda env: ${CONDA_PREFIX:-${QC_ENV}}"

# snakemake must be on PATH (the rest of the tools will be resolved by the
# rules at run time -- whichever way the qc env provides them).
if ! command -v snakemake >/dev/null 2>&1; then
  echo "ERROR: snakemake not on PATH in env '${QC_ENV}'." >&2
  exit 1
fi

CORES="${SLURM_CPUS_PER_TASK:-8}"

if [[ ! -s "${SMK_FILE}" ]]; then
  echo "ERROR: Snakefile not found: ${SMK_FILE}" >&2
  exit 1
fi

if [[ "${SNAKEMAKE_SKIP_UNLOCK:-0}" != "1" ]]; then
  echo "[$(date -Is)] snakemake --unlock"
  snakemake --directory "${ROOT_DIR}" -s "${SMK_FILE}" --unlock || true
fi

echo "[$(date -Is)] starting Snakemake (cores=${CORES})"
# shellcheck disable=SC2086
exec snakemake \
  --directory "${ROOT_DIR}" \
  -s "${SMK_FILE}" \
  -c "${CORES}" \
  --rerun-incomplete \
  --printshellcmds \
  ${SNAKEMAKE_EXTRA:-}
