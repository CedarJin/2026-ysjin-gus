#!/usr/bin/env bash
#SBATCH --job-name=gus-smk
#SBATCH --output=logs/snakemake-%j.out
#SBATCH --error=logs/snakemake-%j.err
#SBATCH --time=24:00:00
# SBATCH CPU count must be ≤ your QOS "MaxCpuPerJob" (e.g. publicgrp+high was 8 CPUs/job on this cluster).
# If you see: QOSMaxCpuPerJobLimit — lower this or try another partition/QOS (e.g. low).
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
## cluster-specific (uncomment / edit for your site):
##SBATCH --account=publicgrp
##SBATCH --partition=high
#
# Run the Snakemake workflow (Snakefile in repo root).
#
# From repository root:
#   mkdir -p logs
#   sbatch scripts/submit_snakemake.sbatch
#
# More CPUs only if your account/QOS allows (check `sacctmgr show assoc ...` or cluster docs):
#   sbatch --cpus-per-task=8 scripts/submit_snakemake.sbatch   # safe default for many QOS
#   sbatch --partition=low --cpus-per-task=32 scripts/submit_snakemake.sbatch  # example only if permitted
#
# Extra Snakemake flags (dry run, touch, etc.):
#   sbatch --export=ALL,SNAKEMAKE_EXTRA='-n' scripts/submit_snakemake.sbatch
#
# Locks: if a previous Snakemake run was killed (OOM, time limit), the working directory
# can stay locked. This script runs `snakemake --unlock` by default before the workflow.
# Skip that only if you know no stale lock exists (e.g. another Snakemake is NOT running):
#   sbatch --export=ALL,SNAKEMAKE_SKIP_UNLOCK=1 scripts/submit_snakemake.sbatch

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"
ROOT_DIR="$(pwd)"

mkdir -p logs

GUS_ENV="${ROOT_DIR}/envs/gus-map"

# Non-interactive batch jobs do not run `conda init`; `conda activate` only works after
# loading conda’s shell hook (and often `module load conda` on HPC).
module load conda/latest 2>/dev/null || true

if command -v conda >/dev/null 2>&1; then
  _conda_base="$(conda info --base 2>/dev/null || true)"
  if [[ -n "${_conda_base}" && -f "${_conda_base}/etc/profile.d/conda.sh" ]]; then
    # shellcheck disable=SC1091
    source "${_conda_base}/etc/profile.d/conda.sh"
  else
    eval "$(conda shell.bash hook 2>/dev/null)" || true
  fi
  if ! conda activate "${GUS_ENV}"; then
    echo "WARN: conda activate ${GUS_ENV} failed; prepending env bin to PATH." >&2
    export PATH="${GUS_ENV}/bin:${PATH}"
  fi
else
  echo "WARN: conda not on PATH; prepending ${GUS_ENV}/bin (no activate)." >&2
  export PATH="${GUS_ENV}/bin:${PATH}"
fi

# Use cluster modules only if a binary is still missing from the env
command -v prodigal >/dev/null || module load prodigal/2.6.3 2>/dev/null || true
command -v diamond >/dev/null || module load diamond/2.1.7 2>/dev/null || true

CORES="${SLURM_CPUS_PER_TASK:-8}"
CONFIG="${ROOT_DIR}/config/config.yaml"

command -v snakemake >/dev/null || {
  echo "ERROR: snakemake not on PATH. Install in gus-map env, e.g. conda env update -f environment.yml" >&2
  exit 1
}

if [[ "${SNAKEMAKE_SKIP_UNLOCK:-0}" != "1" ]]; then
  echo "Removing stale Snakemake lock (if any) in ${ROOT_DIR}/.snakemake ..."
  snakemake --directory "${ROOT_DIR}" --unlock
fi

# Merge CLI --config over config.yaml so DIAMOND threads track the allocation
# shellcheck disable=SC2086
exec snakemake \
  --directory "${ROOT_DIR}" \
  -c "${CORES}" \
  --configfile "${CONFIG}" \
  --config "threads_diamond=${CORES}" \
  --rerun-incomplete \
  --printshellcmds \
  ${SNAKEMAKE_EXTRA:-}
