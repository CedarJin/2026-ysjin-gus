#!/usr/bin/env bash
#SBATCH --job-name=gus-map
#SBATCH --output=logs/gus-map-%j.out
#SBATCH --error=logs/gus-map-%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#
# Run from repo root after: mkdir -p logs
#
# This batch file runs the monolithic bash driver. For Snakemake (resume per rule), after
# activating the same env + modules, run:
#   snakemake -c "${SLURM_CPUS_PER_TASK:-8}" --configfile config/config.yaml
#
# Override bash pipeline threads / min length for this sbatch script:
#   sbatch --export=ALL,THREADS=32,MIN_CONTIG_LENGTH=1000 scripts/run_gus_mapping.sbatch

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"
ROOT_DIR="$(pwd)"

mkdir -p logs

module load prodigal/2.6.3 diamond/2.1.7 2>/dev/null || true
# PDF extract needs pypdf; use project env if present:
if [[ -x "${ROOT_DIR}/envs/gus-map/bin/python" ]]; then
  export PATH="${ROOT_DIR}/envs/gus-map/bin:${PATH}"
fi

export THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-8}}"
export MIN_CONTIG_LENGTH="${MIN_CONTIG_LENGTH:-1000}"

bash "${ROOT_DIR}/scripts/run_gus_mapping.sh"
