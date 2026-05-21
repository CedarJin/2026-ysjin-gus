#!/bin/bash
#SBATCH --job-name=gus-abund-1
#SBATCH --output=logs/gus_abundance_1sample-%j.out
#SBATCH --error=logs/gus_abundance_1sample-%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=high
#SBATCH --account=publicgrp

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL

# Run gus_abundance.smk for a single sample (smoke test before full 18-sample run).
#
# From repository root:
#   sbatch scripts/submit_gus_abundance_one_sample.sbatch
#
# Pick another sample:
#   sbatch --export=ALL,SAMPLE=SRR18354894 scripts/submit_gus_abundance_one_sample.sbatch
#
# Dry run:
#   sbatch --export=ALL,SNAKEMAKE_EXTRA='-n' scripts/submit_gus_abundance_one_sample.sbatch

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"
ROOT_DIR="$(pwd)"
mkdir -p logs

# Override at submit time: sbatch --export=ALL,SAMPLE=SRR18354894 ...
SAMPLE="${SAMPLE:-SRR18354890}"

GUS_ENV="${ROOT_DIR}/envs/gus-map"
SMK_FILE="${ROOT_DIR}/scripts/gus_abundance.smk"
CONFIG="${ROOT_DIR}/config/gus_abundance.yaml"

for f in \
  "${ROOT_DIR}/cleandata/${SAMPLE}/${SAMPLE}_1.fastq" \
  "${ROOT_DIR}/cleandata/${SAMPLE}/${SAMPLE}_2.fastq" \
  "${ROOT_DIR}/assembly/megahit/${SAMPLE}/final.contigs.fa"
do
  [[ -f "${f}" ]] || {
    echo "ERROR: missing input for sample ${SAMPLE}: ${f}" >&2
    exit 1
  }
done

module load conda/base/latest
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${GUS_ENV}"

for cmd in snakemake bowtie2-build bowtie2 samtools featureCounts prodigal diamond; do
  command -v "${cmd}" >/dev/null || {
    echo "ERROR: ${cmd} not on PATH in gus-map env." >&2
    exit 1
  }
done

echo "[$(date -Is)] sample=${SAMPLE} cores=${SLURM_CPUS_PER_TASK:-8}"

if [[ "${SNAKEMAKE_SKIP_UNLOCK:-0}" != "1" ]]; then
  snakemake --directory "${ROOT_DIR}" --unlock
fi

exec snakemake \
  --directory "${ROOT_DIR}" \
  -s "${SMK_FILE}" \
  -c "${SLURM_CPUS_PER_TASK:-8}" \
  --configfile "${CONFIG}" \
  --config "sample=${SAMPLE}" \
  --rerun-incomplete \
  --printshellcmds \
  ${SNAKEMAKE_EXTRA:-}
