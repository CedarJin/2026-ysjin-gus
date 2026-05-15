#!/bin/bash

#SBATCH --job-name=rclone_raw_fastqgz
#SBATCH --output=logs/rclone_raw_fastqgz_%j.out
#SBATCH --error=logs/rclone_raw_fastqgz_%j.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL

set -euo pipefail

if command -v module >/dev/null 2>&1; then
  module load rclone
fi

# ── edit if needed ────────────────────────────────────────────────────
PROJECT_DIR="${PROJECT_DIR:-/home/jys0914/2026-ysjin-gus}"
SRC="${PROJECT_DIR}/rawdata"
REMOTE_ROOT="${RCLONE_REMOTE_ROOT:-mybox:/2026-ysjin-gus-backup}"
DEST="${REMOTE_ROOT%/}/rawdata"

mkdir -p "${PROJECT_DIR}/logs"

if [[ ! -d "${SRC}" ]]; then
  echo "ERROR: source dir missing: ${SRC}" >&2
  exit 1
fi

echo "[$(date)] SRC=${SRC}"
echo "[$(date)] DEST=${DEST}"

if ! command -v rclone >/dev/null 2>&1; then
  echo "ERROR: rclone not in PATH (module load rclone / conda install rclone?)" >&2
  exit 1
fi

rclone copy "${SRC}" "${DEST}" \
  --include '*.fastq.gz' \
  --stats 1m \
  --stats-log-level NOTICE \
  --transfers 8 \
  --checkers 16 \
  --buffer-size 64M

echo "[$(date)] DONE"
