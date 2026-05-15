#!/usr/bin/env bash
#SBATCH --job-name=sra-dl-PRJNA816918
#SBATCH --output=logs/sra-dl-%j.out
#SBATCH --error=logs/sra-dl-%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#
# Download all 18 SRA runs for PRJNA816918 using NCBI sra-tools
# (prefetch -> fasterq-dump -> pigz). Files land flat in rawdata/ as:
#     rawdata/SRR18354890_1.fastq.gz
#     rawdata/SRR18354890_2.fastq.gz
#     ...
#
# Submit from repo root:
#   sbatch scripts/download_PRJNA816918_sratools.sbatch
#
# Optional overrides:
#   sbatch --export=ALL,RAWDATA=/scratch/me/rawdata scripts/download_PRJNA816918_sratools.sbatch
#   sbatch --export=ALL,SRA_CACHE=/scratch/me/sracache scripts/download_PRJNA816918_sratools.sbatch
#   sbatch --export=ALL,ACCESSIONS=config/PRJNA816918_accessions.txt scripts/download_PRJNA816918_sratools.sbatch
#
# Notes:
# - Intermediate .sra files are cached in SRA_CACHE and removed per-run after fastq
#   extraction succeeds. To keep them, set KEEP_SRA=1.
# - prefetch can use up to ~50 GB of temp space per accession; SRA_CACHE should be on
#   a filesystem with sufficient free space (the .sra files total ~50-70 GB for 18 runs).

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"
ROOT_DIR="$(pwd)"

mkdir -p logs

RAWDATA="${RAWDATA:-${ROOT_DIR}/rawdata}"
SRA_CACHE="${SRA_CACHE:-${ROOT_DIR}/rawdata/.sracache}"
ACCESSIONS="${ACCESSIONS:-${ROOT_DIR}/config/PRJNA816918_accessions.txt}"
KEEP_SRA="${KEEP_SRA:-0}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${RAWDATA}" "${SRA_CACHE}"

if [[ ! -s "${ACCESSIONS}" ]]; then
  echo "ERROR: accession list not found: ${ACCESSIONS}" >&2
  exit 1
fi

GUS_ENV="${ROOT_DIR}/envs/gus-map"
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
  export PATH="${GUS_ENV}/bin:${PATH}"
fi

command -v prefetch     >/dev/null || { echo "ERROR: prefetch not on PATH (install sra-tools)" >&2; exit 1; }
command -v fasterq-dump >/dev/null || { echo "ERROR: fasterq-dump not on PATH (install sra-tools)" >&2; exit 1; }
command -v pigz         >/dev/null || { echo "ERROR: pigz not on PATH" >&2; exit 1; }

echo "[$(date -Is)] sra-tools download starting"
echo "[$(date -Is)] RAWDATA=${RAWDATA}"
echo "[$(date -Is)] SRA_CACHE=${SRA_CACHE}"
echo "[$(date -Is)] THREADS=${THREADS}"

n_ok=0
n_fail=0
failed=()

while IFS= read -r acc; do
  [[ -z "${acc}" || "${acc}" =~ ^# ]] && continue
  echo "[$(date -Is)] === ${acc} ==="

  fq1="${RAWDATA}/${acc}_1.fastq.gz"
  fq2="${RAWDATA}/${acc}_2.fastq.gz"
  if [[ -s "${fq1}" && -s "${fq2}" ]]; then
    echo "  [skip] ${acc} already has both fastq.gz files"
    n_ok=$((n_ok + 1))
    continue
  fi

  # 1) prefetch the .sra (writes to ${SRA_CACHE}/${acc}/${acc}.sra)
  if ! prefetch --max-size 200g -O "${SRA_CACHE}" "${acc}"; then
    echo "  [FAIL] prefetch ${acc}" >&2
    n_fail=$((n_fail + 1))
    failed+=("${acc}")
    continue
  fi

  # 2) extract paired fastq to RAWDATA (split-3 -> _1/_2/.fastq for unpaired)
  if ! fasterq-dump --split-3 --skip-technical --threads "${THREADS}" \
        --temp "${SRA_CACHE}" -O "${RAWDATA}" "${SRA_CACHE}/${acc}/${acc}.sra"; then
    echo "  [FAIL] fasterq-dump ${acc}" >&2
    n_fail=$((n_fail + 1))
    failed+=("${acc}")
    continue
  fi

  # 3) gzip with pigz (parallel)
  for f in "${RAWDATA}/${acc}_1.fastq" "${RAWDATA}/${acc}_2.fastq" "${RAWDATA}/${acc}.fastq"; do
    if [[ -s "${f}" ]]; then
      pigz -p "${THREADS}" -f "${f}"
    fi
  done

  # 4) clean .sra cache for this run unless KEEP_SRA=1
  if [[ "${KEEP_SRA}" != "1" ]]; then
    rm -rf "${SRA_CACHE}/${acc}"
  fi

  if [[ -s "${fq1}" && -s "${fq2}" ]]; then
    echo "  [ok  ] ${acc}: ${fq1##*/}, ${fq2##*/}"
    n_ok=$((n_ok + 1))
  else
    echo "  [FAIL] ${acc}: paired fastq missing after extract" >&2
    n_fail=$((n_fail + 1))
    failed+=("${acc}")
  fi
done < "${ACCESSIONS}"

echo "[$(date -Is)] DONE. ok=${n_ok}  failed=${n_fail}"
if (( n_fail > 0 )); then
  echo "Failed accessions: ${failed[*]}" >&2
  exit 1
fi
