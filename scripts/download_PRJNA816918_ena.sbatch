#!/usr/bin/env bash
#SBATCH --job-name=ena-dl-PRJNA816918
#SBATCH --output=logs/ena-dl-%j.out
#SBATCH --error=logs/ena-dl-%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL


# Download all 18 paired-end FASTQ.gz files for PRJNA816918 directly from ENA
# (no sra-tools required). Files land flat in rawdata/ as:
#     rawdata/SRR18354890_1.fastq.gz
#     rawdata/SRR18354890_2.fastq.gz
#     ...
#
# Total size ~78 GB compressed. Job uses 1 sample at a time but each sample's
# _1 and _2 download in parallel; MD5 is verified against ENA's reported md5.
#
# Submit from repo root:
#   sbatch scripts/download_PRJNA816918_ena.sbatch
#
# Optional overrides:
#   sbatch --export=ALL,RAWDATA=/scratch/me/rawdata scripts/download_PRJNA816918_ena.sbatch
#   sbatch --export=ALL,RUNS_TSV=config/PRJNA816918_runs.tsv scripts/download_PRJNA816918_ena.sbatch

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"
ROOT_DIR="$(pwd)"

mkdir -p logs

RAWDATA="${RAWDATA:-${ROOT_DIR}/rawdata}"
RUNS_TSV="${RUNS_TSV:-${ROOT_DIR}/config/PRJNA816918_runs.tsv}"

mkdir -p "${RAWDATA}"

if [[ ! -s "${RUNS_TSV}" ]]; then
  echo "ERROR: runs metadata TSV not found: ${RUNS_TSV}" >&2
  exit 1
fi

command -v curl >/dev/null || { echo "ERROR: curl not on PATH" >&2; exit 1; }
command -v md5sum >/dev/null || { echo "ERROR: md5sum not on PATH" >&2; exit 1; }

echo "[$(date -Is)] ENA download starting -> ${RAWDATA}"
echo "[$(date -Is)] Reading runs from: ${RUNS_TSV}"

# Per-file retry with curl. Re-uses partial downloads (-C -).
download_one() {
  local url="$1" out="$2" want_md5="$3"
  local attempt=0 max_attempts=5

  if [[ -s "${out}" ]]; then
    local got
    got="$(md5sum "${out}" | awk '{print $1}')"
    if [[ "${got}" == "${want_md5}" ]]; then
      echo "  [skip] ${out##*/} (md5 ok)"
      return 0
    fi
    echo "  [redo] ${out##*/} md5 mismatch (got=${got} want=${want_md5}); re-downloading"
    rm -f "${out}"
  fi

  while (( attempt < max_attempts )); do
    attempt=$((attempt + 1))
    echo "  [get ] attempt ${attempt}/${max_attempts}: ${url}"
    if curl -fL --retry 3 --retry-delay 10 --connect-timeout 30 \
        --speed-time 60 --speed-limit 102400 \
        -C - -o "${out}" "https://${url}"; then
      local got
      got="$(md5sum "${out}" | awk '{print $1}')"
      if [[ "${got}" == "${want_md5}" ]]; then
        echo "  [ok  ] ${out##*/} md5 verified"
        return 0
      fi
      echo "  [warn] ${out##*/} md5 mismatch after download (got=${got} want=${want_md5}); retrying"
      rm -f "${out}"
    else
      echo "  [warn] curl failed for ${url} (attempt ${attempt})"
    fi
    sleep 15
  done

  echo "  [FAIL] ${out##*/} after ${max_attempts} attempts" >&2
  return 1
}

n_ok=0
n_fail=0
failed_files=()

# Skip header; use process substitution so counters live in this shell
while IFS=$'\t' read -r run experiment sample lib strat src reads bases ftps md5s sizes; do
  echo "[$(date -Is)] === ${run} (lib=${lib}) ==="
  IFS=';' read -r -a urls   <<< "${ftps}"
  IFS=';' read -r -a md5s_a <<< "${md5s}"
  if [[ "${#urls[@]}" -ne "${#md5s_a[@]}" ]]; then
    echo "  [skip] ${run}: url/md5 count mismatch (urls=${#urls[@]}, md5=${#md5s_a[@]})" >&2
    n_fail=$((n_fail + 1))
    failed_files+=("${run}")
    continue
  fi
  pids=()
  for i in "${!urls[@]}"; do
    url="${urls[$i]}"
    want_md5="${md5s_a[$i]}"
    out="${RAWDATA}/$(basename "${url}")"
    download_one "${url}" "${out}" "${want_md5}" &
    pids+=($!)
  done
  ok=1
  for pid in "${pids[@]}"; do
    if ! wait "${pid}"; then
      ok=0
    fi
  done
  if (( ok == 1 )); then
    n_ok=$((n_ok + 1))
  else
    n_fail=$((n_fail + 1))
    failed_files+=("${run}")
  fi
done < <(tail -n +2 "${RUNS_TSV}")

echo "[$(date -Is)] DONE. ok=${n_ok}  failed=${n_fail}"
if (( n_fail > 0 )); then
  echo "Failed runs: ${failed_files[*]}" >&2
  exit 1
fi
