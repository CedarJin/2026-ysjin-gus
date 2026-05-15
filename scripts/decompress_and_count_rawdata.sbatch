#!/usr/bin/env bash
#SBATCH --job-name=fq-decompress
#SBATCH --output=logs/fq-decompress-%j.out
#SBATCH --error=logs/fq-decompress-%j.err
#SBATCH --partition=high
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=END,FAIL

# Decompress every *.fastq.gz in rawdata/ to *.fastq (kept side-by-side)
# and count reads (lines / 4) per file using gunzip + wc.
#
# Summary TSV is written to logs/read_counts_<jobid>.tsv with columns:
#   file        reads      size_gz_bytes   size_fq_bytes
#
# Submit from repo root:
#   sbatch scripts/decompress_and_count_rawdata.sbatch
#
# Optional overrides:
#   sbatch --export=ALL,RAWDATA=/path/to/rawdata scripts/decompress_and_count_rawdata.sbatch
#   sbatch --export=ALL,REMOVE_GZ=1 scripts/decompress_and_count_rawdata.sbatch  # delete .gz after success
#   sbatch --export=ALL,FORCE=1     scripts/decompress_and_count_rawdata.sbatch  # re-decompress even if .fastq exists

set -eo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"
ROOT_DIR="$(pwd)"
mkdir -p logs

RAWDATA="${RAWDATA:-${ROOT_DIR}/rawdata}"
REMOVE_GZ="${REMOVE_GZ:-0}"
FORCE="${FORCE:-0}"

if [[ ! -d "${RAWDATA}" ]]; then
  echo "ERROR: RAWDATA not found: ${RAWDATA}" >&2
  exit 1
fi

command -v gunzip >/dev/null || { echo "ERROR: gunzip not on PATH" >&2; exit 1; }

SUMMARY="logs/read_counts_${SLURM_JOB_ID:-local}.tsv"
printf 'file\treads\tsize_gz_bytes\tsize_fq_bytes\n' > "${SUMMARY}"
echo "[$(date -Is)] decompressing with gunzip"
echo "[$(date -Is)] summary will be written to: ${SUMMARY}"

shopt -s nullglob
gz_files=("${RAWDATA}"/*.fastq.gz)
if (( ${#gz_files[@]} == 0 )); then
  echo "ERROR: no *.fastq.gz files in ${RAWDATA}" >&2
  exit 1
fi
echo "[$(date -Is)] found ${#gz_files[@]} .fastq.gz files in ${RAWDATA}"

n_ok=0
n_skip=0
n_fail=0
total_reads=0
i=0

for gz in "${gz_files[@]}"; do
  i=$((i + 1))
  fq="${gz%.gz}"
  base="$(basename "${gz}")"
  size_gz="$(stat -c '%s' "${gz}")"

  echo "[$(date -Is)] (${i}/${#gz_files[@]}) --- ${base} ---"

  # Reuse existing .fastq if present (and not forced).
  if [[ -s "${fq}" && "${FORCE}" != "1" ]]; then
    lines=$(wc -l < "${fq}")
    reads=$(( lines / 4 ))
    size_fq="$(stat -c '%s' "${fq}")"
    echo "  [skip] $(basename "${fq}") already exists; counted ${reads} reads"
    printf '%s\t%d\t%d\t%d\n' "$(basename "${fq}")" "${reads}" "${size_gz}" "${size_fq}" >> "${SUMMARY}"
    n_skip=$((n_skip + 1))
    total_reads=$(( total_reads + reads ))
    continue
  fi

  tmp="${fq}.tmp"
  rm -f "${tmp}"

  # Step 1: decompress
  if ! gunzip -c "${gz}" > "${tmp}"; then
    echo "  [FAIL] gunzip failed for ${base}" >&2
    rm -f "${tmp}"
    n_fail=$((n_fail + 1))
    continue
  fi

  # Step 2: count
  lines=$(wc -l < "${tmp}")
  if (( lines % 4 != 0 )); then
    echo "  [warn] ${base}: line count ${lines} not divisible by 4 (truncated file?)" >&2
  fi
  reads=$(( lines / 4 ))

  mv -f "${tmp}" "${fq}"
  size_fq="$(stat -c '%s' "${fq}")"
  echo "  [ok  ] $(basename "${fq}")  reads=${reads}  ($(numfmt --to=iec "${size_gz}") -> $(numfmt --to=iec "${size_fq}"))"
  printf '%s\t%d\t%d\t%d\n' "$(basename "${fq}")" "${reads}" "${size_gz}" "${size_fq}" >> "${SUMMARY}"

  if [[ "${REMOVE_GZ}" == "1" ]]; then
    rm -f "${gz}"
    echo "  [rm  ] ${base}"
  fi

  n_ok=$((n_ok + 1))
  total_reads=$(( total_reads + reads ))
done

echo
echo "[$(date -Is)] DONE.  ok=${n_ok}  skipped=${n_skip}  failed=${n_fail}  total_reads=${total_reads}"
echo "[$(date -Is)] summary TSV: ${SUMMARY}"

# Per-sample (SRR) read totals across mates (_1 + _2).
echo
echo "Per-run paired totals (reads in _1 + _2):"
awk -F'\t' 'NR>1 {
  split($1, a, "_");
  run=a[1];
  reads[run]+=$2;
  count[run]++;
}
END {
  for (r in reads) printf "  %s\t%d reads across %d file(s)\n", r, reads[r], count[r];
}' "${SUMMARY}" | sort

if (( n_fail > 0 )); then
  exit 1
fi
