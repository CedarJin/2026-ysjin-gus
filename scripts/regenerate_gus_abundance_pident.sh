#!/usr/bin/env bash
# Regenerate GUS summary/abundance at alternate filter profiles.
# Requires existing search/, mapping BAM, and prodigal outputs from gus_abundance.smk.
#
# Profiles: OUT_SUBDIR MIN_IDENTITY MIN_QUERY_COV MIN_SUBJECT_COV
# Default runs pident90, pident95, and pident95_cov90.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

MINLEN="${MIN_CONTIG_LENGTH:-1000}"
RESULTS="${RESULTS_DIR:-results}"
ASSEMBLY_DIR="${ASSEMBLY_DIR:-assembly/megahit}"
CONTIGS_BASENAME="${CONTIGS_BASENAME:-final.contigs.fa}"
THREADS_FC="${THREADS_FEATURECOUNTS:-8}"
MAX_EVALUE="${MAX_EVALUE:-1e-5}"

DATABASES=(total candidate)

# out_subdir min_identity min_query_cov min_subject_cov
DEFAULT_PROFILES=(
  "pident90 90 0.5 0.5"
  "pident95 95 0.5 0.5"
  "pident95_cov90 95 0.9 0.9"
)

if [[ $# -gt 0 ]]; then
  PROFILES=("$@")
else
  PROFILES=("${DEFAULT_PROFILES[@]}")
fi

mapfile -t SAMPLES < <(
  for f in "${RESULTS}/search/"*.min"${MINLEN}".total_gus.tsv; do
    [[ -f "${f}" ]] || continue
    basename "${f}" ".min${MINLEN}.total_gus.tsv"
  done | sort -u
)

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "No samples found under ${RESULTS}/search/*.min${MINLEN}.total_gus.tsv" >&2
  exit 1
fi

mkdir -p "${RESULTS}/counts"
for profile in "${PROFILES[@]}"; do
  read -r out_subdir _ _ _ <<< "${profile}"
  mkdir -p "${RESULTS}/summary/${out_subdir}"
done

for sample in "${SAMPLES[@]}"; do
  counts="${RESULTS}/counts/${sample}.min${MINLEN}.featurecounts.txt"
  bam="${RESULTS}/mapping/${sample}.min${MINLEN}.sorted.bam"
  gff="${RESULTS}/prodigal/${sample}.min${MINLEN}.gff"
  faa="${RESULTS}/prodigal/${sample}.min${MINLEN}.faa"
  contigs="${ASSEMBLY_DIR}/${sample}/${CONTIGS_BASENAME}"

  for path in "${bam}" "${gff}" "${faa}" "${contigs}"; do
    if [[ ! -f "${path}" ]]; then
      echo "Missing required file for ${sample}: ${path}" >&2
      exit 1
    fi
  done

  if [[ ! -f "${counts}" ]]; then
    echo "[${sample}] featureCounts"
    featureCounts -T "${THREADS_FC}" -p -B -C -t CDS -g ID -a "${gff}" -o "${counts}" "${bam}"
  fi

  for db in "${DATABASES[@]}"; do
    hits="${RESULTS}/search/${sample}.min${MINLEN}.${db}_gus.tsv"
    if [[ "${db}" == "total" ]]; then
      metadata="${RESULTS}/reference/gus_279_plus_fmn_metadata.tsv"
    else
      metadata="${RESULTS}/reference/candidate_gus_metadata.tsv"
    fi

    for profile in "${PROFILES[@]}"; do
      read -r out_subdir min_identity min_query_cov min_subject_cov <<< "${profile}"
      out_dir="${RESULTS}/summary/${out_subdir}"
      best="${out_dir}/${sample}.min${MINLEN}.${db}_gus_best_hits.tsv"
      clade="${out_dir}/${sample}.min${MINLEN}.${db}_gus_clade_abundance.tsv"
      abundance="${out_dir}/${sample}.min${MINLEN}.${db}_gus_abundance.tsv"

      echo "[${sample}] ${db} -> ${out_subdir} (pident>=${min_identity}, qcov>=${min_query_cov}, scov>=${min_subject_cov})"
      python3 scripts/summarize_gus_hits.py \
        --hits "${hits}" \
        --contigs "${contigs}" \
        --metadata "${metadata}" \
        --best-hits "${best}" \
        --clade-summary "${clade}" \
        --min-identity "${min_identity}" \
        --min-query-cov "${min_query_cov}" \
        --min-subject-cov "${min_subject_cov}" \
        --max-evalue "${MAX_EVALUE}"

      python3 scripts/calculate_gus_abundance.py \
        --hits "${best}" \
        --counts "${counts}" \
        --faa "${faa}" \
        --gff "${gff}" \
        --output "${abundance}"
    done
  done
done

echo "Done. Profiles written under ${RESULTS}/summary/ for: ${PROFILES[*]}"
