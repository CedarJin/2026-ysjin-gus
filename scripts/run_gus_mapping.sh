#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONTIGS="${ROOT_DIR}/assembly/SRR10692699_50M_seed11/final.contigs.fa"
MIN_CONTIG_LENGTH="${MIN_CONTIG_LENGTH:-1000}"
FILTERED_DIR="${ROOT_DIR}/results/filtered_contigs"
FILTERED_CONTIGS="${FILTERED_DIR}/SRR10692699_50M_seed11.min${MIN_CONTIG_LENGTH}.fa"
RUN_LABEL="SRR10692699_50M_seed11.min${MIN_CONTIG_LENGTH}"
PDF="${ROOT_DIR}/reference_from_paper/279_gus_seq.pdf"
FMN_XLSX="${ROOT_DIR}/reference_from_paper/Supplemental_Table_2.xlsx"
REF_DIR="${ROOT_DIR}/results/reference"
PRODIGAL_DIR="${ROOT_DIR}/results/prodigal"
SEARCH_DIR="${ROOT_DIR}/results/search"
SUMMARY_DIR="${ROOT_DIR}/results/summary"

mkdir -p "${FILTERED_DIR}" "${REF_DIR}" "${PRODIGAL_DIR}" "${SEARCH_DIR}" "${SUMMARY_DIR}"

python3 "${ROOT_DIR}/scripts/build_gus_database.py" \
  --paper-fasta "${REF_DIR}/gus_279.faa" \
  --fmn-xlsx "${FMN_XLSX}" \
  --fasta "${REF_DIR}/gus_279_plus_fmn.faa" \
  --metadata "${REF_DIR}/gus_279_plus_fmn_metadata.tsv" \
  --expected-paper-count 279

python3 "${ROOT_DIR}/scripts/filter_contigs_by_length.py" \
  --input "${CONTIGS}" \
  --output "${FILTERED_CONTIGS}" \
  --min-length "${MIN_CONTIG_LENGTH}"

prodigal \
  -i "${FILTERED_CONTIGS}" \
  -p meta \
  -q \
  -a "${PRODIGAL_DIR}/${RUN_LABEL}.faa" \
  -d "${PRODIGAL_DIR}/${RUN_LABEL}.ffn" \
  -o "${PRODIGAL_DIR}/${RUN_LABEL}.gff" \
  -f gff

diamond makedb \
  --in "${REF_DIR}/gus_279_plus_fmn.faa" \
  --db "${REF_DIR}/gus_279_plus_fmn"

diamond blastp \
  --query "${PRODIGAL_DIR}/${RUN_LABEL}.faa" \
  --db "${REF_DIR}/gus_279_plus_fmn" \
  --out "${SEARCH_DIR}/${RUN_LABEL}.predicted_proteins_vs_gus.tsv" \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
  --evalue 1e-5 \
  --max-target-seqs 25 \
  --threads "${THREADS:-4}"

python3 "${ROOT_DIR}/scripts/summarize_gus_hits.py" \
  --hits "${SEARCH_DIR}/${RUN_LABEL}.predicted_proteins_vs_gus.tsv" \
  --contigs "${FILTERED_CONTIGS}" \
  --metadata "${REF_DIR}/gus_279_plus_fmn_metadata.tsv" \
  --best-hits "${SUMMARY_DIR}/${RUN_LABEL}.gus_best_hits.tsv" \
  --clade-summary "${SUMMARY_DIR}/${RUN_LABEL}.gus_clade_abundance.tsv"
