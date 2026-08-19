#!/usr/bin/env python3
"""Compute 'active GUS' abundance, replicating the total/candidate GUS pipeline.

Active GUS = the first 11 entries of Table S3 (candidate GUS), in order:
  E. coli, C. perfringens, R. gnavus, Eubacterium (UNC 361-28), S. agalactiae,
  E. eligens, B. fragilis, F. saccharivorans, B. mass, B. dorei, P. copri.

Because the candidate DIAMOND database holds 18 sequences and was searched with
--max-target-seqs 25 (>= 18), every significant hit (e<=1e-5) to any candidate
reference is already present in results/search/*.candidate_gus.tsv. Filtering
those rows to the 11 active reference IDs is therefore identical to running
DIAMOND against an 11-sequence 'active' database. We reuse the exact same
summarize_gus_hits.py + calculate_gus_abundance.py steps as total/candidate.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MINLEN = 1000
REF_DIR = ROOT / "results/reference"
SEARCH_DIR = ROOT / "results/search"
SUMMARY_DIR = ROOT / "results/summary"
COUNTS_DIR = ROOT / "results/counts"
PRODIGAL_DIR = ROOT / "results/prodigal"
ASSEMBLY_DIR = ROOT / "assembly/megahit"
CONTIGS_BASENAME = "final.contigs.fa"

# Order matches Table S3 rows 1-11 (== first 11 rows of candidate metadata/faa).
ACTIVE_IDS = [
    "NP_416134.1",     # E. coli
    "WP_003467686.1",  # C. perfringens
    "6MVG_A",          # R. gnavus
    "WP_022501306",    # Eubacterium (UNC 361-28)
    "WP_000966715",    # S. agalactiae
    "WP_012740861",    # E. eligens
    "3CMG_A",          # B. fragilis
    "6NCY_A",          # F. saccharivorans
    "WP_270534885",    # B. mass
    "WP_007841259",    # B. dorei
    "WP_217761416",    # P. copri
]
ACTIVE_ID_SET = set(ACTIVE_IDS)

# (name, min_identity, min_query_cov, min_subject_cov, abundance_dir, hits_dir, keep_clade)
# hits_dir=None -> best_hits/clade written to a throwaway temp dir (not kept),
# matching the profile where candidate keeps only the abundance file.
PROFILES = [
    ("default40", 40.0, 0.5, 0.5, SUMMARY_DIR, None, False),
    ("pident90", 90.0, 0.5, 0.5, SUMMARY_DIR / "pident90", SUMMARY_DIR / "pident90", True),
    ("pident95", 95.0, 0.5, 0.5, SUMMARY_DIR / "pident95/abundance", SUMMARY_DIR / "pident95", True),
    ("pident95_cov90", 95.0, 0.9, 0.9, SUMMARY_DIR / "pident95_cov90/abundance", SUMMARY_DIR / "pident95_cov90/hits", False),
]

TMP_DIR = SUMMARY_DIR / ".active_tmp"


def build_active_reference() -> None:
    # active_gus.faa: keep the 11 records from candidate_gus.faa (order preserved).
    cand_faa = (REF_DIR / "candidate_gus.faa").read_text().splitlines()
    out_lines: list[str] = []
    keep = False
    for line in cand_faa:
        if line.startswith(">"):
            rid = line[1:].split()[0]
            keep = rid in ACTIVE_ID_SET
        if keep:
            out_lines.append(line)
    (REF_DIR / "active_gus.faa").write_text("\n".join(out_lines) + "\n")

    # metadata + taxonomy: header + rows whose reference_id is active.
    for src, dst in [
        ("candidate_gus_metadata.tsv", "active_gus_metadata.tsv"),
        ("candidate_gus_metadata.taxonomy.tsv", "active_gus_metadata.taxonomy.tsv"),
    ]:
        src_path = REF_DIR / src
        if not src_path.exists():
            continue
        rows = src_path.read_text().splitlines()
        header = rows[0]
        kept = [r for r in rows[1:] if r.split("\t", 1)[0] in ACTIVE_ID_SET]
        (REF_DIR / dst).write_text("\n".join([header, *kept]) + "\n")
    print(f"Built active reference ({len(ACTIVE_IDS)} sequences).")


def samples() -> list[str]:
    return sorted(
        p.name[: -len(f".min{MINLEN}.candidate_gus.tsv")]
        for p in SEARCH_DIR.glob(f"*.min{MINLEN}.candidate_gus.tsv")
    )


def build_active_search(sample: str) -> Path:
    cand = SEARCH_DIR / f"{sample}.min{MINLEN}.candidate_gus.tsv"
    active = SEARCH_DIR / f"{sample}.min{MINLEN}.active_gus.tsv"
    with cand.open() as fin, active.open("w") as fout:
        for line in fin:
            # column 2 (0-based index 1) is sseqid
            parts = line.split("\t")
            if len(parts) > 1 and parts[1] in ACTIVE_ID_SET:
                fout.write(line)
    return active


def run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True, cwd=ROOT)


def main() -> None:
    build_active_reference()
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    metadata = REF_DIR / "active_gus_metadata.tsv"

    sample_list = samples()
    if not sample_list:
        sys.exit("No candidate search files found under results/search/")

    for sample in sample_list:
        contigs = ASSEMBLY_DIR / sample / CONTIGS_BASENAME
        counts = COUNTS_DIR / f"{sample}.min{MINLEN}.featurecounts.txt"
        faa = PRODIGAL_DIR / f"{sample}.min{MINLEN}.faa"
        gff = PRODIGAL_DIR / f"{sample}.min{MINLEN}.gff"
        for path in (contigs, counts, faa, gff):
            if not path.exists():
                sys.exit(f"Missing required file for {sample}: {path}")

        active_search = build_active_search(sample)

        for name, min_id, qcov, scov, abundance_dir, hits_dir, keep_clade in PROFILES:
            abundance_dir.mkdir(parents=True, exist_ok=True)
            best_dir = hits_dir if hits_dir is not None else TMP_DIR
            best_dir.mkdir(parents=True, exist_ok=True)
            clade_dir = best_dir if keep_clade else TMP_DIR

            best = best_dir / f"{sample}.min{MINLEN}.active_gus_best_hits.tsv"
            clade = clade_dir / f"{sample}.min{MINLEN}.active_gus_clade_abundance.tsv"
            abundance = abundance_dir / f"{sample}.min{MINLEN}.active_gus_abundance.tsv"

            print(f"[{sample}] active -> {name} (id>={min_id}, qcov>={qcov}, scov>={scov})")
            run([
                "python3", "scripts/summarize_gus_hits.py",
                "--hits", str(active_search),
                "--contigs", str(contigs),
                "--metadata", str(metadata),
                "--best-hits", str(best),
                "--clade-summary", str(clade),
                "--min-identity", str(min_id),
                "--min-query-cov", str(qcov),
                "--min-subject-cov", str(scov),
                "--max-evalue", "1e-5",
            ])
            run([
                "python3", "scripts/calculate_gus_abundance.py",
                "--hits", str(best),
                "--counts", str(counts),
                "--faa", str(faa),
                "--gff", str(gff),
                "--output", str(abundance),
            ])

    print("Done. Active GUS abundance generated for all profiles.")


if __name__ == "__main__":
    main()
