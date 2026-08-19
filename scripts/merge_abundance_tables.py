#!/usr/bin/env python3
"""Merge per-sample abundance tables into cross-subject matrices."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


TABLE_SPECS = {
    "total_gus_abundance": "total_gus_abundance.tsv",
    "candidate_gus_abundance": "candidate_gus_abundance.tsv",
}

VALUE_COLS = {
    "gus_db": "tpm_gus_db",
    "gus_genome": "tpm_genome_wide",
}


def subject_id_from_path(path: Path) -> str:
    return path.name.split(".min1000.", 1)[0]


def load_sample_table(path: Path, value_col: str) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    if "sseqid" not in df.columns:
        raise SystemExit(f"Missing sseqid column in {path}")
    if value_col not in df.columns:
        raise SystemExit(f"Missing {value_col} column in {path}")

    series = df.groupby("sseqid", sort=False)[value_col].sum()
    series.name = subject_id_from_path(path)
    return series


def build_matrix(input_dir: Path, suffix: str, value_col: str) -> pd.DataFrame:
    paths = sorted(input_dir.glob(f"*.{suffix}"))
    if not paths:
        raise SystemExit(f"No files matching *.{suffix} found in {input_dir}")

    series_list = [load_sample_table(path, value_col) for path in paths]
    matrix = pd.concat(series_list, axis=1).fillna(0.0)
    matrix.index.name = "sseqid"
    matrix = matrix.reindex(sorted(matrix.index))
    matrix = matrix.reindex(sorted(matrix.columns), axis=1)
    return matrix


def to_relative_abundance_percent_within_subset(matrix: pd.DataFrame) -> pd.DataFrame:
    column_sums = matrix.sum(axis=0).replace(0, pd.NA)
    percent = matrix.divide(column_sums, axis=1) * 100.0
    return percent.fillna(0.0)


def to_relative_abundance_percent_of_genome(matrix: pd.DataFrame) -> pd.DataFrame:
    # tpm_genome_wide sums to 1e6 across all genes in a sample.
    # Convert TPM to percent of the full genome by dividing by 10,000.
    return matrix / 10_000.0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("results/summary/pident95_cov90/abundance"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/summary/pident95_cov90/abundance_all_subject"),
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    outputs = [
        (
            "total_gus_abundance_within_gus_db_percent_all_subject.tsv",
            "total_gus_abundance.tsv",
            "tpm_gus_db",
        ),
        (
            "candidate_gus_abundance_within_gus_db_percent_all_subject.tsv",
            "candidate_gus_abundance.tsv",
            "tpm_gus_db",
        ),
        (
            "total_gus_abundance_of_genome_percent_all_subject.tsv",
            "total_gus_abundance.tsv",
            "tpm_genome_wide",
        ),
        (
            "candidate_gus_abundance_of_genome_percent_all_subject.tsv",
            "candidate_gus_abundance.tsv",
            "tpm_genome_wide",
        ),
        (
            "active_gus_abundance_within_gus_db_percent_all_subject.tsv",
            "active_gus_abundance.tsv",
            "tpm_gus_db",
        ),
        (
            "active_gus_abundance_of_genome_percent_all_subject.tsv",
            "active_gus_abundance.tsv",
            "tpm_genome_wide",
        ),
    ]

    tpm_outputs = [
        (
            "active_gus_abundance_gus_genome_all_subject.tsv",
            "active_gus_abundance.tsv",
            "tpm_genome_wide",
        ),
    ]

    for output_name, suffix, value_col in outputs:
        matrix = build_matrix(args.input_dir, suffix, value_col)
        if value_col == "tpm_gus_db":
            matrix = to_relative_abundance_percent_within_subset(matrix)
        else:
            matrix = to_relative_abundance_percent_of_genome(matrix)
        output_path = args.output_dir / output_name
        matrix.to_csv(output_path, sep="\t", float_format="%.6f")
        print(f"Wrote {output_path} with shape {matrix.shape}")

    for output_name, suffix, value_col in tpm_outputs:
        matrix = build_matrix(args.input_dir, suffix, value_col)
        output_path = args.output_dir / output_name
        matrix.to_csv(output_path, sep="\t", float_format="%.6f")
        print(f"Wrote {output_path} with shape {matrix.shape}")


if __name__ == "__main__":
    main()
