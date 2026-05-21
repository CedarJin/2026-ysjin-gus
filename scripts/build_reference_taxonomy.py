#!/usr/bin/env python3
"""Build taxonomy annotation tables for GUS reference metadata.

This script uses metadata source names, queries TaxonKit for lineage,
and writes annotated TSVs.

Expected TaxonKit setup:
  - `taxonkit` executable available on PATH, or use `--taxonkit`
  - NCBI taxdump installed and reachable via `--data-dir` or `TAXONKIT_DB`

Outputs include the original metadata columns plus:
  source_name, matched_name, matched_level, taxid, kingdom, phylum, class, order, family, genus, species
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd


RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

PHYLUM_SYNONYMS = {
    "Firmicutes": "Bacillota",
    "Bacteroidetes": "Bacteroidota",
}

SPECIES_SYNONYMS = {
    "R.gnavus": "Ruminococcus gnavus",
    "B.mass": "Phocaeicola massiliensis",
    "B.dorei": "Phocaeicola",
    "P.copri": "Segatella copri",
    "6MVG_A": "Mediterraneibacter gnavus",
    "WP_012740861": "Lachnospira eligens",
    "E.eligens": "Lachnospira eligens",
}

def extract_source_name(raw: str) -> str:
    raw = str(raw).strip()
    if "source=" in raw:
        raw = raw.split("source=", 1)[1].strip()
    raw = raw.replace("_", " ")
    raw = re.sub(r"\s+", " ", raw).strip()
    return raw


def normalize_candidate(raw: str) -> str:
    raw = raw.strip()
    raw = re.sub(r"\bCAG:\s*\d+\b", "", raw).strip()
    raw = re.sub(r"\bATCC\s+\d+\b", "", raw).strip()
    raw = re.sub(r"\bSL\d+(?:/\d+)?\b", "", raw).strip()
    raw = re.sub(r"\s+", " ", raw).strip()
    return raw


def candidate_priority(level: str) -> int:
    order = {
        "source_name": 0,
        "abbr_expansion": 1,
    }
    if level.startswith("prefix_"):
        try:
            n = int(level.split("_", 1)[1])
        except ValueError:
            n = 0
        # Larger prefix length means more specific, so it should sort first.
        return 1000 - n
    return order.get(level, 500)


def has_taxid(value) -> bool:
    return pd.notna(value) and str(value).strip() != ""


def candidate_names(source_name: str) -> list[tuple[str, str]]:
    """Return (candidate_name, matched_level) ordered from most to least specific."""
    source_name = extract_source_name(source_name)
    candidates: list[tuple[str, str]] = []

    def add(name: str, level: str):
        name = normalize_candidate(name)
        if name and name not in {c for c, _ in candidates}:
            candidates.append((name, level))

    add(source_name, "source_name")

    if source_name in SPECIES_SYNONYMS:
        add(SPECIES_SYNONYMS[source_name], "species_synonym")

    if source_name in PHYLUM_SYNONYMS:
        add(PHYLUM_SYNONYMS[source_name], "synonym")
    elif source_name.startswith("Firmicutes "):
        add("Bacillota", "synonym")
    elif source_name.startswith("Bacteroidetes "):
        add("Bacteroidota", "synonym")

    parts = source_name.split()
    # Full prefix ladder: species -> genus -> family -> ... by chopping trailing words.
    for i in range(len(parts) - 1, 0, -1):
        add(" ".join(parts[:i]), f"prefix_{i}")

    # If the string starts with an abbreviated genus like B.fragilis, expand a few common ones.
    if "." in source_name:
        first = source_name.split()[0]
        if first in {
            "E.coli",
            "C.perfringens",
            "R.gnavus",
            "E.eligens",
            "B.fragilis",
            "F.saccharvivorans",
            "B.mass",
            "B.dorei",
            "P.copri",
            "R.hominis",
            "F.prausnitzii",
            "P.merdae",
            "B.ovatus",
            "S.agalactiae",
        }:
            add(
                {
                    "E.coli": "Escherichia coli",
                    "C.perfringens": "Clostridium perfringens",
                    "R.gnavus": "Ruminococcus gnavus",
                    "E.eligens": "Eubacterium eligens",
                    "B.fragilis": "Bacteroides fragilis",
                    "F.saccharvivorans": "Fusicatenibacter saccharivorans",
                    "B.mass": "Bacteroides massiliensis",
                    "B.dorei": "Bacteroides dorei",
                    "P.copri": "Prevotella copri",
                    "R.hominis": "Roseburia hominis",
                    "F.prausnitzii": "Faecalibacterium prausnitzii",
                    "P.merdae": "Parabacteroides merdae",
                    "B.ovatus": "Bacteroides ovatus",
                    "S.agalactiae": "Streptococcus agalactiae",
                }[first],
                "abbr_expansion",
            )

    return candidates


def unique_sources(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["source_name"] = out["source"].map(extract_source_name)
    return out


def run_taxonkit(name_file: Path, taxonkit: str, data_dir: str | None) -> pd.DataFrame:
    cmd_name2taxid = [taxonkit]
    if data_dir:
        cmd_name2taxid += ["--data-dir", data_dir]
    cmd_name2taxid += ["name2taxid", "--sci-name", "--show-rank", "-i", "1"]

    name2taxid = subprocess.run(
        cmd_name2taxid,
        stdin=name_file.open("r", encoding="utf-8"),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if name2taxid.returncode != 0:
        raise RuntimeError(
            "taxonkit name2taxid failed:\n"
            + name2taxid.stderr
            + "\nTip: make sure the NCBI taxdump is installed and pass --data-dir."
        )

    with tempfile.NamedTemporaryFile("w+", delete=False, encoding="utf-8") as tf:
        tf.write(name2taxid.stdout)
        tmp_path = Path(tf.name)

    cmd_lineage = [taxonkit]
    if data_dir:
        cmd_lineage += ["--data-dir", data_dir]
    cmd_lineage += [
        "reformat2",
        "-I",
        "2",
        "-f",
        "{domain|superkingdom};{phylum};{class};{order};{family};{genus};{species}",
    ]
    lineage = subprocess.run(
        cmd_lineage,
        stdin=tmp_path.open("r", encoding="utf-8"),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    tmp_path.unlink(missing_ok=True)
    if lineage.returncode != 0:
        raise RuntimeError(
            "taxonkit lineage failed:\n"
            + lineage.stderr
            + "\nTip: make sure the NCBI taxdump is installed and pass --data-dir."
        )

    rows = []
    for line in lineage.stdout.strip().splitlines():
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        # candidate_name, taxid, rank, reformatted lineage
        matched_name, taxid = parts[0], parts[1]
        lineage_names = parts[-1].split(";") if parts[-1] else []
        row = {"matched_name": matched_name, "taxid": taxid}
        for idx, rank in enumerate(RANKS):
            row[rank] = lineage_names[idx] if idx < len(lineage_names) else ""
        rows.append(row)
    return pd.DataFrame(rows)


def annotate_metadata(metadata_path: Path, taxonkit: str, data_dir: str | None) -> pd.DataFrame:
    df = pd.read_csv(metadata_path, sep="\t")
    df = unique_sources(df)
    source_to_candidates = {
        source_name: candidate_names(source_name)
        for source_name in df["source_name"].drop_duplicates().sort_values()
    }
    lookup_rows = [
        {
            "source_name": source_name,
            "matched_name": candidate,
            "matched_level": level,
        }
        for source_name, candidates in source_to_candidates.items()
        for candidate, level in candidates
    ]
    lookup = pd.DataFrame(lookup_rows)

    with tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8") as fh:
        for name in lookup["matched_name"].drop_duplicates():
            fh.write(f"{name}\n")
        name_file = Path(fh.name)

    try:
        tax_df = run_taxonkit(name_file, taxonkit=taxonkit, data_dir=data_dir)
    finally:
        name_file.unlink(missing_ok=True)

    tax_hits = {}
    for _, row in tax_df.iterrows():
        tax_hits[row["matched_name"]] = row.to_dict()
    selected_rows = []
    for source_name, candidates in source_to_candidates.items():
        best = None
        for candidate, level in candidates:
            hit = tax_hits.get(candidate)
            if hit and has_taxid(hit.get("taxid")):
                best = {
                    "source_name": source_name,
                    "matched_name": candidate,
                    "matched_level": level,
                    **hit,
                }
                break
        if best is None:
            candidate, level = candidates[0]
            hit = tax_hits.get(candidate, {})
            best = {
                "source_name": source_name,
                "matched_name": candidate,
                "matched_level": level,
                "taxid": hit.get("taxid"),
            }
            for rank in RANKS:
                best[rank] = hit.get(rank, "")
        selected_rows.append(best)

    tax_df = pd.DataFrame(selected_rows)
    annotated = df.merge(tax_df, on="source_name", how="left", suffixes=("", "_tax"))
    return annotated


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--metadata",
        nargs="+",
        required=True,
        help="Input metadata TSVs to annotate.",
    )
    parser.add_argument(
        "--out-dir",
        default="results/reference",
        help="Directory to write annotated TSVs.",
    )
    parser.add_argument(
        "--taxonkit",
        default="taxonkit",
        help="Path to the taxonkit executable.",
    )
    parser.add_argument(
        "--data-dir",
        default=os.environ.get("TAXONKIT_DB"),
        help="TaxonKit data directory containing names.dmp and nodes.dmp.",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for metadata in args.metadata:
        in_path = Path(metadata)
        annotated = annotate_metadata(in_path, taxonkit=args.taxonkit, data_dir=args.data_dir)
        out_path = out_dir / f"{in_path.stem}.taxonomy.tsv"
        annotated.to_csv(out_path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
        print(f"Wrote {out_path}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
