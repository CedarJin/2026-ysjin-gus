#!/usr/bin/env python3
"""Summarize GUS similarity hits and contig-level abundance weights."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


CONTIG_HEADER = re.compile(r"^>(?P<contig>\S+).*?\bmulti=(?P<multi>[0-9.]+).*?\blen=(?P<length>\d+)")


def load_contigs(path: Path) -> dict[str, dict[str, float]]:
    contigs: dict[str, dict[str, float]] = {}
    with path.open() as handle:
        for line in handle:
            match = CONTIG_HEADER.match(line)
            if not match:
                continue
            length = int(match.group("length"))
            multi = float(match.group("multi"))
            contigs[match.group("contig")] = {
                "contig_length": length,
                "contig_multi": multi,
                "contig_weight": length * multi,
            }
    return contigs


def load_metadata(path: Path) -> dict[str, dict[str, str]]:
    with path.open() as handle:
        return {row["reference_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def contig_from_gene(gene_id: str) -> str:
    # Prodigal appends _1, _2, ... to the source contig id.
    return gene_id.rsplit("_", 1)[0]


def parse_hits(
    hits_path: Path,
    contigs: dict[str, dict[str, float]],
    metadata: dict[str, dict[str, str]],
    min_identity: float,
    min_query_cov: float,
    min_subject_cov: float,
    max_evalue: float,
) -> list[dict[str, str]]:
    best_by_gene: dict[str, dict[str, str]] = {}
    fields = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qlen",
        "slen",
    ]
    with hits_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t", fieldnames=fields):
            pident = float(row["pident"])
            evalue = float(row["evalue"])
            aln_len = int(row["length"])
            qcov = aln_len / int(row["qlen"])
            scov = aln_len / int(row["slen"])
            if (
                pident < min_identity
                or evalue > max_evalue
                or qcov < min_query_cov
                or scov < min_subject_cov
            ):
                continue
            previous = best_by_gene.get(row["qseqid"])
            if previous is None or float(row["bitscore"]) > float(previous["bitscore"]):
                contig_id = contig_from_gene(row["qseqid"])
                contig = contigs.get(contig_id, {})
                ref = metadata.get(row["sseqid"], {})
                best_by_gene[row["qseqid"]] = {
                    **row,
                    "query_coverage": f"{qcov:.4f}",
                    "subject_coverage": f"{scov:.4f}",
                    "contig_id": contig_id,
                    "contig_length": str(int(contig.get("contig_length", 0))),
                    "contig_multi": f"{contig.get('contig_multi', 0):.4f}",
                    "contig_weight": f"{contig.get('contig_weight', 0):.4f}",
                    "reference_clade": ref.get("clade", ""),
                    "reference_source": ref.get("source", ""),
                }
    return sorted(best_by_gene.values(), key=lambda row: float(row["bitscore"]), reverse=True)


def write_best_hits(rows: list[dict[str, str]], path: Path) -> None:
    fields = [
        "qseqid",
        "contig_id",
        "sseqid",
        "reference_clade",
        "reference_source",
        "pident",
        "length",
        "query_coverage",
        "subject_coverage",
        "evalue",
        "bitscore",
        "qlen",
        "slen",
        "contig_length",
        "contig_multi",
        "contig_weight",
    ]
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_clade_summary(rows: list[dict[str, str]], path: Path) -> None:
    summary: dict[str, dict[str, float]] = {}
    for row in rows:
        clade = row["reference_clade"] or "unknown"
        item = summary.setdefault(clade, {"genes": 0, "contig_weight": 0.0})
        item["genes"] += 1
        item["contig_weight"] += float(row["contig_weight"])

    with path.open("w") as handle:
        handle.write("clade\tgenes\tcontig_weight\n")
        for clade, item in sorted(summary.items(), key=lambda kv: kv[1]["contig_weight"], reverse=True):
            handle.write(f"{clade}\t{int(item['genes'])}\t{item['contig_weight']:.4f}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits", required=True, type=Path)
    parser.add_argument("--contigs", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--best-hits", required=True, type=Path)
    parser.add_argument("--clade-summary", required=True, type=Path)
    parser.add_argument("--min-identity", type=float, default=40.0)
    parser.add_argument("--min-query-cov", type=float, default=0.5)
    parser.add_argument("--min-subject-cov", type=float, default=0.5)
    parser.add_argument("--max-evalue", type=float, default=1e-5)
    args = parser.parse_args()

    rows = parse_hits(
        args.hits,
        load_contigs(args.contigs),
        load_metadata(args.metadata),
        args.min_identity,
        args.min_query_cov,
        args.min_subject_cov,
        args.max_evalue,
    )
    args.best_hits.parent.mkdir(parents=True, exist_ok=True)
    write_best_hits(rows, args.best_hits)
    write_clade_summary(rows, args.clade_summary)
    print(f"Kept {len(rows)} best gene-level GUS hits")


if __name__ == "__main__":
    main()
