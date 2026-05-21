#!/usr/bin/env python3
"""Calculate length- and library-size-normalized abundance for GUS hits."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_featurecounts(path: Path) -> tuple[dict[str, int], int]:
    counts: dict[str, int] = {}
    with path.open() as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    header_index = next((i for i, line in enumerate(lines) if line.startswith("Geneid\t")), None)
    if header_index is None:
        raise SystemExit(f"Could not find featureCounts header in {path}")
    for line in lines[header_index + 1 :]:
        parts = line.split("\t")
        gene_id = parts[0]
        count = int(float(parts[-1]))
        counts[gene_id] = count
    library_size = sum(counts.values())
    return counts, library_size


def load_best_hits(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_qseqid_to_featurecounts_id(faa_path: Path) -> dict[str, str]:
    """Map Prodigal FAA headers (k141_*) to featureCounts gene IDs (numeric_*)."""
    mapping: dict[str, str] = {}
    with faa_path.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            body = line[1:].strip()
            qseqid, _, rest = body.partition(" ")
            if "ID=" not in rest:
                raise SystemExit(
                    f"Prodigal FAA header missing ID= for {qseqid!r} in {faa_path}"
                )
            featurecounts_id = rest.split("ID=", 1)[1].split(";", 1)[0]
            mapping[qseqid] = featurecounts_id
    if not mapping:
        raise SystemExit(f"No sequences found in {faa_path}")
    return mapping


def load_gene_lengths_from_gff(gff_path: Path) -> dict[str, int]:
    """CDS length in bp from Prodigal GFF (matches featureCounts Geneid)."""
    lengths: dict[str, int] = {}
    with gff_path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            start, end = int(parts[3]), int(parts[4])
            attrs = parts[8]
            if "ID=" not in attrs:
                continue
            gene_id = attrs.split("ID=", 1)[1].split(";", 1)[0]
            lengths[gene_id] = end - start + 1
    if not lengths:
        raise SystemExit(f"No CDS features found in {gff_path}")
    return lengths


def rpk(count: int, length_bp: int) -> float:
    if length_bp <= 0:
        return 0.0
    return count / (length_bp / 1000.0)


def genome_rpk_sum(counts: dict[str, int], gene_lengths_bp: dict[str, int]) -> float:
    total = 0.0
    for gene_id, count in counts.items():
        if count <= 0:
            continue
        length_bp = gene_lengths_bp.get(gene_id)
        if length_bp is None:
            continue
        total += rpk(count, length_bp)
    return total


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits", required=True, type=Path)
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument(
        "--faa",
        required=True,
        type=Path,
        help="Prodigal FAA used for DIAMOND; maps qseqid to featureCounts Geneid",
    )
    parser.add_argument(
        "--gff",
        required=True,
        type=Path,
        help="Prodigal GFF for CDS lengths (genome-wide TPM denominator)",
    )
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    counts, library_size = load_featurecounts(args.counts)
    if library_size == 0:
        raise SystemExit(f"No counts found in {args.counts}")

    qseqid_to_fc = load_qseqid_to_featurecounts_id(args.faa)
    gene_lengths_bp = load_gene_lengths_from_gff(args.gff)
    genome_rpk_total = genome_rpk_sum(counts, gene_lengths_bp)
    if genome_rpk_total == 0:
        raise SystemExit(f"Genome-wide RPK sum is zero (check {args.gff} vs {args.counts})")

    rows = load_best_hits(args.hits)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) + [
        "featurecounts_gene_id",
        "cds_length_bp",
        "raw_count",
        "gene_length_bp",
        "library_size",
        "rpk",
        "cpm",
        "tpm_gus_db",
        "tpm_genome_wide",
    ]
    with args.output.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        gus_rpk_values: list[float] = []
        enriched: list[dict[str, str | int | float]] = []
        for row in rows:
            qseqid = row["qseqid"]
            featurecounts_gene_id = qseqid_to_fc.get(qseqid, qseqid)
            raw_count = counts.get(featurecounts_gene_id, 0)
            query_length_aa = int(row["qlen"])
            query_length_kb = query_length_aa / 1000.0
            gus_rpk = raw_count / query_length_kb if query_length_kb > 0 else 0.0
            cpm = raw_count / library_size * 1_000_000.0
            cds_length_bp = gene_lengths_bp.get(featurecounts_gene_id, 0)
            genome_wide_rpk = rpk(raw_count, cds_length_bp)
            item: dict[str, str | int | float] = dict(row)
            item.update(
                {
                    "featurecounts_gene_id": featurecounts_gene_id,
                    "cds_length_bp": cds_length_bp,
                    "raw_count": raw_count,
                    "gene_length_bp": query_length_aa,
                    "library_size": library_size,
                    "rpk": f"{gus_rpk:.6f}",
                    "cpm": f"{cpm:.6f}",
                    "tpm_genome_wide": (
                        f"{(genome_wide_rpk / genome_rpk_total * 1_000_000.0):.6f}"
                    ),
                }
            )
            enriched.append(item)
            gus_rpk_values.append(gus_rpk)

        gus_rpk_sum = sum(gus_rpk_values)
        for item, gus_rpk in zip(enriched, gus_rpk_values):
            item["tpm_gus_db"] = (
                f"{(gus_rpk / gus_rpk_sum * 1_000_000.0) if gus_rpk_sum > 0 else 0.0:.6f}"
            )
            writer.writerow(item)


if __name__ == "__main__":
    main()
