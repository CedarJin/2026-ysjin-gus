#!/usr/bin/env python3
"""Extract GUS protein FASTA records and metadata from the paper PDF/text."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


PAGE_MARKER = re.compile(r"^--\s+\d+\s+of\s+\d+\s+--$")
HEADER = re.compile(r"^>(?P<protein_id>[^-\s]+)-(?P<clade>[A-Za-z0-9]+)\s*(?P<source>.*)$")
AA_LINE = re.compile(r"^[A-Za-zxX]+$")


def read_input(path: Path) -> str:
    if path.suffix.lower() != ".pdf":
        return path.read_text()

    try:
        from pypdf import PdfReader
    except ImportError as exc:
        raise SystemExit(
            "pypdf is required for PDF extraction. Install it in the analysis "
            "environment, or pass a text file converted from the PDF."
        ) from exc

    reader = PdfReader(str(path))
    return "\n".join(page.extract_text() or "" for page in reader.pages)


def parse_records(text: str) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    seq_parts: list[str] = []

    def flush() -> None:
        nonlocal current, seq_parts
        if current is None:
            return
        sequence = "".join(seq_parts).upper().replace("X", "X")
        if not sequence:
            raise ValueError(f"Empty sequence for {current['reference_id']}")
        current["sequence"] = sequence
        current["length"] = str(len(sequence))
        records.append(current)
        current = None
        seq_parts = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or PAGE_MARKER.match(line):
            continue
        if line.startswith(">"):
            flush()
            match = HEADER.match(line)
            if not match:
                raise ValueError(f"Could not parse FASTA header: {line}")
            protein_id = match.group("protein_id")
            clade = match.group("clade")
            source = " ".join(match.group("source").split())
            current = {
                "reference_id": f"{protein_id}-{clade}",
                "protein_id": protein_id,
                "clade": clade,
                "source": source,
            }
            continue
        if current is None:
            continue
        cleaned = re.sub(r"[^A-Za-zxX]", "", line)
        if cleaned and AA_LINE.match(cleaned):
            seq_parts.append(cleaned)

    flush()
    return records


def write_outputs(records: list[dict[str, str]], fasta_path: Path, metadata_path: Path) -> None:
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)

    with fasta_path.open("w") as fasta:
        for record in records:
            fasta.write(
                f">{record['reference_id']} clade={record['clade']} "
                f"source={record['source']}\n"
            )
            seq = record["sequence"]
            for i in range(0, len(seq), 80):
                fasta.write(seq[i : i + 80] + "\n")

    fields = ["reference_id", "protein_id", "clade", "source", "length"]
    with metadata_path.open("w") as metadata:
        metadata.write("\t".join(fields) + "\n")
        for record in records:
            metadata.write("\t".join(record[field] for field in fields) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--expected-count", type=int, default=279)
    args = parser.parse_args()

    records = parse_records(read_input(args.input))
    if len(records) != args.expected_count:
        raise SystemExit(
            f"Extracted {len(records)} records, expected {args.expected_count}. "
            "Inspect the PDF text extraction before continuing."
        )
    write_outputs(records, args.fasta, args.metadata)
    print(f"Extracted {len(records)} records to {args.fasta}")


if __name__ == "__main__":
    main()
