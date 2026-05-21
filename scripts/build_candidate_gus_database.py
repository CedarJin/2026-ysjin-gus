#!/usr/bin/env python3
"""Build the candidate GUS database from Table S3 in the paper supplement."""

from __future__ import annotations

import argparse
import csv
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path


NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
AA_RE = re.compile(r"^[A-Z]+$")


def _paragraph_texts(docx_path: Path) -> list[str]:
    with zipfile.ZipFile(docx_path) as zf:
        root = ET.fromstring(zf.read("word/document.xml"))
    texts: list[str] = []
    for para in root.findall(".//w:p", NS):
        text = "".join(node.text or "" for node in para.findall(".//w:t", NS))
        text = re.sub(r"\s+", " ", text).strip()
        if text:
            texts.append(text)
    return texts


def _table_rows(docx_path: Path) -> list[list[str]]:
    with zipfile.ZipFile(docx_path) as zf:
        root = ET.fromstring(zf.read("word/document.xml"))
    tables = root.findall(".//w:tbl", NS)
    if not tables:
        raise SystemExit(f"No tables found in {docx_path}")
    rows: list[list[str]] = []
    for row in tables[0].findall("./w:tr", NS):
        cells = []
        for cell in row.findall("./w:tc", NS):
            texts = [t.text or "" for t in cell.findall(".//w:t", NS)]
            cells.append(" ".join("".join(texts).split()))
        if cells:
            rows.append(cells)
    return rows


def _normalize_source(name: str) -> str:
    name = name.replace(". ", ".")
    name = re.sub(r"\s*\([^)]*\)", "", name)
    name = re.sub(r"\s+GUS protein$", "", name).strip()
    return name


def parse_docx(docx_path: Path) -> list[dict[str, str]]:
    table_rows = _table_rows(docx_path)
    if len(table_rows) != 19:
        raise SystemExit(f"Expected 19 table rows including header, found {len(table_rows)}")

    paragraphs = _paragraph_texts(docx_path)
    try:
        seq_start = paragraphs.index("* Amino acid sequences of the GUS isoforms are shown below:")
    except ValueError as exc:
        raise SystemExit("Could not find amino-acid sequence section in Table S3.docx") from exc

    header_rows = table_rows[1:]
    seq_paragraphs = paragraphs[seq_start + 1 :]

    records: list[dict[str, str]] = []
    i = 0
    while i < len(seq_paragraphs):
        header = seq_paragraphs[i]
        if not header.endswith("GUS protein"):
            i += 1
            continue
        if i + 1 >= len(seq_paragraphs):
            raise SystemExit(f"Missing sequence after header: {header}")
        sequence = seq_paragraphs[i + 1].replace(" ", "").upper()
        if not AA_RE.match(sequence):
            raise SystemExit(f"Invalid sequence found after header {header}")

        meta = header_rows[len(records)]
        clade = meta[1]
        protein_id = meta[3]
        reference_id = meta[4]
        source = _normalize_source(header)
        if _normalize_source(meta[0]) != source:
            raise SystemExit(
                f"Header mismatch at record {len(records)+1}: table={meta[0]!r}, text={header!r}"
            )
        records.append(
            {
                "reference_id": reference_id,
                "protein_id": protein_id,
                "clade": clade,
                "source": source,
                "sequence": sequence,
                "length": str(len(sequence)),
            }
        )
        i += 2

    if len(records) != 18:
        raise SystemExit(f"Expected 18 candidate GUS records, found {len(records)}")
    return records


def write_outputs(records: list[dict[str, str]], fasta_path: Path, metadata_path: Path) -> None:
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)

    with fasta_path.open("w") as fasta:
        for record in records:
            fasta.write(
                f">{record['reference_id']} clade={record['clade']} source={record['source']}\n"
            )
            seq = record["sequence"]
            for i in range(0, len(seq), 80):
                fasta.write(seq[i : i + 80] + "\n")

    fields = ["reference_id", "protein_id", "clade", "source", "length"]
    with metadata_path.open("w") as metadata:
        writer = csv.DictWriter(metadata, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows({k: r[k] for k in fields} for r in records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    args = parser.parse_args()

    records = parse_docx(args.input)
    write_outputs(records, args.fasta, args.metadata)
    print(f"Extracted {len(records)} candidate records to {args.fasta}")


if __name__ == "__main__":
    main()
