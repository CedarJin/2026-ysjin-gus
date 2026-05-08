#!/usr/bin/env python3
"""Filter FASTA contigs by the len= value in the header or sequence length."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


HEADER_LEN = re.compile(r"\blen=(\d+)")


def keep_record(header: str, sequence: str, min_length: int) -> bool:
    match = HEADER_LEN.search(header)
    length = int(match.group(1)) if match else len(sequence)
    return length >= min_length


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--min-length", type=int, default=1000)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    total = 0
    with args.input.open() as in_handle, args.output.open("w") as out_handle:
        header: str | None = None
        seq_parts: list[str] = []

        def flush() -> None:
            nonlocal kept, total, header, seq_parts
            if header is None:
                return
            total += 1
            sequence = "".join(seq_parts)
            if keep_record(header, sequence, args.min_length):
                kept += 1
                out_handle.write(header + "\n")
                for i in range(0, len(sequence), 80):
                    out_handle.write(sequence[i : i + 80] + "\n")
            header = None
            seq_parts = []

        for raw_line in in_handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line
            else:
                seq_parts.append(line)
        flush()

    print(f"Kept {kept} of {total} contigs with length >= {args.min_length}")


if __name__ == "__main__":
    main()
