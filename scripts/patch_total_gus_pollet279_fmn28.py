#!/usr/bin/env python3
"""Restrict total GUS merged tables to Pollet 279 + Simpson Loop_Class==FMN (28).

Does not re-run DIAMOND. Filters existing all-subject matrices.

- tpm_genome_wide / of_genome_percent: drop non-keep rows (denominator is the
  full gene catalog, so values stay as-is).
- tpm_gus_db / within_gus_db_percent: drop non-keep rows, then renormalize
  each sample to 1e6 or 100% over the remaining references.
"""

from __future__ import annotations

import shutil
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ABUND_DIR = ROOT / "results/summary/pident95_cov90/abundance_all_subject"
BACKUP_DIR = ABUND_DIR / "_backup_before_pollet279_fmn28"
POLLET_META = ROOT / "results/reference/gus_279_metadata.tsv"
FMN_XLSX = ROOT / "reference_from_paper/Supplemental_Table_2.xlsx"
SRA = ROOT / "reference_from_paper/SraRunTable.csv"

TABLES = [
    "total_gus_abundance_gus_db_all_subject.tsv",
    "total_gus_abundance_gus_genome_all_subject.tsv",
    "total_gus_abundance_gus_genome_all_subject_sample_id.tsv",
    "total_gus_abundance_of_genome_percent_all_subject.tsv",
    "total_gus_abundance_within_gus_db_percent_all_subject.tsv",
]


def simpson_fmn_ids(xlsx: Path) -> list[str]:
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    with zipfile.ZipFile(xlsx) as zf:
        ss: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in root.findall(f"{ns}si"):
                ss.append("".join(t.text or "" for t in si.iter(f"{ns}t")))
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rid_to_target = {r.attrib["Id"]: r.attrib["Target"] for r in rels}
        wb = ET.fromstring(zf.read("xl/workbook.xml"))
        target = None
        for sheet in wb.find(f"{ns}sheets"):
            if sheet.attrib.get("name") == "GUS_Summary":
                rid = sheet.attrib[
                    "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
                ]
                target = "xl/" + rid_to_target[rid].lstrip("/")
                break
        if target is None:
            raise SystemExit("GUS_Summary sheet not found")
        grid: dict[tuple[int, int], str] = {}

        def col_row(ref: str) -> tuple[int, int]:
            col = "".join(c for c in ref if c.isalpha())
            row = int("".join(c for c in ref if c.isdigit()))
            n = 0
            for ch in col:
                n = n * 26 + (ord(ch) - 64)
            return n, row

        def cell_val(c: ET.Element) -> str | None:
            t = c.attrib.get("t")
            v = c.find(f"{ns}v")
            if v is None or v.text is None:
                return None
            if t == "s":
                return ss[int(v.text)]
            return v.text

        root = ET.fromstring(zf.read(target))
        max_col = 0
        for row in root.findall(f".//{ns}sheetData/{ns}row"):
            for c in row.findall(f"{ns}c"):
                col, ridx = col_row(c.attrib["r"])
                grid[(col, ridx)] = cell_val(c) or ""
                max_col = max(max_col, col)

    header = {grid.get((c, 1)): c for c in range(1, max_col + 1) if grid.get((c, 1))}
    id_c = header["Gene_ID"]
    loop_c = header["Loop_Class"]
    ids: list[str] = []
    r = 2
    while (id_c, r) in grid:
        gid = grid.get((id_c, r), "")
        loop = grid.get((loop_c, r), "")
        if gid and loop == "FMN":
            ids.append(f"{gid}-FMN")
        r += 1
    return ids


def run_to_sample() -> dict[str, str]:
    sra = pd.read_csv(SRA)
    return {
        str(run): f"S{int(name):02d}"
        for run, name in zip(sra["Run"], sra["Sample Name"])
    }


def keep_ids() -> set[str]:
    pollet = pd.read_csv(POLLET_META, sep="\t")
    fmn28 = simpson_fmn_ids(FMN_XLSX)
    if len(fmn28) != 28:
        raise SystemExit(f"Expected 28 Simpson FMN IDs, got {len(fmn28)}")
    return set(pollet["reference_id"]) | set(fmn28)


def read_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if df.columns[0] != "sseqid":
        raise SystemExit(f"Unexpected first column in {path}: {df.columns[0]}")
    return df.set_index("sseqid")


def write_matrix(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", float_format="%.6f")


def renormalize(df: pd.DataFrame, target_sum: float) -> pd.DataFrame:
    col_sums = df.sum(axis=0).replace(0, pd.NA)
    return df.divide(col_sums, axis=1).fillna(0.0) * target_sum


def main() -> None:
    keep = keep_ids()
    BACKUP_DIR.mkdir(parents=True, exist_ok=True)

    for name in TABLES:
        src = ABUND_DIR / name
        if not src.exists():
            raise SystemExit(f"Missing {src}")
        shutil.copy2(src, BACKUP_DIR / name)

    genome = read_matrix(ABUND_DIR / "total_gus_abundance_gus_genome_all_subject.tsv")
    gus_db = read_matrix(ABUND_DIR / "total_gus_abundance_gus_db_all_subject.tsv")
    of_genome = read_matrix(ABUND_DIR / "total_gus_abundance_of_genome_percent_all_subject.tsv")
    within = read_matrix(ABUND_DIR / "total_gus_abundance_within_gus_db_percent_all_subject.tsv")
    genome_sid = read_matrix(
        ABUND_DIR / "total_gus_abundance_gus_genome_all_subject_sample_id.tsv"
    )

    n_before = genome.shape[0]
    keep_index = [i for i in genome.index if i in keep]
    n_dropped = n_before - len(keep_index)
    n_fmn_kept = sum(i.endswith("-FMN") for i in keep_index)
    n_pollet_kept = len(keep_index) - n_fmn_kept

    genome = genome.loc[keep_index]
    gus_db = gus_db.loc[keep_index]
    of_genome = of_genome.loc[keep_index]
    within = within.loc[keep_index]
    genome_sid = genome_sid.loc[keep_index]

    gus_db = renormalize(gus_db, 1_000_000.0)
    within = renormalize(within, 100.0)

    mapping = run_to_sample()
    expected_sid = [mapping[c] for c in genome.columns]
    if list(genome_sid.columns) != expected_sid:
        genome_sid = genome.copy()
        genome_sid.columns = expected_sid

    write_matrix(gus_db, ABUND_DIR / "total_gus_abundance_gus_db_all_subject.tsv")
    write_matrix(genome, ABUND_DIR / "total_gus_abundance_gus_genome_all_subject.tsv")
    write_matrix(
        genome_sid,
        ABUND_DIR / "total_gus_abundance_gus_genome_all_subject_sample_id.tsv",
    )
    write_matrix(
        of_genome,
        ABUND_DIR / "total_gus_abundance_of_genome_percent_all_subject.tsv",
    )
    write_matrix(
        within,
        ABUND_DIR / "total_gus_abundance_within_gus_db_percent_all_subject.tsv",
    )

    print(f"Keep set: {len(keep)} (279 Pollet + 28 Simpson FMN)")
    print(
        f"Matrix rows: {n_before} -> {len(keep_index)} "
        f"(dropped {n_dropped}; kept {n_pollet_kept} Pollet + {n_fmn_kept} FMN)"
    )
    print("tpm_gus_db column sums:", gus_db.sum().round(3).tolist()[:3], "...")
    print("within-db percent column sums:", within.sum().round(6).tolist()[:3], "...")
    print("Backups in", BACKUP_DIR)


if __name__ == "__main__":
    main()
