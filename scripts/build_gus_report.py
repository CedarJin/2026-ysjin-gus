#!/usr/bin/env python3
"""Build an English report for the strict pident95/cov90 GUS analysis.

The current container does not provide a usable R runtime, so this script
generates the report artifacts directly from the existing TSV outputs using
only the Python standard library plus pandas/numpy.
"""

from __future__ import annotations

import argparse
import html
import math
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


STANDARD_CLADE_ORDER = [
    "FMN",
    "L1",
    "L2",
    "NL",
    "NC",
    "hybrid",
    "m1,2",
    "mL1",
    "mL2",
]

STANDARD_CLADE_MAP = {
    "mL1mL2": "m1,2",
    "mL2Bacteroides": "mL2",
    "mL2Parabacteroides": "mL2",
}

CLADE_COLORS = {
    "FMN": "#1b9e77",
    "L1": "#d95f02",
    "L2": "#7570b3",
    "NL": "#e7298a",
    "NC": "#66a61e",
    "hybrid": "#e6ab02",
    "m1,2": "#a6761d",
    "mL1": "#1f78b4",
    "mL2": "#b2df8a",
    "Other": "#999999",
    "unclear": "#666666",
}

BG = "#ffffff"
TEXT = "#1f2933"
GRID = "#d9e2ec"
AXIS = "#52606d"
BOX = "#cbd2d9"
BOX_FILL = "#f1f5f9"
HIST = "#4c78a8"
HIST_2 = "#f58518"
HEAT_LOW = "#f7fbff"
HEAT_HIGH = "#08306b"


def clean_str(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def standardize_clade(value: object) -> str:
    value = clean_str(value)
    return STANDARD_CLADE_MAP.get(value, value)


def taxonomy_label(row: pd.Series) -> str:
    genus = clean_str(row.get("genus"))
    phylum = clean_str(row.get("phylum"))
    if genus:
        return genus
    if phylum:
        return phylum
    return "unclear"


def safe_float(value: object) -> float:
    if pd.isna(value):
        return 0.0
    return float(value)


def subject_id_from_path(path: Path) -> str:
    return path.name.split(".min1000.", 1)[0]


def load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index.name = "sseqid"
    return df


def load_metadata(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "reference_id" not in df.columns:
        raise SystemExit(f"Missing reference_id in {path}")
    df = df.copy()
    df["clade_std"] = df["clade"].map(standardize_clade)
    df["taxon_label"] = df.apply(taxonomy_label, axis=1)
    return df


def detect_presence(matrix: pd.DataFrame) -> pd.Series:
    return (matrix > 0).sum(axis=1)


def series_from_row_sums(df: pd.DataFrame) -> pd.Series:
    return df.sum(axis=1)


def subject_order_from_total_burden(total_burden: pd.Series) -> list[str]:
    return total_burden.sort_values(ascending=False).index.tolist()


def percentile(values: Iterable[float], q: float) -> float:
    arr = np.asarray(list(values), dtype=float)
    if arr.size == 0:
        return 0.0
    return float(np.quantile(arr, q))


def compute_prevalence_summary(matrix: pd.DataFrame) -> dict[str, object]:
    prevalence = detect_presence(matrix)
    return {
        "unique_genes": int(matrix.shape[0]),
        "prevalence": prevalence,
        "q25": percentile(prevalence, 0.25),
        "median": percentile(prevalence, 0.50),
        "q75": percentile(prevalence, 0.75),
        "core_genes": int((prevalence >= 9).sum()),
        "universal_genes": int((prevalence == 18).sum()),
    }


def subject_burden_table(
    total_of_genome: pd.DataFrame,
    candidate_of_genome: pd.DataFrame,
    total_within: pd.DataFrame,
    candidate_within: pd.DataFrame,
) -> pd.DataFrame:
    total_order = subject_order_from_total_burden(total_of_genome.sum(axis=0))
    data = []
    for subject in total_order:
        data.append(
            {
                "subject": subject,
                "total_gus_of_genome_percent": float(total_of_genome[subject].sum()),
                "candidate_gus_of_genome_percent": float(candidate_of_genome[subject].sum()),
                "total_gus_detected_count": int((total_within[subject] > 0).sum()),
                "candidate_gus_detected_count": int((candidate_within[subject] > 0).sum()),
            }
        )
    return pd.DataFrame(data)


def clade_abundance_tables(
    abundance_within: pd.DataFrame, meta: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta2 = meta.set_index("reference_id").loc[abundance_within.index].copy()
    clades = meta2["clade_std"].tolist()
    clade_table = abundance_within.groupby(clades).sum().T
    clade_table = clade_table.reindex(columns=STANDARD_CLADE_ORDER, fill_value=0.0)
    clade_table.index.name = "subject"
    clade_table.columns.name = "clade"
    clade_table = clade_table.fillna(0.0)
    subject_totals = clade_table.sum(axis=1)
    pct_table = clade_table.divide(subject_totals.replace(0, np.nan), axis=0) * 100.0
    pct_table = pct_table.fillna(0.0)
    return clade_table, pct_table


def taxonomy_resolution_summary(meta: pd.DataFrame, detected_ids: list[str]) -> pd.DataFrame:
    subset = meta.set_index("reference_id").loc[detected_ids].copy()
    categories = []
    for _, row in subset.iterrows():
        genus = clean_str(row.get("genus"))
        phylum = clean_str(row.get("phylum"))
        if genus:
            categories.append("genus")
        elif phylum:
            categories.append("phylum")
        else:
            categories.append("unclear")
    counts = Counter(categories)
    total = len(detected_ids)
    rows = []
    for category in ["genus", "phylum", "unclear"]:
        count = counts.get(category, 0)
        rows.append(
            {
                "category": category,
                "count": count,
                "percent": (count / total * 100.0) if total else 0.0,
            }
        )
    return pd.DataFrame(rows)


def add_grouped_bar(
    parts: list[str],
    left: float,
    top: float,
    width: float,
    height: float,
    labels: list[str],
    series: list[tuple[str, pd.Series, str]],
    x_max: float,
    x_label: str,
    legend_x: float | None = None,
) -> None:
    """Append a grouped horizontal bar chart to parts."""

    n = len(labels)
    band = height / max(n, 1)
    label_w = min(max(len(max(labels, key=len, default="")) * 6.2, 120), 200)
    plot_left = left + label_w
    plot_width = width - label_w - 24
    tick_count = 5
    tick_vals = np.linspace(0, x_max, tick_count)

    parts.append(
        f'<rect x="{left}" y="{top}" width="{width}" height="{height}" fill="{BG}" />'
    )
    # grid and ticks
    for t in tick_vals:
        x = plot_left + (t / x_max * plot_width if x_max else 0)
        parts.append(
            f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + height}" stroke="{GRID}" stroke-width="1" />'
        )
        parts.append(
            f'<line x1="{x:.2f}" y1="{top + height}" x2="{x:.2f}" y2="{top + height + 4}" stroke="{AXIS}" stroke-width="1" />'
        )
        parts.append(
            f'<text x="{x:.2f}" y="{top + height + 18}" font-size="11" fill="{AXIS}" text-anchor="middle">{t:.2f}</text>'
        )
    # axes
    parts.append(
        f'<line x1="{plot_left}" y1="{top}" x2="{plot_left}" y2="{top + height}" stroke="{AXIS}" stroke-width="1.2" />'
    )
    parts.append(
        f'<line x1="{plot_left}" y1="{top + height}" x2="{plot_left + plot_width}" y2="{top + height}" stroke="{AXIS}" stroke-width="1.2" />'
    )
    parts.append(
        f'<text x="{plot_left + plot_width / 2:.2f}" y="{top + height + 38}" font-size="12" fill="{TEXT}" text-anchor="middle">{html.escape(x_label)}</text>'
    )

    bar_series = len(series)
    bar_h = min(10.0, band * 0.32)
    gap = 2.0 if bar_series > 1 else 0.0
    total_h = bar_series * bar_h + (bar_series - 1) * gap

    for i, label in enumerate(labels):
        cy = top + band * (i + 0.5)
        parts.append(
            f'<text x="{left + label_w - 8}" y="{cy + 4:.2f}" font-size="12" fill="{TEXT}" text-anchor="end">{html.escape(label)}</text>'
        )
        y0 = cy - total_h / 2
        for j, (name, values, color) in enumerate(series):
            v = float(values.get(label, 0.0))
            bar_w = plot_width * (v / x_max if x_max else 0)
            y = y0 + j * (bar_h + gap)
            parts.append(
                f'<rect x="{plot_left}" y="{y:.2f}" width="{bar_w:.2f}" height="{bar_h:.2f}" rx="2" ry="2" fill="{color}" />'
            )
            if bar_w > 22:
                parts.append(
                    f'<text x="{plot_left + bar_w - 4:.2f}" y="{y + bar_h - 2:.2f}" font-size="10" fill="#ffffff" text-anchor="end">{v:.3f}</text>'
                )
            elif bar_w > 0:
                parts.append(
                    f'<text x="{plot_left + bar_w + 4:.2f}" y="{y + bar_h - 2:.2f}" font-size="10" fill="{AXIS}" text-anchor="start">{v:.3f}</text>'
                )

    if legend_x is not None:
        ly = top
        parts.append(f'<g transform="translate({legend_x},{ly})">')
        for idx, (name, _, color) in enumerate(series):
            y = idx * 18
            parts.append(
                f'<rect x="0" y="{y}" width="12" height="12" fill="{color}" rx="2" ry="2" />'
            )
            parts.append(
                f'<text x="18" y="{y + 10}" font-size="11" fill="{TEXT}">{html.escape(name)}</text>'
            )
        parts.append("</g>")


def make_bar_chart_horizontal(
    labels: list[str],
    series: list[tuple[str, pd.Series, str]],
    x_label: str,
    width: int,
    height: int,
    title: str,
    subtitle: str | None = None,
) -> str:
    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="100%" height="100%" fill="{BG}" />',
    ]
    parts.append(
        f'<text x="{width/2:.1f}" y="34" text-anchor="middle" font-size="20" font-weight="700" fill="{TEXT}">{html.escape(title)}</text>'
    )
    if subtitle:
        parts.append(
            f'<text x="{width/2:.1f}" y="56" text-anchor="middle" font-size="12" fill="{AXIS}">{html.escape(subtitle)}</text>'
        )
    top = 72
    left = 20
    chart_height = height - top - 36
    x_max = max(float(s.max()) for _, s, _ in series)
    x_max = x_max * 1.12 if x_max > 0 else 1.0
    add_grouped_bar(
        parts,
        left,
        top,
        width - 40,
        chart_height,
        labels,
        series,
        x_max,
        x_label,
        legend_x=width - 150 if len(series) > 1 else None,
    )
    parts.append("</svg>")
    return "\n".join(parts)


def make_histogram_svg(
    counts: pd.Series,
    width: int,
    height: int,
    title: str,
    color: str,
    subtitle: str | None = None,
) -> str:
    labels = [str(i) for i in range(1, 19)]
    max_count = int(counts.max()) if len(counts) else 0
    left = 54
    right = 18
    top = 72
    bottom = 54
    plot_w = width - left - right
    plot_h = height - top - bottom
    bar_w = plot_w / 18.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="100%" height="100%" fill="{BG}" />',
        f'<text x="{width/2:.1f}" y="34" text-anchor="middle" font-size="18" font-weight="700" fill="{TEXT}">{html.escape(title)}</text>',
    ]
    if subtitle:
        parts.append(
            f'<text x="{width/2:.1f}" y="56" text-anchor="middle" font-size="12" fill="{AXIS}">{html.escape(subtitle)}</text>'
        )
    parts.append(
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="{AXIS}" stroke-width="1.2" />'
    )
    parts.append(
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="{AXIS}" stroke-width="1.2" />'
    )
    tick_max = max(max_count, 1)
    for t in np.linspace(0, tick_max, 5):
        y = top + plot_h - (t / tick_max) * plot_h
        parts.append(
            f'<line x1="{left - 3}" y1="{y:.2f}" x2="{left + plot_w}" y2="{y:.2f}" stroke="{GRID}" stroke-width="1" />'
        )
        parts.append(
            f'<text x="{left - 8}" y="{y + 4:.2f}" text-anchor="end" font-size="11" fill="{AXIS}">{int(round(t))}</text>'
        )
    for i, lab in enumerate(labels):
        x = left + i * bar_w + 1
        value = int(counts.get(i + 1, 0))
        h = (value / tick_max) * plot_h if tick_max else 0
        parts.append(
            f'<rect x="{x:.2f}" y="{top + plot_h - h:.2f}" width="{bar_w - 3:.2f}" height="{h:.2f}" fill="{color}" rx="2" ry="2" />'
        )
        parts.append(
            f'<text x="{x + (bar_w - 3)/2:.2f}" y="{top + plot_h + 16}" text-anchor="middle" font-size="10" fill="{AXIS}">{lab}</text>'
        )
        if value > 0:
            parts.append(
                f'<text x="{x + (bar_w - 3)/2:.2f}" y="{top + plot_h - h - 6:.2f}" text-anchor="middle" font-size="10" fill="{TEXT}">{value}</text>'
            )
    parts.append(
        f'<text x="{left + plot_w / 2:.2f}" y="{height - 12}" text-anchor="middle" font-size="12" fill="{TEXT}">Number of subjects with detection</text>'
    )
    parts.append("</svg>")
    return "\n".join(parts)


def hex_to_rgb(hex_color: str) -> tuple[int, int, int]:
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))


def blend_hex(c1: str, c2: str, t: float) -> str:
    r1, g1, b1 = hex_to_rgb(c1)
    r2, g2, b2 = hex_to_rgb(c2)
    r = int(round(r1 + (r2 - r1) * t))
    g = int(round(g1 + (g2 - g1) * t))
    b = int(round(b1 + (b2 - b1) * t))
    return f"#{r:02x}{g:02x}{b:02x}"


def color_scale(value: float, vmin: float, vmax: float) -> str:
    if vmax <= vmin:
        return HEAT_HIGH
    t = min(max((value - vmin) / (vmax - vmin), 0.0), 1.0)
    return blend_hex(HEAT_LOW, HEAT_HIGH, t)


def make_heatmap_svg(
    data: pd.DataFrame,
    row_colors: pd.Series,
    width: int,
    height: int,
    title: str,
    subtitle: str | None = None,
    value_label: str = "Genome-wide TPM",
) -> str:
    rows = list(data.index)
    cols = list(data.columns)
    left = 220
    top = 84
    right = 20
    bottom = 42
    plot_w = width - left - right
    plot_h = height - top - bottom
    cell_w = plot_w / max(len(cols), 1)
    cell_h = plot_h / max(len(rows), 1)
    vmin = float(np.nanmin(data.values)) if data.size else 0.0
    vmax = float(np.nanmax(data.values)) if data.size else 1.0
    if vmax <= vmin:
        vmax = vmin + 1.0

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="100%" height="100%" fill="{BG}" />',
        f'<text x="{width/2:.1f}" y="34" text-anchor="middle" font-size="20" font-weight="700" fill="{TEXT}">{html.escape(title)}</text>',
    ]
    if subtitle:
        parts.append(
            f'<text x="{width/2:.1f}" y="56" text-anchor="middle" font-size="12" fill="{AXIS}">{html.escape(subtitle)}</text>'
        )
    parts.append(
        f'<text x="{left}" y="72" font-size="12" fill="{AXIS}">{html.escape(value_label)}</text>'
    )
    # column labels
    for j, col in enumerate(cols):
        x = left + j * cell_w
        parts.append(
            f'<text x="{x + cell_w/2:.2f}" y="{top - 8}" font-size="10" fill="{AXIS}" text-anchor="end" transform="rotate(-45 {x + cell_w/2:.2f},{top - 8})">{html.escape(col)}</text>'
        )
    # row labels and cells
    for i, row in enumerate(rows):
        y = top + i * cell_h
        clade = clean_str(row_colors.get(row, ""))
        color = CLADE_COLORS.get(clade, "#cccccc")
        parts.append(
            f'<rect x="{left - 18}" y="{y + 1:.2f}" width="10" height="{cell_h - 2:.2f}" fill="{color}" />'
        )
        parts.append(
            f'<text x="{left - 24}" y="{y + cell_h*0.75:.2f}" font-size="10" fill="{TEXT}" text-anchor="end">{html.escape(str(row))}</text>'
        )
        for j, col in enumerate(cols):
            val = float(data.iloc[i, j])
            x = left + j * cell_w
            fill = color_scale(val, vmin, vmax)
            parts.append(
                f'<rect x="{x:.2f}" y="{y:.2f}" width="{cell_w:.2f}" height="{cell_h:.2f}" fill="{fill}" stroke="#ffffff" stroke-width="0.5" />'
            )
            if cell_w >= 22 and cell_h >= 18:
                txt = f"{val:.2f}"
                txt_color = "#ffffff" if val > (vmax * 0.55) else TEXT
                parts.append(
                    f'<text x="{x + cell_w/2:.2f}" y="{y + cell_h*0.65:.2f}" font-size="8.5" fill="{txt_color}" text-anchor="middle">{txt}</text>'
                )
    # color bar
    bar_x = width - 40
    bar_y = top
    bar_h = max(plot_h, 120)
    bar_w = 12
    n_grad = 80
    for k in range(n_grad):
        frac = k / (n_grad - 1)
        y = bar_y + (1 - frac) * bar_h
        val = vmin + frac * (vmax - vmin)
        parts.append(
            f'<rect x="{bar_x}" y="{y:.2f}" width="{bar_w}" height="{bar_h / n_grad + 0.5:.2f}" fill="{color_scale(val, vmin, vmax)}" stroke="none" />'
        )
    parts.append(
        f'<rect x="{bar_x}" y="{bar_y}" width="{bar_w}" height="{bar_h}" fill="none" stroke="{AXIS}" stroke-width="0.6" />'
    )
    for frac, lab in [(0, vmin), (0.5, (vmin + vmax) / 2), (1, vmax)]:
        y = bar_y + (1 - frac) * bar_h
        parts.append(
            f'<line x1="{bar_x + bar_w}" y1="{y:.2f}" x2="{bar_x + bar_w + 6}" y2="{y:.2f}" stroke="{AXIS}" stroke-width="1" />'
        )
        parts.append(
            f'<text x="{bar_x + bar_w + 8}" y="{y + 3:.2f}" font-size="10" fill="{AXIS}">{lab:.1f}</text>'
        )
    parts.append("</svg>")
    return "\n".join(parts)


def boxplot_stats(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {"q1": 0, "median": 0, "q3": 0, "low": 0, "high": 0}
    q1 = float(np.quantile(values, 0.25))
    med = float(np.quantile(values, 0.50))
    q3 = float(np.quantile(values, 0.75))
    iqr = q3 - q1
    low = float(np.min(values[values >= q1 - 1.5 * iqr])) if values.size else 0.0
    high = float(np.max(values[values <= q3 + 1.5 * iqr])) if values.size else 0.0
    return {"q1": q1, "median": med, "q3": q3, "low": low, "high": high}


def make_boxplot_svg(
    table: pd.DataFrame,
    subject_colors: dict[str, str],
    width: int,
    height: int,
    title: str,
    subtitle: str | None = None,
) -> str:
    labels = list(table.columns)
    left = 60
    right = 40
    top = 84
    bottom = 70
    plot_w = width - left - right
    plot_h = height - top - bottom
    n = len(labels)
    band = plot_w / max(n, 1)
    y_max = max(float(table.max().max()) if table.size else 0.0, 1.0)
    y_max = y_max * 1.1

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="100%" height="100%" fill="{BG}" />',
        f'<text x="{width/2:.1f}" y="34" text-anchor="middle" font-size="20" font-weight="700" fill="{TEXT}">{html.escape(title)}</text>',
    ]
    if subtitle:
        parts.append(
            f'<text x="{width/2:.1f}" y="56" text-anchor="middle" font-size="12" fill="{AXIS}">{html.escape(subtitle)}</text>'
        )
    # axes and grid
    parts.append(
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="{AXIS}" stroke-width="1.2" />'
    )
    parts.append(
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="{AXIS}" stroke-width="1.2" />'
    )
    for t in np.linspace(0, y_max, 5):
        y = top + plot_h - (t / y_max * plot_h if y_max else 0)
        parts.append(
            f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" y2="{y:.2f}" stroke="{GRID}" stroke-width="1" />'
        )
        parts.append(
            f'<text x="{left - 8}" y="{y + 4:.2f}" text-anchor="end" font-size="10" fill="{AXIS}">{t:.1f}</text>'
        )
    parts.append(
        f'<text x="{20}" y="{top + plot_h / 2:.2f}" font-size="12" fill="{TEXT}" transform="rotate(-90 20,{top + plot_h / 2:.2f})">Abundance (% of GUS database)</text>'
    )

    # boxes and points
    for i, label in enumerate(labels):
        vals = table[label].to_numpy(dtype=float)
        stats = boxplot_stats(vals)
        cx = left + band * (i + 0.5)
        box_w = min(26.0, band * 0.48)
        # whiskers
        y_low = top + plot_h - (stats["low"] / y_max * plot_h if y_max else 0)
        y_high = top + plot_h - (stats["high"] / y_max * plot_h if y_max else 0)
        y_q1 = top + plot_h - (stats["q1"] / y_max * plot_h if y_max else 0)
        y_q3 = top + plot_h - (stats["q3"] / y_max * plot_h if y_max else 0)
        y_med = top + plot_h - (stats["median"] / y_max * plot_h if y_max else 0)
        parts.append(
            f'<line x1="{cx:.2f}" y1="{y_low:.2f}" x2="{cx:.2f}" y2="{y_high:.2f}" stroke="{AXIS}" stroke-width="1.2" />'
        )
        parts.append(
            f'<line x1="{cx - box_w/3:.2f}" y1="{y_low:.2f}" x2="{cx + box_w/3:.2f}" y2="{y_low:.2f}" stroke="{AXIS}" stroke-width="1.2" />'
        )
        parts.append(
            f'<line x1="{cx - box_w/3:.2f}" y1="{y_high:.2f}" x2="{cx + box_w/3:.2f}" y2="{y_high:.2f}" stroke="{AXIS}" stroke-width="1.2" />'
        )
        parts.append(
            f'<rect x="{cx - box_w/2:.2f}" y="{y_q3:.2f}" width="{box_w:.2f}" height="{max(y_q1 - y_q3, 0.2):.2f}" fill="{BOX_FILL}" stroke="{BOX}" stroke-width="1" rx="2" ry="2" />'
        )
        parts.append(
            f'<line x1="{cx - box_w/2:.2f}" y1="{y_med:.2f}" x2="{cx + box_w/2:.2f}" y2="{y_med:.2f}" stroke="{TEXT}" stroke-width="1.5" />'
        )
        # points
        for j, subject in enumerate(table.index):
            v = float(table.loc[subject, label])
            x = cx + math.sin((i + 1) * (j + 1) * 1.7) * (box_w * 0.35)
            y = top + plot_h - (v / y_max * plot_h if y_max else 0)
            color = subject_colors.get(subject, "#444444")
            parts.append(
                f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.2" fill="{color}" fill-opacity="0.9" stroke="#ffffff" stroke-width="0.5" />'
            )
        parts.append(
            f'<text x="{cx:.2f}" y="{top + plot_h + 18:.2f}" text-anchor="middle" font-size="11" fill="{TEXT}">{html.escape(label)}</text>'
        )
    # subject legend
    legend_x = width - 250
    legend_y = top + 10
    parts.append(f'<g transform="translate({legend_x},{legend_y})">')
    parts.append(
        f'<text x="0" y="-6" font-size="12" fill="{TEXT}" font-weight="700">Subjects</text>'
    )
    for idx, subject in enumerate(table.index):
        row = idx // 2
        col = idx % 2
        x = col * 120
        y = row * 16
        parts.append(
            f'<rect x="{x}" y="{y}" width="10" height="10" fill="{subject_colors.get(subject, "#444444")}" rx="2" ry="2" />'
        )
        parts.append(
            f'<text x="{x + 14}" y="{y + 9}" font-size="10" fill="{TEXT}">{html.escape(subject)}</text>'
        )
    parts.append("</g>")
    parts.append("</svg>")
    return "\n".join(parts)


def make_stacked_bar_svg(
    table: pd.DataFrame,
    category_colors: dict[str, str],
    width: int,
    height: int,
    title: str,
    subtitle: str | None = None,
    y_label: str = "Abundance (% of GUS database)",
    legend_title: str = "Clade",
    right: int = 35,
    legend_x: float | None = None,
) -> str:
    subjects = list(table.index)
    categories = list(table.columns)
    left = 70
    top = 84
    bottom = 92
    plot_w = width - left - right
    plot_h = height - top - bottom
    band = plot_w / max(len(subjects), 1)
    y_max = 100.0

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="100%" height="100%" fill="{BG}" />',
        f'<text x="{width/2:.1f}" y="34" text-anchor="middle" font-size="20" font-weight="700" fill="{TEXT}">{html.escape(title)}</text>',
    ]
    if subtitle:
        parts.append(
            f'<text x="{width/2:.1f}" y="56" text-anchor="middle" font-size="12" fill="{AXIS}">{html.escape(subtitle)}</text>'
        )
    parts.append(
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="{AXIS}" stroke-width="1.2" />'
    )
    parts.append(
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="{AXIS}" stroke-width="1.2" />'
    )
    for t in [0, 25, 50, 75, 100]:
        y = top + plot_h - (t / y_max * plot_h)
        parts.append(
            f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" y2="{y:.2f}" stroke="{GRID}" stroke-width="1" />'
        )
        parts.append(
            f'<text x="{left - 8}" y="{y + 4:.2f}" text-anchor="end" font-size="10" fill="{AXIS}">{t}</text>'
        )
    parts.append(
        f'<text x="18" y="{top + plot_h / 2:.2f}" font-size="12" fill="{TEXT}" transform="rotate(-90 18,{top + plot_h / 2:.2f})">{html.escape(y_label)}</text>'
    )
    for i, subject in enumerate(subjects):
        x0 = left + band * i + band * 0.12
        bar_w = band * 0.76
        baseline = top + plot_h
        running = 0.0
        for cat in categories:
            v = float(table.loc[subject, cat])
            h = (v / y_max) * plot_h
            y = baseline - running - h
            color = category_colors.get(cat, "#999999")
            if v > 0:
                parts.append(
                    f'<rect x="{x0:.2f}" y="{y:.2f}" width="{bar_w:.2f}" height="{h:.2f}" fill="{color}" />'
                )
            running += h
        parts.append(
            f'<text x="{x0 + bar_w/2:.2f}" y="{top + plot_h + 18}" text-anchor="end" font-size="10" fill="{AXIS}" transform="rotate(-45 {x0 + bar_w/2:.2f},{top + plot_h + 18})">{html.escape(subject)}</text>'
        )
    # legend
    if legend_x is None:
        legend_x = width - 150
    legend_y = top
    parts.append(f'<g transform="translate({legend_x},{legend_y})">')
    parts.append(
        f'<text x="0" y="-6" font-size="12" fill="{TEXT}" font-weight="700">{html.escape(legend_title)}</text>'
    )
    for idx, cat in enumerate(categories):
        y = idx * 18
        parts.append(
            f'<rect x="0" y="{y}" width="12" height="12" fill="{category_colors.get(cat, "#999999")}" rx="2" ry="2" />'
        )
        parts.append(
            f'<text x="18" y="{y + 10}" font-size="11" fill="{TEXT}">{html.escape(cat)}</text>'
        )
    parts.append("</g>")
    parts.append("</svg>")
    return "\n".join(parts)


def make_taxonomy_resolution_svg(
    summary: pd.DataFrame,
    width: int,
    height: int,
    title: str,
) -> str:
    labels = summary["category"].tolist()
    values = summary["count"].tolist()
    max_v = max(values) if values else 1
    left = 90
    right = 180
    top = 82
    bottom = 48
    plot_w = width - left - right
    plot_h = height - top - bottom
    band = plot_h / max(len(labels), 1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="100%" height="100%" fill="{BG}" />',
        f'<text x="{width/2:.1f}" y="34" text-anchor="middle" font-size="20" font-weight="700" fill="{TEXT}">{html.escape(title)}</text>',
        f'<text x="{width/2:.1f}" y="56" text-anchor="middle" font-size="12" fill="{AXIS}">Detected unique total GUS only</text>',
    ]
    parts.append(
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="{AXIS}" stroke-width="1.2" />'
    )
    parts.append(
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="{AXIS}" stroke-width="1.2" />'
    )
    for t in np.linspace(0, max_v, 5):
        x = left + (t / max_v * plot_w if max_v else 0)
        parts.append(
            f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_h}" stroke="{GRID}" stroke-width="1" />'
        )
        parts.append(
            f'<text x="{x:.2f}" y="{top + plot_h + 16}" text-anchor="middle" font-size="10" fill="{AXIS}">{int(round(t))}</text>'
        )
    for i, row in summary.iterrows():
        y = top + band * i
        label = row["category"]
        count = int(row["count"])
        pct = float(row["percent"])
        parts.append(
            f'<text x="{left - 8}" y="{y + band*0.65:.2f}" text-anchor="end" font-size="12" fill="{TEXT}">{html.escape(label)}</text>'
        )
        bar_w = (count / max_v) * plot_w if max_v else 0
        color = {
            "genus": "#1f78b4",
            "phylum": "#33a02c",
            "unclear": "#e31a1c",
        }.get(label, "#666666")
        parts.append(
            f'<rect x="{left}" y="{y + 3:.2f}" width="{bar_w:.2f}" height="{band - 8:.2f}" fill="{color}" rx="2" ry="2" />'
        )
        parts.append(
            f'<text x="{left + bar_w + 8:.2f}" y="{y + band*0.65:.2f}" font-size="12" fill="{TEXT}">{count} ({pct:.1f}%)</text>'
        )
    parts.append(
        f'<text x="{left + plot_w / 2:.2f}" y="{height - 12}" text-anchor="middle" font-size="12" fill="{TEXT}">Number of detected genes</text>'
    )
    parts.append("</svg>")
    return "\n".join(parts)


def make_stacked_taxonomy_svg(
    table: pd.DataFrame,
    colors: dict[str, str],
    width: int,
    height: int,
    title: str,
    subtitle: str,
    legend_title: str = "Taxonomy label",
) -> str:
    clades = list(table.index)
    labels = list(table.columns)
    left = 72
    right = 260
    top = 84
    bottom = 56
    plot_w = width - left - right
    plot_h = height - top - bottom
    band = plot_h / max(len(clades), 1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="100%" height="100%" fill="{BG}" />',
        f'<text x="{width/2:.1f}" y="34" text-anchor="middle" font-size="20" font-weight="700" fill="{TEXT}">{html.escape(title)}</text>',
        f'<text x="{width/2:.1f}" y="56" text-anchor="middle" font-size="12" fill="{AXIS}">{html.escape(subtitle)}</text>',
    ]
    parts.append(
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="{AXIS}" stroke-width="1.2" />'
    )
    parts.append(
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="{AXIS}" stroke-width="1.2" />'
    )
    for t in [0, 25, 50, 75, 100]:
        y = top + plot_h - (t / 100.0 * plot_h)
        parts.append(
            f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" y2="{y:.2f}" stroke="{GRID}" stroke-width="1" />'
        )
        parts.append(
            f'<text x="{left - 8}" y="{y + 4:.2f}" text-anchor="end" font-size="10" fill="{AXIS}">{t}</text>'
        )
    parts.append(
        f'<text x="18" y="{top + plot_h / 2:.2f}" font-size="12" fill="{TEXT}" transform="rotate(-90 18,{top + plot_h / 2:.2f})">Percent of clade abundance</text>'
    )

    for i, clade in enumerate(clades):
        y0 = top + band * i
        parts.append(
            f'<text x="{left - 10}" y="{y0 + band*0.65:.2f}" text-anchor="end" font-size="12" fill="{TEXT}">{html.escape(clade)}</text>'
        )
        x0 = left
        running = 0.0
        for label in labels:
            v = float(table.loc[clade, label])
            if v <= 0:
                continue
            w = (v / 100.0) * plot_w
            color = colors.get(label, "#999999")
            parts.append(
                f'<rect x="{x0 + running:.2f}" y="{y0 + 3:.2f}" width="{w:.2f}" height="{band - 8:.2f}" fill="{color}" />'
            )
            if w > 26:
                parts.append(
                    f'<text x="{x0 + running + w/2:.2f}" y="{y0 + band*0.65:.2f}" text-anchor="middle" font-size="9" fill="#ffffff">{html.escape(label)}</text>'
                )
            running += w

    legend_x = width - right + 20
    legend_y = top
    parts.append(f'<g transform="translate({legend_x},{legend_y})">')
    parts.append(f'<text x="0" y="-6" font-size="12" fill="{TEXT}" font-weight="700">{html.escape(legend_title)}</text>')
    for idx, label in enumerate(labels):
        y = idx * 18
        parts.append(
            f'<rect x="0" y="{y}" width="12" height="12" fill="{colors.get(label, "#999999")}" rx="2" ry="2" />'
        )
        parts.append(
            f'<text x="18" y="{y + 10}" font-size="11" fill="{TEXT}">{html.escape(label)}</text>'
        )
    parts.append("</g>")
    parts.append("</svg>")
    return "\n".join(parts)


def make_sankey_svg(
    links: pd.DataFrame,
    left_nodes: list[str],
    right_nodes: list[str],
    node_colors: dict[str, str],
    width: int,
    height: int,
    title: str,
    subtitle: str,
) -> str:
    left = 52
    right = 52
    top = 74
    bottom = 42
    plot_h = height - top - bottom
    node_w = 18
    gap = 10
    x_left = left
    x_right = width - right - node_w
    left_total = links.groupby("left")["weight"].sum().reindex(left_nodes).fillna(0.0)
    right_total = links.groupby("right")["weight"].sum().reindex(right_nodes).fillna(0.0)
    total_weight = float(links["weight"].sum())
    usable_h = plot_h - gap * (max(len(left_nodes), len(right_nodes)) - 1)
    scale = usable_h / total_weight if total_weight > 0 else 1.0
    # heights per node
    left_h = {n: max(8.0, float(left_total[n]) * scale) for n in left_nodes}
    right_h = {n: max(8.0, float(right_total[n]) * scale) for n in right_nodes}
    # positions
    left_pos = {}
    y = top
    for n in left_nodes:
        left_pos[n] = y
        y += left_h[n] + gap
    right_pos = {}
    y = top
    for n in right_nodes:
        right_pos[n] = y
        y += right_h[n] + gap

    # link offsets
    left_offsets = {n: 0.0 for n in left_nodes}
    right_offsets = {n: 0.0 for n in right_nodes}

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="100%" height="100%" fill="{BG}" />',
        f'<text x="{width/2:.1f}" y="34" text-anchor="middle" font-size="20" font-weight="700" fill="{TEXT}">{html.escape(title)}</text>',
        f'<text x="{width/2:.1f}" y="56" text-anchor="middle" font-size="12" fill="{AXIS}">{html.escape(subtitle)}</text>',
    ]

    # links
    for _, row in links.sort_values(["left", "right"]).iterrows():
        l = row["left"]
        r = row["right"]
        w = float(row["weight"])
        if w <= 0:
            continue
        sy = left_pos[l] + left_offsets[l] * scale + (w * scale) / 2
        ty = right_pos[r] + right_offsets[r] * scale + (w * scale) / 2
        left_offsets[l] += w
        right_offsets[r] += w
        thickness = max(1.0, w * scale)
        color = node_colors.get(l, "#999999")
        parts.append(
            f'<path d="M {x_left + node_w} {sy:.2f} C {x_left + 160} {sy:.2f}, {x_right - 160} {ty:.2f}, {x_right} {ty:.2f}" stroke="{color}" stroke-opacity="0.34" stroke-width="{thickness:.2f}" fill="none" stroke-linecap="round" />'
        )

    # nodes
    for n in left_nodes:
        h = left_h[n]
        y0 = left_pos[n]
        parts.append(
            f'<rect x="{x_left}" y="{y0:.2f}" width="{node_w}" height="{h:.2f}" fill="{node_colors.get(n, "#999999")}" stroke="#ffffff" stroke-width="1" rx="2" ry="2" />'
        )
        parts.append(
            f'<text x="{x_left - 10}" y="{y0 + h/2 + 4:.2f}" text-anchor="end" font-size="11" fill="{TEXT}">{html.escape(n)}</text>'
        )
    for n in right_nodes:
        h = right_h[n]
        y0 = right_pos[n]
        parts.append(
            f'<rect x="{x_right}" y="{y0:.2f}" width="{node_w}" height="{h:.2f}" fill="{node_colors.get(n, "#bbbbbb")}" stroke="#ffffff" stroke-width="1" rx="2" ry="2" />'
        )
        parts.append(
            f'<text x="{x_right + node_w + 10}" y="{y0 + h/2 + 4:.2f}" text-anchor="start" font-size="11" fill="{TEXT}">{html.escape(n)}</text>'
        )

    parts.append("</svg>")
    return "\n".join(parts)


def html_table(df: pd.DataFrame, max_rows: int | None = None, classes: str = "data-table") -> str:
    if max_rows is not None and len(df) > max_rows:
        df = df.head(max_rows)
    return df.to_html(index=False, classes=classes, border=0, escape=True, justify="left")


def markdown_table(df: pd.DataFrame) -> str:
    """Render a simple GitHub-flavored markdown table without tabulate."""

    def fmt(value: object) -> str:
        if pd.isna(value):
            return ""
        if isinstance(value, float):
            if abs(value) >= 1000:
                return f"{value:.2f}"
            return f"{value:.4f}".rstrip("0").rstrip(".")
        return str(value)

    cols = list(df.columns)
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join(["---"] * len(cols)) + " |",
    ]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(fmt(row[c]) for c in cols) + " |")
    return "\n".join(lines)


def write_svg(path: Path, svg: str) -> None:
    path.write_text(svg, encoding="utf-8")


def write_tsv(path: Path, df: pd.DataFrame) -> None:
    df.to_csv(path, sep="\t", index=False, float_format="%.6f")


def collapse_taxonomy_labels(labels: list[str], top_n: int = 9) -> tuple[list[str], set[str]]:
    counts = Counter(labels)
    ordered = [lab for lab, _ in counts.most_common() if lab != "unclear"]
    top = ordered[:top_n]
    keep = set(top)
    keep.add("unclear")
    collapsed = [lab if lab in keep else "Other" for lab in labels]
    return collapsed, keep


def build_report(args: argparse.Namespace) -> None:
    out = args.output_dir
    tables_dir = out / "tables"
    figures_dir = out / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    total_of_genome = load_matrix(args.abundance_dir / "total_gus_abundance_of_genome_percent_all_subject.tsv")
    candidate_of_genome = load_matrix(args.abundance_dir / "candidate_gus_abundance_of_genome_percent_all_subject.tsv")
    total_within = load_matrix(args.abundance_dir / "total_gus_abundance_within_gus_db_percent_all_subject.tsv")
    candidate_within = load_matrix(args.abundance_dir / "candidate_gus_abundance_within_gus_db_percent_all_subject.tsv")
    total_meta = load_metadata(args.total_meta)
    candidate_meta = load_metadata(args.candidate_meta)

    # Optional active GUS set (first 11 of Table S3). When its merged tables and
    # metadata are present we emit active-only versions of figures 3, 4 and 6.
    active_within_path = args.abundance_dir / "active_gus_abundance_within_gus_db_percent_all_subject.tsv"
    active_within = None
    active_meta = None
    if active_within_path.exists() and args.active_meta.exists():
        active_within = load_matrix(active_within_path)
        active_meta = load_metadata(args.active_meta)

    # subject summaries
    subjects = subject_burden_table(total_of_genome, candidate_of_genome, total_within, candidate_within)
    write_tsv(tables_dir / "subject_burden.tsv", subjects)
    subject_order = subjects["subject"].tolist()

    # prevalence summaries
    total_prev = compute_prevalence_summary(total_of_genome)
    candidate_prev = compute_prevalence_summary(candidate_of_genome)
    prevalence_summary = pd.DataFrame(
        [
            {
                "dataset": "total",
                **{k: v for k, v in total_prev.items() if k != "prevalence"},
            },
            {
                "dataset": "candidate",
                **{k: v for k, v in candidate_prev.items() if k != "prevalence"},
            },
        ]
    )
    write_tsv(tables_dir / "prevalence_summary.tsv", prevalence_summary)

    total_prev_counts = total_prev["prevalence"].sort_values(ascending=True)
    candidate_prev_counts = candidate_prev["prevalence"].sort_values(ascending=True)
    core_ids = total_prev_counts[total_prev_counts >= 9].index.tolist()
    core_table = (total_of_genome.loc[core_ids].T * 10000.0).reindex(subject_order)
    core_table = core_table.loc[:, core_ids]
    core_meta = total_meta.set_index("reference_id").loc[core_ids, ["clade_std"]]
    core_heatmap_data = core_table.T
    core_row_colors = core_meta["clade_std"]
    core_gene_table = pd.DataFrame(
        {
            "reference_id": core_ids,
            "prevalence_subjects": [int(total_prev_counts.loc[i]) for i in core_ids],
            "clade": [core_meta.loc[i, "clade_std"] for i in core_ids],
            "taxon_label": [total_meta.set_index("reference_id").loc[i, "taxon_label"] for i in core_ids],
        }
    ).sort_values(["prevalence_subjects", "reference_id"], ascending=[False, True])
    write_tsv(tables_dir / "core_genes_total.tsv", core_gene_table)

    # taxonomy resolution
    total_detected_ids = total_prev_counts.index.tolist()
    tax_res = taxonomy_resolution_summary(total_meta, total_detected_ids)
    write_tsv(tables_dir / "taxonomy_resolution_total.tsv", tax_res)

    # clade abundance
    total_clade_counts, total_clade_percent = clade_abundance_tables(total_within, total_meta)
    candidate_clade_counts, candidate_clade_percent = clade_abundance_tables(candidate_within, candidate_meta)
    write_tsv(tables_dir / "total_clade_abundance_within_gus_db_percent.tsv", total_clade_percent.reset_index())
    write_tsv(tables_dir / "candidate_clade_abundance_within_gus_db_percent.tsv", candidate_clade_percent.reset_index())

    # taxonomy composition, abundance-weighted
    total_meta_indexed = total_meta.set_index("reference_id").loc[total_within.index].copy()
    total_mean_abundance = total_within.mean(axis=1)
    tax_labels = total_meta_indexed["taxon_label"].tolist()
    collapsed_labels, keep_labels = collapse_taxonomy_labels(tax_labels, top_n=9)
    total_meta_indexed["taxon_collapsed"] = collapsed_labels
    clade_tax_weight = (
        pd.DataFrame(
            {
                "clade": total_meta_indexed["clade_std"].tolist(),
                "taxon": total_meta_indexed["taxon_collapsed"].tolist(),
                "weight": total_mean_abundance.values,
            }
        )
        .groupby(["clade", "taxon"], as_index=False)["weight"]
        .sum()
    )
    clade_tax_pct = clade_tax_weight.copy()
    clade_tax_pct["pct"] = clade_tax_pct.groupby("clade")["weight"].transform(lambda x: x / x.sum() * 100.0 if x.sum() else 0.0)
    write_tsv(tables_dir / "clade_taxonomy_abundance.tsv", clade_tax_pct)

    # sankey links using the same collapsed taxonomy labels
    sankey_links = clade_tax_weight.rename(columns={"clade": "left", "taxon": "right", "weight": "weight"})
    write_tsv(tables_dir / "sankey_links.tsv", sankey_links)

    # figures
    fig1 = make_bar_chart_horizontal(
        subject_order,
        [
            ("total GUS of genome %", subjects.set_index("subject")["total_gus_of_genome_percent"], "#4c78a8"),
            ("candidate GUS of genome %", subjects.set_index("subject")["candidate_gus_of_genome_percent"], "#f58518"),
        ],
        "Genome percent",
        1500,
        720,
        "Subject-level GUS abundance",
        "Per-subject burden across the 18 subjects",
    )
    write_svg(figures_dir / "figure1_subject_burden.svg", fig1)

    fig2_parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1500" height="980" viewBox="0 0 1500 980">',
        f'<rect width="100%" height="100%" fill="{BG}" />',
        f'<text x="750" y="34" text-anchor="middle" font-size="20" font-weight="700" fill="{TEXT}">Detection prevalence and core GUS heatmap</text>',
        f'<text x="750" y="56" text-anchor="middle" font-size="12" fill="{AXIS}">Top row: total GUS; middle row: candidate GUS; bottom: core total GUS with prevalence >= 9/18</text>',
    ]
    fig2_parts.append(make_histogram_svg(total_prev_counts.value_counts().sort_index(), 700, 280, "Total GUS prevalence", HIST, "Number of subjects detected per GUS").split("\n", 1)[1].rsplit("</svg>", 1)[0])
    fig2_parts.append(
        f'<g transform="translate(770,72)">{make_histogram_svg(candidate_prev_counts.value_counts().sort_index(), 690, 280, "Candidate GUS prevalence", HIST_2, "Number of subjects detected per GUS").split("<svg", 1)[1].split("</svg>", 1)[0]}</g>'
    )
    # core heatmap
    heat_svg = make_heatmap_svg(
        core_heatmap_data,
        core_row_colors,
        1460,
        560,
        "Core total GUS heatmap",
        "Rows are total GUS detected in at least 9/18 subjects; values are genome-wide TPM.",
        value_label="Genome-wide TPM",
    )
    fig2_parts.append(f'<g transform="translate(0,340)">{heat_svg.split("<svg",1)[1].split("</svg>",1)[0]}</g>')
    fig2_parts.append("</svg>")
    write_svg(figures_dir / "figure2_prevalence_and_heatmap.svg", "\n".join(fig2_parts))

    box_svg = make_boxplot_svg(
        total_clade_percent.reindex(subject_order),
        {s: f"hsl({int(360 * i / max(len(subject_order),1))},65%,50%)" for i, s in enumerate(subject_order)},
        1540,
        760,
        "Clade abundance distribution",
        "Boxplots show the distribution across subjects; points are colored by subject.",
    )
    write_svg(figures_dir / "figure3_clade_abundance_boxplot.svg", box_svg)

    stacked_svg = make_stacked_bar_svg(
        total_clade_percent.reindex(subject_order),
        CLADE_COLORS,
        1540,
        760,
        "Subject-by-subject clade composition",
        "Stacked bars use within-GUS-database percent and sum to 100 for each subject.",
    )
    write_svg(figures_dir / "figure4_clade_stacked_bars.svg", stacked_svg)

    tax_res_svg = make_taxonomy_resolution_svg(
        tax_res,
        1100,
        420,
        "Taxonomy resolution among detected total GUS",
    )
    write_svg(figures_dir / "figure5_taxonomy_resolution.svg", tax_res_svg)

    clade_tax_table = clade_tax_pct.pivot(index="clade", columns="taxon", values="pct").fillna(0.0)
    clade_tax_table = clade_tax_table.reindex([c for c in STANDARD_CLADE_ORDER if c in clade_tax_table.index], fill_value=0.0)
    tax_labels_sorted = [c for c in clade_tax_table.columns if c != "Other"] + (["Other"] if "Other" in clade_tax_table.columns else [])
    clade_tax_table = clade_tax_table.reindex(columns=tax_labels_sorted, fill_value=0.0)
    tax_color_labels = {lab: f"hsl({int(360 * i / max(len(tax_labels_sorted),1))},60%,52%)" for i, lab in enumerate(tax_labels_sorted)}
    tax_color_labels["Other"] = "#999999"
    tax_color_labels["unclear"] = "#666666"
    clade_tax_svg = make_stacked_taxonomy_svg(
        clade_tax_table,
        tax_color_labels,
        1540,
        860,
        "Abundance-weighted taxonomy composition by clade",
        "Each clade is shown as a 100% stacked bar. Rare taxonomy labels are grouped into Other.",
    )
    write_svg(figures_dir / "figure6_clade_taxonomy_composition.svg", clade_tax_svg)

    # ---- active GUS (first 11 of Table S3): parallel figures 3/4/6 ----
    active_figures: list[tuple[str, str]] = []
    if active_within is not None and active_meta is not None:
        active_clade_counts, active_clade_percent = clade_abundance_tables(
            active_within, active_meta
        )
        write_tsv(
            tables_dir / "active_clade_abundance_within_gus_db_percent.tsv",
            active_clade_percent.reset_index(),
        )

        active_box_svg = make_boxplot_svg(
            active_clade_percent.reindex(subject_order),
            {s: f"hsl({int(360 * i / max(len(subject_order),1))},65%,50%)" for i, s in enumerate(subject_order)},
            1540,
            760,
            "Clade abundance distribution (active GUS)",
            "Active GUS = first 11 of Table S3. Boxplots across subjects; points colored by subject.",
        )
        write_svg(figures_dir / "figure3_clade_abundance_boxplot_active.svg", active_box_svg)

        active_stacked_svg = make_stacked_bar_svg(
            active_clade_percent.reindex(subject_order),
            CLADE_COLORS,
            1540,
            760,
            "Subject-by-subject clade composition (active GUS)",
            "Active GUS within-database percent; bars sum to 100 for each subject.",
            right=170,
            legend_x=1540 - 150,
        )
        write_svg(figures_dir / "figure4_clade_stacked_bars_active.svg", active_stacked_svg)

        active_meta_indexed = active_meta.set_index("reference_id").loc[active_within.index].copy()
        active_mean_abundance = active_within.mean(axis=1)
        active_tax_labels = active_meta_indexed["taxon_label"].tolist()
        active_collapsed, _ = collapse_taxonomy_labels(active_tax_labels, top_n=9)
        active_meta_indexed["taxon_collapsed"] = active_collapsed
        active_clade_tax_weight = (
            pd.DataFrame(
                {
                    "clade": active_meta_indexed["clade_std"].tolist(),
                    "taxon": active_meta_indexed["taxon_collapsed"].tolist(),
                    "weight": active_mean_abundance.values,
                }
            )
            .groupby(["clade", "taxon"], as_index=False)["weight"]
            .sum()
        )
        active_clade_tax_weight["pct"] = active_clade_tax_weight.groupby("clade")["weight"].transform(
            lambda x: x / x.sum() * 100.0 if x.sum() else 0.0
        )
        active_clade_tax_table = active_clade_tax_weight.pivot(
            index="clade", columns="taxon", values="pct"
        ).fillna(0.0)
        active_clade_tax_table = active_clade_tax_table.reindex(
            [c for c in STANDARD_CLADE_ORDER if c in active_clade_tax_table.index], fill_value=0.0
        )
        active_tax_sorted = [c for c in active_clade_tax_table.columns if c != "Other"] + (
            ["Other"] if "Other" in active_clade_tax_table.columns else []
        )
        active_clade_tax_table = active_clade_tax_table.reindex(columns=active_tax_sorted, fill_value=0.0)
        active_tax_colors = {lab: f"hsl({int(360 * i / max(len(active_tax_sorted),1))},60%,52%)" for i, lab in enumerate(active_tax_sorted)}
        active_tax_colors["Other"] = "#999999"
        active_tax_colors["unclear"] = "#666666"
        active_clade_tax_svg = make_stacked_taxonomy_svg(
            active_clade_tax_table,
            active_tax_colors,
            1540,
            860,
            "Abundance-weighted taxonomy composition by clade (active GUS)",
            "Active GUS (first 11 of Table S3). Each clade is a 100% stacked bar.",
        )
        write_svg(figures_dir / "figure6_clade_taxonomy_composition_active.svg", active_clade_tax_svg)

        # Per-subject stacked bars of the individual active GUS proteins
        # (11 references), each subject summing to 100% within the active set.
        active_meta_map = active_meta.set_index("reference_id")
        active_ref_mean = active_within.mean(axis=1)

        def active_ref_sort_key(rid: str) -> tuple[int, float]:
            clade = clean_str(active_meta_map.loc[rid, "clade_std"])
            clade_rank = STANDARD_CLADE_ORDER.index(clade) if clade in STANDARD_CLADE_ORDER else 99
            return (clade_rank, -float(active_ref_mean.get(rid, 0.0)))

        # Include all 11 active references (undetected ones stay at 0 but remain
        # in the legend), ordered by clade then mean abundance.
        active_all_ids = list(active_meta_map.index)
        active_ref_order = sorted(active_all_ids, key=active_ref_sort_key)
        active_ref_labels = [
            clean_str(active_meta_map.loc[rid, "source"])
            for rid in active_ref_order
        ]
        active_composition = (
            active_within.reindex(active_ref_order).T.reindex(subject_order).fillna(0.0)
        )
        active_composition.columns = active_ref_labels
        active_ref_colors = {
            lab: f"hsl({int(360 * i / max(len(active_ref_labels), 1))},62%,52%)"
            for i, lab in enumerate(active_ref_labels)
        }
        active_composition_svg = make_stacked_bar_svg(
            active_composition,
            active_ref_colors,
            1540,
            760,
            "Active GUS composition per subject",
            "Each bar sums to 100% across the 11 active GUS proteins (first 11 of Table S3).",
            y_label="Abundance (% of active GUS set)",
            legend_title="Active GUS",
            right=280,
            legend_x=1540 - 268,
        )
        write_svg(figures_dir / "figure4b_active_gus_composition.svg", active_composition_svg)

        active_figures = [
            (
                "figure3_clade_abundance_boxplot_active.svg",
                "Active GUS clade abundance distribution across subjects.",
            ),
            (
                "figure4_clade_stacked_bars_active.svg",
                "Active GUS subject-by-subject clade composition.",
            ),
            (
                "figure4b_active_gus_composition.svg",
                "Active GUS per-subject composition by individual protein (sums to 100%).",
            ),
            (
                "figure6_clade_taxonomy_composition_active.svg",
                "Active GUS abundance-weighted taxonomy composition by clade.",
            ),
        ]

    sankey_nodes_right = tax_labels_sorted
    sankey_svg = make_sankey_svg(
        sankey_links,
        [c for c in STANDARD_CLADE_ORDER if c in sankey_links["left"].unique()],
        sankey_nodes_right,
        {**CLADE_COLORS, **tax_color_labels},
        1540,
        860,
        "Clade-to-taxonomy Sankey",
        "Link thickness reflects mean within-GUS-database abundance across the 18 subjects.",
    )
    write_svg(figures_dir / "figure7_sankey.svg", sankey_svg)

    # report HTML
    summary_numbers = pd.DataFrame(
        [
            {
                "dataset": "total",
                "unique_genes": total_prev["unique_genes"],
                "q25_prevalence": total_prev["q25"],
                "median_prevalence": total_prev["median"],
                "q75_prevalence": total_prev["q75"],
                "core_genes_ge_9_of_18": total_prev["core_genes"],
                "universal_genes_18_of_18": total_prev["universal_genes"],
            },
            {
                "dataset": "candidate",
                "unique_genes": candidate_prev["unique_genes"],
                "q25_prevalence": candidate_prev["q25"],
                "median_prevalence": candidate_prev["median"],
                "q75_prevalence": candidate_prev["q75"],
                "core_genes_ge_9_of_18": candidate_prev["core_genes"],
                "universal_genes_18_of_18": candidate_prev["universal_genes"],
            },
        ]
    )
    write_tsv(tables_dir / "summary_numbers.tsv", summary_numbers)

    # markdown report
    md = []
    md.append("# GUS report")
    md.append("")
    md.append("This report summarizes strict pident95/cov90 GUS detection across 18 subjects.")
    md.append("")
    md.append("## Key metrics")
    md.append("")
    md.append(markdown_table(summary_numbers))
    md.append("")
    md.append("## Subject-level burden")
    md.append("")
    md.append(markdown_table(subjects))
    md.append("")
    md.append("## Figures")
    for idx, (fn, cap) in enumerate(
        [
            ("figure1_subject_burden.svg", "Subject-level GUS burden and detected gene counts."),
            ("figure2_prevalence_and_heatmap.svg", "Detection prevalence and the core total GUS heatmap."),
            ("figure3_clade_abundance_boxplot.svg", "Clade abundance distribution across subjects."),
            ("figure4_clade_stacked_bars.svg", "Subject-by-subject clade composition."),
            ("figure5_taxonomy_resolution.svg", "Taxonomy resolution among detected total GUS."),
            ("figure6_clade_taxonomy_composition.svg", "Abundance-weighted taxonomy composition by clade."),
            ("figure7_sankey.svg", "Clade-to-taxonomy Sankey summary."),
        ],
        start=1,
    ):
        md.append(f"### Figure {idx}")
        md.append("")
        md.append(f"![{cap}](figures/{fn})")
        md.append("")
        md.append(cap)
        md.append("")
    if active_figures:
        md.append("## Active GUS figures")
        md.append("")
        md.append("Active GUS = the first 11 entries of Table S3. These figures mirror figures 3, 4 and 6 but are restricted to the active set.")
        md.append("")
        for fn, cap in active_figures:
            md.append(f"![{cap}](figures/{fn})")
            md.append("")
            md.append(cap)
            md.append("")
    (out / "report.md").write_text("\n".join(md), encoding="utf-8")

    html_parts = [
        "<!doctype html>",
        "<html lang=\"en\">",
        "<head>",
        '<meta charset="utf-8" />',
        '<meta name="viewport" content="width=device-width, initial-scale=1" />',
        "<title>GUS report</title>",
        "<style>",
        "body{font-family:Arial,Helvetica,sans-serif;color:#1f2933;background:#f8fafc;margin:0;padding:0;}",
        ".wrap{max-width:1400px;margin:0 auto;padding:32px 28px 64px;}",
        "h1,h2,h3{color:#102a43;}",
        "p,li{line-height:1.55;}",
        ".card{background:#fff;border:1px solid #e3e8ee;border-radius:14px;padding:20px 22px;margin:18px 0;box-shadow:0 1px 3px rgba(16,24,40,.05);}",
        "img.figure{width:100%;height:auto;border:1px solid #e3e8ee;border-radius:10px;background:#fff;}",
        "table.dataframe{border-collapse:collapse;width:100%;font-size:13px;}",
        "table.dataframe th, table.dataframe td{border:1px solid #d9e2ec;padding:6px 8px;text-align:left;vertical-align:top;}",
        "table.dataframe th{background:#eef2f7;}",
        ".small{font-size:12px;color:#52606d;}",
        "</style>",
        "</head>",
        "<body>",
        '<div class="wrap">',
        "<h1>GUS report</h1>",
        "<p>This report summarizes strict pident95/cov90 GUS detection across 18 subjects. The main analysis uses the total GUS set; candidate GUS are retained for comparison where relevant.</p>",
        "<div class=\"card\">",
        "<h2>Key metrics</h2>",
        html_table(summary_numbers),
        "</div>",
        "<div class=\"card\">",
        "<h2>Subject-level burden</h2>",
        "<p>Subject order is sorted by descending total GUS genome percent.</p>",
        html_table(subjects),
        "</div>",
    ]
    for idx, (fn, cap) in enumerate(
        [
            ("figure1_subject_burden.svg", "Subject-level GUS burden and detected gene counts."),
            ("figure2_prevalence_and_heatmap.svg", "Detection prevalence and the core total GUS heatmap."),
            ("figure3_clade_abundance_boxplot.svg", "Clade abundance distribution across subjects."),
            ("figure4_clade_stacked_bars.svg", "Subject-by-subject clade composition."),
            ("figure5_taxonomy_resolution.svg", "Taxonomy resolution among detected total GUS."),
            ("figure6_clade_taxonomy_composition.svg", "Abundance-weighted taxonomy composition by clade."),
            ("figure7_sankey.svg", "Clade-to-taxonomy Sankey summary."),
        ],
        start=1,
    ):
        html_parts.extend(
            [
                '<div class="card">',
                f"<h2>Figure {idx}</h2>",
                f'<img class="figure" src="figures/{fn}" alt="{html.escape(cap)}" />',
                f"<p class=\"small\">{html.escape(cap)}</p>",
                "</div>",
            ]
        )
    if active_figures:
        html_parts.append('<div class="card">')
        html_parts.append("<h2>Active GUS figures</h2>")
        html_parts.append(
            "<p class=\"small\">Active GUS = the first 11 entries of Table S3. These mirror figures 3, 4 and 6 but are restricted to the active set.</p>"
        )
        html_parts.append("</div>")
        for fn, cap in active_figures:
            html_parts.extend(
                [
                    '<div class="card">',
                    f'<img class="figure" src="figures/{fn}" alt="{html.escape(cap)}" />',
                    f"<p class=\"small\">{html.escape(cap)}</p>",
                    "</div>",
                ]
            )
    html_parts.extend(["</div>", "</body>", "</html>"])
    (out / "index.html").write_text("\n".join(html_parts), encoding="utf-8")

    # short README
    readme = [
        "GUS report outputs",
        "",
        f"Open {out / 'index.html'} in a browser to view the report.",
        "",
        "Tables are in report/tables and figures are in report/figures.",
    ]
    (out / "README.txt").write_text("\n".join(readme), encoding="utf-8")

    rscript = shutil.which("Rscript")
    figure_script = Path(__file__).with_name("build_gus_report_figures.R")
    if rscript and figure_script.exists():
        subprocess.run(
            [rscript, str(figure_script), "--report-dir", str(out)],
            check=True,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--abundance-dir",
        type=Path,
        default=Path("results/summary/pident95_cov90/abundance_all_subject"),
        help="Directory containing the merged all-subject abundance tables.",
    )
    parser.add_argument(
        "--total-meta",
        type=Path,
        default=Path("results/reference/gus_279_plus_fmn_metadata.taxonomy.tsv"),
    )
    parser.add_argument(
        "--candidate-meta",
        type=Path,
        default=Path("results/reference/candidate_gus_metadata.taxonomy.tsv"),
    )
    parser.add_argument(
        "--active-meta",
        type=Path,
        default=Path("results/reference/active_gus_metadata.taxonomy.tsv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/summary/pident95_cov90/report"),
    )
    args = parser.parse_args()
    build_report(args)


if __name__ == "__main__":
    main()
