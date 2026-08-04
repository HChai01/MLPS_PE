#!/usr/bin/env python3
"""
Script Name:  draw_qc_plot.py
Description:  Draw QC plots from seqtk fqchk TSV output.
Author:       Dr. Haiting Chai
Last update:  2026-APR-11
Version:      ver.2026-04-11

DEPENDENCIES:
    - Python 3.10
    - argparse (standard library)
    - gzip (standard library)
    - logging (standard library)
    - io (standard library)
    - math (standard library)
    - pathlib (standard library)
    - re (standard library)
    - sys (standard library)
    - time (standard library)
    - numpy (external package)
    - pandas (external package)
    - matplotlib (external package)
    - assets.utils_banner (project module)
"""

# ========================================================================================
# Imports
# ========================================================================================

import argparse
import gzip
import logging
import re
import sys
import time
from io import StringIO
from math import ceil, floor
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator, NullFormatter, ScalarFormatter

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from assets.utils_banner import print_header, print_summary, setup_logging


# ========================================================================================
# Configuration
# ========================================================================================

logger = setup_logging()


# ========================================================================================
# Helper Functions
# ========================================================================================

def get_parser():
    """Build the CLI parser for seqtk `fqchk` QC plotting."""
    parser = argparse.ArgumentParser(
        description="Draw QC plots from seqtk fqchk TSV output.",
    )
    parser.add_argument("--input", required=True, help="TSV input from seqtk fqchk TSV output.")
    parser.add_argument("--qthres", default="20", help="Phred Q threshold used in seqtk fqchk.")
    parser.add_argument("--output", default="out.pdf", help="QC plots.")

    return parser

def choose_col(df, candidates):
    """Pick the first matching column name, allowing for minor label variations."""
    # Prefer exact column names before falling back to normalized comparisons.
    for c in candidates:
        if c in df.columns:
            return c
    norm = {x.lower().replace("%", "").replace("#", "").strip(): x for x in df.columns}
    for c in candidates:
        key = c.lower().replace("%", "").replace("#", "").strip()
        if key in norm:
            return norm[key]
    return None

def read_file_find_table(path):
    """Locate the `POS` table inside `fqchk` output and return its summary header."""
    p = Path(path)
    opener = gzip.open if p.suffix.lower() in (".gz", ".gzip") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    if not lines:
        raise ValueError(f"Empty file: {path}")

    pos_idx = None
    for i, L in enumerate(lines):
        if L.strip().upper().startswith("POS"):
            pos_idx = i
            break
    if pos_idx is None:
        raise ValueError("Could not find POS header line in file.")

    header_lines = lines[:pos_idx]
    summary_line = header_lines[0].strip() if header_lines else ""
    table_text = "".join(lines[pos_idx:])
    df = pd.read_csv(StringIO(table_text), sep=r"\s+", engine="python")
    return summary_line, df

def parse_header(line):
    """Extract high-level read-length statistics from the free-text fqchk header."""
    out = {}
    m = re.search(r"min_len:\s*(\d+)", line)
    if m:
        out["min_len"] = int(m.group(1))
    m = re.search(r"max_len:\s*(\d+)", line)
    if m:
        out["max_len"] = int(m.group(1))
    m = re.search(r"avg_len:\s*([\d.]+)", line)
    if m:
        out["avg_len"] = float(m.group(1))
    return out

def vaf_redline_from_avgQ(avgQ):
    """Estimate a conservative variant-frequency floor from average base quality."""
    multiplier = 5
    min_vaf = 0.005
    max_vaf = 0.100

    base_error = 10 ** (-avgQ / 10.0)
    vaf = multiplier * base_error
    # Clamp to practical reporting bounds so extreme quality values do not produce
    # implausibly small or large thresholds in the figure title.
    vaf = max(min_vaf, min(vaf, max_vaf))
    return vaf

def classify_quality(avgQ, errQ):
    """Classify overall sequencing quality from avgQ and error-derived quality."""
    d = avgQ - errQ

    # Failed or unusable data: globally bad + heterogeneous, suggest to be excluded from analysis
    # Q20 means 1 error/100 bases or 99.0% base accuracy
    # Q17 means 1 error/50 bases or 98.0% base accuracy
    # avgQ becomes misleading
    # Variant calls become unstable
    # Adjustment for quality-related drop-off parameters is required if these data must be processed
    if avgQ < 20.0 or errQ < 17.0:
        return f"Failed: Sample with poor accuracy; trimming may provide limited improvement, but resequencing is recommended."

    # Clean data: high accuracy, low heterogeneity, error ≤ ~0.2% and no strong low-Q tail
    # Q30 means 1 error/1000 bases or 99.9% base accuracy
    # Q27 means 1 error/500 bases or 99.8% base accuracy
    # d3.0 means true error rate is ~2× the uniform-quality expectation 
    # Error are spread evenly, thus trimming has little effect
    # Low-quality bases may exist but do not dominate
    # Variant callers can behave as expected
    if avgQ >= 30.0 and errQ >= 27.0 and d <= 3.0:
        return f"Clean_data: High-quality sample with low heterogeneity; estimated error rate ≤0.2% and no pronounced low-quality tail."

    # d3.0 means true error rate is ~2× the uniform-quality expectation
    # avgQ may from 20 ~ 30
    # Aggressive trimming will only reduce depth without improving errQ much
    if d <= 3.0:
        return f"Acceptable_data: Sample with low heterogeneity; moderate absolute accuracy and no pronounced low-quality tail."

    # QC required data: high overall quality but moderate heterogeneity (2–6× inflation)
    # Q20 means 1 error/100 bases or 99.0% base accuracy
    # d3.0 means true error rate is ~2× the uniform-quality expectation
    # d7.8 means true error rate is ~6× the uniform-quality expectation
    # Need standard trimming
    if avgQ >= 20.0 and d > 3.0 and d <= 7.8:
        return f"QC_required_data: Sample with moderate heterogeneity; elevated error driven by low-quality tails, expected to improve after standard trimming."

    # Otherwise may still be usable but need to be treated with caution
    return f"Risky_data: Sample with high heterogeneity; error dominated by low-quality outliers and downstream analyses should be interpreted with caution."

def font_sizes(n_positions):
    """Scale font sizes so plots remain readable across different read lengths."""
    if n_positions < 50:
        return 12, 10, 10, 12
    if n_positions < 200:
        return 10, 9, 9, 11
    if n_positions < 800:
        return 8, 7, 7, 9
    return 7, 6, 6, 8

def add_bands(ax, y_min, y_max, step=10, alpha=0.06):
    """Add alternating background bands to help track positions across panels."""
    start = int(y_min)
    end = int(y_max)
    for y in range(start, end + 1, step * 2):
        ax.axhspan(y - 0.5, min(y + step - 0.5, end + 0.5), color="#777777", alpha=alpha, zorder=0, linewidth=0)

def style_shared_y_axis(ax, y_min, y_max, show_y_labels=False, tick_fs=8):
    """Apply the shared vertical layout used by all panels in the composite figure."""
    ax.invert_yaxis()
    ax.set_ylim(y_min - 0.5, y_max + 0.5)

    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_minor_formatter(NullFormatter())

    if not show_y_labels:
        ax.yaxis.set_major_formatter(NullFormatter())

    # x labels ON everywhere (major); minor ticks kept but unlabeled
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis="x", which="major", labelsize=tick_fs)
    ax.tick_params(axis="y", which="major", labelsize=tick_fs)

    # frame
    for s in ax.spines.values():
        s.set_visible(True)

    add_bands(ax, y_min, y_max)

def nice_round_step(x: float) -> int:
    """Round a step size to the nearest plotting-friendly 1/2/5 multiple."""
    if x <= 0:
        return 1
    exp = int(np.floor(np.log10(x)))
    base = x / (10 ** exp)
    if base <= 1:
        nice = 1
    elif base <= 2:
        nice = 2
    elif base <= 5:
        nice = 5
    else:
        nice = 10
    return int(nice * (10 ** exp))

def force_y_labels(ax_cov, y_values, step=10, fontsize=8):
    """Force y-axis labels onto the master coverage panel at regular intervals."""
    y_values = np.asarray(y_values, dtype=int)
    ticks = [int(v) for v in y_values if int(v) % step == 0]
    ax_cov.set_yticks(ticks)
    ax_cov.set_yticklabels([str(v) for v in ticks], fontsize=fontsize)


# ========================================================================================
# Main Entry Point
# ========================================================================================

def main():
    """Render the multi-panel QC figure from a seqtk `fqchk` report."""
    parser = get_parser()
    args = parser.parse_args()

    summary_line, df = read_file_find_table(args.input)
    header = parse_header(summary_line)

    if "POS" not in df.columns:
        raise KeyError("Input table missing 'POS' column. Columns: " + ", ".join(df.columns))

    # Split ALL vs per-position
    all_mask = df["POS"].astype(str).str.upper() == "ALL"
    all_row = df[all_mask].iloc[0] if all_mask.any() else None
    pos_df = df[~all_mask].copy()

    # Coerce POS to int, sort
    pos_df["POS"] = pd.to_numeric(pos_df["POS"], errors="coerce").astype(int)
    pos_df = pos_df.sort_values("POS")

    y = pos_df["POS"].values
    y_min, y_max = int(y.min()), int(y.max())
    yrange = float(y_max) - float(y_min)
    if yrange <= 0:
        step = 1
    else:
        approx = max(1.0, yrange / 10.0)
        step = nice_round_step(ceil(approx))

    # Extract ALL stats for title
    reads = None
    avgQ_all = np.nan
    errQ_all = np.nan
    if all_row is not None:
        try:
            avgQ_all = float(all_row["avgQ"])
        except Exception:
            pass
        try:
            errQ_all = float(all_row["errQ"])
        except Exception:
            pass

        # Read count = #bases at POS=1 (each read contributes one base at cycle 1)
        cov_col = choose_col(pos_df, ["#bases", "bases", "n_bases"])
        if cov_col is None:
            raise KeyError("Coverage column not found. Available columns: " + ", ".join(pos_df.columns))

        row1 = pos_df.loc[pos_df["POS"] == 1, cov_col]
        if row1.empty:
            raise ValueError("Cannot infer read count: POS=1 row not found.")
        reads = int(pd.to_numeric(row1.iloc[0], errors="coerce"))

    # Fall back to mean across positions if needed
    if np.isnan(avgQ_all) and "avgQ" in pos_df.columns:
        avgQ_all = float(np.nanmean(pd.to_numeric(pos_df["avgQ"], errors="coerce")))
    if np.isnan(errQ_all) and "errQ" in pos_df.columns:
        errQ_all = float(np.nanmean(pd.to_numeric(pos_df["errQ"], errors="coerce")))

    # Auto fonts
    base_fs, tick_fs, leg_fs, sup_fs = font_sizes(len(y))
    mpl.rcParams.update({"font.size": base_fs})

    # Coverage
    cov_col = choose_col(pos_df, ["#bases", "bases", "n_bases"])
    if cov_col is None:
        raise KeyError("Coverage column not found. Available columns: " + ", ".join(pos_df.columns))
    coverage = pd.to_numeric(pos_df[cov_col], errors="coerce").fillna(0).astype(int).values

    # Frequencies
    bases = ["A", "C", "G", "T", "N"]
    freqs = {}
    for b in bases:
        col = choose_col(pos_df, [f"%{b}", b, b.lower(), f"{b}%"])
        if col is None:
            continue
        freqs[b] = pd.to_numeric(pos_df[col], errors="coerce").fillna(0).values

    # Quality (as lines)
    if "avgQ" not in pos_df.columns or "errQ" not in pos_df.columns:
        raise KeyError("Missing avgQ/errQ columns. Available columns: " + ", ".join(pos_df.columns))
    avgQ = pd.to_numeric(pos_df["avgQ"], errors="coerce").fillna(0).values
    errQ = pd.to_numeric(pos_df["errQ"], errors="coerce").fillna(0).values

    # Dropoff bases
    dropoff = pd.to_numeric(pos_df["%low"], errors="coerce").fillna(0).values

    # Keep all metrics on a shared y-axis so cycle-specific artefacts line up visually
    # across coverage, base composition, quality, and low-quality-tail panels.
    fig = plt.figure(figsize=(15, 8))
    gs = fig.add_gridspec(1, 8, width_ratios=[1.4, 1, 1, 1, 1, 1, 1.4, 1], wspace=0.08)

    # Coverage plot (master sharey)
    ax_cov = fig.add_subplot(gs[0, 0])
    ax_cov.plot(coverage, y, linewidth=1.3, color="#ce181e", zorder=2)
    style_shared_y_axis(ax_cov, y_min, y_max, show_y_labels=True, tick_fs=tick_fs)
    ax_cov.set_xlabel("Covered reads", fontsize=base_fs)
    ax_cov.set_ylabel("Read position (bp)", fontsize=base_fs)

    # Frequency plots
    colors = {"A": "#1f77b4", "C": "#ff7f0e", "G": "#2ca02c", "T": "#d62728", "N": "#9467bd"}
    freq_axes = []
    for i, b in enumerate(bases):
        axf = fig.add_subplot(gs[0, i + 1], sharey=ax_cov)
        if b in freqs:
            axf.plot(freqs[b], y, color=colors.get(b, "k"), linewidth=1.3, zorder=2)
        style_shared_y_axis(axf, y_min, y_max, show_y_labels=False, tick_fs=tick_fs)
        axf.set_xlabel(f"{b}%", fontsize=base_fs)
        freq_axes.append(axf)

    # Quality plot
    ax_q = fig.add_subplot(gs[0, 6], sharey=ax_cov)
    # Light grey dashed connectors between avgQ and errQ
    for xi, xj, yi in zip(avgQ, errQ, y):
        ax_q.plot(
            [xi, xj], [yi, yi],
            color="#bdbdbd",
            linewidth=0.8,
            linestyle="--",
            zorder=1
        )
    ax_q.plot(avgQ, y, linewidth=1.3, color="#ce181e", label="avgQ", zorder=2)
    ax_q.plot(errQ, y, linewidth=1.3, color="#007cc0", label="errQ", zorder=2)
    style_shared_y_axis(ax_q, y_min, y_max, show_y_labels=False, tick_fs=tick_fs)
    ax_q.set_xlabel("Phred", fontsize=base_fs)
    ax_q.legend(fontsize=leg_fs, loc="best")

    # Dropoff plot, only include bases with Phred Q <= args.qthres
    ax_drop = fig.add_subplot(gs[0, 7], sharey=ax_cov)
    ax_drop.plot(dropoff, y, linewidth=1.2, color="#49c0b6", zorder=3)
    style_shared_y_axis(ax_drop, y_min, y_max, show_y_labels=False, tick_fs=tick_fs)
    ax_drop.set_xlabel(f"bases% below Q{args.qthres}", fontsize=base_fs)

    # Force y labels ONLY on coverage, then hide on others
    force_y_labels(ax_cov, y, step, fontsize=tick_fs)
    for a in freq_axes + [ax_q] + [ax_drop]:
        a.tick_params(axis="y", which="both", labelleft=False)

    # Estimate convincing VAF
    vaf = vaf_redline_from_avgQ(float(avgQ_all))
    vaf_pct = vaf * 100.0

    # Title
    qc = classify_quality(float(avgQ_all), float(errQ_all)) if (not np.isnan(avgQ_all) and not np.isnan(errQ_all)) else "unknown"
    reads_text = f"{reads} reads" if reads is not None else ""
    title = (
        f"FOUND {reads_text} from {header.get('min_len','?')}–{header.get('max_len','?')} bp "
        f"(ave={header.get('avg_len','?')} bp), aveQ={avgQ_all:.1f}, errQ={errQ_all:.1f}. Convincing VAF ≥ {vaf_pct:.2f}%\n"
        f"{qc}"
    )
    fig.suptitle(title, fontsize=sup_fs, y=0.98)

    plt.subplots_adjust(left=0.07, right=0.99, top=0.92)

    outp = Path(args.output)
    if outp.suffix.lower() == ".pdf":
        fig.savefig(args.output, bbox_inches="tight")
    else:
        fig.savefig(args.output, dpi=300, bbox_inches="tight")

    logger.info("QC plots saved to %s", args.output)

if __name__ == "__main__":
    start_time = time.time()
    print_header("draw_qc_plot.py")
    main()
    print_summary(start_time, "draw_qc_plot.py")
