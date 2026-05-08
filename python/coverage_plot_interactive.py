#!/usr/bin/env python3
"""
Script Name:  coverage_plot_interactive.py
Description:  Generate an interactive combined coverage plot from depth tables.
Author:       Dr. Haiting Chai
Last update:  2026-APR-11
Version:      ver.2026-04-11

DEPENDENCIES:
    - Python 3.10
    - __future__ (standard library)
    - argparse (standard library)
    - logging (standard library)
    - pathlib (standard library)
    - sys (standard library)
    - time (standard library)
    - numpy (external package)
    - pandas (external package)
    - plotly (external package)
    - assets.utils_banner (project module)
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from assets.utils_banner import print_header, print_summary, setup_logging


logger = setup_logging()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate an interactive combined coverage plot from samtools depth files.",
    )
    parser.add_argument("--indir", required=True, help="Directory containing depth `.txt` files")
    parser.add_argument("--output", required=True, help="Output HTML file path")
    parser.add_argument("--yaxis", required=True, choices=["raw", "log10"], help="Y-axis scaling mode")
    return parser


def transform_depths(depths: np.ndarray, yaxis_mode: str) -> tuple[np.ndarray, str]:
    """Transform coverage values for the requested y-axis mode."""
    if yaxis_mode == "raw":
        return depths, "Depth"
    if yaxis_mode == "log10":
        return np.log10(np.where(depths == 0, 0.1, depths)), "log10(Depth)"
    raise ValueError(f"Unsupported y-axis mode: {yaxis_mode}")


def add_depth_trace(figure: go.Figure, depth_file: Path, yaxis_mode: str) -> bool:
    """Load one depth file and add it as a trace to the combined figure."""
    sample_name = depth_file.stem
    dataframe = pd.read_csv(
        depth_file,
        sep="\t",
        header=None,
        names=["reference", "position", "depth"],
    )

    if dataframe.empty:
        logger.warning("Skipping empty depth file: %s", depth_file.name)
        return False

    y_values, y_label = transform_depths(dataframe["depth"].to_numpy(), yaxis_mode)
    figure.add_trace(
        go.Scatter(
            x=dataframe["position"].to_numpy(),
            y=y_values,
            mode="lines",
            name=sample_name,
            line=dict(width=1),
            hovertemplate=(
                "Sample: %{fullData.name}<br>"
                "Position: %{x}<br>"
                "Depth: %{y}<extra></extra>"
            ),
        )
    )
    figure.update_yaxes(title_text=y_label)
    return True


def plot_combined(input_dir: Path, output_path: Path, yaxis_mode: str) -> Path | None:
    """Create a combined interactive plot from all depth files in the input directory."""
    depth_files = sorted(input_dir.glob("*.txt"))
    if not depth_files:
        logger.warning("No .txt depth files were found in %s", input_dir)
        return None

    figure = go.Figure()
    added_traces = 0

    for depth_file in depth_files:
        try:
            if add_depth_trace(figure, depth_file, yaxis_mode):
                added_traces += 1
        except Exception as error:
            logger.error("Failed to process %s: %s", depth_file.name, error)

    if added_traces == 0:
        logger.warning("No valid depth traces were added to the figure")
        return None

    figure.update_layout(
        title="Interactive Combined Coverage Plot",
        xaxis_title="Position",
        template="plotly_white",
        hovermode="closest",
        legend=dict(orientation="v", x=1.02, y=1),
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.write_html(output_path)
    logger.info("Interactive coverage plot written to %s", output_path)
    return output_path


def main() -> None:
    args = build_parser().parse_args()
    input_dir = Path(args.indir)
    if not input_dir.is_dir():
        raise SystemExit(f"Input directory does not exist: {input_dir}")

    output_path = Path(args.output)
    plot_combined(input_dir, output_path, args.yaxis)


if __name__ == "__main__":
    start_time = time.time()
    print_header("coverage_plot_interactive.py")
    main()
    print_summary(start_time, "coverage_plot_interactive.py")
