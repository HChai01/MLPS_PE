#!/usr/bin/env python3
"""
Script Name:  built_nt_freq.py
Description:  Expand `formatpileup.py` nucleotide-frequency output into the non-aligned nt summary table.
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
    - pandas (external library)
    - assets.utils_banner (project module)
"""

from __future__ import annotations

# ========================================================================================
# Imports
# ========================================================================================

import argparse
import logging
import sys
import time
from pathlib import Path

import pandas as pd

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

def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for nt-frequency table expansion."""
    parser = argparse.ArgumentParser(
        description="Expand nt-frequency CSV rows into long-form nucleotide count records.",
    )
    parser.add_argument("input", help="Input nt-frequency CSV from formatpileup.py")
    parser.add_argument("output", help="Output CSV path")
    return parser


def expand_nt_counts(input_path: Path, output_path: Path) -> None:
    """Convert wide nucleotide count rows into one record per observed base."""
    dataframe = pd.read_csv(input_path)
    records: list[dict[str, object]] = []

    for _, row in dataframe.iterrows():
        counts = {base: int(row[base]) for base in ("A", "C", "G", "T")}
        total_acgt = sum(counts.values())
        depth = int(row["Depth"])

        # Preserve positions with zero aligned A/C/G/T support so downstream reports
        # can distinguish empty sites from missing rows.
        if total_acgt == 0:
            records.append(
                {
                    "non-aligned_pos": row["Pos"],
                    "consensus_nt": row["CalledConsensus"],
                    "nt": "NaN",
                    "nt_count": 0,
                    "nt_percent": "0.00",
                }
            )
            continue

        for nucleotide, count in counts.items():
            if count <= 0:
                continue
            percent = (count / depth) * 100 if depth > 0 else 0.0
            records.append(
                {
                    "non-aligned_pos": row["Pos"],
                    "consensus_nt": row["CalledConsensus"],
                    "nt": nucleotide,
                    "nt_count": count,
                    "nt_percent": f"{percent:.2f}",
                }
            )

    pd.DataFrame(records).to_csv(output_path, index=False)
    logger.info("Expanded nucleotide counts written to %s", output_path)


# ========================================================================================
# Main Entry Point
# ========================================================================================

def main() -> None:
    """Parse CLI arguments and expand the nucleotide-frequency table."""
    args = build_parser().parse_args()
    expand_nt_counts(Path(args.input), Path(args.output))


if __name__ == "__main__":
    start_time = time.time()
    print_header("built_nt_freq.py")
    main()
    print_summary(start_time, "built_nt_freq.py")
