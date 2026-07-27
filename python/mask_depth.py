#!/usr/bin/env python3
"""
Script Name:  mask_depth.py
Description:  Convert low-depth positions from a samtools depth file into BED intervals.
Author:       Dr. Haiting Chai
Last update:  2026-APR-12
Version:      ver.2026-04-12

DEPENDENCIES:
    - Python 3.10
    - __future__ (standard library)
    - argparse (standard library)
    - logging (standard library)
    - pathlib (standard library)
    - sys (standard library)
    - time (standard library)
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
    """Build the CLI parser for low-depth BED masking."""
    parser = argparse.ArgumentParser(
        description="Mask positions below a depth threshold and export them as BED intervals.",
    )
    parser.add_argument("--input", required=True, help="Input depth file in samtools depth format")
    parser.add_argument("--bq", type=float, required=True, help="Depth threshold")
    parser.add_argument("--output", required=True, help="Output BED file")
    return parser


def mask_depth(input_path: Path, threshold: float, output_path: Path) -> int:
    """Write BED intervals for positions with depth values below the requested threshold."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    masked_sites = 0

    with input_path.open("r", encoding="utf-8") as input_handle, output_path.open("w", encoding="utf-8") as output_handle:
        for line in input_handle:
            if not line.strip():
                continue

            parts = line.strip().split()
            if len(parts) < 3:
                logger.warning("Skipping malformed depth line: %s", line.rstrip("\n"))
                continue

            chrom = parts[0]
            position = int(parts[1])
            depth = float(parts[2])

            # BED format uses 0-based starts and 1-based ends.
            if depth < threshold:
                output_handle.write(f"{chrom}\t{position - 1}\t{position}\n")
                masked_sites += 1

    logger.info("Wrote %s masked BED intervals to %s", masked_sites, output_path)
    return masked_sites


# ========================================================================================
# Main Entry Point
# ========================================================================================

def main() -> None:
    """Parse CLI arguments and convert low-depth positions into BED intervals."""
    args = build_parser().parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.is_file():
        raise SystemExit(f"Input depth file was not found: {input_path}")

    mask_depth(input_path, args.bq, output_path)


if __name__ == "__main__":
    start_time = time.time()
    print_header("mask_depth.py")
    main()
    print_summary(start_time, "mask_depth.py")
