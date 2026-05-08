#!/usr/bin/env python3
"""
Script Name:  get_coverage.py
Description:  Calculate average and median coverage per reference from a depth table.
Author:       Dr. Haiting Chai
Last update:  2026-APR-11
Version:      ver.2026-04-11

DEPENDENCIES:
    - Python 3.10
    - __future__ (standard library)
    - argparse (standard library)
    - collections (standard library)
    - csv (standard library)
    - logging (standard library)
    - pathlib (standard library)
    - statistics (standard library)
    - sys (standard library)
    - time (standard library)
    - assets.utils_banner (project module)
"""

from __future__ import annotations

import argparse
import csv
import logging
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from assets.utils_banner import print_header, print_summary, setup_logging


logger = setup_logging()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Calculate average and median coverage per reference.",
    )
    parser.add_argument("--input", required=True, help="Input depth table")
    parser.add_argument("--prefix", required=True, help="Sample prefix written to the output summary")
    parser.add_argument("--output", default="coverage_summary.txt", help="Output summary file")
    return parser


def load_depths(input_path: Path) -> dict[str, list[int]]:
    """Read per-position depth values and group them by reference name."""
    depths: dict[str, list[int]] = defaultdict(list)
    with input_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            reference, _position, depth = line.strip().split()
            depths[reference].append(int(depth))

    if not depths:
        raise SystemExit(f"No depth records were found in: {input_path}")
    return depths


def summarise_depths(depths: dict[str, list[int]]) -> tuple[dict[str, float], dict[str, float]]:
    """Calculate average and median depth for each reference."""
    average_coverage = {reference: sum(values) / len(values) for reference, values in depths.items()}
    median_coverage = {reference: statistics.median(values) for reference, values in depths.items()}
    return average_coverage, median_coverage


def write_summary(
    prefix: str,
    output_path: Path,
    average_coverage: dict[str, float],
    median_coverage: dict[str, float],
) -> None:
    """Append the coverage summary as a single TSV row, creating a header when needed."""
    sorted_references = sorted(average_coverage, key=lambda reference: average_coverage[reference], reverse=True)
    top_reference = sorted_references[0]

    average_string = ";".join(f"{reference}:{average_coverage[reference]:.6f}" for reference in sorted_references)
    median_string = ";".join(f"{reference}:{median_coverage[reference]:.6f}" for reference in sorted_references)

    with output_path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
            quoting=csv.QUOTE_NONE,
            escapechar="\\",
        )
        if output_path.stat().st_size == 0:
            writer.writerow(["sample", "ave_covs", "median_covs", "highest_ave_cov", "highest_median_cov"])
        writer.writerow(
            [
                prefix,
                average_string,
                median_string,
                f"{average_coverage[top_reference]:.6f}",
                f"{median_coverage[top_reference]:.6f}",
            ]
        )

    logger.info("Coverage summary appended to %s", output_path)


def main() -> None:
    args = build_parser().parse_args()
    input_path = Path(args.input)
    if not input_path.is_file():
        raise SystemExit(f"Input depth table was not found: {input_path}")

    depths = load_depths(input_path)
    average_coverage, median_coverage = summarise_depths(depths)
    write_summary(args.prefix, Path(args.output), average_coverage, median_coverage)


if __name__ == "__main__":
    start_time = time.time()
    print_header("get_coverage.py")
    main()
    print_summary(start_time, "get_coverage.py")
