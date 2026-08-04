#!/usr/bin/env python3
"""
Script Name:  cigar_summary.py
Description:  Summarize key CIGAR operation counts from a text file containing one CIGAR string per line.
Author:       Dr. Haiting Chai
Last update:  2026-APR-13
Version:      ver.2026-04-13

DEPENDENCIES:
    - Python 3.10
    - __future__ (standard library)
    - argparse (standard library)
    - logging (standard library)
    - pathlib (standard library)
    - re (standard library)
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
import re
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
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


# ========================================================================================
# Helper Functions
# ========================================================================================

def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for CIGAR summary export."""
    parser = argparse.ArgumentParser(
        description="Summarize CIGAR operations and append intelligently to a TSV output file.",
    )
    parser.add_argument("--input", required=True, help="Input CIGAR file")
    parser.add_argument("--outputtsv", required=True, help="Output TSV file")
    parser.add_argument("--prefix", required=True, help="Prefix used for the summary column names")
    parser.add_argument("--samplename", default=None, help="Optional sample name")
    return parser


def parse_cigar_file(path: Path) -> tuple[int, int, int, int, int]:
    """Accumulate the clipping, match, insertion, and deletion counts from the input CIGAR file."""
    hardclip_bases = 0
    softclip_bases = 0
    match_bases = 0
    insertion_bases = 0
    deletion_bases = 0

    with path.open("r", encoding="utf-8") as input_handle:
        for line_number, line in enumerate(input_handle, start=1):
            cigar = line.strip()
            if not cigar or cigar == "*":
                continue

            operations = CIGAR_RE.findall(cigar)
            if not operations:
                raise ValueError(f"Invalid CIGAR at line {line_number}: {cigar}")

            # Count only the operations requested for the exported summary.
            for length_string, operation in operations:
                length = int(length_string)
                if operation == "H":
                    hardclip_bases += length
                elif operation == "S":
                    softclip_bases += length
                elif operation in {"M", "=", "X"}:
                    match_bases += length
                elif operation == "I":
                    insertion_bases += length
                elif operation == "D":
                    deletion_bases += length

    return hardclip_bases, softclip_bases, match_bases, insertion_bases, deletion_bases


def infer_sample_name(input_path: Path, user_sample: str | None) -> str:
    """Return user-specified sample name or infer from filename."""
    if user_sample:
        return user_sample
    return input_path.stem


def expected_header(prefix: str) -> str:
    """Build the expected TSV header for the requested prefix."""
    return (
        f"sample\t{prefix}_hardclip_bases\t{prefix}_softclip_bases\t"
        f"{prefix}_match_bases\t{prefix}_insertion_bases\t{prefix}_deletion_bases"
    )


def header_matches(output_path: Path, prefix: str) -> bool:
    """Return True when the existing output file already has the expected header."""
    if not output_path.exists():
        return False

    with output_path.open("r", encoding="utf-8") as input_handle:
        first_line = input_handle.readline().strip()
    return first_line == expected_header(prefix)


def write_summary(
    output_path: Path,
    prefix: str,
    sample_name: str,
    hardclip_bases: int,
    softclip_bases: int,
    match_bases: int,
    insertion_bases: int,
    deletion_bases: int,
) -> None:
    """Append to a matching TSV file or overwrite when the header layout changes."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    append_mode = header_matches(output_path, prefix)
    mode = "a" if append_mode else "w"

    with output_path.open(mode, encoding="utf-8") as output_handle:
        if not append_mode:
            output_handle.write(expected_header(prefix) + "\n")

        output_handle.write(
            f"{sample_name}\t{hardclip_bases}\t{softclip_bases}\t{match_bases}\t"
            f"{insertion_bases}\t{deletion_bases}\n"
        )

    logger.info("CIGAR summary written to %s", output_path)


# ========================================================================================
# Main Entry Point
# ========================================================================================

def main() -> None:
    """Parse CLI arguments and export the CIGAR operation summary."""
    args = build_parser().parse_args()
    input_path = Path(args.input)
    output_path = Path(args.outputtsv)

    if not input_path.is_file():
        raise SystemExit(f"Input CIGAR file was not found: {input_path}")

    sample_name = infer_sample_name(input_path, args.samplename)
    hardclip_bases, softclip_bases, match_bases, insertion_bases, deletion_bases = parse_cigar_file(input_path)
    write_summary(
        output_path,
        args.prefix,
        sample_name,
        hardclip_bases,
        softclip_bases,
        match_bases,
        insertion_bases,
        deletion_bases,
    )


if __name__ == "__main__":
    start_time = time.time()
    print_header("cigar_summary.py")
    main()
    print_summary(start_time, "cigar_summary.py")
