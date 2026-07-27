#!/usr/bin/env python3
"""
Script Name:  dedouble_depth.py
Description:  Collapse doubled-reference depth values back onto the original reference coordinates.
Author:       Dr. Haiting Chai
Last update:  2026-APR-12
Version:      ver.2026-04-12

DEPENDENCIES:
    - Python 3.10
    - __future__ (standard library)
    - argparse (standard library)
    - collections (standard library)
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
from collections import OrderedDict
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
    """Build the CLI parser for doubled-reference depth collapse."""
    parser = argparse.ArgumentParser(
        description=(
            "Dedouble a samtools depth -a file produced from a reference where each contig "
            "was concatenated with itself."
        ),
    )
    parser.add_argument("--input", required=True, help="Input depth file from samtools depth -a")
    parser.add_argument("--originalref", required=True, help="Original FASTA reference")
    parser.add_argument("--output", required=True, help="Output dedoubled depth file")
    return parser


def parse_fasta_lengths(fasta_path: Path) -> OrderedDict[str, int]:
    """Read FASTA records and return contig lengths while preserving input order."""
    lengths: OrderedDict[str, int] = OrderedDict()
    current_name: str | None = None
    current_length = 0

    with fasta_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    lengths[current_name] = current_length
                current_name = line[1:].split()[0]
                current_length = 0
            else:
                current_length += len(line)

    if current_name is not None:
        lengths[current_name] = current_length
    if not lengths:
        raise ValueError(f"No FASTA records were found in: {fasta_path}")
    return lengths


def flush_contig(output_handle, contig: str, combined_depths: list[int], length: int) -> None:
    """Write one dedoubled contig block to the output depth file."""
    for position in range(1, length + 1):
        output_handle.write(f"{contig}\t{position}\t{combined_depths[position]}\n")


def dedouble_depth_file(input_path: Path, original_ref: Path, output_path: Path) -> None:
    """Collapse both halves of a doubled-reference depth table onto original positions."""
    reference_lengths = parse_fasta_lengths(original_ref)
    current_contig: str | None = None
    current_length: int | None = None
    combined_depths: list[int] | None = None

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with input_path.open("r", encoding="utf-8") as input_handle, output_path.open("w", encoding="utf-8") as output_handle:
        for line_number, line in enumerate(input_handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                raise ValueError(f"Malformed depth line at {line_number}: {line}")

            contig, pos_string, depth_string = parts[0], parts[1], parts[2]
            if contig not in reference_lengths:
                raise ValueError(f"Contig '{contig}' at line {line_number} is not present in --originalref")

            position = int(pos_string)
            depth = int(depth_string)
            contig_length = reference_lengths[contig]

            if current_contig != contig:
                if current_contig is not None and combined_depths is not None and current_length is not None:
                    flush_contig(output_handle, current_contig, combined_depths, current_length)

                current_contig = contig
                current_length = contig_length
                combined_depths = [0] * (contig_length + 1)

            # Fold the duplicated coordinate space back onto the original contig length
            # so positions from the first and second halves accumulate together.
            target_position = ((position - 1) % contig_length) + 1
            assert combined_depths is not None
            combined_depths[target_position] += depth

        if current_contig is not None and combined_depths is not None and current_length is not None:
            flush_contig(output_handle, current_contig, combined_depths, current_length)

    logger.info("Dedoubled depth file written to %s", output_path)


# ========================================================================================
# Main Entry Point
# ========================================================================================

def main() -> None:
    """Parse CLI arguments and collapse the doubled depth table."""
    args = build_parser().parse_args()
    input_path = Path(args.input)
    original_ref = Path(args.originalref)
    output_path = Path(args.output)

    if not input_path.is_file():
        raise SystemExit(f"Input depth file was not found: {input_path}")
    if not original_ref.is_file():
        raise SystemExit(f"Original reference FASTA was not found: {original_ref}")

    dedouble_depth_file(input_path, original_ref, output_path)


if __name__ == "__main__":
    start_time = time.time()
    print_header("dedouble_depth.py")
    main()
    print_summary(start_time, "dedouble_depth.py")
