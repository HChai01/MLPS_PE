#!/usr/bin/env python3
"""
Script Name:  get_fasta_summary.py
Description:  Summarize FASTA sequence composition into high-confidence, low-confidence, ambiguous, and gap counts.
Author:       Dr. Haiting Chai
Last update:  2026-APR-13
Version:      ver.2026-04-13

DEPENDENCIES:
    - Python 3.10
    - __future__ (standard library)
    - argparse (standard library)
    - csv (standard library)
    - gzip (standard library)
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
import csv
import gzip
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
EXPECTED_HEADER_TEMPLATE = [
    "sample",
    "{prefix}_high_cof_pos",
    "{prefix}_low_cof_pos",
    "{prefix}_ambiguous_pos",
    "{prefix}_gap_pos",
]


# ========================================================================================
# Helper Functions
# ========================================================================================

def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for FASTA composition summarization."""
    parser = argparse.ArgumentParser(
        description="Summarize FASTA sequence composition.",
    )
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--prefix", required=True, help="Prefix for output column names")
    parser.add_argument("--outputtsv", required=True, help="Output TSV file")
    parser.add_argument("--samplename", default=None, help="Optional sample name override")
    parser.add_argument(
        "--seqtype",
        required=True,
        choices=["DNA", "RNA", "NT", "AA", "dna", "rna", "nt", "aa"],
        help="Sequence type",
    )
    return parser


def open_maybe_gzip(path: Path, mode: str = "rt"):
    """Open a plain FASTA or gzipped FASTA transparently."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode, encoding="utf-8")


def get_sample_name(fasta_path: Path, user_sample: str | None) -> str:
    """Return the user-specified sample name or infer it from the FASTA filename."""
    if user_sample:
        return user_sample

    sample_name = fasta_path.name
    for suffix in [
        ".fasta",
        ".fa",
        ".fna",
        ".fas",
        ".faa",
        ".fasta.gz",
        ".fa.gz",
        ".fna.gz",
        ".fas.gz",
        ".faa.gz",
    ]:
        if sample_name.endswith(suffix):
            return sample_name[: -len(suffix)]
    return sample_name


def parse_fasta_sequence(fasta_path: Path) -> str:
    """Read all non-header FASTA lines into a single sequence string."""
    sequence_parts: list[str] = []
    with open_maybe_gzip(fasta_path, "rt") as input_handle:
        for line in input_handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            sequence_parts.append(line)
    return "".join(sequence_parts)


def get_normal_bases(sequence_type: str) -> set[str]:
    """Return the canonical alphabet used to separate normal and ambiguous positions."""
    sequence_type = sequence_type.upper()
    if sequence_type == "DNA":
        return set("ATCG")
    if sequence_type == "RNA":
        return set("AUCG")
    if sequence_type == "NT":
        return set("ATUCG")
    if sequence_type == "AA":
        return set("ARNDCQEGHILKMFPSTWYV")
    raise ValueError(f"Unsupported --seqtype: {sequence_type}")


def summarize_sequence(sequence: str, normal_bases: set[str]) -> tuple[int, int, int, int]:
    """Count uppercase normal bases, lowercase normal bases, ambiguous symbols, and gaps."""
    high_confidence = 0
    low_confidence = 0
    ambiguous = 0
    gaps = 0

    for character in sequence:
        if character == "-":
            gaps += 1
        elif character.upper() in normal_bases:
            if character.isupper():
                high_confidence += 1
            else:
                low_confidence += 1
        else:
            ambiguous += 1

    return high_confidence, low_confidence, ambiguous, gaps


def expected_header(prefix: str) -> str:
    """Build the expected TSV header for the requested prefix."""
    return "\t".join(column.format(prefix=prefix) for column in EXPECTED_HEADER_TEMPLATE)


def header_matches(output_path: Path, prefix: str) -> bool:
    """Return True when the existing output file already has the expected header."""
    if not output_path.exists() or output_path.stat().st_size == 0:
        return False

    with output_path.open("r", encoding="utf-8") as input_handle:
        first_line = input_handle.readline().rstrip("\n")
    return first_line == expected_header(prefix)


def write_summary(output_path: Path, prefix: str, row: list[str]) -> None:
    """Append to a matching TSV file or overwrite when the header layout changes."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    append_mode = header_matches(output_path, prefix)
    mode = "a" if append_mode else "w"

    with output_path.open(mode, newline="", encoding="utf-8") as output_handle:
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        if not append_mode:
            writer.writerow(expected_header(prefix).split("\t"))
        writer.writerow(row)

    logger.info("FASTA summary written to %s", output_path)


# ========================================================================================
# Main Entry Point
# ========================================================================================

def main() -> None:
    """Parse CLI arguments and export the FASTA composition summary."""
    args = build_parser().parse_args()
    input_path = Path(args.input)
    output_path = Path(args.outputtsv)

    if not input_path.is_file():
        raise SystemExit(f"Input FASTA file was not found: {input_path}")

    sample_name = get_sample_name(input_path, args.samplename)
    sequence = parse_fasta_sequence(input_path)
    normal_bases = get_normal_bases(args.seqtype)
    high_confidence, low_confidence, ambiguous, gaps = summarize_sequence(sequence, normal_bases)

    # Export one summary row using the same append-if-header-matches pattern as other summary utilities.
    write_summary(
        output_path,
        args.prefix,
        [
            sample_name,
            str(high_confidence),
            str(low_confidence),
            str(ambiguous),
            str(gaps),
        ],
    )


if __name__ == "__main__":
    start_time = time.time()
    print_header("get_fasta_summary.py")
    main()
    print_summary(start_time, "get_fasta_summary.py")
