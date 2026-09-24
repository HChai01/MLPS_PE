#!/usr/bin/env python3
"""
Script Name:  format_fasta.py
Description:  Convert a FASTA file into a two-line normalized FASTA with optional sequence edits.
Author:       Dr. Haiting Chai
Last update:  2026-SEP-24
Version:      ver.2026-09-24

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

from assets.utils_banner import parse_yes_no, print_header, print_summary, setup_logging


# ========================================================================================
# Configuration
# ========================================================================================

logger = setup_logging()


# ========================================================================================
# Helper Functions
# ========================================================================================

def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for FASTA normalization."""
    parser = argparse.ArgumentParser(
        description="Convert multi-line FASTA to a normalized two-line FASTA.",
    )
    parser.add_argument("--input", help="Input FASTA file")
    parser.add_argument("--output", help="Output FASTA file")
    parser.add_argument("--split", help="Delimiter used to trim FASTA headers")
    parser.add_argument("--uppercase", default="no", help="Convert sequence to uppercase: yes/no, default: no")
    parser.add_argument("--removegap", default="no", help="Remove '-' and '.' gap characters: yes/no, , default: no")
    parser.add_argument("--rmanchor", default="yes", help="Remove '^' and '$' anchor characters: yes/no, default: yes")
    parser.add_argument("--revcomp", default="no", type=str.lower, choices=("no", "dna", "rna"), help="Reverse-complement using DNA or RNA bases: no/dna/rna, , default: no")
    parser.add_argument("--double", default="no", help="Duplicate the sequence body: yes/no, default: no")
    parser.add_argument("--5pedit", dest="edit_5p", default="", help="Sequence to add at the 5' end, default:")
    parser.add_argument("--3pedit", dest="edit_3p", default="", help="Sequence to add at the 3' end, default:")
    return parser


def split_header(header: str, user_delimiter: str | None = None) -> str:
    """Trim FASTA headers to the first useful token for downstream tools."""
    clean_header = header[1:]
    if user_delimiter:
        return clean_header.split(user_delimiter)[0]
    for delimiter in [",", " ", "|", "@"]:
        if delimiter in clean_header:
            return clean_header.split(delimiter)[0]
    return clean_header


DNA_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")
RNA_COMPLEMENT = str.maketrans("ACGURYSWKMBDHVNacguryswkmbdhvn", "UGCAYRSWMKVHDBNugcayrswmkvhdbn")


def reverse_complement(sequence: str, sequence_type: str) -> str:
    """Return a DNA or RNA reverse complement, including IUPAC ambiguity codes."""
    complement = DNA_COMPLEMENT if sequence_type == "dna" else RNA_COMPLEMENT
    return sequence.translate(complement)[::-1]


def process_sequence(
    sequence: str,
    *,
    uppercase: bool = False,
    remove_gap: bool = False,
    remove_anchor: bool = True,
    revcomp: str = "no",
    edit_5p: str = "",
    edit_3p: str = "",
    double_sequence: bool = False,
) -> str:
    """Apply the optional sequence normalization steps requested by the CLI."""
    if remove_gap:
        sequence = sequence.replace("-", "").replace(".", "")
    if remove_anchor:
        sequence = sequence.replace("^", "").replace("$", "")
    if uppercase:
        sequence = sequence.upper()
    if revcomp != "no":
        sequence = reverse_complement(sequence, revcomp)
    if double_sequence:
        sequence = sequence + sequence
    sequence = edit_5p + sequence + edit_3p
    return sequence


def convert_fasta(
    input_path: Path,
    *,
    output_path: Path | None = None,
    split_char: str | None = None,
    uppercase: bool = False,
    remove_gap: bool = False,
    remove_anchor: bool = False,
    revcomp: str = "no",
    edit_5p: str = "",
    edit_3p: str = "",
    double_sequence: bool = False,
) -> Path:
    """Convert a FASTA file into a normalized two-line FASTA output."""
    if output_path is None:
        if input_path.suffix == ".fasta":
            output_path = input_path.with_name(input_path.stem + ".format.fasta")
        else:
            output_path = input_path.with_suffix(input_path.suffix + ".format.fasta")

    header: str | None = None
    sequence_lines: list[str] = []

    with input_path.open(encoding="utf-8") as input_handle, output_path.open("w", encoding="utf-8") as output_handle:
        for raw_line in input_handle:
            line = raw_line.strip()
            if line.startswith(">"):
                if header is not None:
                    # Emit each completed record as soon as the next header is seen so
                    # large FASTA files can be streamed without loading them fully.
                    sequence = process_sequence(
                        "".join(sequence_lines),
                        uppercase=uppercase,
                        remove_gap=remove_gap,
                        remove_anchor=remove_anchor,
                        revcomp=revcomp,
                        edit_5p=edit_5p,
                        edit_3p=edit_3p,
                        double_sequence=double_sequence,
                    )
                    formatted_header = split_header(header, split_char)
                    if revcomp != "no":
                        formatted_header += "_RC"
                    output_handle.write(f">{formatted_header}\n{sequence}\n")
                header = line
                sequence_lines = []
            else:
                sequence_lines.append(line)

        if header is not None:
            sequence = process_sequence(
                "".join(sequence_lines),
                uppercase=uppercase,
                remove_gap=remove_gap,
                remove_anchor=remove_anchor,
                revcomp=revcomp,
                edit_5p=edit_5p,
                edit_3p=edit_3p,
                double_sequence=double_sequence,
            )
            formatted_header = split_header(header, split_char)
            if revcomp != "no":
                formatted_header += "_RC"
            output_handle.write(f">{formatted_header}\n{sequence}\n")

    logger.info("Formatted FASTA written to %s", output_path)
    return output_path


# ========================================================================================
# Main Entry Point
# ========================================================================================

def main() -> None:
    """Parse CLI arguments and normalize the requested FASTA file."""
    args = build_parser().parse_args()
    input_path = Path(args.input)
    if not input_path.is_file():
        raise SystemExit(f"Input FASTA file was not found: {input_path}")

    convert_fasta(
        input_path,
        output_path=Path(args.output) if args.output else None,
        split_char=args.split,
        uppercase=parse_yes_no(args.uppercase, param_name="--uppercase"),
        remove_gap=parse_yes_no(args.removegap, param_name="--removegap"),
        remove_anchor=parse_yes_no(args.rmanchor, param_name="--rmanchor"),
        revcomp=args.revcomp,
        edit_5p=args.edit_5p,
        edit_3p=args.edit_3p,
        double_sequence=parse_yes_no(args.double, param_name="--double"),
    )


if __name__ == "__main__":
    start_time = time.time()
    print_header("format_fasta.py")
    main()
    print_summary(start_time, "format_fasta.py")
