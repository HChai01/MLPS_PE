#!/usr/bin/env python3
"""
Script Name:  format_fasta.py
Description:  Convert a FASTA file into a two-line normalized FASTA with optional sequence cleanup.
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
    - assets.utils_banner (project module)
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from assets.utils_banner import parse_yes_no, print_header, print_summary, setup_logging


logger = setup_logging()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert multi-line FASTA to a normalized two-line FASTA.",
    )
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("--output", help="Output FASTA file")
    parser.add_argument("--split", help="Delimiter used to trim FASTA headers")
    parser.add_argument("--uppercase", default="no", help="Convert sequence to uppercase: yes/no")
    parser.add_argument("--removegap", default="no", help="Remove '-' and '.' gap characters: yes/no")
    parser.add_argument("--double", default="no", help="Duplicate the sequence body: yes/no")
    return parser


def split_header(header: str, user_delimiter: str | None = None) -> str:
    clean_header = header[1:]
    if user_delimiter:
        return clean_header.split(user_delimiter)[0]
    for delimiter in [",", " ", "|", "@"]:
        if delimiter in clean_header:
            return clean_header.split(delimiter)[0]
    return clean_header


def process_sequence(
    sequence: str,
    *,
    uppercase: bool = False,
    remove_gap: bool = False,
    double_sequence: bool = False,
) -> str:
    if remove_gap:
        sequence = sequence.replace("-", "").replace(".", "")
    if uppercase:
        sequence = sequence.upper()
    if double_sequence:
        sequence = sequence + sequence
    return sequence


def convert_fasta(
    input_path: Path,
    *,
    output_path: Path | None = None,
    split_char: str | None = None,
    uppercase: bool = False,
    remove_gap: bool = False,
    double_sequence: bool = False,
) -> Path:
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
                    sequence = process_sequence(
                        "".join(sequence_lines),
                        uppercase=uppercase,
                        remove_gap=remove_gap,
                        double_sequence=double_sequence,
                    )
                    output_handle.write(f">{split_header(header, split_char)}\n{sequence}\n")
                header = line
                sequence_lines = []
            else:
                sequence_lines.append(line)

        if header is not None:
            sequence = process_sequence(
                "".join(sequence_lines),
                uppercase=uppercase,
                remove_gap=remove_gap,
                double_sequence=double_sequence,
            )
            output_handle.write(f">{split_header(header, split_char)}\n{sequence}\n")

    logger.info("Formatted FASTA written to %s", output_path)
    return output_path


def main() -> None:
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
        double_sequence=parse_yes_no(args.double, param_name="--double"),
    )


if __name__ == "__main__":
    start_time = time.time()
    print_header("format_fasta.py")
    main()
    print_summary(start_time, "format_fasta.py")
