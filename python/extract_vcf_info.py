#!/usr/bin/env python3
"""
Script Name:  extract_vcf_info.py
Description:  Extract key VCF columns and INFO-derived metrics into a TSV summary table.
Author:       Dr. Haiting Chai
Last update:  2026-APR-12
Version:      ver.2026-04-12

DEPENDENCIES:
    - Python 3.10
    - __future__ (standard library)
    - argparse (standard library)
    - csv (standard library)
    - gzip (standard library)
    - logging (standard library)
    - math (standard library)
    - pathlib (standard library)
    - sys (standard library)
    - time (standard library)
    - typing (standard library)
    - assets.utils_banner (project module)
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import math
import sys
import time
from pathlib import Path
from typing import TextIO

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from assets.utils_banner import print_header, print_summary, setup_logging


logger = setup_logging()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract CHROM, POS, REF, ALT, QUAL, DP, DP4, allele frequency, and error rate from a VCF.",
    )
    parser.add_argument("--input", required=True, help="Input .vcf or .vcf.gz file")
    parser.add_argument("--outputtsv", required=True, help="Output TSV file path")
    return parser


def parse_dp4(dp4_value: str | None) -> list[int | None]:
    """Parse the DP4 INFO value into ref/alt forward and reverse counts."""
    if not dp4_value or dp4_value == ".":
        return [None, None, None, None]

    parts = dp4_value.split(",")
    if len(parts) != 4:
        return [None, None, None, None]

    try:
        return [int(value) for value in parts]
    except ValueError:
        return [None, None, None, None]


def qual_to_error_rate(qual_value: str | None) -> str:
    """Convert QUAL to an estimated sequencing error rate."""
    if qual_value is None or qual_value == ".":
        return "NA"

    try:
        quality = float(qual_value)
    except ValueError:
        return "NA"

    return f"{math.pow(10.0, -quality / 10.0):.6g}"


def format_float(value: float | None) -> str:
    """Format floating-point outputs consistently for TSV export."""
    if value is None:
        return "NA"
    return f"{value:.6f}"


def open_text_maybe_gzip(path: Path) -> TextIO:
    """Open either a plain-text or gzipped VCF as a text stream."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def parse_info_field(info_field: str) -> dict[str, str | bool]:
    """Convert the VCF INFO column into a simple dictionary."""
    info_values: dict[str, str | bool] = {}
    for item in info_field.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info_values[key] = value
        else:
            info_values[item] = True
    return info_values


def compute_allele_frequency(dp4_values: list[int | None]) -> str:
    """Calculate alternate allele frequency from DP4 counts when available."""
    if not all(value is not None for value in dp4_values):
        return "NA"

    ref_total = dp4_values[0] + dp4_values[1]
    alt_total = dp4_values[2] + dp4_values[3]
    total = ref_total + alt_total
    if total == 0:
        return "NA"
    return format_float(alt_total / total)


def extract_vcf_rows(input_path: Path, output_path: Path) -> None:
    """Read the input VCF and write the selected fields to a TSV summary file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open_text_maybe_gzip(input_path) as input_handle, output_path.open("w", encoding="utf-8", newline="") as output_handle:
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "CHROM",
                "POS",
                "REF",
                "ALT",
                "QUAL",
                "DP",
                "Q13_ref_fwd",
                "Q13_ref_rev",
                "Q13_alt_fwd",
                "Q13_alt_rev",
                "allele_freq",
                "error_rate",
            ]
        )

        # Skip VCF metadata/header lines and transform each variant line into one TSV row.
        for line in input_handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                logger.warning("Skipping malformed VCF line: %s", line.rstrip("\n"))
                continue

            info_dict = parse_info_field(fields[7])
            dp4_values = parse_dp4(info_dict.get("DP4") if isinstance(info_dict.get("DP4"), str) else None)

            writer.writerow(
                [
                    fields[0],
                    fields[1],
                    fields[3],
                    fields[4],
                    fields[5],
                    info_dict.get("DP", "NA"),
                    dp4_values[0] if dp4_values[0] is not None else "NA",
                    dp4_values[1] if dp4_values[1] is not None else "NA",
                    dp4_values[2] if dp4_values[2] is not None else "NA",
                    dp4_values[3] if dp4_values[3] is not None else "NA",
                    compute_allele_frequency(dp4_values),
                    qual_to_error_rate(fields[5]),
                ]
            )

    logger.info("VCF summary written to %s", output_path)


def main() -> None:
    args = build_parser().parse_args()
    input_path = Path(args.input)
    output_path = Path(args.outputtsv)

    if not input_path.is_file():
        raise SystemExit(f"Input VCF file was not found: {input_path}")
    if output_path.suffix != ".tsv":
        raise SystemExit("--outputtsv must end with .tsv")

    extract_vcf_rows(input_path, output_path)


if __name__ == "__main__":
    start_time = time.time()
    print_header("extract_vcf_info.py")
    main()
    print_summary(start_time, "extract_vcf_info.py")
