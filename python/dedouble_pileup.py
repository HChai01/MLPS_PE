#!/usr/bin/env python3
"""
Script Name:  dedouble_pileup.py
Description:  Collapse doubled-reference pileup coordinates back onto the original reference length.
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

from assets.utils_banner import print_header, print_summary, setup_logging


logger = setup_logging()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Collapse pileup positions from a doubled reference back to original coordinates.",
    )
    parser.add_argument("--input", required=True, help="Pileup generated from the doubled reference")
    parser.add_argument("--originalref", required=True, help="Original non-doubled FASTA")
    parser.add_argument("--output", required=True, help="Collapsed pileup output")
    return parser


def get_reference_length(fasta_path: Path) -> int:
    length = 0
    with fasta_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                length += len(line.strip())
    return length


def parse_pileup_line(line: str) -> tuple[str, int, str, int, str, str] | None:
    fields = line.rstrip().split("\t")
    if len(fields) < 6:
        return None
    reference, position, ref_base, depth, reads, qualities = fields[:6]
    return reference, int(position), ref_base, int(depth), reads, qualities


def merge_reads_with_star_logic(reads: str, qualities: str) -> tuple[str, str]:
    pairs = list(zip(reads, qualities))
    real_pairs = [(read, quality) for read, quality in pairs if read != "*"]
    star_pairs = [(read, quality) for read, quality in pairs if read == "*"]

    if not real_pairs:
        return ("*", star_pairs[0][1]) if star_pairs else ("", "")

    return (
        "".join(read for read, _ in real_pairs),
        "".join(quality for _, quality in real_pairs),
    )


def collapse_pileup(input_path: Path, original_ref: Path, output_path: Path) -> None:
    reference_length = get_reference_length(original_ref)
    pileup_by_position: dict[tuple[str, int], dict[str, list[str] | str]] = {}

    with input_path.open(encoding="utf-8") as handle:
        for line in handle:
            parsed = parse_pileup_line(line)
            if parsed is None:
                continue

            reference, position, ref_base, _depth, reads, qualities = parsed
            collapsed_position = position if position <= reference_length else position - reference_length
            key = (reference, collapsed_position)

            if key not in pileup_by_position:
                pileup_by_position[key] = {"ref_base": ref_base, "reads": [], "quals": []}

            pileup_by_position[key]["reads"].append(reads)
            pileup_by_position[key]["quals"].append(qualities)

    with output_path.open("w", encoding="utf-8") as output_handle:
        for reference, position in sorted(pileup_by_position, key=lambda item: item[1]):
            data = pileup_by_position[(reference, position)]
            merged_reads = "".join(data["reads"])
            merged_qualities = "".join(data["quals"])
            final_reads, final_qualities = merge_reads_with_star_logic(merged_reads, merged_qualities)
            depth = sum(1 for base in final_reads if base not in {"^", "$"})
            output_handle.write(
                f"{reference}\t{position}\t{data['ref_base']}\t{depth}\t{final_reads}\t{final_qualities}\n"
            )

    logger.info("Collapsed pileup written to %s", output_path)


def main() -> None:
    args = build_parser().parse_args()
    collapse_pileup(Path(args.input), Path(args.originalref), Path(args.output))


if __name__ == "__main__":
    start_time = time.time()
    print_header("dedouble_pileup.py")
    main()
    print_summary(start_time, "dedouble_pileup.py")
