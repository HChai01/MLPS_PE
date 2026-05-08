#!/usr/bin/env python3
"""
Script Name:  extract_top_references.py
Description:  Extract the top-N FASTA records named in a reference summary TSV.
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
        description="Extract the top-N FASTA records listed in a TSV summary.",
    )
    parser.add_argument("--tsv", required=True, help="Summary TSV with reference ids in column 1")
    parser.add_argument("--fasta", required=True, help="Source FASTA file")
    parser.add_argument("--out", required=True, help="Output FASTA file")
    parser.add_argument("--topn", type=int, default=1, help="Number of top rows to extract")
    return parser


def get_target_ids(tsv_path: Path, top_n: int) -> list[str]:
    target_ids: list[str] = []
    with tsv_path.open(encoding="utf-8") as handle:
        header = next(handle, None)
        if header is None:
            raise SystemExit(f"TSV file is empty: {tsv_path}")
        for index, line in enumerate(handle, start=1):
            if index > top_n:
                break
            parts = line.strip().split()
            if parts:
                target_ids.append(parts[0])

    if not target_ids:
        raise SystemExit(f"No target ids were found in: {tsv_path}")
    return target_ids


def extract_fastas(fasta_path: Path, target_ids: list[str], output_path: Path) -> None:
    targets = set(target_ids)
    found_ids: set[str] = set()
    write_sequence = False

    with fasta_path.open(encoding="utf-8") as fasta_handle, output_path.open("w", encoding="utf-8") as output_handle:
        for line in fasta_handle:
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                write_sequence = current_id in targets
                if write_sequence:
                    found_ids.add(current_id)
            if write_sequence:
                output_handle.write(line)

    missing_ids = sorted(targets - found_ids)
    for missing_id in missing_ids:
        logger.warning("Reference '%s' was not found in %s", missing_id, fasta_path)

    logger.info("Selected FASTA records written to %s", output_path)


def main() -> None:
    args = build_parser().parse_args()
    if args.topn < 1:
        raise SystemExit("--topn must be >= 1")

    target_ids = get_target_ids(Path(args.tsv), args.topn)
    extract_fastas(Path(args.fasta), target_ids, Path(args.out))


if __name__ == "__main__":
    start_time = time.time()
    print_header("extract_top_references.py")
    main()
    print_summary(start_time, "extract_top_references.py")
