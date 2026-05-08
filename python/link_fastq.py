#!/usr/bin/env python3
"""
Script Name:  link_fastq.py
Description:  Create canonical symlinks for single-end FASTQ inputs.
Author:       Dr. Haiting Chai
Last update:  2026-APR-11
Version:      ver.2026-04-11

DEPENDENCIES:
    - Python 3.10
    - __future__ (standard library)
    - argparse (standard library)
    - logging (standard library)
    - os (standard library)
    - pathlib (standard library)
    - sys (standard library)
    - time (standard library)
    - assets.utils_banner (project module)
"""

from __future__ import annotations

import argparse
import logging
import os
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
        description="Create symlinks for .fastq.gz and .fq.gz inputs using canonical FASTQ names.",
    )
    parser.add_argument("--inDir", required=True, help="Input directory containing raw FASTQ files")
    parser.add_argument("--outDir", required=True, help="Output directory where symlinks will be created")
    return parser


def canonical_name(path: Path) -> str:
    """Return the canonical `.fastq.gz` output name for an input FASTQ file."""
    if path.name.endswith(".fastq.gz"):
        return path.name
    if path.name.endswith(".fq.gz"):
        return path.name[:-6] + ".fastq.gz"
    raise ValueError(f"Unsupported FASTQ filename: {path.name}")


def create_symlink(source_path: Path, output_dir: Path) -> None:
    """Create one canonical symlink unless the destination already exists."""
    target_path = source_path.resolve()
    link_path = output_dir / canonical_name(source_path)

    if link_path.exists() or link_path.is_symlink():
        logger.info("Skipping existing link target: %s", link_path)
        return

    os.symlink(target_path, link_path)
    logger.info("Linked %s -> %s", link_path, target_path)


def create_symlinks(input_dir: Path, output_dir: Path) -> None:
    """Create symlinks for all supported single-end FASTQ files in the input directory."""
    input_path = input_dir.resolve()
    output_path = output_dir.resolve()

    if not input_path.is_dir():
        raise SystemExit(f"Input directory does not exist: {input_path}")

    output_path.mkdir(parents=True, exist_ok=True)

    # Preserve existing `.fastq.gz` names and normalize `.fq.gz` names to `.fastq.gz`.
    fastq_files = sorted(input_path.glob("*.fastq.gz")) + sorted(input_path.glob("*.fq.gz"))
    if not fastq_files:
        logger.warning("No .fastq.gz or .fq.gz files were found in %s", input_path)
        return

    for fastq_file in fastq_files:
        create_symlink(fastq_file, output_path)


def main() -> None:
    args = build_parser().parse_args()
    create_symlinks(Path(args.inDir), Path(args.outDir))


if __name__ == "__main__":
    start_time = time.time()
    print_header("link_fastq.py")
    main()
    print_summary(start_time, "link_fastq.py")
