#!/usr/bin/env python3
"""
Script Name:  merge_tsv.py
Description:  Merge multiple TSV files by a shared identifier column.
Author:       Dr. Haiting Chai
Last update:  2026-APR-12
Version:      ver.2026-04-12

DEPENDENCIES:
    - Python 3.10
    - __future__ (standard library)
    - argparse (standard library)
    - csv (standard library)
    - logging (standard library)
    - pathlib (standard library)
    - sys (standard library)
    - time (standard library)
    - assets.utils_banner (project module)
"""

from __future__ import annotations

import argparse
import csv
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
        description="Merge all TSV files in a folder by a shared ID column.",
    )
    parser.add_argument("--indir", required=True, type=Path, help="Input directory containing TSV files")
    parser.add_argument("--id", dest="id_col", required=True, help="Shared identifier column name")
    parser.add_argument("--outputtsv", required=True, type=Path, help="Output TSV file path")
    return parser


def merge_tsv_files(input_dir: Path, id_column: str, output_tsv: Path) -> None:
    """Merge TSV rows by identifier, preserving first-seen row and column order."""
    tsv_files = sorted(
        path for path in input_dir.glob("*.tsv") if path.resolve() != output_tsv.resolve()
    )
    if not tsv_files:
        logger.warning("No .tsv files were found in %s", input_dir)
        return

    merged_rows: dict[str, dict[str, str]] = {}
    sample_order: list[str] = []
    ordered_columns: list[str] = []
    seen_columns: set[str] = set()

    for tsv_file in tsv_files:
        with tsv_file.open("r", newline="", encoding="utf-8") as input_handle:
            reader = csv.DictReader(input_handle, delimiter="\t")
            if not reader.fieldnames:
                logger.warning("Skipping %s: empty file or missing header", tsv_file.name)
                continue
            if id_column not in reader.fieldnames:
                logger.warning("Skipping %s: missing identifier column '%s'", tsv_file.name, id_column)
                continue

            file_columns = [column for column in reader.fieldnames if column != id_column]

            # Preserve output column order by first appearance across the input files.
            for column in file_columns:
                if column not in seen_columns:
                    seen_columns.add(column)
                    ordered_columns.append(column)

            # Later files can overwrite earlier values for the same sample/column pair.
            for row in reader:
                sample_id = row.get(id_column, "")
                if sample_id == "":
                    continue

                if sample_id not in merged_rows:
                    merged_rows[sample_id] = {}
                    sample_order.append(sample_id)

                for column in file_columns:
                    value = row.get(column, "")
                    if value != "":
                        merged_rows[sample_id][column] = value

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", newline="", encoding="utf-8") as output_handle:
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow([id_column] + ordered_columns)

        for sample_id in sample_order:
            record = merged_rows.get(sample_id, {})
            writer.writerow([sample_id] + [record.get(column, "NA") for column in ordered_columns])

    logger.info("Merged TSV written to %s", output_tsv)


def main() -> None:
    args = build_parser().parse_args()
    if not args.indir.is_dir():
        raise SystemExit(f"--indir is not a directory: {args.indir}")

    merge_tsv_files(args.indir, args.id_col, args.outputtsv)


if __name__ == "__main__":
    start_time = time.time()
    print_header("merge_tsv.py")
    main()
    print_summary(start_time, "merge_tsv.py")
