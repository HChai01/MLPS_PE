#!/usr/bin/env python3
"""
Script Name:  cleanup_tmp_data.py
Description:  Truncate oversized regular files in traced Nextflow work directories.
Author:       Dr. Haiting Chai
Last update:  2026-APR-12
Version:      ver.2026-04-12
"""

from __future__ import annotations

import argparse
import csv
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
        description="Truncate large regular files in traced Nextflow task work directories.",
    )
    parser.add_argument("--trace", type=Path, required=True, help="Nextflow trace.txt path")
    parser.add_argument("--workdir", type=Path, required=True, help="Nextflow work directory")
    parser.add_argument("--max-bytes", type=int, required=True, help="Maximum file size to keep")
    return parser


def read_traced_task_dirs(trace_path: Path, workdir: Path) -> list[Path]:
    if not trace_path.exists():
        logger.warning("Trace file not found: %s", trace_path)
        return []

    task_dirs: list[Path] = []
    seen: set[Path] = set()
    workdir_resolved = workdir.resolve()

    with trace_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            task_hash = (row.get("hash") or "").strip()
            if "/" not in task_hash:
                continue

            shard, task_id = task_hash.split("/", 1)
            task_dir = (workdir / shard / task_id).resolve()
            if not task_dir.exists() or not task_dir.is_dir():
                continue

            try:
                task_dir.relative_to(workdir_resolved)
            except ValueError:
                logger.warning("Skipping task directory outside workdir: %s", task_dir)
                continue

            if task_dir not in seen:
                seen.add(task_dir)
                task_dirs.append(task_dir)

    return task_dirs


def truncate_large_files(task_dirs: list[Path], workdir: Path, max_bytes: int) -> tuple[int, int]:
    truncated_files = 0
    bytes_reclaimed = 0
    workdir_resolved = workdir.resolve()

    for task_dir in task_dirs:
        for path in task_dir.rglob("*"):
            if path.is_symlink() or not path.is_file():
                continue

            resolved_path = path.resolve()
            try:
                resolved_path.relative_to(workdir_resolved)
            except ValueError:
                logger.warning("Skipping file outside workdir: %s", resolved_path)
                continue

            size = path.stat().st_size
            if size <= max_bytes:
                continue

            path.write_bytes(b"")
            truncated_files += 1
            bytes_reclaimed += size
            logger.info("Truncated %s (%d bytes)", path, size)

    return truncated_files, bytes_reclaimed


def main() -> None:
    args = build_parser().parse_args()
    workdir = args.workdir.resolve()

    if args.max_bytes < 0:
        raise SystemExit("--max-bytes must be >= 0")
    if not workdir.exists() or not workdir.is_dir():
        raise SystemExit(f"--workdir is not a directory: {args.workdir}")

    task_dirs = read_traced_task_dirs(args.trace, workdir)
    if not task_dirs:
        logger.warning("No traced task directories were available for cleanup")
        return

    truncated_files, bytes_reclaimed = truncate_large_files(task_dirs, workdir, args.max_bytes)
    logger.info(
        "Cleanup complete across %d traced task directories; truncated %d file(s), reclaimed %d bytes",
        len(task_dirs),
        truncated_files,
        bytes_reclaimed,
    )


if __name__ == "__main__":
    start_time = time.time()
    print_header("cleanup_tmp_data.py")
    main()
    print_summary(start_time, "cleanup_tmp_data.py")
