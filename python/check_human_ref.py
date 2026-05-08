#!/usr/bin/env python3
"""
Script Name:  check_human_ref.py
Description:  Validate a Bowtie2 human reference prefix and optionally rebuild it.
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
    - subprocess (standard library)
    - sys (standard library)
    - tempfile (standard library)
    - time (standard library)
    - assets.utils_banner (project module)
"""

from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from assets.utils_banner import parse_yes_no, print_header, print_summary, setup_logging


logger = setup_logging()

EXPECTED_SUFFIXES = [
    ".1.bt2",
    ".2.bt2",
    ".3.bt2",
    ".4.bt2",
    ".rev.1.bt2",
    ".rev.2.bt2",
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Check and optionally rebuild a Bowtie2 human reference index.",
    )
    parser.add_argument("--humanref", required=True, help="Bowtie2 index prefix")
    parser.add_argument("--rebuild", default="no", help="Rebuild on failure: yes/no")
    parser.add_argument("--newrefdir", default="./", help="Directory for rebuilt reference")
    parser.add_argument(
        "--reflink",
        default=(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/"
            "GCF_009914755.1_T2T-CHM13v2.0/"
            "GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"
        ),
        help="Reference download link used when rebuilding",
    )
    parser.add_argument("--refprefix", default="t2t", help="Prefix name for rebuilt index")
    return parser


def run_command(command: list[str], *, cwd: Path | None = None) -> None:
    logger.info("Running command: %s", " ".join(command))
    subprocess.run(command, cwd=cwd, check=True)


def check_index_files(prefix: Path) -> bool:
    ok = True
    for suffix in EXPECTED_SUFFIXES:
        index_path = Path(f"{prefix}{suffix}")
        if not index_path.exists():
            logger.error("Missing index file: %s", index_path)
            ok = False
            continue
        if index_path.stat().st_size == 0:
            logger.error("Empty index file: %s", index_path)
            ok = False
    return ok


def check_bowtie2(prefix: Path) -> bool:
    try:
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as fastq_handle:
            fastq_handle.write("@r1\nACGT\n+\nFFFF\n")
            fastq_path = Path(fastq_handle.name)

        command = ["bowtie2", "-x", str(prefix), "-U", str(fastq_path), "-S", "/dev/null"]
        result = subprocess.run(
            command,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        fastq_path.unlink(missing_ok=True)
        if result.returncode != 0:
            logger.error("Bowtie2 alignment smoke test failed for prefix: %s", prefix)
            return False
        return True
    except FileNotFoundError:
        logger.error("bowtie2 is not available on PATH")
        return False


def rebuild_reference(ref_dir: Path, reflink: str, prefix: str) -> None:
    logger.info("Rebuilding reference under %s", ref_dir)
    ref_dir.mkdir(parents=True, exist_ok=True)

    gz_path = ref_dir / Path(reflink).name
    fasta_path = ref_dir / f"{prefix}.fa"

    run_command(["wget", "-c", "-O", str(gz_path), reflink])

    with fasta_path.open("w", encoding="utf-8") as output_handle:
        decompress = subprocess.Popen(["pigz", "-dc", str(gz_path)], stdout=subprocess.PIPE)
        header_cleanup = subprocess.Popen(
            ["awk", r"{if($0 ~ /^>/) print $1; else print}"],
            stdin=decompress.stdout,
            stdout=output_handle,
        )
        assert decompress.stdout is not None
        decompress.stdout.close()
        header_cleanup.communicate()
        if header_cleanup.returncode != 0:
            raise RuntimeError("FASTA header cleanup failed during rebuild")

    gz_path.unlink(missing_ok=True)

    threads = max(1, (os.cpu_count() or 1) - 1)
    run_command(
        [
            "bowtie2-build",
            "--threads",
            str(threads),
            str(fasta_path),
            str(ref_dir / prefix),
        ]
    )
    logger.info("Reference rebuild completed")


def main() -> None:
    args = build_parser().parse_args()
    prefix = Path(args.humanref)
    rebuild = parse_yes_no(args.rebuild, param_name="--rebuild")

    logger.info("Checking human reference prefix: %s", prefix)
    is_valid = check_index_files(prefix) and check_bowtie2(prefix)

    if is_valid:
        logger.info("Reference index is valid")
        return

    logger.warning("Reference validation failed")
    if not rebuild:
        raise SystemExit("Reference validation failed and rebuild is disabled")

    rebuild_reference(Path(args.newrefdir), args.reflink, args.refprefix)
    logger.info("Reference rebuilt. Re-run this script or the workflow to validate the new index.")


if __name__ == "__main__":
    start_time = time.time()
    print_header("check_human_ref.py")
    main()
    print_summary(start_time, "check_human_ref.py")
