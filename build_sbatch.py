#!/usr/bin/env python3
"""
Script Name:  build_sbatch.py
Description:  Generate per-sample SBATCH wrappers and split input CSV files.
Author:       Dr. Haiting Chai
Last update:  2026-MAY-08
Version:      ver.2026-05-08

DEPENDENCIES:
    - Python 3.10
    - __future__ (standard library)
    - argparse (standard library)
    - csv (standard library)
    - logging (standard library)
    - os (standard library)
    - pathlib (standard library)
    - shlex (standard library)
    - sys (standard library)
    - time (standard library)
    - typing (standard library)
    - assets.utils_banner (project module)
"""

from __future__ import annotations

# ========================================================================================
# Imports
# ========================================================================================

import argparse
import csv
import logging
import os
import shlex
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from assets.utils_banner import print_header, print_summary, setup_logging


# ========================================================================================
# Configuration
# ========================================================================================

logger = setup_logging()

REQUIRED_SAMPLE_HEADERS = ["job_id", "sample_name", "r1_path", "r2_path"]


# ========================================================================================
# Helper Functions
# ========================================================================================

def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for SBATCH wrapper generation."""
    parser = argparse.ArgumentParser(
        description="Generate SBATCH scripts and per-sample input CSV files.",
    )
    parser.add_argument("--inputcsv", type=Path, required=True, help="Sample list CSV file")
    parser.add_argument("--rundir", type=Path, required=True, help="Run directory")
    parser.add_argument("--human_db", required=True, help="Human database path")
    parser.add_argument("--target_reference",  type=Path, required=True, help="Target reference path")
    parser.add_argument("--fwd_primer_seq",  type=Path, required=False, help="forward primer sequence to trim in R1 file")
    parser.add_argument("--rev_primer_seq",  type=Path, required=False, help="reverse primer sequence to trim in R2 file")
    parser.add_argument(
        "--target_genome_type",
        choices=["linear", "circular"],
        default="linear",
        help="Target genome type",
    )
    parser.add_argument(
        "--dedup",
        choices=["yes", "no"],
        default="yes",
        help="Deduplication flag",
    )
    parser.add_argument("--sbatcha", required=True, help="SBATCH account")
    parser.add_argument(
        "--sbatchp",
        choices=["short", "long"],
        default="short",
        help="SBATCH partition",
    )
    parser.add_argument("--pipelinepath", type=Path, required=True, help="Pipeline directory")
    return parser


def shell_quote(value: str | Path) -> str:
    """Quote filesystem values before embedding them in generated shell scripts."""
    return shlex.quote(str(value))


def load_samples(input_csv: Path) -> tuple[list[str], list[dict[str, str]]]:
    """
    Load and validate the sample sheet used to generate per-sample SBATCH jobs.

    Parameters
    ----------
    input_csv : Path
        Input CSV containing paired-end sample metadata.

    Returns
    -------
    tuple[list[str], list[dict[str, str]]]
        CSV header order and validated sample rows.
    """
    with input_csv.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{input_csv} is missing a header row")

        missing_headers = [name for name in REQUIRED_SAMPLE_HEADERS if name not in reader.fieldnames]
        if missing_headers:
            raise ValueError(
                f"{input_csv} must contain these header columns: {', '.join(REQUIRED_SAMPLE_HEADERS)}. "
                f"Missing: {', '.join(missing_headers)}"
            )

        rows = [row for row in reader if row]
        if not rows:
            raise ValueError(f"{input_csv} contains no sample rows")

        for row_number, row in enumerate(rows, start=2):
            job_id = (row.get("job_id") or "").strip()
            if not job_id:
                raise ValueError(f"Missing job_id in row {row_number}")

        return list(reader.fieldnames), rows


def write_one_row_csv(output_csv: Path, fieldnames: list[str], row: dict[str, str]) -> None:
    """Write a one-sample CSV that can be passed directly to the workflow."""
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(row)


def make_sbatch_script(
    *,
    job_id: str,
    rundir: Path,
    sbatch_account: str,
    sbatch_partition: str,
    pipeline_path: Path,
    human_db: str,
    fwd_primer_seq: str,
    rev_primer_seq: str,
    target_reference: str,
    target_genome_type: str,
    dedup: str,
) -> str:
    """Render the SBATCH wrapper that launches Nextflow for one sample."""
    backslash = "\\"
    lines = [
        "#!/usr/bin/env bash",
        "",
        f"#SBATCH -A {sbatch_account}",
        f"#SBATCH -J {job_id[-6:]}",
        f"#SBATCH -o {shell_quote(rundir / 'sbatch' / f'{job_id}.sbatch.out')}",
        f"#SBATCH -e {shell_quote(rundir / 'sbatch' / f'{job_id}.sbatch.err')}",
        f"#SBATCH -p {sbatch_partition}",
        "",
        f"cd {shell_quote(rundir / job_id)}",
        "",
        f"{shell_quote(pipeline_path / 'nextflow')} -log {shell_quote(rundir / job_id / f'{job_id}.log')} {backslash}",
        f"  run {shell_quote(pipeline_path / 'main.nf')} {backslash}",
        f"  --inputcsv {shell_quote(rundir / job_id / f'{job_id}.input.csv')} {backslash}",
        f"  --outdir {shell_quote(rundir / job_id)} {backslash}",
        f"  --human_db {shell_quote(human_db)} {backslash}",
        f"  --fwd_primer_seq {shell_quote(fwd_primer_seq)} {backslash}",
        f"  --rev_primer_seq {shell_quote(rev_primer_seq)} {backslash}",
        f"  --target_reference {shell_quote(target_reference)} {backslash}",
        f"  --target_genome_type {shell_quote(target_genome_type)} {backslash}",
        f"  --dedup {shell_quote(dedup)} {backslash}",
        f"  --bq_thr 13 {backslash}",
        f"  --depth_thr 5 {backslash}",
        "  --qual_thr 30"
    ]
    return "\n".join(lines) + "\n"


# ========================================================================================
# Main Entry Point
# ========================================================================================

def main() -> None:
    """Generate per-sample workflow inputs and the matching SBATCH helpers."""
    args = build_parser().parse_args()

    input_csv = args.inputcsv.expanduser().resolve()
    rundir = args.rundir.expanduser().resolve()
    pipeline_path = args.pipelinepath.expanduser().resolve()

    logger.info("Loading sample sheet: %s", input_csv)
    fieldnames, rows = load_samples(input_csv)

    rundir.mkdir(parents=True, exist_ok=True)
    sbatch_dir = rundir / "sbatch"
    sbatch_dir.mkdir(parents=True, exist_ok=True)

    submit_lines: list[str] = []

    for row in rows:
        job_id = row["job_id"].strip()
        job_dir = rundir / job_id
        job_dir.mkdir(parents=True, exist_ok=True)

        per_sample_csv = job_dir / f"{job_id}.input.csv"
        write_one_row_csv(per_sample_csv, fieldnames, row)

        sbatch_script = make_sbatch_script(
            job_id=job_id,
            rundir=rundir,
            sbatch_account=args.sbatcha,
            sbatch_partition=args.sbatchp,
            pipeline_path=pipeline_path,
            human_db=args.human_db,
            target_reference=args.target_reference,
            fwd_primer_seq=args.fwd_primer_seq,
            rev_primer_seq=args.rev_primer_seq,
            target_genome_type=args.target_genome_type,
            dedup=args.dedup,
        )
        sbatch_path = sbatch_dir / f"{job_id}.sbatch.sh"
        sbatch_path.write_text(sbatch_script, encoding="utf-8")
        os.chmod(sbatch_path, 0o755)
        logger.info("Wrote SBATCH wrapper: %s", sbatch_path)

        submit_lines.append(f"sbatch --cpus-per-task=1 --mem=16G {shell_quote(sbatch_path)}")

    submit_path = rundir / "submit_sbatch.sh"
    submit_path.write_text("\n".join(submit_lines) + "\n", encoding="utf-8")
    os.chmod(submit_path, 0o755)
    logger.info("Wrote batch submit helper: %s", submit_path)


if __name__ == "__main__":
    start_time = time.time()
    print_header("build_sbatch.py")
    main()
    print_summary(start_time, "build_sbatch.py")
