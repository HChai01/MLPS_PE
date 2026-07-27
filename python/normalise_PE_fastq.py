#!/usr/bin/env python3
"""
Script Name:  normalise_PE_fastq.py
Description:  Paired-end digital normalization based on canonical k-mer median coverage.
Author:       Dr. Haiting Chai
Last update:  2026-MAY-08
Version:      ver.2026-05-08

DEPENDENCIES:
    - Python 3.10
    - __future__ (standard library)
    - argparse (standard library)
    - collections (standard library)
    - gzip (standard library)
    - itertools (standard library)
    - logging (standard library)
    - pathlib (standard library)
    - statistics (standard library)
    - sys (standard library)
    - time (standard library)
    - typing (standard library)
    - biopython / Bio.SeqIO (external library)
    - assets.utils_banner (project module)
"""

from __future__ import annotations

# ========================================================================================
# Imports
# ========================================================================================

import argparse
import gzip
import logging
import statistics
import sys
import time
from collections import defaultdict
from itertools import zip_longest
from pathlib import Path
from typing import Iterator

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from assets.utils_banner import print_header, print_summary, setup_logging


# ========================================================================================
# Configuration
# ========================================================================================

logger = setup_logging()


# ========================================================================================
# Helper Functions
# ========================================================================================

def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for paired-end digital normalization."""
    parser = argparse.ArgumentParser(
        description="Paired-end digital normalization for FASTQ / FASTQ.gz input.",
    )
    parser.add_argument("-1", "--fastq_r1", type=Path, required=True, help="Input FASTQ for read 1")
    parser.add_argument("-2", "--fastq_r2", type=Path, required=True, help="Input FASTQ for read 2")
    parser.add_argument("-o", "--out_r1", type=Path, required=True, help="Output normalized FASTQ for read 1")
    parser.add_argument("-O", "--out_r2", type=Path, required=True, help="Output normalized FASTQ for read 2")
    parser.add_argument("-k", "--kmer", type=int, default=21, help="k-mer size")
    parser.add_argument("-C", "--coverage", type=int, default=20, help="Target median k-mer coverage")
    return parser


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a nucleotide sequence."""
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def get_kmers(sequence: str, kmer_size: int) -> Iterator[str]:
    """Yield canonicalized k-mers so each strand contributes to the same count key."""
    if len(sequence) < kmer_size:
        return
    for start in range(len(sequence) - kmer_size + 1):
        kmer = sequence[start : start + kmer_size]
        reverse = reverse_complement(kmer)
        yield kmer if kmer <= reverse else reverse


def open_fastq(path: Path, mode: str = "rt"):
    """Open plain-text or gzipped FASTQ files transparently."""
    return gzip.open(path, mode) if path.suffix == ".gz" else path.open(mode, encoding="utf-8")


def normalise_read_id(read_id: str) -> str:
    """Strip `/1` and `/2` suffixes before validating paired read names."""
    if read_id.endswith("/1") or read_id.endswith("/2"):
        return read_id[:-2]
    return read_id


def scale_coverage(median_coverage: float, read_length: int, *, k_from: int, k_to: int, kept: bool) -> float:
    """Project the observed median coverage onto alternate k-mer sizes for reporting."""
    if not kept or median_coverage == 0 or read_length < k_to:
        return 0.0
    denominator = read_length - k_from + 1
    if denominator <= 0:
        return 0.0
    return median_coverage * (read_length - k_to + 1) / denominator


def write_kmer_summary(kmer_counts: dict[str, int], kmer_size: int, coverage_target: int) -> None:
    """Write the final canonical k-mer count table for audit and tuning purposes."""
    output_path = Path(f"kmer_summary_k{kmer_size}_C{coverage_target}.tsv.gz")
    with gzip.open(output_path, "wt", encoding="utf-8") as handle:
        handle.write(f"# k={kmer_size}\tC={coverage_target}\n")
        handle.write("kmer\tcount\n")
        for kmer, count in sorted(kmer_counts.items()):
            handle.write(f"{kmer}\t{count}\n")
    logger.info("K-mer summary written to %s", output_path)


def write_median_summary(
    records: list[tuple[str, float, str, int, int, float, float]],
    kmer_size: int,
    coverage_target: int,
) -> None:
    """Write one per-pair median-coverage record for downstream QC review."""
    output_path = Path(f"median_kmer_coverage_k{kmer_size}_C{coverage_target}.tsv.gz")
    with gzip.open(output_path, "wt", encoding="utf-8") as handle:
        handle.write(
            f"read_id\tmedian_kmer_cov_k{kmer_size}\tstatus\tlen_r1\tlen_r2\t"
            "estimated_median_cov_k33\testimated_median_cov_k55\n"
        )
        for read_id, median_cov, status, len_r1, len_r2, med33, med55 in records:
            handle.write(f"{read_id}\t{median_cov}\t{status}\t{len_r1}\t{len_r2}\t{med33}\t{med55}\n")
    logger.info("Per-pair median coverage written to %s", output_path)


def iter_paired_records(r1_handle, r2_handle) -> Iterator[tuple[int, SeqRecord, SeqRecord]]:
    """Iterate through R1/R2 reads in lockstep and fail fast on pairing mismatches."""
    for pair_index, (record_r1, record_r2) in enumerate(
        zip_longest(SeqIO.parse(r1_handle, "fastq"), SeqIO.parse(r2_handle, "fastq")),
        start=1,
    ):
        if record_r1 is None:
            raise ValueError(f"R1 ended before R2 at pair {pair_index}; extra R2 read: {record_r2.id}")
        if record_r2 is None:
            raise ValueError(f"R2 ended before R1 at pair {pair_index}; extra R1 read: {record_r1.id}")
        if normalise_read_id(record_r1.id) != normalise_read_id(record_r2.id):
            raise ValueError(
                f"Mismatched paired read IDs at pair {pair_index}: R1={record_r1.id} R2={record_r2.id}"
            )
        yield pair_index, record_r1, record_r2


def normalize_paired_fastq(
    input_r1: Path,
    input_r2: Path,
    output_r1: Path,
    output_r2: Path,
    *,
    kmer_size: int,
    coverage_target: int,
) -> None:
    """Perform digital normalization across paired reads using median canonical k-mer coverage."""
    kmer_counts: dict[str, int] = defaultdict(int)
    median_records: list[tuple[str, float, str, int, int, float, float]] = []
    kept_pairs = 0
    discarded_pairs = 0

    output_r1.parent.mkdir(parents=True, exist_ok=True)
    output_r2.parent.mkdir(parents=True, exist_ok=True)

    with open_fastq(input_r1) as input_handle_r1, open_fastq(input_r2) as input_handle_r2:
        with open_fastq(output_r1, "wt") as output_handle_r1, open_fastq(output_r2, "wt") as output_handle_r2:
            for _, record_r1, record_r2 in iter_paired_records(input_handle_r1, input_handle_r2):
                sequence_r1 = str(record_r1.seq)
                sequence_r2 = str(record_r2.seq)
                kmers = list(get_kmers(sequence_r1, kmer_size)) + list(get_kmers(sequence_r2, kmer_size))

                # Pairs shorter than the requested k-mer size cannot contribute counts,
                # so they are preserved unchanged and recorded separately.
                if not kmers:
                    SeqIO.write(record_r1, output_handle_r1, "fastq")
                    SeqIO.write(record_r2, output_handle_r2, "fastq")
                    kept_pairs += 1
                    median_records.append(
                        (record_r1.id, 0.0, "kept_too_short", len(sequence_r1), len(sequence_r2), 0.0, 0.0)
                    )
                    continue

                counts = [kmer_counts[kmer] for kmer in kmers]
                median_coverage = float(statistics.median(counts))
                contains_new_kmer = any(count == 0 for count in counts)
                keep_pair = median_coverage <= coverage_target or contains_new_kmer

                if keep_pair:
                    # Update the k-mer table only for retained pairs so future decisions
                    # reflect the coverage of the normalized read set, not discarded reads.
                    SeqIO.write(record_r1, output_handle_r1, "fastq")
                    SeqIO.write(record_r2, output_handle_r2, "fastq")
                    kept_pairs += 1
                    for kmer in kmers:
                        kmer_counts[kmer] += 1
                    status = "kept" if median_coverage <= coverage_target else "kept_for_new_kmer"
                else:
                    discarded_pairs += 1
                    status = "discarded"

                shortest_read = min(len(sequence_r1), len(sequence_r2))
                median_records.append(
                    (
                        record_r1.id,
                        median_coverage,
                        status,
                        len(sequence_r1),
                        len(sequence_r2),
                        round(
                            scale_coverage(
                                median_coverage,
                                shortest_read,
                                k_from=kmer_size,
                                k_to=33,
                                kept=keep_pair,
                            ),
                            2,
                        ),
                        round(
                            scale_coverage(
                                median_coverage,
                                shortest_read,
                                k_from=kmer_size,
                                k_to=55,
                                kept=keep_pair,
                            ),
                            2,
                        ),
                    )
                )

    logger.info("Normalization finished: kept=%s discarded=%s", kept_pairs, discarded_pairs)
    write_kmer_summary(kmer_counts, kmer_size, coverage_target)
    write_median_summary(median_records, kmer_size, coverage_target)


# ========================================================================================
# Main Entry Point
# ========================================================================================

def main() -> None:
    """Parse CLI arguments and normalize the paired FASTQ inputs."""
    args = build_parser().parse_args()
    normalize_paired_fastq(
        args.fastq_r1,
        args.fastq_r2,
        args.out_r1,
        args.out_r2,
        kmer_size=args.kmer,
        coverage_target=args.coverage,
    )


if __name__ == "__main__":
    start_time = time.time()
    print_header("normalise_PE_fastq.py")
    main()
    print_summary(start_time, "normalise_PE_fastq.py")
