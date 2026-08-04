#!/usr/bin/env python3
"""
Script Name:  em_real_location.py
Description:  Expectation-maximization reallocation for multi-mapped reads using BAM input.
Author:       Dr. Haiting Chai
Last update:  2026-APR-11
Version:      ver.2026-04-11

DEPENDENCIES:
    - Python 3.10
    - __future__ (standard library)
    - argparse (standard library)
    - collections (standard library)
    - logging (standard library)
    - math (standard library)
    - pathlib (standard library)
    - sys (standard library)
    - time (standard library)
    - pysam (external library)
    - assets.utils_banner (project module)
"""

from __future__ import annotations

# ========================================================================================
# Imports
# ========================================================================================

import argparse
import logging
import math
import sys
import time
from collections import defaultdict
from pathlib import Path

import pysam

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
    """Build the CLI parser for EM-based multi-mapper reallocation."""
    parser = argparse.ArgumentParser(
        description="Reallocate multi-mapped reads across references with EM.",
    )
    parser.add_argument("--bam", required=True, help="Input BAM containing primary and secondary alignments")
    parser.add_argument("--power_x", type=float, default=1.0, help="Exponent applied to current abundance estimates")
    parser.add_argument("--beta", type=float, default=0.1, help="Scaling factor for score-derived probabilities")
    parser.add_argument("--gamma", type=float, default=0.02, help="Clipping penalty for the HYBRID model")
    parser.add_argument(
        "--model",
        default="AS+NM",
        choices=["AS", "NM", "AS+NM", "MAPQ", "HYBRID"],
        help="Likelihood model to use: AS/NM/AS+NM/MAPQ/HYBRID",
    )
    parser.add_argument("--max_iter", type=int, default=50, help="Maximum number of EM iterations")
    parser.add_argument("--converge", type=float, default=1e-4, help="Convergence threshold")
    parser.add_argument("--out", default="em_assignments.tsv", help="Per-read posterior output TSV")
    return parser


def load_bam_alignments(bam_path: Path) -> tuple[dict[str, list[dict[str, int | str | None]]], list[str]]:
    """Load mapped BAM alignments and capture the fields needed by the EM model."""
    reads: dict[str, list[dict[str, int | str | None]]] = defaultdict(list)
    reference_names: set[str] = set()

    with pysam.AlignmentFile(bam_path, "rb") as bam_handle:
        for alignment in bam_handle.fetch(until_eof=True):
            if alignment.is_unmapped:
                continue

            reference_name = bam_handle.get_reference_name(alignment.reference_id)
            clipped_bases = 0
            if alignment.cigartuples:
                for operation, length in alignment.cigartuples:
                    if operation in (4, 5):
                        clipped_bases += length

            reads[alignment.query_name].append(
                {
                    "ref": reference_name,
                    "AS": alignment.get_tag("AS") if alignment.has_tag("AS") else None,
                    "NM": alignment.get_tag("NM") if alignment.has_tag("NM") else None,
                    "MAPQ": alignment.mapping_quality,
                    "CLIP": clipped_bases,
                }
            )
            reference_names.add(reference_name)

    if not reads:
        raise SystemExit("No mapped reads were found in the BAM file")

    sorted_references = sorted(reference_names)
    logger.info("Loaded %s reads across %s references", len(reads), len(sorted_references))
    return reads, sorted_references


def score_alignment_info(alignment: dict[str, int | str | None]) -> tuple[float, float, float, float]:
    """Return a sortable tuple that prefers stronger alignments for each reference."""
    as_score = float(alignment["AS"]) if alignment["AS"] is not None else -1e9
    nm_score = -float(alignment["NM"]) if alignment["NM"] is not None else -1e9
    mapq_score = float(alignment["MAPQ"]) if alignment["MAPQ"] is not None else 0.0
    clip_score = -float(alignment["CLIP"])
    return as_score, nm_score, mapq_score, clip_score


def compute_likelihood(
    abundance: float,
    alignment: dict[str, int | str | None],
    *,
    model: str,
    beta: float,
    gamma: float,
    power_x: float,
) -> float:
    """Compute one alignment likelihood under the selected scoring model."""
    as_score = float(alignment["AS"]) if alignment["AS"] is not None else 0.0
    nm_score = float(alignment["NM"]) if alignment["NM"] is not None else 0.0
    mapq_score = float(alignment["MAPQ"]) if alignment["MAPQ"] is not None else 0.0
    clipped_bases = float(alignment["CLIP"])

    if model == "AS":
        score_term = math.exp(beta * as_score)
    elif model == "NM":
        score_term = math.exp(-beta * nm_score)
    elif model == "MAPQ":
        score_term = math.exp(beta * mapq_score)
    elif model == "HYBRID":
        score_term = math.exp(beta * as_score) * math.exp(-gamma * clipped_bases)
    else:
        score_term = math.exp(beta * as_score) * math.exp(-beta * nm_score)

    return (abundance**power_x) * score_term


def run_em(
    reads: dict[str, list[dict[str, int | str | None]]],
    reference_names: list[str],
    *,
    beta: float,
    gamma: float,
    max_iter: int,
    converge: float,
    model: str,
    power_x: float,
) -> tuple[dict[str, float], dict[str, dict[str, float]]]:
    """Run expectation-maximization to estimate per-reference read abundances."""
    abundances = {reference: 1.0 / len(reference_names) for reference in reference_names}
    posteriors: dict[str, dict[str, float]] = defaultdict(dict)

    for iteration in range(1, max_iter + 1):
        posteriors = defaultdict(dict)

        for read_id, alignments in reads.items():
            grouped: dict[str, list[dict[str, int | str | None]]] = defaultdict(list)
            for alignment in alignments:
                grouped[str(alignment["ref"])].append(alignment)

            likelihoods: list[tuple[str, float]] = []
            for reference, entries in grouped.items():
                # Multiple alignments to the same reference are reduced to the strongest
                # representative so the EM step weights references rather than duplicates.
                best_alignment = max(entries, key=score_alignment_info)
                likelihood = compute_likelihood(
                    abundances[reference],
                    best_alignment,
                    model=model,
                    beta=beta,
                    gamma=gamma,
                    power_x=power_x,
                )
                if math.isnan(likelihood) or likelihood <= 0.0:
                    likelihood = 0.0
                likelihoods.append((reference, likelihood))

            total_likelihood = sum(value for _, value in likelihoods)
            if total_likelihood == 0:
                continue

            for reference, likelihood in likelihoods:
                posteriors[reference][read_id] = likelihood / total_likelihood

        total_reads = len(reads)
        new_abundances = {
            reference: sum(posteriors[reference].values()) / total_reads for reference in reference_names
        }
        delta = sum(abs(new_abundances[reference] - abundances[reference]) for reference in reference_names)
        logger.info("Iteration %s delta=%.6f threshold=%.6f", iteration, delta, converge)
        abundances = new_abundances
        if delta < converge:
            break

    return abundances, posteriors


def save_outputs(
    abundances: dict[str, float],
    posteriors: dict[str, dict[str, float]],
    reference_names: list[str],
    reads: dict[str, list[dict[str, int | str | None]]],
    output_path: Path,
) -> None:
    """Write per-read posteriors and the reference-level abundance summary."""
    total_reads = len(reads)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("read_id\tref\tposterior\n")
        for reference in reference_names:
            for read_id, posterior in posteriors[reference].items():
                handle.write(f"{read_id}\t{reference}\t{posterior:.6f}\n")

    summary_path = output_path.with_name(output_path.stem + "_ref_summary.tsv")
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("Reference\tTotal_posterior_mass\tRelative_abundance\n")
        for reference, abundance in sorted(abundances.items(), key=lambda item: item[1], reverse=True):
            total_mass = abundance * total_reads
            handle.write(f"{reference}\t{total_mass:.6f}\t{abundance:.6f}\n")

    logger.info("Saved per-read assignments to %s", output_path)
    logger.info("Saved reference summary to %s", summary_path)


# ========================================================================================
# Main Entry Point
# ========================================================================================

def main() -> None:
    """Parse CLI arguments, run EM, and export the resulting assignments."""
    args = build_parser().parse_args()
    bam_path = Path(args.bam)
    if not bam_path.is_file():
        raise SystemExit(f"BAM file was not found: {bam_path}")

    reads, reference_names = load_bam_alignments(bam_path)
    abundances, posteriors = run_em(
        reads,
        reference_names,
        beta=args.beta,
        gamma=args.gamma,
        max_iter=args.max_iter,
        converge=args.converge,
        model=args.model,
        power_x=args.power_x,
    )
    save_outputs(abundances, posteriors, reference_names, reads, Path(args.out))


if __name__ == "__main__":
    start_time = time.time()
    print_header("em_real_location.py")
    main()
    print_summary(start_time, "em_real_location.py")
