#!/usr/bin/env python3
"""
Script Name:  format_pileup.py
Description:  Process a pileup file into consensus, patched-reference, and nucleotide-frequency outputs.
Author:       Dr. Haiting Chai
Last update:  2026-APR-11
Version:      ver.2026-04-11

DEPENDENCIES:
    - Python 3.10
    - __future__ (standard library)
    - argparse (standard library)
    - collections (standard library)
    - logging (standard library)
    - pathlib (standard library)
    - re (standard library)
    - sys (standard library)
    - time (standard library)
    - assets.utils_banner (project module)
"""

from __future__ import annotations

# ========================================================================================
# Imports
# ========================================================================================

import argparse
import logging
import re
import sys
import time
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from assets.utils_banner import print_header, print_summary, setup_logging


# ========================================================================================
# Configuration
# ========================================================================================

logger = setup_logging()

IUPAC = {
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["A", "C"]): "M",
    frozenset(["G", "T"]): "K",
    frozenset(["A", "T"]): "W",
    frozenset(["C", "G"]): "S",
    frozenset(["C", "G", "T"]): "B",
    frozenset(["A", "G", "T"]): "D",
    frozenset(["A", "C", "T"]): "H",
    frozenset(["A", "C", "G"]): "V",
    frozenset(["A", "C", "G", "T"]): "N",
}


# ========================================================================================
# Helper Functions
# ========================================================================================

def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for pileup formatting."""
    parser = argparse.ArgumentParser(
        description="Generate pileup-derived consensus and summary outputs.",
    )
    parser.add_argument("--input", help="Input pileup file")
    return parser


def clean_reads(reads: str) -> str:
    """Remove pileup control markers so only observed base calls remain."""
    reads = re.sub(r"\^.", "", reads)
    reads = re.sub(r"\$", "", reads)
    reads = re.sub(r"[+-]\d+[ACGTNacgtn]+", "", reads)
    return reads


def base_counts(counts: Counter[str]) -> dict[str, int]:
    """Return counts for the canonical DNA bases in fixed column order."""
    return {base: counts.get(base, 0) for base in "ACGT"}


def consensus_from_top_counts(counts: Counter[str]) -> str:
    """Return the top base or IUPAC ambiguity code for a pileup position."""
    filtered = base_counts(counts)
    max_count = max(filtered.values())
    if max_count == 0:
        return "?"
    top_bases = [base for base, value in filtered.items() if value == max_count]
    if len(top_bases) == 1:
        return top_bases[0]
    return IUPAC.get(frozenset(top_bases), "N")


def get_case_consensus(counts: Counter[str]) -> str:
    """Encode confidence by letter case and low-depth masking symbols."""
    filtered = base_counts(counts)
    max_count = max(filtered.values())
    if max_count == 0:
        return "?"
    if max_count < 5:
        return "*"
    consensus = consensus_from_top_counts(counts)
    if 5 <= max_count <= 30:
        return consensus.lower()
    return consensus


def get_patch_consensus(counts: Counter[str], reference_base: str) -> str:
    """Fall back to the reference base when a pileup position has no calls."""
    filtered = base_counts(counts)
    if max(filtered.values()) == 0:
        return reference_base
    return consensus_from_top_counts(counts)


def get_variant_signal(counts: Counter[str]) -> str:
    """Summarize the set of observed alleles as a single IUPAC symbol."""
    observed = {base for base in "ACGT" if counts.get(base, 0) > 0}
    if not observed:
        return "?"
    if len(observed) == 1:
        return next(iter(observed))
    return IUPAC.get(frozenset(observed), "N")


def parse_pileup(input_path: Path) -> tuple[list[dict[str, str | int]], dict[str, str]]:
    """Build per-position summary rows and the FASTA strings derived from them."""
    rows: list[dict[str, str | int]] = []
    sequences = {
        "case": "",
        "called": "",
        "patched": "",
        "variant": "",
    }

    with input_path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 5:
                continue

            chrom, position, reference_base, depth, reads = fields[:5]
            cleaned_reads = clean_reads(reads)
            normalized_reads = cleaned_reads.upper().replace(".", reference_base.upper()).replace(",", reference_base.upper())
            counts = Counter(normalized_reads)

            # The workflow emits several FASTA variants from the same pileup pass, so
            # consensus strings are assembled here while position-level counts are fresh.
            case_consensus = get_case_consensus(counts)
            called_consensus = consensus_from_top_counts(counts)
            patched_consensus = get_patch_consensus(counts, reference_base.upper())
            variant_signal = get_variant_signal(counts)

            sequences["case"] += case_consensus
            sequences["called"] += called_consensus
            sequences["patched"] += patched_consensus
            sequences["variant"] += variant_signal

            rows.append(
                {
                    "Ref": chrom,
                    "Pos": int(position),
                    "A": counts.get("A", 0),
                    "C": counts.get("C", 0),
                    "G": counts.get("G", 0),
                    "T": counts.get("T", 0),
                    "N": counts.get("N", 0),
                    "Depth": int(depth),
                    "CasedConsensus": case_consensus,
                    "CalledConsensus": called_consensus,
                    "PatchedWithRef": patched_consensus,
                    "VariantSignal": variant_signal,
                }
            )

    return rows, sequences


def write_outputs(input_path: Path, rows: list[dict[str, str | int]], sequences: dict[str, str]) -> None:
    """Write the CSV and FASTA products derived from one pileup input."""
    sample_name = input_path.stem
    ntfreq_path = input_path.with_suffix(".ntfreq.csv")
    consensus_path = input_path.with_suffix(".consensus.fasta")
    bases_path = input_path.with_suffix(".bases.fasta")
    patched_path = input_path.with_suffix(".patchedref.fasta")
    variant_path = input_path.with_suffix(".variantsignals.fasta")

    with ntfreq_path.open("w", encoding="utf-8") as handle:
        handle.write("Ref,Pos,A,C,G,T,N,Depth,CasedConsensus,CalledConsensus,PatchedWithRef,VariantSignal\n")
        for row in rows:
            handle.write(
                "{Ref},{Pos},{A},{C},{G},{T},{N},{Depth},{CasedConsensus},{CalledConsensus},{PatchedWithRef},{VariantSignal}\n".format(
                    **row
                )
            )

    with consensus_path.open("w", encoding="utf-8") as handle:
        handle.write(f">{sample_name}_caseconsensus\n{sequences['case']}\n")
    with bases_path.open("w", encoding="utf-8") as handle:
        handle.write(f">{sample_name}_calledconsensus\n{sequences['called']}\n")
    with patched_path.open("w", encoding="utf-8") as handle:
        handle.write(f">{sample_name}_patchedref\n{sequences['patched']}\n")
    with variant_path.open("w", encoding="utf-8") as handle:
        handle.write(f">{sample_name}_variantsignals\n{sequences['variant']}\n")

    logger.info("NT frequency written to %s", ntfreq_path)
    logger.info("Cased consensus written to %s", consensus_path)
    logger.info("Called consensus written to %s", bases_path)
    logger.info("Patched reference written to %s", patched_path)
    logger.info("Variant signals written to %s", variant_path)


# ========================================================================================
# Main Entry Point
# ========================================================================================

def main() -> None:
    """Parse CLI arguments and export all pileup-derived outputs."""
    args = build_parser().parse_args()
    input_path = Path(args.input)
    rows, sequences = parse_pileup(input_path)
    write_outputs(input_path, rows, sequences)


if __name__ == "__main__":
    start_time = time.time()
    print_header("format_pileup.py")
    main()
    print_summary(start_time, "format_pileup.py")
