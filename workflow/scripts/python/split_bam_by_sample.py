#!/usr/bin/env python3
"""
split_bam_by_sample.py

Split a BAM file into multiple BAM files, one per sample, based on a sample
identifier embedded in each read name.

Read names are expected to look like:
    NYStgBot_1-A1.OddBot_70-F10.EvenBot_46-D10.OddBot_33-C9.EvenBot_11-A11.OddBot_17-B5.DPM_ACTD

The sample name can appear ANYWHERE within the read name, as a plain
substring with no required delimiters. For example, given the sample name
"ACTD", the read name "NY_1-A1.ODD_70-ACTD.EVEN_46-C9" also belongs to
sample "ACTD" (the substring "ACTD" appears inside the field "70-ACTD").
The script searches the full read name for each known sample name as a
substring.

Reads whose read name does not contain exactly one of the supplied sample
names as a substring (no match, or more than one sample name found as a
substring) are written to a separate "_unassigned.bam" file.

Usage:
    python split_bam_by_sample.py --bam input.bam --samples ACTD FVP GFP \
        --outdir /path/to/outdir

Output file naming:
    <input_bam_basename>_<sample>.bam

For example, given input.bam and sample "ACTD", the output file would be:
    input_ACTD.bam

Requires: pysam (pip install pysam --break-system-packages)
"""

import argparse
import os
import sys
from dataclasses import dataclass, field

import pysam


# ---------------------------------------------------------------------------
# Sample matching
# ---------------------------------------------------------------------------

def get_sample_from_read_name(read_name, sample_names):
    """
    Find which sample name(s) occur as a substring anywhere within the read
    name. No delimiters are required around the sample name - it can be
    embedded directly inside a larger field (e.g. "70-ACTD" contains the
    sample name "ACTD").

    Returns a tuple (sample, ambiguous):
        - (sample_name, False) if exactly one known sample name was found
          as a substring of the read name.
        - (None, False) if no sample name was found.
        - (None, True) if more than one distinct sample name was found
          as a substring (ambiguous).
    """
    matched_samples = {sample for sample in sample_names if sample in read_name}

    if len(matched_samples) == 1:
        return next(iter(matched_samples)), False
    elif len(matched_samples) == 0:
        return None, False
    else:
        return None, True


# ---------------------------------------------------------------------------
# CLI / setup
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Split a BAM file into per-sample BAM files based on "
        "sample names embedded in read names."
    )
    parser.add_argument("--bam", required=True, help="Input BAM file path.")
    parser.add_argument(
        "--samples",
        required=True,
        nargs="+",
        help="List of sample names to split by (space-separated).",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Directory to write output BAM files to (default: current directory).",
    )
    parser.add_argument(
        "--write-unassigned",
        action="store_true",
        help="Also write reads that don't match any sample to an "
        "'_unassigned.bam' file (default: such reads are skipped).",
    )
    return parser.parse_args()


def get_bam_stem(bam_path):
    """Return the input BAM's basename with any '.bam' extension stripped."""
    basename = os.path.basename(bam_path)
    if basename.lower().endswith(".bam"):
        return basename[:-4]
    return basename


def build_output_path(outdir, bam_stem, suffix):
    """Build an output BAM path of the form <bam_stem>_<suffix>.bam."""
    return os.path.join(outdir, f"{bam_stem}_{suffix}.bam")


# ---------------------------------------------------------------------------
# BAM I/O
# ---------------------------------------------------------------------------

def open_output_bams(template_bam, outdir, bam_stem, samples, write_unassigned):
    """
    Open one output BAM file per sample (plus an optional unassigned file),
    using template_bam for the header.

    Returns (out_bams, unassigned_bam), where out_bams is a dict mapping
    sample name -> open AlignmentFile, and unassigned_bam is an open
    AlignmentFile or None if write_unassigned is False.
    """
    out_bams = {
        sample: pysam.AlignmentFile(
            build_output_path(outdir, bam_stem, sample), "wb", template=template_bam
        )
        for sample in samples
    }

    unassigned_bam = None
    if write_unassigned:
        unassigned_path = build_output_path(outdir, bam_stem, "unassigned")
        unassigned_bam = pysam.AlignmentFile(unassigned_path, "wb", template=template_bam)

    return out_bams, unassigned_bam


def close_bams(in_bam, out_bams, unassigned_bam):
    in_bam.close()
    for bam in out_bams.values():
        bam.close()
    if unassigned_bam is not None:
        unassigned_bam.close()


def index_bam(bam_path):
    """Index a BAM file, printing a warning to stderr on failure."""
    try:
        pysam.index(bam_path)
    except Exception as e:
        print(f"Warning: could not index {bam_path}: {e}", file=sys.stderr)


def index_output_bams(outdir, bam_stem, samples, write_unassigned):
    for sample in samples:
        index_bam(build_output_path(outdir, bam_stem, sample))
    if write_unassigned:
        index_bam(build_output_path(outdir, bam_stem, "unassigned"))


# ---------------------------------------------------------------------------
# Read processing
# ---------------------------------------------------------------------------

@dataclass
class SplitStats:
    """Tracks per-sample and overall read counts for the run."""
    counts: dict = field(default_factory=dict)
    unassigned_count: int = 0
    ambiguous_count: int = 0
    total: int = 0


def process_reads(in_bam, out_bams, unassigned_bam, samples):
    """
    Iterate over every read in in_bam, route it to the correct output BAM
    based on its read name, and tally statistics as we go.
    """
    stats = SplitStats(counts={sample: 0 for sample in samples})

    for read in in_bam:
        stats.total += 1
        sample, ambiguous = get_sample_from_read_name(read.query_name, samples)

        if sample is not None:
            out_bams[sample].write(read)
            stats.counts[sample] += 1
        else:
            stats.unassigned_count += 1
            if ambiguous:
                stats.ambiguous_count += 1
            if unassigned_bam is not None:
                unassigned_bam.write(read)

    return stats


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def print_summary(stats, write_unassigned):
    print(f"Total reads processed: {stats.total}")
    for sample, count in stats.counts.items():
        print(f"  {sample}: {count} reads")

    unassigned_note = "" if write_unassigned else " (not written, use --write-unassigned to save)"
    print(f"  unassigned: {stats.unassigned_count} reads{unassigned_note}")

    if stats.ambiguous_count:
        print(
            f"    of which {stats.ambiguous_count} were ambiguous "
            "(matched more than one sample name)"
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    if not os.path.isfile(args.bam):
        sys.exit(f"Error: input BAM file not found: {args.bam}")

    os.makedirs(args.outdir, exist_ok=True)
    bam_stem = get_bam_stem(args.bam)

    in_bam = pysam.AlignmentFile(args.bam, "rb")
    out_bams, unassigned_bam = open_output_bams(
        in_bam, args.outdir, bam_stem, args.samples, args.write_unassigned
    )

    stats = process_reads(in_bam, out_bams, unassigned_bam, args.samples)

    close_bams(in_bam, out_bams, unassigned_bam)
    index_output_bams(args.outdir, bam_stem, args.samples, args.write_unassigned)

    print_summary(stats, args.write_unassigned)


if __name__ == "__main__":
    main()
