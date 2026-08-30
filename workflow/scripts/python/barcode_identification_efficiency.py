"""
Calculate the efficiency of barcode identification: how many reads were fully
assigned a barcode at every tag position, broken down by tag and by number of
tags found per read.
"""

import gzip
import argparse
import re
from collections import Counter
from helpers import fastq_parse
import pandas as pd


def parse_tags(config_file, format_file):
    """
    Build the tag lookup tables needed to interpret barcode fields in read names.

    Args
    - config_file: str
        Path to the barcode config file (used by splitcode). Must have "ids" and
        "groups" columns.
    - format_file: str
        Path to the splitcode barcode format file. Each line is a comma-separated
        list of tag names in read2 barcode order (e.g. "DPM,R5_ODD,R4_EVEN,...");
        the leading "DPM"/"BPM" token is stripped since read1 tags are handled
        separately.

    Returns
    - tag_map: dict (str -> str)
        Map from individual barcode ID (e.g. "R5_ODD_A1") to its tag group
        (e.g. "R5_ODD").
    - tag_index: dict (int -> str)
        Map from read2 barcode position (0-indexed) to tag group name.
    """
    df = pd.read_csv(config_file, comment="#", sep="\t", header=0)
    tag_map = df.set_index('ids')['groups'].to_dict()
    tag_index = {}
    with open(format_file, 'r') as f:
        for line in f:
            #line=DPM,R5_ODD,R4_EVEN,R3_ODD,R2_EVEN,R1_ODD
            #R2 tags only
            line = line.replace('DPM,','').replace('BPM,','')
            tags = line.strip().split(',')
            for index, t in enumerate(tags):
                tag_index[index] = t
    return tag_map, tag_index


def parse_reads(filelist, tag_map):
    """
    Count how many barcode tags were identified per read, both per-tag-group and
    per-read (i.e. how many tags total were found on each read).

    Args
    - filelist: list of str
        Paths to gzip-compressed FASTQ files to scan.
    - tag_map: dict (str -> str)
        Map from barcode ID to tag group, as returned by parse_tags().

    Returns
    - num_reads: int
        Total number of reads scanned across all files.
    - individual_tag_counts: Counter (str -> int)
        Number of reads containing each tag group.
    - total_counts: Counter (int -> int)
        Number of reads that had exactly N tags identified, keyed by N.
    """
    total_counts = Counter()
    individual_tag_counts = Counter()
    num_reads = 0

    pattern = re.compile(r"\[([a-zA-Z0-9_\-]+)\]")
    for read_path in filelist:
        with gzip.open(read_path, 'rt') as f:
            for qname, seq, thrd, qual in fastq_parse(f):
                barcodes = pattern.findall(qname.split('::')[1])
                groups = [tag_map[b] for b in barcodes] 
                individual_tag_counts.update(groups)
                total_counts[len(groups)] +=1 
                num_reads += 1

    return num_reads, individual_tag_counts, total_counts

def percentage(numerator, denominator):
    """Format numerator/denominator as a percentage string with 1 decimal place;
    returns "0.0" instead of raising if denominator is 0."""
    return "0.0" if denominator == 0 else f"{(100 * numerator / denominator):.1f}"

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description = 'Calculate barcode identification efficiency'
    )
    parser.add_argument(
        "--assigned",
        metavar="FILE",
        action="store",
        required=True,
        help = "The fully assigned fastq file",
    )
    parser.add_argument(
        "--unassigned",
        metavar="FILE",
        action="store",
        required=True,
        help = "The unassigned fastq file",
    )
    parser.add_argument(
        "--config",
        metavar="FILE",
        action="store",
        required=True,
        help = "The config file used by splitcode",
    )
    parser.add_argument(
        "--format",
        metavar="FILE",
        action="store",
        required=True,
        help = "The format file used by splitcode",
    )
    return parser.parse_args()



def main():
    """Compute and print barcode identification efficiency statistics for one
    sample from its assigned/unassigned FASTQ files."""
    args = parse_arguments()


    tag_map, tag_index = parse_tags(args.config, args.format)
    num_tags = max(list(tag_index.keys()))+1 #split-pool tags, plus RPM

    num_full, tag_counts_full, total_counts_full = parse_reads(args.assigned.split(' '), tag_map)
    num_unassigned, tag_counts_unassigned, total_counts_unassigned = parse_reads(args.unassigned.split(' '), tag_map)
    total_reads = num_full + num_unassigned
    tag_counts_all = tag_counts_full + tag_counts_unassigned
    total_counts_all = total_counts_full + total_counts_unassigned

    pct_full = percentage(num_full, total_reads)
    print()

    #Print tags for each group
    print("Read 1 tags")
    for r1 in ['DPM', 'BPM']:
        count = tag_counts_all.get(r1,0)
        pct = percentage(count, total_reads)
        print(f"{count} ({pct}%) barcodes found with tag {r1} (read1, position 1).")
    print("Read 2 tags")
    for i in range(num_tags):
        group = tag_index[i]
        count = tag_counts_all.get(group, 0)
        pct = percentage(count, total_reads)
        print(f"{count} ({pct}%) barcodes found with tag {group} (read2, position {i}).")

    #Print reads with total number of tags
    print()
    for i in range(num_tags + 1):
        count = total_counts_all.get(i+1, 0)
        pct = percentage(count, total_reads)
        print(f"{count} ({pct}%) reads found with {i+1} barcodes.")

    print()
    print(f"{num_full} ({pct_full}%) reads found with all barcodes and matching pattern.\n")


if __name__ == "__main__":
    main()

