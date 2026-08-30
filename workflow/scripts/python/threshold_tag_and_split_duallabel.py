import argparse
from collections import defaultdict, Counter
from itertools import combinations
import os
from pathlib import Path
import re
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import cluster as c
import pysam
from tqdm import tqdm

"""
Assign reads to antibodies and generate individual BAM files for each antibody.
Reads are deduplicated as individual files are generated: only one read per position (start and end) per cluster is
kept.
"""


def parse_args():
    parser = argparse.ArgumentParser(
        description="Assign reads to antibodies and generate individual BAM files for each antibody"
    )
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bam", help="Output labeled BAM file with antibody assignment added via read group (RG) tag")
    parser.add_argument("clusters", help="Clusters file to assign antibody labels")
    parser.add_argument("output_dir", help="Directory to write individual antibody BAMs to")
    parser.add_argument(
        "--num_tags",
        type=int,
        required=True,
        help="Number of tags in barcode",
    )
    parser.add_argument(
        "--min_oligos",
        type=int,
        required=True,
        help="The minimum number of oligos needed to call a cluster",
    )
    parser.add_argument(
        "--proportion",
        type=float,
        required=True,
        help="The minimum representation proportion of oligos needed to call a cluster",
    )
    parser.add_argument(
        "--max_size",
        type=int,
        required=True,
        help="The maximum cluster size to keep",
    )
    parser.add_argument(
        "--progress",
        action="store_true",
        help="Display a progress bar",
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    print("Using min_oligos: ", args.min_oligos, file=sys.stderr)
    print("Using proportion: ", args.proportion, file=sys.stderr)
    print("Using max_size: ", args.max_size, file=sys.stderr)
    print("Writing labeled BAM to: ", args.output_bam, file=sys.stderr)
    print("Writing splitbams to: ", args.output_dir, file=sys.stderr)

    labels = assign_labels(args.clusters, args.min_oligos, args.proportion, args.max_size, progress=args.progress)
    label_bam_file(args.input_bam, args.output_bam, labels, args.num_tags)
    split_bam_by_RG(args.output_bam, args.output_dir)


def label_cluster_reads(reads, min_oligos, threshold, max_size):
    """
    Assign a label to a cluster based on oligo reads

    Args:
        reads(list): list of cluster formated reads
        threshold(float): fraction of oligo reads that are of one type needed to assign the cluster
        min_oligos(int): minimum number of oligos of one type in a cluster in order to assign the cluster (i.e., greater than or equal to)
        max_size(int): maximum DNA (DPM) reads allowed per cluster
    """
    bead_reads = [read for read in reads if read.startswith("BPM")]
    if len(bead_reads) == 0:
        return "none"
    cluster_size = len(reads) - len(bead_reads)
    if int(cluster_size) > int(max_size):
        return "filtered"
    bead_labels = Counter([read.split(":")[0].split("_", 1)[1] for read in bead_reads])
    candidate, count = bead_labels.most_common()[0]
    total = sum(bead_labels.values())
    if count < min_oligos:
        dual_label = label_dual_cluster_reads(bead_labels, min_oligos, threshold, total)
        return dual_label if dual_label is not None else "uncertain"
    elif count / total < threshold:
        dual_label = label_dual_cluster_reads(bead_labels, min_oligos, threshold, total)
        return dual_label if dual_label is not None else "ambiguous"
    else:
        return candidate


def label_dual_cluster_reads(bead_labels, min_oligos, threshold, total):
    """
    Check whether a cluster can be assigned to a pair of antibody labels,
    used as a fallback when no single label meets the min_oligos/threshold
    criteria on its own.

    Every pair of labels is checked against three criteria:
    - Both labels individually have a count >= min_oligos
    - The two labels' counts are within 1 of each other
    - Their combined count meets the proportion threshold (relative to the
      total number of bead reads in the cluster)

    It's fine if more than two labels each separately satisfy one or two of
    these criteria (e.g. several labels individually pass min_oligos, or
    several different pairs happen to have counts within 1 of each other).
    What matters is how many labels pass all three criteria *together* as a
    pair: if exactly one pair of labels passes all three criteria jointly
    (i.e. exactly two labels qualify in total), that pair is returned;
    otherwise (no qualifying pair, or labels spread across more than one
    qualifying pair) None is returned.

    Args
    - bead_labels: Counter(str -> int)
        Counts of bead reads per antibody label in the cluster
    - min_oligos: int
        The minimum number of oligos needed to call a cluster
    - threshold: float
        The minimum representation proportion of oligos needed to call a cluster
    - total: int
        Total number of bead reads in the cluster

    Returns
    - str or None
        Hyphenated combination of the two qualifying labels in alphabetical
        order, or None if no such unique pair of labels exists
    """
    qualifying_pairs = []
    for (label1, count1), (label2, count2) in combinations(bead_labels.items(), 2):
        if count1 < min_oligos or count2 < min_oligos:
            continue
        if (count1 + count2) / total < threshold:
            continue
        qualifying_pairs.append((label1, label2))
    qualifying_labels = {label for pair in qualifying_pairs for label in pair}
    if len(qualifying_labels) != 2:
        return None
    ordered = sorted(qualifying_labels)
    return f"{ordered[0]}-{ordered[1]}"


def assign_labels(path_clusters, min_oligos, proportion, max_size, progress=False):
    """
    Assign labels to all clusters based on oligo reads

    Args
    - path_clusters: str
        Path to clusters file
    - min_oligos: int
        The minimum number of oligos needed to call a cluster
    - proportion: float
        The minimum representation proportion of oligos needed to call a cluster
    - max_size: int
        The maximum cluster size to keep
    - progress: bool
        Whether to display a progress bar

    Returns
    - labels: dict(str -> str)
        Map from barcode string to antibody label
    """
    labels = {}
    with open(path_clusters, "r") as clusters:
        for line in tqdm(clusters, disable=not progress):
            barcode, *reads = line.rstrip("\n").split("\t")
            label = label_cluster_reads(reads, min_oligos, proportion, max_size)
            labels[barcode] = label
    return labels


def label_bam_file(input_bam, output_bam, labels, num_tags, progress=False):
    """
    Add antibody label to individual reads in a BAM file.
    - Antibody labels are added as read group (RG) tags
    - Reads are deduplicated based on cluster assignment, start position, and end position (last aligned residue)

    Args
    - input_bam: str
        Path to input BAM file
    - output_bam: str or file object
        Path to write labeled BAM file
    - labels: dict(str -> str)
        Map from barcode string to antibody label
    - num_tags: int
        Number of tags in barcode
    - progress: bool
        Whether to display a progress bar
    """
    count, duplicates, skipped = 0, 0, 0

    written = defaultdict(int) #Counter of reads per bead type written
    pattern = re.compile("::" + num_tags * r"\[([a-zA-Z0-9_\-]+)\]")
    found = defaultdict(set) #Tracker for deduplicating reads
    
    library = os.path.basename(input_bam).split(".")[0]
    header = construct_read_group_header(input_bam, labels)
    with pysam.AlignmentFile(input_bam, "rb") as in_bam, \
         pysam.AlignmentFile(output_bam, "wb", header=header) as out_bam:
        for read in tqdm(in_bam.fetch(until_eof=True), disable=not progress):
            count += 1
            name = read.query_name
            match = pattern.search(name)
            barcode = list(match.groups())
            if "DPM" in name or "BEAD" in name:
                barcode_drop = barcode[1:]
            elif "RPM" in name:
                barcode_drop = barcode[:-1]
            barcode_drop.append(library)
            full_barcode = ".".join(barcode)
            position = read.reference_name + ":" + str(read.reference_start) + '-' + str(read.reference_end)
            if position in found[full_barcode]:
                duplicates += 1
            else:
                try:
                    readlabel = labels[full_barcode]
                    found[full_barcode].add(position)
                    read.set_tag("RG", readlabel, replace=True)
                    out_bam.write(read)
                    written[readlabel] += 1
                except KeyError:
                    skipped += 1

    print("Total reads:", count)
    print("Reads written:", written)
    print("Duplicate reads:", duplicates)
    print("Reads with an error not written out:", skipped)


def construct_read_group_header(input_bam, labels):
    """
    Add read group tags for each antibody to the BAM file header

    Args
    - input_bam: str
        Path to input BAM file
    - labels: dict(str -> str)
        Map from barcode string to antibody label

    Returns
    - bam_header: dict
        Updated BAM header
    """
    proteins = set(labels.values())
    sample_name = input_bam.split(".", 1)[0]
    with pysam.AlignmentFile(input_bam, "rb") as input_file:
        bam_header = input_file.header.to_dict()
        read_group_dict = [{"ID": name, "SM": sample_name} for name in list(proteins)]
        bam_header["RG"] = read_group_dict
        return bam_header


def split_bam_by_RG(output_bam, output_dir):
    """
    Split a BAM file based on the read group

    Args
    - output_bam: str
        Path to BAM file containing reads with read group assignments
    - output_dir: str
        Directory to write read group specific BAM files to
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    format_string = output_dir + "/%*_%!.%."
    pysam.split("-f", format_string, output_bam)


if __name__ == "__main__":
    main()
