"""
Core data structures and functions for building, reading, and merging
barcoding "cluster" files, which group together all BPM (bead oligo), DPM
(DNA proximity marker), and RPM (RNA proximity marker) reads that share the
same split-pool barcode.
"""

from collections import Counter
from dataclasses import dataclass, field
import re
import os
import sys
import pysam


@dataclass(frozen=True)
class Position:
    """A single genomic position for one read (DPM, RPM, or BPM).

    Note: `_feature` (e.g. strand) is deliberately excluded from equality/hashing
    (via `field(compare=False)`) so that two reads differing only in feature at the
    same (type, chromosome, start, end) are treated as duplicates. This matches the
    position-only deduplication used elsewhere in the pipeline (see
    threshold_tag_and_split.py's label_bam_file()).
    """

    _type: str
    _feature: str = field(compare=False)
    _chromosome: str
    # start/end are int when built from BAM records (get_clusters) but str when
    # parsed back out of a cluster file's text (merge_clusters) - not enforced at
    # runtime, just documenting the actual mixed usage.
    _start_coordinate: int | str
    _end_coordinate: int | str

    def to_string(self):
        """Serialize to the cluster-file format: "TYPE[feature]_chrom:start-end",
        e.g. "DPM[+]_chr1:1000-1050"."""
        try:
            out = (
                self._type
                + "["
                + self._feature
                + "]"
                + "_"
                + self._chromosome
                + ":"
                + str(self._start_coordinate)
                + "-"
                + str(self._end_coordinate)
            )
        except Exception:
            print(self._type, self._feature, self._chromosome)
            print("Elements are not as expect!")
            sys.exit()
        return out

    def score(self):
        """Order: BPM, RPM(chromosome), RPM(custom genome), DPM"""
        if self._type == "RPM":
            if self._chromosome.startswith("chr"):
                return 2
            else:
                return 3
        elif self._type == "BPM":
            return 1
        else:
            return 4


class Cluster:
    """A barcoding cluster: a deduplicated collection of genomic `Position`s that
    share the same split-pool barcode.

    The underlying data structure is a set, so duplicate positions (per
    `Position`'s equality rules) are discarded automatically.
    """

    def __init__(self):
        self._positions = set()

    def add_position(self, position):
        """Add a `Position` to this cluster (no-op if already present)."""
        self._positions.add(position)

    def to_string(self):
        """Serialize this cluster's positions to a tab-delimited string, sorted
        by `Position.score()` (BPM, then RPM, then DPM)."""
        positions_sorted = sorted(list(self._positions), key=lambda x: x.score())
        position_strings = [position.to_string() for position in positions_sorted]
        return "\t".join(position_strings)


class Clusters:
    """A collection of barcoding `Cluster`s, keyed by barcode string."""

    def __init__(self):
        self._clusters = {}

    def get_cluster(self, barcode):
        """Return the `Cluster` for `barcode`, creating an empty one first if
        it doesn't exist yet."""
        if barcode not in self._clusters:
            self._clusters[barcode] = Cluster()
        return self._clusters[barcode]

    def add_position(self, barcode, position):
        """Add `position` to the cluster identified by `barcode` (creating
        that cluster if needed)."""
        self.get_cluster(barcode).add_position(position)

    def to_strings(self):
        """Yield '<barcode>\\t<serialized positions>' for every cluster."""
        for barcode, cluster in self._clusters.items():
            yield barcode + "\t" + cluster.to_string()


##############################################################################################
# FUNCTIONS
##############################################################################################


def get_clusters(filelist, num_tags):
    """
    Generate a cluster file

    Args:
        filelist(list) = list of BAM files with barcoded reads
        num_tags(int) = number of tags in barcode
    """

    clusters = Clusters()
    pattern = re.compile("::" + num_tags * r"\[([a-zA-Z0-9_\-]+)\]")
    dpm_counts = 0
    bpm_counts = 0
    rpm_counts = 0
    for sample in filelist:
        file_name = os.path.basename(sample)
        sample_name = file_name.split(".")[0]
        try:
            with pysam.AlignmentFile(sample, "rb") as f:
                for read in f.fetch(until_eof=True):
                    name = read.query_name
                    match = pattern.search(name)
                    barcode = list(match.groups())
                    if "DPM" in name:
                        barcode_drop = barcode[1:]
                        dpm_counts += 1
                        strand = "+" if not read.is_reverse else "-"
                        position = Position(
                            "DPM",
                            strand,
                            read.reference_name,
                            read.reference_start,
                            read.reference_end,
                        )
                    elif "RPM" in name:
                        barcode_drop = barcode[1:]
                        rpm_counts +=1
                        strand = "+" if not read.is_reverse else "-"
                        position = Position(
                            "RPM",
                            strand,
                            read.reference_name,
                            read.reference_start,
                            read.reference_end
                        )
                    elif "BEAD" in name:
                        barcode_drop = barcode[1:]
                        bpm_counts += 1
                        UMI = read.reference_start
                        position = Position("BPM", 
                            "", 
                            read.reference_name, 
                            UMI, 
                            0
                        )
                    barcode_drop.append(sample_name)
                    barcode_str = ".".join(barcode_drop)
                    clusters.add_position(barcode_str, position)
        except ValueError:
            print("File provided has issues")
    print("Total BPM: ", bpm_counts)
    print("Total DPM: ", dpm_counts)
    print("Total RPM: ", rpm_counts)
    return clusters


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
    if count < min_oligos:
        return "uncertain"
    elif count / sum(bead_labels.values()) < threshold:
        return "ambiguous"
    else:
        return candidate


def write_clusters_to_file(clusters, outfile):
    """
    Writes a Clusters object to a file

    Args:
        clusters(obj): cluster object
        outfile(str): path to save clusterfile
    """

    count = 0
    with open(outfile, "w") as f:
        for cluster_string in clusters.to_strings():
            f.write(cluster_string)
            f.write("\n")
            count += 1
    print("Number of clusters written: ", count)


def write_single_cluster(barcode, reads, out):
    """
    Write a single cluster to a file in its string representation

    Arg:
        barcode(str): barocde
        reads(list): list of position objects representing reads
        out(file): open file to which cluster should be written
    """
    positions_sorted = sorted(list(reads), key=lambda x: x.score())
    position_strings = [position.to_string() for position in positions_sorted]
    out_string = "\t".join([barcode] + position_strings)
    out.write(out_string)
    out.write("\n")


def merge_clusters(in_file, out_file):
    """
    Merge clusters that contain the same barcode

    Notes:
        clusters are in alphabetical order by barcode
        reads are deduplicated during merging

    Args:
        in_file(str): filepath of cluster file, potentially with multiple clusters sharing the same barcodes
        out_file(str): filepath to write deduplicated, merged clusters to
    """
    current_barcode = ""
    current_reads = set()
    count = 0
    pattern = re.compile(r"([a-zA-Z0-9]+)\[(.*)\]_(.+):([0-9]+)\-([0-9]+)")
    with open(in_file, "r") as in_clusters, \
         open(out_file, "w") as out_clusters:
        for line in in_clusters:
            barcode, *reads = line.rstrip("\n").split("\t")
            if barcode != current_barcode:
                if current_barcode != "":
                    write_single_cluster(current_barcode, current_reads, out_clusters)
                    count += 1
                current_barcode = barcode
                current_reads = set()
            for read in reads:
                try:
                    match = pattern.search(read)
                    read_type, feature, chrom, start, end = match.groups()
                    position = Position(read_type, feature, chrom, start, end)
                    current_reads.add(position)
                except Exception:
                    print(read)
                    raise Exception("Pattern did not match above printed string")
        write_single_cluster(current_barcode, current_reads, out_clusters)
        count += 1
    print("Total clusters written: ", count)
