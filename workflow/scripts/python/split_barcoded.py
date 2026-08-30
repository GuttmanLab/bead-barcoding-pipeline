"""
Program to split barcoded reads into two files based on adaptor (DPM, BPM) tags
"""

import argparse
import gzip
import os
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from helpers import fastq_parse, file_open

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Split fastq based on DPM, RPM or BPM barcode"
    )
    parser.add_argument(
        "--input", dest="input", type=str, required=True, help="Path to input fastq"
    )
    parser.add_argument(
        "--type", type=str, help="Experiment type (DNA, RNA)"
    )
    opts = parser.parse_args()

    return opts

def main():
    """Split the input R1 FASTQ into per-adaptor-type output files: DPM/BPM for
    DNA experiments, or RPM (R1+R2)/BPM for RNA experiments."""
    opts = parse_args()
    base_path = os.path.splitext(os.path.splitext(opts.input)[0])[0]
    if opts.type == 'DNA':
        split_dpm_reads(opts.input, base_path)
    elif opts.type == 'RNA':
        split_rpm_reads(opts.input, base_path)

def split_dpm_reads(read_1_path, base_path):
    """
    Split a barcoded R1 FASTQ into separate DPM and BPM (bead oligo) output
    FASTQ files, based on which tag appears in each read's name. Raises if a
    read has neither tag.

    Args
    - read_1_path: str
        Path to the input barcoded R1 FASTQ (gzip-compressed or plain text).
    - base_path: str
        Output path prefix; writes "<base_path>_dpm.fastq.gz" and
        "<base_path>_bpm.fastq.gz".
    """
    dpm_out_path = base_path + "_dpm.fastq.gz"
    bpm_out_path = base_path + "_bpm.fastq.gz"
    dpm_count = 0
    bpm_count = 0
    counter = 0
    with file_open(read_1_path, mode="rt") as read_1, \
         gzip.open(dpm_out_path, "wt") as dpm_out, \
         gzip.open(bpm_out_path, "wt") as bpm_out:
        for qname, seq, thrd, qual in fastq_parse(read_1):
            counter += 1
            if counter % 10000 == 0:
                print(counter)
            if "DPM" in qname:
                dpm_count += 1
                dpm_out.write(qname + "\n" + seq + "\n" + thrd + "\n" + qual + "\n")
            elif "BEAD" in qname:
                bpm_count += 1
                bpm_out.write(qname + "\n" + seq + "\n" + thrd + "\n" + qual + "\n")
            else:
                print("Read with unexpected barcode", file=sys.stderr)
                print(qname + "\n" + seq + "\n" + thrd + "\n" + qual + "\n", file=sys.stderr)
                raise Exception
    print("DPM reads out:", dpm_count)
    print("BPM reads out:", bpm_count)

def split_rpm_reads(read_1_path, base_path):
    """
    Split barcoded R1/R2 FASTQs into RPM and BPM (bead oligo) output FASTQ
    files. R1 reads are split into RPM or BPM by tag; R2 reads are written to
    the RPM output only if their mate (same read) was tagged RPM in R1.

    For RPM reads, also relocates the "[RPM]" tag from wherever splitcode put
    it to immediately after "::", matching where "[DPM]"/"[BEAD_...]" tags sit
    for DPM/BPM reads. This keeps the tag position consistent across read
    types for downstream parsing in cluster.py and threshold_tag_and_split.py,
    both of which assume the read-type tag is the first bracketed group after
    "::".

    Args
    - read_1_path: str
        Path to the input barcoded R1 FASTQ. The matching R2 path is derived
        by replacing "_R1." with "_R2." in this path.
    - base_path: str
        Output path prefix; writes "<base_path>_rpm.fastq.gz" (R1),
        "<base_path>_rpm.fastq.gz" with "_R1."->"_R2." (R2), and
        "<base_path>_bpm.fastq.gz" (R1 BPM reads only).
    """
    rpm_out_path = base_path + "_rpm.fastq.gz"
    bpm_out_path = base_path + "_bpm.fastq.gz"
    # Paths for R2
    read_2_path = read_1_path.replace('_R1.', '_R2.')
    rpm_out_path_r2 = rpm_out_path.replace('_R1.', '_R2.')
    rpm_count = 0
    rpm_count2 = 0
    bpm_count = 0
    counter = 0
    counter2 = 0
    with file_open(read_1_path, mode="rt") as read_1, \
         gzip.open(rpm_out_path, "wt") as rpm_out, \
         gzip.open(bpm_out_path, "wt") as bpm_out: 
        for qname, seq, thrd, qual in fastq_parse(read_1):
            counter += 1
            if counter % 10000 == 0:
                print(counter)
            if "RPM" in qname:
                rpm_count += 1
                qname = qname.replace("[RPM]", "").replace("::", "::[RPM]")
                rpm_out.write(qname + "\n" + seq + "\n" + thrd + "\n" + qual + "\n")
            elif "BEAD" in qname:
                bpm_count += 1
                bpm_out.write(qname + "\n" + seq + "\n" + thrd + "\n" + qual + "\n")
            else:
                print("Read with unexpected barcode", file=sys.stderr)
                print(qname + "\n" + seq + "\n" + thrd + "\n" + qual + "\n", file=sys.stderr)
                raise Exception

    with file_open(read_2_path, mode="rt") as read_2, \
            gzip.open(rpm_out_path_r2, "wt") as rpm_out2:
        for qname, seq, thrd, qual in fastq_parse(read_2):
            counter2 += 1
            if counter2 % 10000 == 0:
                print(counter2)
            if "RPM" in qname:
                rpm_count2 += 1
                qname = qname.replace("[RPM]", "").replace("::", "::[RPM]")
                rpm_out2.write(qname + "\n" + seq + "\n" + thrd + "\n" + qual + "\n")
    print("RPM reads out R1:", rpm_count)
    print("RPM reads out R2:", rpm_count2)
    print("BPM reads out:", bpm_count)

if __name__ == "__main__":
    main()
