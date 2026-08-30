"""
Count 1) number of clusters, 2) number of DPM reads (aligned), 3) number of RPM reads (aligned)  and 4) number of  BPM reads within each clusterfile for a directory of clusterfiles.
"""

import argparse
import glob
import tqdm

def main():
    """Print cluster/DPM/BPM/RPM statistics for every clusterfile matching
    <directory>/*<pattern>."""
    args = parse_arguments()
    search = args.directory + "/*" + args.pattern
    files = glob.glob(search)
    for f in files:
        count_statistics(f)

def count_statistics(clusterfile):
    """
    Loop through all clusters within a clusterfile, counting DPM, RPM  and BPM reads

    Args:
        clusterfile(str): Path to clusterfile
    """
    cluster = 0
    dpm = 0
    bpm = 0
    rpm = 0
    with open(clusterfile, "r") as clusters:
        for line in tqdm.tqdm(clusters):
            cluster += 1
            dpm += line.count("DPM[")
            bpm += line.count("BPM[")
            rpm += line.count("RPM[")
    print("For clusterfile ", clusterfile)
    print("Total number of clusters: ", cluster)
    print("Total number of BPM: ", bpm)
    print("Total number of DPM: ", dpm)
    print("Total number of RPM: ", rpm)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate the statistics for all clusterfiles in a directory"
    )
    parser.add_argument(
        "--directory",
        metavar="FILE",
        action="store",
        help="The directory of clusters file",
    )
    parser.add_argument(
        "--pattern", action="store", help="The pattern of cluster file names"
    )
    return parser.parse_args()

if __name__ == "__main__":
    main()
