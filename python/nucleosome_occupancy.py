# Write some python boiler plate code
# Import the necessary modules
import numpy as np
import argparse
import csv
from collections import defaultdict
import matplotlib.pyplot as plt

def plot_coverage(coverage):
    """
    Plot the coverage of the nucleosome occupancy in a line plot
    """
    plt.plot(coverage)
    plt.savefig('nucleosome_occupancy.png')


def aggregate_coverage(coverage, tfbs, window):
    """
    aggregate the coverage vectors across every window in the tfbs file

    :param coverage: a dictonary of chromosome names and coverage vectors
    :paramtfbs: path to the tfbs file
    :param window: the size of the window to aggregate over
    :return: aggregated coverage vector
    """
    aggregate_coverage = np.zeros(2*window+1)
    with open(tfbs, 'r') as tfbs:
        reader = csv.reader(tfbs, delimiter='\t')
        for chr, _, _, _, midpoint in reader:
            start = midpoint - window
            end = midpoint + window
            aggregate_coverage += coverage[chr][start:end]

    return aggregate_coverage

def get_coverage(bed_file):
    """
    Get the coverage at every locus in the interval [start, end]

    :param bed_file: The bed file to get the coverage from
    :return: A list of coverage values
    """
    coverage = {'chr{}'.format(i): np.zeros(249*10**6) for i in range(1, 23)}
    with open(bed_file, 'r') as bed:
        reader = csv.reader(bed, delimiter='\t')
        for chr, start, end in reader:
            for i in range(start, end+1):
                coverage[chr][i] += 1

    return coverage



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input bed file")
    parser.add_argument('-t', '--tfbs', help='TFBS file')
    parser.add_argument('-w', '--window', help='Window size to search around TFBS', type=int)
    args = parser.parse_args()

    coverage = get_coverage(args.input)
    aggregate_coverage = aggregate_coverage(coverage, args.tfbs, args.window)
    plot_coverage(aggregate_coverage)

if __name__ == '__main__':
    main()
