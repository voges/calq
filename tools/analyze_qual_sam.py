#!/usr/bin/env python

import argparse
import os
import sys


def main():
    # set up the command line parser
    parser = argparse.ArgumentParser(
        description='analyze quality scores in a SAM file'
    )
    parser.add_argument('sam_file_name', help='SAM file name')
    args = parser.parse_args()

    # retrieve the command line arguments
    sam_file_name = args.sam_file_name

    # check if the user-supplied input file exists
    if not os.path.isfile(sam_file_name):
        sys.exit("error: this is not a file: {}".format(sam_file_name))

    # check if the user-supplied SAM file can be suspected to be a SAM file
    if not sam_file_name.endswith(".sam"):
        sys.exit("error: SAM file name must end with '.sam'")

    # open the input file
    sam_file = open(sam_file_name, 'r')

    # initialize statistics
    qual_min = sys.maxsize
    qual_max = -sys.maxsize - 1
    qual_size = 0
    total_line_cnt = 0
    header_line_cnt = 0
    alignment_line_cnt = 0
    qual_dist = []
    for x in range(0, 128):
        qual_dist.append(0)

    # parse SAM file
    while 1:
        line = sam_file.readline()
        if not line:
            break
        if not line.startswith('@'):
            fields = line.split('\t')
            qual = fields[10]
            if len(qual) > 0:
                qual_size += len(qual)
                for q in qual:
                    if ord(q) > qual_max:
                        qual_max = ord(q)
                    if ord(q) < qual_min:
                        qual_min = ord(q)
                    qual_dist[ord(q)] += 1
            else:
                sys.exit("error: no quality scores in line {}"
                         .format(total_line_cnt))
            alignment_line_cnt += 1
        else:
            header_line_cnt += 1
        total_line_cnt += 1

    # print statistics
    print("SAM file: {}".format(sam_file_name))
    print("  lines: {}".format(total_line_cnt))
    print("    header lines: {}".format(header_line_cnt))
    print("    alignment lines: {}".format(alignment_line_cnt))
    print("  quality scores: {}".format(qual_size))
    print("    range (inclusive): [{},{}]".format(qual_min, qual_max))
    print("    distribution:")
    for x in range(0, 128):
        if qual_dist[x] != 0:
            print("      {}: {}".format(x, qual_dist[x]))

    # close the input file
    sam_file.close()


if __name__ == "__main__":
    main()
