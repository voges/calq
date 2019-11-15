#!/usr/bin/env python

import argparse
import os
import sys


def main():
    # set up the command line parser
    parser = argparse.ArgumentParser(
        description='analyze quality scores in a FASTQ file'
    )
    parser.add_argument('fastq_file_name', help='FASTQ file name')
    args = parser.parse_args()

    # retrieve the command line arguments
    fastq_file_name = args.fastq_file_name

    # check if the user-supplied input file exists
    if not os.path.isfile(fastq_file_name):
        sys.exit("error: this is not a file: {}".format(fastq_file_name))

    # check if the user-supplied FASTQ file can be suspected to be a FASTQ file
    if (not fastq_file_name.endswith(".fastq")
            and not fastq_file_name.endswith(".fq")):
        sys.exit("error: input file name must end with '.fastq' or '.fq'")

    # open the input file
    fastq_file = open(fastq_file_name, 'r')

    # initialize statistics
    qual_min = sys.maxsize
    qual_max = -sys.maxsize - 1
    qual_size = 0
    line_cnt = 0
    record_cnt = 0
    qual_dist = []
    for x in range(0, 128):
        qual_dist.append(0)

    # parse FASTQ file
    while 1:
        sequence_identifier = fastq_file.readline().rstrip('\n')
        if not sequence_identifier:
            break  # reached end of file, everything ok
        fastq_file.readline().rstrip('\n')  # sequence
        fastq_file.readline().rstrip('\n')  # description
        quality_scores = fastq_file.readline().rstrip('\n')

        line_cnt += 4
        record_cnt += 1

        if len(quality_scores) > 0:
            qual_size += len(quality_scores)
            for q in quality_scores:
                if ord(q) > qual_max:
                    qual_max = ord(q)
                if ord(q) < qual_min:
                    qual_min = ord(q)
                qual_dist[ord(q)] += 1
        else:
            sys.exit("error: no quality score in record {} (line {})"
                     .format(record_cnt - 1, line_cnt - 1))

    print("FASTQ file: {}".format(fastq_file_name))
    print("  lines: {}".format(line_cnt))
    print("  records: {}".format(record_cnt))
    print("  quality scores: {}".format(qual_size))
    print("    range (inclusive): [{},{}]".format(qual_min, qual_max))
    print("    distribution:")
    for x in range(0, 128):
        if qual_dist[x] != 0:
            print("      {}: {}".format(x, qual_dist[x]))

    # close the input file
    fastq_file.close()


if __name__ == "__main__":
    main()
