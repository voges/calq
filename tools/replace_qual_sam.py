#!/usr/bin/env python

import argparse
import os
import sys


def main():
    # set up the command line parser
    parser = argparse.ArgumentParser(
        description='replace quality scores in a SAM file'
    )
    parser.add_argument('sam_file_name', help='SAM file name')
    parser.add_argument('qual_file_name', help='QUAL file name')
    args = parser.parse_args()

    # retrieve the command line arguments
    sam_file_name = args.sam_file_name
    qual_file_name = args.qual_file_name

    # check if the user-supplied input files exists
    if not os.path.isfile(sam_file_name):
        sys.exit("error: this is not a file: {}".format(sam_file_name))
    if not os.path.isfile(qual_file_name):
        sys.exit("error: this is not a file: {}".format(qual_file_name))

    # check if the user-supplied SAM file can be suspected to be a SAM file
    if not sam_file_name.endswith(".sam"):
        sys.exit("error: SAM file name must end with '.sam'")

    # open the input files
    sam_file = open(sam_file_name, 'r')
    qual_file = open(qual_file_name, 'r')

    # write new SAM contents to stdout
    while 1:
        line = sam_file.readline().rstrip('\n')
        if not line:
            break
        if not line.startswith('@'):
            new_qual = qual_file.readline().rstrip('\n')
            fields = line.split('\t')
            new_line = (fields[0] + '\t' +
                        fields[1] + '\t' +
                        fields[2] + '\t' +
                        fields[3] + '\t' +
                        fields[4] + '\t' +
                        fields[5] + '\t' +
                        fields[6] + '\t' +
                        fields[7] + '\t' +
                        fields[8] + '\t' +
                        fields[9] + '\t' +
                        new_qual)
            num_fields = len(fields)
            if num_fields > 11:  # check if we have AUX fields at all
                for i in range(11, num_fields):
                    new_line += '\t' + fields[i]
                new_line += '\n'
            else:  # no AUX fields, so just write a line break
                new_line += '\n'
            sys.stdout.write(new_line)
        else:
            sys.stdout.write(line)

    # close the input files
    sam_file.close()
    qual_file.close()


if __name__ == "__main__":
    main()
