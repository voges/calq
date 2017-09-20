#!/usr/bin/env python

import sys

# Usage
if len(sys.argv) != 2:
    sys.exit("Usage: python {} file.sam 1>file.fastq".format(sys.argv[0]))

# Get SAM file
sam_file_name = sys.argv[1]
if not sam_file_name.endswith(".sam"):
    sys.exit("SAM file name must end with '.sam'")
sam_file = open(sam_file_name, 'r')

while 1:
    line = sam_file.readline()
    if not line:
        break
    if line[0] != '@':
        fields = line.split('\t')
        sys.stdout.write("@" + fields[0] + "\n")  # rname
        sys.stdout.write(fields[9] + "\n")  # seq
        sys.stdout.write("+" + "\n")
        sys.stdout.write(fields[10] + "\n")  # qual

sys.exit()
