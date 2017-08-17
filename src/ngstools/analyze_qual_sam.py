#!/usr/bin/env python

import sys

# Usage
if len(sys.argv) != 2:
    sys.exit("Usage: python {} file.sam".format(sys.argv[0]))

# Get SAM file
sam_file_name = sys.argv[1]
if not sam_file_name.endswith(".sam"):
    sys.exit("SAM file name must end with '.sam'")
sam_file = open(sam_file_name, 'r')
print("SAM file: {}".format(sam_file_name))

# Initialize statistics
qual_min = sys.maxint
qual_max = -sys.maxint - 1
qual_size = 0
total_line_cnt = 0
header_line_cnt = 0
alignment_line_cnt = 0

# Parse SAM file
while 1:
    line = sam_file.readline()
    if not line:
        break
    if line[0] != '@':
        fields = line.split('\t')
        qual = fields[10]

        if len(qual) > 0:
            qual_size += len(qual)
            for q in qual:
                if ord(q) > qual_max:
                    qual_max = ord(q)
                if ord(q) < qual_min:
                    qual_min = ord(q)
        else:
            sys.exit("Error: No quality scores in line {}".format(total_line_cnt))

        alignment_line_cnt += 1
    else:
        header_line_cnt += 1
    total_line_cnt += 1

# Print statistics
print("Lines: {}".format(total_line_cnt))
print("  Header lines: {}".format(header_line_cnt))
print("  Alignment lines: {}".format(alignment_line_cnt))
print("Quality score range: [{}:{}]".format(qual_min, qual_max))
print("Quality score size: {} bytes".format(qual_size))

sys.exit()
