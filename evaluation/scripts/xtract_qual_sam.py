#!/usr/bin/env python

###############################################################################
#                 Extract the quality values from a SAM file                  #
###############################################################################

import sys

if len(sys.argv) != 2:
    sys.exit("Usage: python {} file.sam 2> file.qual".format(sys.argv[0]))

sam_file_name = sys.argv[1]
sam_file = open(sam_file_name, 'r')

print("SAM file: {}".format(sam_file_name))

line_cnt = 0
header_line_cnt = 0
alignment_line_cnt = 0

while 1:
    print "\rLine: {}".format(line_cnt),
    sys.stdout.flush()

    line = sam_file.readline()
    if not line:
        break
    if line[0] != '@':
        fields = line.split('\t')
        qual = fields[10]
        sys.stderr.write(qual)
        sys.stderr.write('\n')
        alignment_line_cnt += 1
    else:
        header_line_cnt += 1
    line_cnt += 1
print ""
print "Processed {} header + {} alignment = {} lines".format(header_line_cnt, alignment_line_cnt, line_cnt)
