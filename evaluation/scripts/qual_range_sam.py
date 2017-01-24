#!/usr/bin/env python

###############################################################################
#      Get the minimum and maximum quality values from a SAM file             #
###############################################################################

import os
import sys


if len(sys.argv) != 2:
    sys.exit("Usage: python {} file.sam".format(sys.argv[0]))

sam_file_name = sys.argv[1]
sam_file = open(sam_file_name, 'r')

print "SAM file: {}".format(sam_file_name)

sam_file.seek(0, os.SEEK_END)
sam_file_size = sam_file.tell()
sam_file.seek(0, os.SEEK_SET)

quality_score_min = sys.maxint
quality_score_max = -sys.maxint - 1

line_cnt = 0
header_lines = 0
alignment_lines = 0

while 1:
    print "\r{0:.2f}%".format(100*(float(sam_file.tell())/float(sam_file_size))),
    sys.stdout.flush()

    line = sam_file.readline()
    if not line:
        break
    if line[0] != '@':
        fields = line.split('\t')
        qual = fields[10]

        if len(qual) > 0:
            for q in qual:
                if ord(q) > quality_score_max:
                    quality_score_max = ord(q)
                if ord(q) < quality_score_min:
                    quality_score_min = ord(q)
        else:
            print "Error: No quality scores in line {}".format(line_cnt)

        alignment_lines += 1
    else:
        header_lines += 1
    line_cnt += 1
print ""
print "Processed {} header + {} alignment = {} lines".format(header_lines, alignment_lines, line_cnt)
print "Quality score range: [{}:{}]".format(quality_score_min, quality_score_max)

sys.exit()
