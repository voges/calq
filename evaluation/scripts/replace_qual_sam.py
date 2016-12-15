#!/usr/bin/env python

###############################################################################
#                    Replace quality values in a SAM file                     #
###############################################################################

import os
import sys

if len(sys.argv) != 3:
    sys.exit("Usage: python {} file.sam new.qual".format(sys.argv[0]))

sam_file_name = sys.argv[1]
sam_file = open(sam_file_name, 'r')
new_qual_file_name = sys.argv[2]
new_qual_file = open(new_qual_file_name, 'r')
new_sam_file_name = sam_file_name + '.new_qual.sam'
new_sam_file = open(new_sam_file_name, 'w')

print "SAM file: {}".format(sam_file_name)
print "QUAL file: {}".format(new_qual_file_name)
print "New SAM file: {}".format(new_sam_file_name)

sam_file.seek(0, os.SEEK_END)
sam_file_size = sam_file.tell()
sam_file.seek(0, os.SEEK_SET)

idx = 0
header_lines = 0
alignment_lines = 0

print "Replacing quality values"
while 1:
    print "\r{0:.2f}%".format(100*(float(sam_file.tell())/float(sam_file_size))),
    sys.stdout.flush()

    line = sam_file.readline()
    if not line:
        break
    if line[0] != '@':
        new_qual = new_qual_file.readline()
        a = line.split('\t')
        new_line = a[0] + '\t' + a[1] + '\t' + a[2] + '\t' + a[3] + '\t' + a[4] + '\t' + a[5] + '\t'
        new_line += a[6] + '\t' + a[7] + '\t' + a[8] + '\t' + a[9] + '\t'
        new_qual = new_qual.rstrip('\n')
        new_line += new_qual
        if len(a) > 11:
            for i in range(11, len(a)):
                new_line += '\t' + a[i]
        else:
            new_line += '\n'
        new_sam_file.write(new_line)
        alignment_lines += 1
    else:
        new_sam_file.write(line)
        header_lines += 1
    idx += 1
print ""
print "Processed {} lines ({} header lines, {} alignment lines)".format(idx, header_lines, alignment_lines)

sys.exit()
