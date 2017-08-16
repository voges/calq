#!/usr/bin/env python

###############################################################################
#                            Convert SAM to FASTQ                             #
###############################################################################

import sys


if len(sys.argv) != 3:
    sys.exit("Usage: python {} file.sam file.fastq".format(sys.argv[0]))

sam_file_name = sys.argv[1]
sam_file = open(sam_file_name, 'r')
print("SAM file: {}".format(sam_file_name))

fastq_file_name = sys.argv[2]
fastq_file = open(fastq_file_name, 'w')
print("FASTQ file: {}".format(fastq_file_name))

line_cnt = 0
header_line_cnt = 0
alignment_line_cnt = 0

while 1:
    print("\rLine: {}".format(line_cnt)),
    sys.stdout.flush()

    line = sam_file.readline()
    if not line:
        break
    if line[0] != '@':
        fields = line.split('\t')
        rname = fields[0]
        seq = fields[9]
        qual = fields[10]
        fastq_file.write(rname)
        fastq_file.write(seq)
        fastq_file.write("+")
        fastq_file.write(qual)
        fastq_file.write('\n')
        alignment_line_cnt += 1
    else:
        header_line_cnt += 1
    line_cnt += 1
print("")
print("Processed {} header + {} alignment = {} lines".format(header_line_cnt, alignment_line_cnt, line_cnt))

sys.exit()

