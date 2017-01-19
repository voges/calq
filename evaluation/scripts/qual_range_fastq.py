#!/usr/bin/env python

###############################################################################
#      Get the minimum and maximum quality values from a FASTQ file           #
###############################################################################

import os
import sys


if len(sys.argv) != 2:
    sys.exit("Usage: python {} file.[fastq|fq]".format(sys.argv[0]))

# Get FASTQ file
fastq_file_name = sys.argv[1]
if not fastq_file_name.endswith(".fastq") and not fastq_file_name.endswith(".fq"):
    sys.exit("FASTQ file name must end with either '.fastq' or '.fq'")
fastq_file = open(fastq_file_name, 'r')
print "FASTQ file: {}".format(fastq_file_name)

fastq_file.seek(0, os.SEEK_END)
fastq_file_size = fastq_file.tell()
fastq_file.seek(0, os.SEEK_SET)

quality_score_min = sys.maxint
quality_score_max = -sys.maxint - 1

record_cnt = 0
line_cnt = 0

while 1:
    print "\r{0:.2f}%".format(100 * (float(fastq_file.tell()) / float(fastq_file_size))),
    sys.stdout.flush()

    # Try to read a whole record (i.e. 4 lines)
    sequence_identifier = fastq_file.readline().rstrip('\n')
    if not sequence_identifier:
        break  # reached end of file, everything ok
    line_cnt += 1

    raw_sequence = fastq_file.readline().rstrip('\n')
    if not raw_sequence:
        print "Parse error (record {}, line {}): Record is not complete".format(record_cnt+1, line_cnt+1)
        break
    line_cnt += 1

    description = fastq_file.readline().rstrip('\n')
    if not description:
        print "Parse error (record {}, line {}): Record is not complete".format(record_cnt + 1, line_cnt + 1)
        break
    line_cnt += 1

    quality_scores = fastq_file.readline().rstrip('\n')
    if not quality_scores:
        print "Parse error (record {}, line {}): Record is not complete".format(record_cnt + 1, line_cnt + 1)
        break
    line_cnt += 1

    record_cnt += 1

    # Check quality score line

    # Typical ranges:
    #  Sanger         Phred+33   [0,40]
    #  Solexa         Solexa+64  [-5,40]
    #  Illumina 1.3+  Phred+64   [0,40]
    #  Illumina 1.5+  Phred+64   [0,40] with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator ('B')
    #  Illumina 1.8+  Phred+33   [0,41]

    if len(quality_scores) > 0:
        for q in quality_scores:
            if ord(q) > quality_score_max:
                quality_score_max = ord(q)
            if ord(q) < quality_score_min:
                quality_score_min = ord(q)
    else:
        print "Error (record {}, line {}): Quality score line is empty".format(record_cnt, line_cnt)
print ""
print "Processed {} FASTQ record(s) (i.e. {} lines)".format(record_cnt, line_cnt)
print "Quality score range: [{}:{}]".format(quality_score_min, quality_score_max)

sys.exit()
