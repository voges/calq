#!/usr/bin/env python

###############################################################################
#                  Clip the quality values in a FASTQ file                    #
###############################################################################

import os
import sys


if len(sys.argv) != 2:
    sys.exit("Usage: python {} file.[fastq|fq]".format(sys.argv[0]))

# Get FASTQ file
fastq_file_name = sys.argv[1]
if not fastq_file_name.endswith(".fastq") and not fastq_file_name.endswith(".fq"):
    print("FASTQ file name must end with either '.fastq' or '.fq'")
    sys.exit()
fastq_file = open(fastq_file_name, 'r')
print("FASTQ file: {}".format(fastq_file_name))

# Get FASTQ file size
fastq_file.seek(0, os.SEEK_END)
fastq_file_size = fastq_file.tell()
fastq_file.seek(0, os.SEEK_SET)

# Make new output FASTQ file
clipped_fastq_file_name = fastq_file_name + ".clipped_qual.fastq"
clipped_fastq_file = open(clipped_fastq_file_name, 'w')
print("Clipped FASTQ file: {}".format(clipped_fastq_file_name))

print("Clipping quality values into Illumina 1.8+ format ([0,41]+33)")

record_cnt = 0
line_cnt = 0

qual_total_cnt = 0
qual_raise_cnt = 0
qual_decrease_cnt = 0

while 1:
    print("\r{0:.2f}%".format(100*(float(fastq_file.tell())/float(fastq_file_size)))),
    sys.stdout.flush()

    sequence_identifier = fastq_file.readline().rstrip('\n')
    if not sequence_identifier:
        break  # reached end of file, everything ok
    line_cnt += 1

    raw_sequence = fastq_file.readline().rstrip('\n')
    line_cnt += 1

    description = fastq_file.readline().rstrip('\n')
    line_cnt += 1

    quality_scores = fastq_file.readline().rstrip('\n')
    line_cnt += 1

    record_cnt += 1

    clipped_fastq_file.write(sequence_identifier + "\n")
    clipped_fastq_file.write(raw_sequence + "\n")
    clipped_fastq_file.write(description + "\n")

    # Check and clip quality scores

    # Typical ranges:
    #  Sanger         Phred+33   [0,40]
    #  Solexa         Solexa+64  [-5,40]
    #  Illumina 1.3+  Phred+64   [0,40]
    #  Illumina 1.5+  Phred+64   [0,40] with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator ('B')
    #  Illumina 1.8+  Phred+33   [0,41]

    # Currently using Illumina 1.8+
    for q in quality_scores:
        qual_total_cnt += 1
        if ord(q) < 33:
            q = chr(33)
            qual_raise_cnt += 1
        if ord(q) > 33+41:
            q = chr(33+41)
            qual_decrease_cnt += 1
        clipped_fastq_file.write(q)

    clipped_fastq_file.write("\n")

print("")
print("Processed {} FASTQ record(s) (i.e. {} lines)".format(record_cnt, line_cnt))
print("Clipped {} out of {} quality values ({} raised + {} decreased)".format(qual_raise_cnt+qual_decrease_cnt, qual_total_cnt, qual_raise_cnt, qual_decrease_cnt))

sys.exit()

