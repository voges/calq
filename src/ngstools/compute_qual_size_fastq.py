#!/usr/bin/env python

###############################################################################
#             Compute the quality values size in a FASTQ file                 #
###############################################################################

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

record_cnt = 0
line_cnt = 0
quality_scores_size = 0

while 1:
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

    quality_scores_size += len(quality_scores)

    sys.stdout.flush()

print("Processed {} FASTQ record(s) (i.e. {} lines)".format(record_cnt, line_cnt))
print("Quality values size: {} bytes".format(quality_scores_size))

sys.exit()

