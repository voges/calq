#!/usr/bin/env python

###############################################################################
#               Extract the quality values from a FASTQ file                  #
###############################################################################

import sys


if len(sys.argv) != 2:
    sys.exit("Usage: python {} file.[fastq|fq] 2> file.qual".format(sys.argv[0]))

# Get FASTQ file
fastq_file_name = sys.argv[1]
if not fastq_file_name.endswith(".fastq") and not fastq_file_name.endswith(".fq"):
    print("FASTQ file name must end with either '.fastq' or '.fq'")
    sys.exit()
fastq_file = open(fastq_file_name, 'r')
print("FASTQ file: {}".format(fastq_file_name))

record_cnt = 0
line_cnt = 0

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

    sys.stderr.write(quality_scores + "\n")

    sys.stdout.flush()

print("Processed {} FASTQ record(s) (i.e. {} lines)".format(record_cnt, line_cnt))

sys.exit()

