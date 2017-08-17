#!/usr/bin/env python

import sys

# Usage
if len(sys.argv) != 3:
    sys.exit("Usage: python {} file.[fastq|fq] part_id".format(sys.argv[0]))

# Get FASTQ file
fastq_file_name = sys.argv[1]
if not fastq_file_name.endswith(".fastq") and not fastq_file_name.endswith(".fq"):
    sys.exit("FASTQ file name must end with either '.fastq' or '.fq'")
fastq_file = open(fastq_file_name, 'r')

# Get part ID
part_id = int(sys.argv[2])
if not 0 <= part_id <= 3:
    sys.exit("part_id must be in range [0,3]")

# Write part to stdout
while 1:
    sequence_identifier = fastq_file.readline().rstrip('\n')
    if not sequence_identifier:
        break  # reached end of file, everything ok
    raw_sequence = fastq_file.readline().rstrip('\n')
    description = fastq_file.readline().rstrip('\n')
    quality_scores = fastq_file.readline().rstrip('\n')

    if part_id == 0:
        sys.stdout.write(sequence_identifier + "\n")
    if part_id == 1:
        sys.stdout.write(raw_sequence + "\n")
    if part_id == 2:
        sys.stdout.write(description + "\n")
    if part_id == 3:
        sys.stdout.write(quality_scores + "\n")

sys.exit()
