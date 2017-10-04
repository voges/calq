#!/usr/bin/env python

import sys

# Usage
if len(sys.argv) != 3:
    sys.exit("Usage: python {} file.fastq new.qual 1>new.fastq".format(sys.argv[0]))

# Get FASTQ file
fastq_file_name = sys.argv[1]
if not fastq_file_name.endswith(".fastq") and not fastq_file_name.endswith(".fq"):
    sys.exit("FASTQ file name must end with either '.fastq' or '.fq'")
fastq_file = open(fastq_file_name, 'r')

# Get new QUAL file
new_qual_file_name = sys.argv[2]
if not new_qual_file_name.endswith(".qual"):
    sys.exit("QUAL file name must end with '.qual'")
new_qual_file = open(new_qual_file_name, 'r')

while 1:
    sequence_identifier = fastq_file.readline().rstrip('\n')
    if not sequence_identifier:
        break  # reached end of file, everything ok
    raw_sequence = fastq_file.readline().rstrip('\n')
    description = fastq_file.readline().rstrip('\n')
    quality_scores = fastq_file.readline().rstrip('\n')

    new_quality_scores = new_qual_file.readline().rstrip('\n')

    sys.stdout.write(sequence_identifier + "\n")
    sys.stdout.write(raw_sequence + "\n")
    sys.stdout.write(description + "\n")
    sys.stdout.write(new_quality_scores + "\n")

sys.exit()
