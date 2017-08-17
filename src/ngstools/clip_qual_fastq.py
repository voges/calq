#!/usr/bin/env python

import sys

# Usage
if len(sys.argv) != 2:
    sys.exit("Usage: python {} file.[fastq|fq]".format(sys.argv[0]))

# Get FASTQ file
fastq_file_name = sys.argv[1]
if not fastq_file_name.endswith(".fastq") and not fastq_file_name.endswith(".fq"):
    sys.exit("FASTQ file name must end with either '.fastq' or '.fq'")
fastq_file = open(fastq_file_name, 'r')

sys.stderr.write("Clipping quality values into Illumina 1.8+ format ([0,41]+33)\n")

# Initialize statistics
record_cnt = 0
line_cnt = 0
qual_total_cnt = 0
qual_raise_cnt = 0
qual_decrease_cnt = 0

# Parse FASTQ file and clip the quality scores
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

    sys.stdout.write(sequence_identifier + "\n")
    sys.stdout.write(raw_sequence + "\n")
    sys.stdout.write(description + "\n")

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
            sys.stdout.write(q)

    sys.stdout.write("\n")

sys.stderr.write("Clipped {} out of {} quality values\n".format(qual_raise_cnt+qual_decrease_cnt, qual_total_cnt))
sys.stderr.write("  Raised: {}\n".format(qual_raise_cnt))
sys.stderr.write("  Decreased: {}\n".format(qual_decrease_cnt))

sys.exit()
