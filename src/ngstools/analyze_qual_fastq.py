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
print("FASTQ file: {}".format(fastq_file_name))

# Initialize statistics
qual_min = sys.maxint
qual_max = -sys.maxint - 1
qual_size = 0
record_cnt = 0
line_cnt = 0

# Parse FASTQ file
while 1:
    # Try to read a whole record (i.e. 4 lines)
    sequence_identifier = fastq_file.readline().rstrip('\n')
    if not sequence_identifier:
        break  # reached end of file, everything ok
    line_cnt += 1

    raw_sequence = fastq_file.readline().rstrip('\n')
    if not raw_sequence:
        sys.exit("Error in record {}, line {}: Record is not complete".format(record_cnt+1, line_cnt+1))
    line_cnt += 1

    description = fastq_file.readline().rstrip('\n')
    if not description:
        sys.exit("Error in record {}, line {}: Record is not complete".format(record_cnt+1, line_cnt+1))
    line_cnt += 1

    quality_scores = fastq_file.readline().rstrip('\n')
    if not quality_scores:
        sys.exit("Error in record {}, line {}: Record is not complete".format(record_cnt+1, line_cnt+1))
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
        qual_size += len(quality_scores)
        for q in quality_scores:
            if ord(q) > qual_max:
                qual_max = ord(q)
            if ord(q) < qual_min:
                qual_min = ord(q)
    else:
        sys.exit("Error in record {}, line {}: Quality score line is empty".format(record_cnt, line_cnt))

print("Line(s): {}".format(line_cnt))
print("FASTQ record(s): {}".format(record_cnt))
print("Quality score range: [{}:{}]".format(qual_min, qual_max))
print("Quality score size: {} bytes".format(qual_size))

sys.exit()
